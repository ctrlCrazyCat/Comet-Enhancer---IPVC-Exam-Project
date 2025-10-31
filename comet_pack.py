import os
import numpy as np
import cv2
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RangeSlider,Button
from matplotlib.cm import get_cmap
import tkinter as tk
from tkinter import filedialog as fd
TWOPI = 2*np.pi

class Params:
    def __init__(self):
        return

OPTIONS = [
            'Division by Azimuthal Average',
            'Division by Azimuthal Median',
            'Azimuthal Renormalization',
            'Division by 1/rho profile',
            'Radially Variable Spatial Filtering'
        ]
def get_options():
    return OPTIONS
    
def get_input_data(filepath):
    hdul = fits.open(filepath)
    hdu=hdul[0]
    imold = hdu.data
    imold[np.isnan(imold)] = np.min(imold[~np.isnan(imold)])
    imold=imold.astype(float)
    imold = 100*imold/np.mean(imold)
    
    return hdul,imold



def rho_division(imold,xnuc:float,ynuc:float,xmin:int,xmax:int,ymin:int,ymax:int):
    (NROW,NCOL) = imold.shape
    linspace_vec = 0.1*np.linspace(1,10,10)
    i_coords = np.arange(NROW).reshape(-1, 1) # (lcol, 1) per broadcasting
    ynew_vec_all = np.square(ynuc - i_coords - 0.55 + linspace_vec) # (lcol, 10)
    j_coords = np.arange(NCOL).reshape(-1, 1) # (lrow, 1) per broadcasting
    xnew_vec_all = np.square(xnuc - j_coords - 0.55 + linspace_vec) # (lrow, 10)
    ynew_vec_expanded = ynew_vec_all[:, np.newaxis, :, np.newaxis] # (lcol, 1, 10, 1)
    xnew_vec_expanded = xnew_vec_all[np.newaxis, :, np.newaxis, :] # (1, lrow, 1, 10)
    sum_of_squares_all = xnew_vec_expanded + ynew_vec_expanded
    sqrt_values_all = np.sqrt(sum_of_squares_all)
    rho_all = np.sum(sqrt_values_all, axis=(-2, -1)) # Somma sulle ultime due dimensioni
    rho_all = rho_all * 0.01
    imn = imold * rho_all
    return imn[ymin:ymax-1,xmin:xmax-1]



def azimuthal_median_division(imiun):
    (nrad,ntheta)=imiun.shape
    jjj_all = np.sum(imiun < 0.0, axis=1)
    sorted_bvect_all = np.sort(imiun, axis=1)
    nmid_fortran_style_all = ((ntheta + jjj_all) // 2).astype(int)
    median_idx_all = np.maximum(0, nmid_fortran_style_all - 1)

    median_all = sorted_bvect_all[np.arange(nrad), median_idx_all]
    imien = imiun / (median_all[:, np.newaxis] + 1.0e-6)
    return imien


def azimuthal_average_division(imiun,rejsig):
    mask_pass1 = (imiun >= 0.0)
    jj_all = np.sum(mask_pass1, axis=1)
    sum_pass1 = np.sum(imiun * mask_pass1, axis=1)
    mean_pass1 = np.where(jj_all > 0, sum_pass1 / jj_all, 0.0)
    diff_sq_pass1 = np.where(mask_pass1, np.square(imiun - mean_pass1[:, np.newaxis]), 0.0)
    sigma_sum_pass1 = np.sum(diff_sq_pass1, axis=1)
    sigma_pass1 = np.where(
        jj_all > 1,
        np.sqrt(sigma_sum_pass1 / (jj_all - 1)),
        0.0 
    )
    rej_all = rejsig * sigma_pass1
    rejmin_all = np.maximum(1.0e-5, mean_pass1 - rej_all)
    rejmax_all = mean_pass1 + rej_all
    mask_pass2 = (imiun >= rejmin_all[:, np.newaxis]) & \
                (imiun <= rejmax_all[:, np.newaxis])
    kk_all = np.sum(mask_pass2, axis=1)
    sum_pass2 = np.sum(imiun * mask_pass2, axis=1)
    mean_final = np.where(kk_all > 0, sum_pass2 / kk_all, 0.0)
    imien=imiun / (mean_final[:, np.newaxis] + 1.0e-6)
    return imien


def azimuthal_renormalization(imiun,rejsig,nsig):
    
    (nrad,ntheta)=imiun.shape
    imien = np.zeros((nrad, ntheta))
    # Ciclo esterno su 'i' (nrad) rimane, dato che le statistiche vengono calcolate per ogni riga indipendentemente
    for i in range(nrad):
        avect = imiun[i]

        # --- Primo blocco di calcolo di media e sigma (filtering avect[j] >= 0.0) ---
        # Maschera per i valori validi (avect[j] >= 0.0)
        valid_mask_1 = (avect >= 0.0)
        valid_values_1 = avect[valid_mask_1]

        # Calcolo di jj, mean e sigma in modo vettorizzato
        jj = len(valid_values_1) # Equivalente a np.sum(valid_mask_1)

        if jj > 0:
            mean_1 = np.mean(valid_values_1)
            # np.std di default calcola la deviazione standard con N. Per N-1 (unbiased) usa ddof=1
            sigma_1 = np.std(valid_values_1, ddof=1) if jj > 1 else 0.0
        else:
            mean_1 = 0.0 # Se non ci sono valori validi, la media è 0
            sigma_1 = 0.0 # E la deviazione standard è 0

        # Calcolo di rej, rejmin, rejmax
        rej = rejsig * sigma_1
        rejmin = max(1.0e-5, (mean_1 - rej))
        rejmax = mean_1 + rej

        # --- Secondo blocco di calcolo di mean e sigma (filtering avect[j] >= rejmin and avect[j] <= rejmax) ---
        # Maschera per i valori validi in base a rejmin e rejmax
        valid_mask_2 = (avect >= rejmin) & (avect <= rejmax)
        valid_values_2 = avect[valid_mask_2]

        # Calcolo di kk, mean e sigma in modo vettorizzato
        kk = len(valid_values_2)

        if kk > 0:
            mean_2 = np.mean(valid_values_2)
            sigma_2 = np.std(valid_values_2, ddof=1) if kk > 1 else 0.0
        else:
            mean_2 = 0.0
            sigma_2 = 0.0

        # Calcolo di rej, rowmin, rowmax (secondo set di valori)
        rej_final = nsig * sigma_2
        rowmin = max(1.0e-5, (mean_2 - rej_final))
        rowmax = mean_2 + rej_final

        # Debug print come nell'originale (se necessario)
        # print(i, float(i)*rowmin, float(i)*rowmax, rowmin, rowmax)

        # --- Calcolo di bvect e assegnazione a imien ---
        diff = rowmax - rowmin
        
        # Vettorizzazione del calcolo di bvect
        if diff >= 1.0e-5:
            # Se diff è sufficientemente grande, normalizza
            bvect = (avect - rowmin) * 255.0 / diff
            # Clipa i valori di bvect per essere nel range [0, 255] se necessario
            # (Il tuo codice originale non lo faceva, ma è comune in normalizzazione di immagini)
            #bvect = np.clip(bvect, 0.0, 255.0)
        else:
            # Se diff è troppo piccolo, tutti i valori di bvect dovrebbero essere 0 (o un altro valore di default)
            bvect = np.zeros(ntheta)

        imien[i] = bvect
    return imien



def reconstruct_from_polar(imien,NCOL:int,NROW:int,xnuc:float,ynuc:float,xmin:int,xmax:int,ymin:int,ymax:int):
    
    (nrad,ntheta)=imien.shape

    
    imn = np.zeros((NCOL, NROW))

    
    # Pre-calculate constants
    tenth_ii = np.arange(1, 11) * 0.1
    tenth_jj = np.arange(1, 11) * 0.1

    # Creazione delle griglie di coordinate in modo vettorizzato
    # Invece di creare meshgrid all'interno del ciclo, creiamo qui le griglie complete
    # per tutte le iterazioni i e j.

    # Vettorizzazione di i e j per le coordinate
    i_coords_all = np.arange(NCOL).reshape(-1, 1, 1) # (lcol, 1, 1)
    j_coords_all = np.arange(NCOL).reshape(1, -1, 1) # (1, lrow, 1)

    # Espandiamo tenth_ii e tenth_jj per broadcasting
    tenth_ii_reshaped = tenth_ii.reshape(1, 1, -1) # (1, 1, 10)
    tenth_jj_reshaped = tenth_jj.reshape(1, 1, -1) # (1, 1, 10)

    # Calcolo vettorizzato di xdist e ydist per tutti i pixel e tutte le sotto-posizioni (tenth_ii, tenth_jj)
    # xdist_grid sarà (lcol, lrow, 10)
    xdist_grid_all = (j_coords_all + tenth_jj_reshaped - 0.55 - xnuc)
    # ydist_grid sarà (lcol, lrow, 10)
    ydist_grid_all = (i_coords_all + tenth_ii_reshaped - 0.55 - ynuc)


    # Calcolo vettorizzato di rpix_squared
    rpix_squared_all = np.square(xdist_grid_all) + np.square(ydist_grid_all)
    rpix_all = np.round(np.sqrt(rpix_squared_all)).astype(int)

    # Filtro vettorizzato per valori validi di rpix
    valid_mask_all = (rpix_all >= 1) & (rpix_all <= nrad)

    # Calcolo vettorizzato di angthe
    angthe_all = np.arctan2(-xdist_grid_all, ydist_grid_all)

    # Regolazione vettorizzata di angthe per essere in [0, TWOPI)
    angthe_all = np.where(angthe_all < 0.0, angthe_all + TWOPI, angthe_all)
    angthe_all = np.where(angthe_all >= TWOPI, angthe_all - TWOPI, angthe_all)

    # Calcolo vettorizzato di thpix
    thpix_all = np.round(0.5 + angthe_all * (ntheta / TWOPI)).astype(int)

    # Clipping vettorizzato di thpix
    thpix_all = np.clip(thpix_all, 0, ntheta - 1)

    # Preparazione degli indici per l'accesso a imien
    # Dobbiamo considerare che rpix_all va da 1 a nrad, mentre gli indici di imien sono 0-indexed
    rpix_indices = np.clip(rpix_all - 1, 0, nrad - 1)

    # Utilizzo di advanced indexing per ottenere i valori da imien
    # Creiamo un array temporaneo che conterrà i valori di imien per tutti i punti
    # dove valid_mask_all è True, altrimenti 0.
    # Dobbiamo estrarre i valori da imien usando gli indici calcolati
    # e poi applicare la maschera.

    # Flattening degli indici per l'advanced indexing su imien
    # Creiamo una versione "flattata" di rpix_indices e thpix_all che include solo i punti validi
    # Questo sarà un approccio più efficiente.
    # Per ogni punto (i, j, k) dove k è l'indice delle sotto-posizioni (tenth_ii/jj),
    # vogliamo accedere a imien[rpix_indices[i, j, k], thpix_all[i, j, k]]

    # Prima, applichiamo la maschera ai nostri indici
    rpix_masked = rpix_indices[valid_mask_all]
    thpix_masked = thpix_all[valid_mask_all]

    # Ora, accediamo a imien usando gli indici mascherati
    # Questo creerà un array 1D di tutti i valori di imien che sono validi
    imien_values = imien[rpix_masked, thpix_masked]

    # Creiamo un array di zeri delle stesse dimensioni di rpix_all
    # e ci mettiamo i valori calcolati solo dove la maschera è True
    # Questo ci permette di sommare correttamente per ciascun pixel (i, j)
    temp_sum_array = np.zeros(rpix_all.shape)
    temp_sum_array[valid_mask_all] = imien_values

    # Sommiamo lungo l'asse delle sotto-posizioni (l'ultima dimensione)
    # per ottenere la somma per ogni pixel (i, j)
    cvect_sum_all = np.sum(temp_sum_array, axis=2) # Somma su k (le 10 sotto-posizioni)

    # Assegnazione finale a imn
    # cvect_sum_all ha la forma (lcol, lrow), che è esattamente quella che vogliamo per imn
    imn = cvect_sum_all * 0.01
    return imn[ymin:ymax-1,xmin:xmax-1]




def polarize(imold,nrad,ntheta,xnuc:float,ynuc:float):
    (NROW,NCOL) = imold.shape
    print("row",NROW,"col",NCOL)
    mtheta=ntheta*10
    fmthet=float(mtheta)
    angle = (np.arange(mtheta + 1, dtype=float) - 0.5) * (TWOPI / fmthet)
    angcos=-np.sin(angle)
    angsin=np.cos(angle)
    imiun = np.zeros((nrad, ntheta))

    linspace_vec1 = 0.1 * np.linspace(1, 10, 10) # Array di 10 elementi: [0.1, 0.2, ..., 1.0]
    linspace_vec2 = np.arange(9) # Array di 9 elementi: [0, 1, ..., 8]

    # Vettorizzazione completa
    # Creazione di indici per i e j per broadcasting
    i_coords = np.arange(nrad).reshape(-1, 1) # (nrad, 1) per broadcasting con linspace_vec1
    j_coords_base = np.arange(ntheta).reshape(1, -1) # (1, ntheta) per broadcasting con linspace_vec2

    # Calcolo di ai_v per tutte le 'i' e linspace_vec1
    # Forma desiderata: (nrad, 10)
    ai_v_all = (i_coords - 0.55 + linspace_vec1) # (nrad, 10)

    # Calcolo di jnew_v per tutte le 'j' e linspace_vec2
    # Forma desiderata: (ntheta, 9)
    jnew_v_all_base = (j_coords_base.T * 10 + linspace_vec2) # .T per (ntheta, 1) * (9,) -> (ntheta, 9)

    # Clipa gli indici per evitare errori se jnew_v_all va fuori range di angcos/angsin
    jnew_v_clipped = np.clip(jnew_v_all_base, 0, len(angcos) - 1).astype(int)

    # Estendere angcos e angsin per tutti i jnew_v
    # Forma desiderata: (ntheta, 9)
    angcos_vals = angcos[jnew_v_clipped]
    angsin_vals = angsin[jnew_v_clipped]


    # --- Ora, prepariamo per la moltiplicazione "outer" vettorizzata ---
    # Per la moltiplicazione (ai_v * angcos_vals), vogliamo ottenere una forma (nrad, ntheta, 10, 9)
    # ai_v_all ha forma (nrad, 10)
    # angcos_vals ha forma (ntheta, 9)

    # Espandiamo ai_v_all per essere (nrad, 1, 10, 1)
    ai_v_expanded = ai_v_all[:, np.newaxis, :, np.newaxis] # (nrad, 1, 10, 1)

    # Espandiamo angcos_vals per essere (1, ntheta, 1, 9)
    angcos_expanded_for_mul = angcos_vals[np.newaxis, :, np.newaxis, :] # (1, ntheta, 1, 9)
    angsin_expanded_for_mul = angsin_vals[np.newaxis, :, np.newaxis, :] # (1, ntheta, 1, 9)

    # Calcolo vettorizzato di xpix_v e ypix_v
    # Il risultato sarà (nrad, ntheta, 10, 9)
    xpix_vals_unrounded = xnuc + (ai_v_expanded * angcos_expanded_for_mul)
    ypix_vals_unrounded = ynuc + (ai_v_expanded * angsin_expanded_for_mul)

    # Arrotondamento e conversione a int32
    xpix_all = np.round(xpix_vals_unrounded).astype(np.int32)
    ypix_all = np.round(ypix_vals_unrounded).astype(np.int32)

    # --- Gestione delle condizioni di "out of bounds" in modo vettorizzato ---
    # Crea una maschera per i valori validi (in-bounds)
    valid_coords_mask = (xpix_all >= 0) & (xpix_all < NCOL) & \
                        (ypix_all >= 0) & (ypix_all < NROW)

    # Vogliamo sapere per ogni (i, j) se TUTTI i 90 sotto-punti sono validi.
    all_subpoints_valid_per_ij = np.all(valid_coords_mask, axis=(-2, -1))

    # Inizializza un array per i risultati di `aarray[ypix_v, xpix_v]`
    aarray_values = np.zeros(xpix_all.shape, dtype=imold.dtype)

    # Usa advanced indexing per riempire aarray_values solo dove valid_coords_mask è True
    flat_valid_coords_mask = valid_coords_mask.flatten()
    flat_xpix_all = xpix_all.flatten()
    flat_ypix_all = ypix_all.flatten()

    valid_flat_xpix = flat_xpix_all[flat_valid_coords_mask]
    valid_flat_ypix = flat_ypix_all[flat_valid_coords_mask]

    values_from_aarray = imold[valid_flat_ypix, valid_flat_xpix]

    temp_flat_aarray_values = np.zeros(xpix_all.size, dtype=imold.dtype)
    temp_flat_aarray_values[flat_valid_coords_mask] = values_from_aarray
    aarray_values_reshaped = temp_flat_aarray_values.reshape(xpix_all.shape)

    # Somma i valori di aarray_values_reshaped per ottenere il risultato per ogni (i, j)
    sum_of_aarray_for_ij = np.sum(aarray_values_reshaped, axis=(-2, -1))

    # Applica il moltiplicatore 0.01
    calculated_avect_values = 0.01 * sum_of_aarray_for_ij

    # Inizialmente tutti i valori calcolati
    imiun = calculated_avect_values

    # Se almeno un sotto-punto non è valido, imposta l'intero blocco (i, j) a -1.0
    imiun[~all_subpoints_valid_per_ij] = -1.0

    # Gestione della condizione avect[j] < 1e-5
    small_value_mask = (imiun < 1e-5) & (imiun != -1.0)
    imiun[small_value_mask] = -1.0e-5
    return imiun

def radially_variable_spatial_filtering(imold,A,B,N,NUMLOG,xnuc,ynuc,xmin,xmax,ymin,ymax):
    imavg=np.mean(imold)
    imold=imold*100.0/imavg
    (numrows,numcols)=imold.shape
    if NUMLOG:
        if imold[imold<=1e-15].shape[0] ==0:
            imold=np.log10(imold)
        else:
            NUMLOG=False     #SERVE PER SCRIVERLO NEL FITS
    imn = np.zeros_like(imold)

    # Pre-calcolo delle coordinate sub-pixel m e n
    m_coords = -5.5e-1 + (np.arange(10) + 1) * 1.0e-1
    n_coords = -5.5e-1 + (np.arange(10) + 1) * 1.0e-1

    # Creiamo un meshgrid di am e an per calcolare rho e a0 per tutti i sub-pixel (10x10)
    am_grid_subpixel, an_grid_subpixel = np.meshgrid(m_coords, n_coords)

    # Preparazione per gli offset base per edge e corner. Questi non dipendono da i o j.
    # Offsets per gli edge pixels: (+/-1, 0) o (0, +/-1)
    mult_i_edge = np.array([-1, 1, 0, 0])
    mult_j_edge = np.array([0, 0, -1, 1])

    # Offsets per i corner pixels: (+/-1, +/-1)
    mult_i_crn = np.array([-1, -1, 1, 1])
    mult_j_crn = np.array([-1, 1, -1, 1])
    for i in range(numrows):
        
        # Calcolo di a0 vettorizzato per la riga 'i' corrente
        # J_coords sarà una riga (1, xlim) per tutti i j
        J_coords = np.arange(numcols).reshape(1, numcols, 1, 1) # Reshape per broadcasting con sub-pixel
        
        # Le griglie sub-pixel devono essere espanse per broadcasting con J_coords
        am_expanded = am_grid_subpixel[np.newaxis, np.newaxis, :, :] # (1, 1, 10, 10)
        an_expanded = an_grid_subpixel[np.newaxis, np.newaxis, :, :] # (1, 1, 10, 10)

        # rho avrà forma (1, xlim, 10, 10)
        rho = np.sqrt(np.square(J_coords - xnuc + am_expanded) + np.square(float(i) - ynuc + an_expanded))
        a0 = A + (B * (np.power(rho, N))) # a0 ha forma (1, xlim, 10, 10)
        
        # Arrotonda e appiattisci a0 per gli offset
        offsets_a0_flat = np.round(a0).astype(int).flatten() # Forma (xlim * 10 * 10,)

        # Prepara gli indici originali 'j' per la mappatura dei risultati
        # Ogni 'j' si ripete 10*10 volte (per i sub-pixel m,n)
        j_base_flat = np.repeat(np.arange(numcols), 10 * 10)

        # --- Calcolo di sumedg per l'intera riga 'i' ---
        # Espandiamo offsets_a0_flat e i moltiplicatori per poterli combinare
        repeated_offsets_edge = np.repeat(offsets_a0_flat, len(mult_i_edge))
        
        tiled_mult_i_edge = np.tile(mult_i_edge, len(offsets_a0_flat))
        tiled_mult_j_edge = np.tile(mult_j_edge, len(offsets_a0_flat))

        all_delta_i_edge = repeated_offsets_edge * tiled_mult_i_edge
        all_delta_j_edge = repeated_offsets_edge * tiled_mult_j_edge

        # iii e jjj sono relativi alla riga 'i' corrente
        iii_edge = i + all_delta_i_edge
        
        j_base_flat_repeated = np.repeat(j_base_flat, len(mult_j_edge)) # Deve corrispondere a all_delta_j_edge
        jjj_edge = j_base_flat_repeated + all_delta_j_edge

        # Filtra gli indici che sono all'interno dei limiti dell'immagine
        valid_mask_edge = (iii_edge >= 0) & (iii_edge < numrows) & \
                            (jjj_edge >= 0) & (jjj_edge < numcols) & \
                            (((all_delta_i_edge == 0) & (all_delta_j_edge != 0)) | \
                            ((all_delta_j_edge == 0) & (all_delta_i_edge != 0)))

        valid_imold_values_edge = imold[iii_edge[valid_mask_edge], jjj_edge[valid_mask_edge]]
        
        # Mappa i contributi ai pixel (i,j) originali nella riga corrente
        target_j_flat_edge = j_base_flat_repeated[valid_mask_edge]
        
        # Usiamo bincount per sommare i contributi per ogni 'j' nella riga 'i'
        summed_contributions_edge = np.bincount(target_j_flat_edge, weights=valid_imold_values_edge, minlength=numcols)
        total_sumedg_row = summed_contributions_edge / 100.0

        # --- Calcolo di sumcrn per l'intera riga 'i' ---
        repeated_offsets_crn = np.repeat(offsets_a0_flat, len(mult_i_crn))
        
        tiled_mult_i_crn = np.tile(mult_i_crn, len(offsets_a0_flat))
        tiled_mult_j_crn = np.tile(mult_j_crn, len(offsets_a0_flat))

        all_delta_i_crn = repeated_offsets_crn * tiled_mult_i_crn
        all_delta_j_crn = repeated_offsets_crn * tiled_mult_j_crn

        iii_crn = i + all_delta_i_crn
        jjj_crn = j_base_flat_repeated + all_delta_j_crn # j_base_flat_repeated è lo stesso di prima

        valid_mask_crn = (iii_crn >= 0) & (iii_crn < numrows) & \
                            (jjj_crn >= 0) & (jjj_crn < numcols) & \
                            (all_delta_i_crn != 0) & (all_delta_j_crn != 0)

        valid_imold_values_crn = imold[iii_crn[valid_mask_crn], jjj_crn[valid_mask_crn]]
        
        target_j_flat_crn = j_base_flat_repeated[valid_mask_crn]

        summed_contributions_crn = np.bincount(target_j_flat_crn, weights=valid_imold_values_crn, minlength=numcols)
        total_sumcrn_row = summed_contributions_crn * 1.0e-2

        # Calcolo finale di imn per la riga 'i' corrente
        imn[i, :] = 1024.0 * imold[i, :] - 192.0 * total_sumedg_row - 64.0 * total_sumcrn_row
    return imn[ymin:ymax-1,xmin:xmax-1]





#TODO: GESTIRE CASO IMMAGINE MONOCOLORE (DATAMAX == DATAMIN)
def interactive_image_viewer(p:Params, gamma_step=0.05):
    """
    Mostra un'immagine astronomica in modo interattivo con controlli per lo stretching 
    dei livelli (RangeSlider verticale) e la correzione Gamma (scroll del mouse).

    Args:
        image_data (np.ndarray): Array NumPy (2D) contenente i dati dell'immagine.
        gamma_step (float): Il passo di incremento/decremento della correzione Gamma.
    """
    
    # Rimuovi i valori NaN e assicurati che i dati siano float32
    data = cv2.flip(p.imn.astype(np.float32),0)
    data[np.isnan(data)] = np.min(data[~np.isnan(data)])
    
    # --- Variabili di Stato Locali ---
    # Usiamo un wrapper per tenere traccia dello stato di Gamma
    class State:
        gamma_val = 1.0
        img_data = data
    
    state = State()
    

    # --- Calcolo Range Iniziale ---
    data_min = np.min(data)
    data_max = np.max(data)

    print('data_min',data_min,'data_max',data_max)
    

    # --- Funzioni di Elaborazione (Nested) ---
    if not(data_max>data_min):
        print('AAAA')
    
    def normalizza_immagine( min_val, max_val):
        """Esegue il clipping e lo stretching lineare [min_val, max_val] -> [0.0, 1.0]."""
        
        img_clippata = np.clip(data, min_val, max_val)
        range_attuale = max_val - min_val

        if range_attuale <= 0:
            return np.zeros_like(data, dtype=np.float32)
        
        # Stretching lineare a 0.0-1.0 float
        return (img_clippata - min_val) / range_attuale 
    

    


    def applica_gamma(img_normalized, gamma):
        """Applica la correzione Gamma (I_out = I_in ^ gamma)."""
        gamma = max(0.001, gamma)
        # Applica la formula: I_out = I_in ^ gamma
        return np.power(img_normalized, gamma)

    

    # --- Configurazione della Figura e degli Assi ---
    fig, ax = plt.subplots(figsize=(8, 8))
    # Lascia spazio a destra per lo slider verticale
    plt.subplots_adjust(left=0.05, right=0.85, bottom=0.1) 

    # 1. Visualizzazione Iniziale (Usa dati iniziali, verrà aggiornata dalla prima chiamata a update_image)
    im = ax.imshow(np.zeros_like(state.img_data, dtype=np.float32), cmap='gray', vmin=0, vmax=1)
    ax.set_title(f"Immagine FITS (Gamma: {state.gamma_val:.2f})")
    ax.axis('off')

    # 2. Creazione del Range Slider Verticale (Min/Max)
    ax_range_slider = plt.axes([0.9, 0.1, 0.03, 0.8]) # Posizionato a destra
    ax_range_slider.set_ylim(data_min, data_max)
    cmap_gradient = get_cmap('binary') 
    gradient_data = np.linspace(0, 1, 100).reshape(-1, 1)

    # Aggiunge il gradiente come sfondo della traccia
    ax_range_slider.imshow(
        gradient_data,
        aspect='auto',
        cmap=cmap_gradient,
        origin='lower',
        extent=[0, 1, data_min, data_max]
    )
    ax_range_slider.set_xticks([]); ax_range_slider.set_yticks([])
    for spine in ax_range_slider.spines.values():
        spine.set_visible(False)
    ax_range_slider.set_facecolor('none')
    slider_range = RangeSlider(
        ax_range_slider, 
        'Livelli Min/Max', 
        data_min, 
        data_max,
        valinit=(data_min,data_max),
        valstep=(data_max-data_min)/1000,
        orientation="vertical",
        track_color='none' # Colore della traccia del cursore (non della barra intera)
    )
    print(data_min, data_max)
    print("SL VALS",slider_range.val)
    
    # Nasconde la barra colorata tra i manici
    slider_range.poly.set_facecolor('none')
    slider_range.poly.set_edgecolor('none')


    def update_image(val):
        """Aggiorna l'immagine quando lo slider Min/Max viene mosso."""
        
        # Ottiene i valori Min e Max dal RangeSlider
        min_val, max_val = slider_range.val
        
        # Esegue lo stretching
        img_float_0_1 = normalizza_immagine(min_val, max_val)
        
        # Applica la correzione Gamma
        img_final = applica_gamma(img_float_0_1, state.gamma_val)
        
        # Aggiorna i dati mostrati sull'asse
        im.set_data(img_final)
        
        # Forza il ridisegno della figura
        fig.canvas.draw_idle()

    def scroll_gamma(event):
        """Callback per la rotellina del mouse (scroll) per la correzione Gamma."""
        
        # 'up' o 'down' indica la direzione dello scroll
        if event.button == 'up':
            state.gamma_val = min(5.0, state.gamma_val + gamma_step)
        elif event.button == 'down':
            state.gamma_val = max(0.1, state.gamma_val - gamma_step)
        else:
            return
        
        # Aggiorna l'immagine con il nuovo valore Gamma
        update_image(slider_range.val)
        
        # Aggiorna il titolo per mostrare il valore Gamma corrente
        ax.set_title(f"Gamma: {state.gamma_val:.2f}")
        fig.canvas.draw_idle()

    ax_button = fig.add_axes([0.7, 0.05, 0.1, 0.05])
    button = Button(ax_button, 'SAVE')
    def save(event):
        print('SALVATAGGIO')
        save_enhanced_image(p)

    button.on_clicked(save)

    # 3. Connessione degli Eventi
    slider_range.on_changed(update_image)
    fig.canvas.mpl_connect('scroll_event', scroll_gamma)

    # 4. Avvio Iniziale e Display
    print("\nVisualizzatore FITS Interattivo Avviato:")
    print("  - Trascina i manici dello slider (a destra) per lo Stretching Min/Max.")
    print("  - Scorri la rotellina del mouse sopra l'immagine per regolare il Gamma.")
    update_image(None) # Chiama la funzione una volta per visualizzare l'immagine iniziale
    plt.show()


def save_enhanced_image(p:Params):
    print()
    directory_path = os.path.basename(p.input_path)
    temp = p.input_path.split('.',1)
    print(temp)
    extension = temp[1]
    name=temp[0]
    print(name,extension)
    out_path_recommended = name+"_enhanced"+"_"+p.option
    if p.option == OPTIONS[0]:
        out_path_recommended=out_path_recommended+"_nrad_"+str(p.nrad)+"_ntheta_"+str(p.ntheta)+"_stdtheta_"+str(1/p.rejsig)+"_"
    if p.option == OPTIONS [1]:
        out_path_recommended=out_path_recommended+"_nrad_"+str(p.nrad)+"_ntheta_"+str(p.ntheta)+"_"
    if p.option == OPTIONS [3]:
        out_path_recommended=name+"_enhanced"+"Division by 1_rho profile"+"_nrad_"+str(p.nrad)+"_ntheta_"+str(p.ntheta)+"_"
    if p.option == OPTIONS[2]:
        out_path_recommended=out_path_recommended+"_nrad_"+str(p.nrad)+"_ntheta_"+str(p.ntheta)+"_stdtheta_"+str(1/p.rejsig)+"_nsig_"+str(p.nsig)+"_"
    if p.option == OPTIONS[4]:
        out_path_recommended=out_path_recommended+"_A_"+str(p.A)+"_B_"+str(p.B)+"_N_"+str(p.N)+"_NUMLOG_"+str(p.NUMLOG)+"_"

    

    root = tk.Tk()
    root.withdraw() # Nasconde la finestra principale vuota

    # 2. Definisci i tipi di file filtrabili
    filetypes = [("Immagini FITS", "*.fit *.fits")] 
    
    # 3. Apri la finestra di dialogo "Salva con nome"
    filepath = fd.asksaveasfilename(
        title='Salva il file come...',
        initialdir=directory_path,       # Directory iniziale (es. radice del sistema)
        initialfile=out_path_recommended,
        defaultextension="."+extension, # Estensione predefinita se l'utente non la specifica
        filetypes=filetypes
    )
    if not filepath:
        print("Saving Aborted")
        return
    p.hdul[0].data = p.imn
    p.hdul[0].header['DATAMIN']=np.min(p.imn)
    p.hdul[0].header['AVISUMIN']=np.min(p.imn)
    p.hdul[0].header['DATAMAX']=np.max(p.imn)
    p.hdul[0].header['AVISUMAX']=np.max(p.imn)
    p.hdul.writeto(filepath,overwrite=True)

