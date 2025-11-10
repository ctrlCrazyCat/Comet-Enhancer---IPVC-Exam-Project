import tkinter as tk
from tkinter import filedialog , messagebox, ttk,font as tkFont

from astropy.io import fits
import matplotlib.pyplot as plt
import cv2
import os
import comet_pack
import numpy as np

from matplotlib.cm import get_cmap
import matplotlib
import sys
import matplotlib.pyplot as plt
from matplotlib.widgets import RangeSlider,Button
from matplotlib.cm import get_cmap
import tkinter as tk

matplotlib.use('TkAgg')

from tooltip import Tooltip


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


def conditional_config(entry_widget_info,configuration):
        (condition,text) = configuration
        entry_widget_info.config(fg='green',cursor="question_arrow")
        if not condition:
            entry_widget_info.config(fg='red',cursor="question_arrow")
        entry_widget_info.tooltip.text = text

def _get_info_label(frame,tooltip_enabled = False):
    out = tk.Label(frame, text="l",font=tkFont.Font(family="Wingdings", size=12, weight="bold"),fg="gray")
    out.tooltip= Tooltip(out,"",tooltip_enabled)
    
    return out
def _get_guide_label(frame,text):
    out=tk.Label(frame,text=" ?  ",font=tkFont.Font(size=14,weight="bold"),fg="blue", cursor="question_arrow")
    out.tooltip = Tooltip(out,text,tooltip_enabled=True)
    return out

def _get_custom_entry(frame,text_variable,text_shown,tooltip_enabled = False):
    out = tk.Entry(frame, textvariable=text_variable, width=15, state=tk.DISABLED)
    out.info_label = _get_info_label(frame,tooltip_enabled)
    out.info_label.pack(side=tk.LEFT)
    tk.Label(frame, text=text_shown).pack(side=tk.LEFT)
    out.pack(side=tk.LEFT)
    return out

class ImageProcessingGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Comet Enhancer")

        # --- Variabili di stato dell'applicazione ---
        self.image_path = None
        self.image_width = 0
        self.image_height = 0
        self.log_abeled = True
        
        # Aggiungi una variabile per tenere traccia se il submit è avvenuto con successo
        self.submitted_successfully = False

        # --- Variabili per le caselle di testo principali ---
        self.center_x_var = tk.StringVar()
        
        
        self.center_y_var = tk.StringVar()
        self.xmin_var = tk.StringVar()
        self.xmax_var = tk.StringVar()
        self.ymin_var = tk.StringVar()
        self.ymax_var = tk.StringVar()

        # --- Variabili per TUTTI i campi condizionali ---
        self.rho_pixels_var = tk.StringVar()
        self.rho_max = 0
        self.rho_default = 0
        self.theta_pixels_var = tk.StringVar()
        self.std_dev_theta_var = tk.StringVar()
        self.min_max_std_dev_var = tk.StringVar()
        self.kernel_a_term_var = tk.StringVar()
        self.kernel_b_term_var = tk.StringVar()
        self.kernel_n_term_var = tk.StringVar()
        self.transform_log_var = tk.BooleanVar()

        # --- Variabile per il Combobox ---
        self.combobox_choice = tk.StringVar()
        
        self.combobox_choice.set(OPTIONS[0])

        # --- Tracciamento delle modifiche nelle caselle di testo (tutte le variabili) ---
        self.center_x_var.trace("w", self._validate_all_inputs)
        self.center_y_var.trace("w", self._validate_all_inputs)
        self.xmin_var.trace("w", self._validate_all_inputs)
        self.xmax_var.trace("w", self._validate_all_inputs)
        self.ymin_var.trace("w", self._validate_all_inputs)
        self.ymax_var.trace("w", self._validate_all_inputs)
        
        self.rho_pixels_var.trace("w", self._validate_all_inputs)
        self.theta_pixels_var.trace("w", self._validate_all_inputs)
        self.std_dev_theta_var.trace("w", self._validate_all_inputs)
        self.min_max_std_dev_var.trace("w", self._validate_all_inputs)
        self.kernel_a_term_var.trace("w", self._validate_all_inputs)
        self.kernel_b_term_var.trace("w", self._validate_all_inputs)
        self.kernel_n_term_var.trace("w", self._validate_all_inputs)
        self.transform_log_var.trace("w", self._validate_all_inputs)
        
        # --- Creazione dei Widget ---

        # Pulsante per scegliere l'immagine
        self.select_image_button = tk.Button(root, text="Scegli Immagine FIT/FITS", command=self._select_image,cursor="hand2") # Modificato il testo
        self.select_image_button.pack(pady=10)

        self.image_info_label = tk.Label(root, text="Nessuna immagine selezionata.")
        self.image_info_label.pack(pady=5)

        # Frame per i valori del Centro (Double)
        self.center_frame = ttk.LabelFrame(root, text="Centro Nucleo Cometa (X, Y - double)")
        self.center_frame.pack(pady=10, padx=10, fill=tk.X)

        self.all_entries=[]

        
        
        self.entry_center_x = _get_custom_entry(self.center_frame,self.center_x_var,"Centro X:")   
        self.all_entries.append(self.entry_center_x)

        tk.Label(self.center_frame).pack(side=tk.LEFT,padx=20)

        self.entry_center_y = _get_custom_entry(self.center_frame,self.center_y_var,"Centro Y:")
        self.all_entries.append(self.entry_center_y)

        _get_guide_label(self.center_frame,"Coordinate (x,y) dell'optocentro - indicano dove centrare il miglioramento.\nTipicamente la parte più luminosa dell'immagine.\nNota bene che il pixel all'origine è (1,1) e non (0,0).").pack(side=tk.LEFT)

        # Frame per i Limiti (Int)
        self.limits_frame = ttk.LabelFrame(root, text="Limiti di Elaborazione (X Min/Max, Y Min/Max - int)")
        self.limits_frame.pack(pady=10, padx=10, fill=tk.X)

        # Riga X Min/Max
        self.xminmax_frame = tk.Frame(self.limits_frame)
        self.xminmax_frame.pack(fill=tk.X, pady=2)

             
        
        self.entry_xmin = _get_custom_entry(self.xminmax_frame,self.xmin_var,"    X Min:")
        self.all_entries.append(self.entry_xmin)

        tk.Label(self.xminmax_frame).pack(side=tk.LEFT,padx=20)

        self.entry_xmax = _get_custom_entry(self.xminmax_frame,self.xmax_var,"    X Max:")

        self.all_entries.append(self.entry_xmax)
        
        _get_guide_label(self.xminmax_frame,"Modificare questi parametri se si desidera elaborare solo una sezione dell'immagine").pack(side=tk.LEFT)

        # Riga Y Min/Max
        self.yminmax_frame = tk.Frame(self.limits_frame)
        self.yminmax_frame.pack(fill=tk.X, pady=2)
        
        self.entry_ymin = _get_custom_entry(self.yminmax_frame,self.ymin_var,"    Y Min:")
        self.all_entries.append(self.entry_ymin)

        tk.Label(self.yminmax_frame).pack(side=tk.LEFT,padx=20)

        self.entry_ymax = _get_custom_entry(self.yminmax_frame,self.ymax_var,"    Y Max:")  
        self.all_entries.append(self.entry_ymax)

        # Combobox
        self.style = ttk.Style()
        self.style.configure('TCombobox', fieldbackground='white', background='lightgrey')
        self.combobox = ttk.Combobox(root, textvariable=self.combobox_choice, values=OPTIONS, state='readonly', width=40,cursor="hand2")
        self.combobox.bind("<<ComboboxSelected>>", self._on_combobox_selected_event_for_bind)
        self.combobox.pack(pady=10)

        # --- Singolo Frame per tutti i campi condizionali ---
        self.conditional_params_frame = ttk.LabelFrame(root, text="Parametri Specifici dell'Opzione")
        self.conditional_params_frame.pack(pady=10, padx=10, fill=tk.X) 

        # Inizializzazione di TUTTI i widget condizionali e memorizzazione per un facile accesso
        self.conditional_widgets_map = {
            "rho_pixels": {
                "label": tk.Label(self.conditional_params_frame, text="Numero di pixel asse ρ (rho):"), 
                "info_label":_get_info_label(self.conditional_params_frame,tooltip_enabled=True),
                "widget": tk.Entry(self.conditional_params_frame, textvariable=self.rho_pixels_var, width=20),
                "variable": self.rho_pixels_var,
                "guide":_get_guide_label(self.conditional_params_frame,"Questa routine converte l'immagine di ingresso in coordinate polari (ρ, θ).\nQuesto parametro permette di scegliere la quantità di pixel da impostare sull'asse ρ.")
            },
            "theta_pixels": {
                "label": tk.Label(self.conditional_params_frame, text="Numero di pixel asse θ:"), 
                "info_label":_get_info_label(self.conditional_params_frame,tooltip_enabled=True),
                "widget": tk.Entry(self.conditional_params_frame, textvariable=self.theta_pixels_var, width=20),
                "variable": self.theta_pixels_var,
                "guide":_get_guide_label(self.conditional_params_frame,"Questa routine converte l'immagine di ingresso in coordinate polari (ρ, θ).\nQuesto parametro permette di scegliere la quantità di pixel da impostare sull'asse θ.")
            },
            "std_dev_theta": {
                "label": tk.Label(self.conditional_params_frame, text="Quante deviazioni standard dovrebbero essere accettate per i pixel asse θ?"), 
                "info_label":_get_info_label(self.conditional_params_frame,tooltip_enabled=True),
                "widget": tk.Entry(self.conditional_params_frame, textvariable=self.std_dev_theta_var, width=20),
                "variable": self.std_dev_theta_var,
                "guide":_get_guide_label(self.conditional_params_frame,'Gli algoritmi che migliorano le immagini possono essere molto sensibi ai pixel "morti", corpi luminosi e altre anomalie.\nPer mitigare questo problema si ignorano i pixel anomali rispetto a quelli vicini ad essi.\nQuesto parametro permette di impostare quante deviazioni standard dal pixel medio sono accettabili. Una scelta comune è 3.')
            },
            "min_max_std_dev": {
                "label": tk.Label(self.conditional_params_frame, text="How many standard deviations from the mean should\nthe minimum and maximum pixel values be?"), 
                "info_label":_get_info_label(self.conditional_params_frame,tooltip_enabled=True),
                "widget": tk.Entry(self.conditional_params_frame, textvariable=self.min_max_std_dev_var, width=20),
                "variable": self.min_max_std_dev_var,
                "guide":_get_guide_label(self.conditional_params_frame,"Nella rinormalizzazione azimutale i valori dei pixel sono scalati ad un nuovo range determinato dagli altri pixel nell'azimut.\nQuesto parametro permette di impostare quante deviazioni standard dal pixel medio i valore dovrebbero essere all'interno dell'immagine migliorata.\nUna scelta comune è 3.")
            },
            "kernel_a_term": {
                "label": tk.Label(self.conditional_params_frame, text="valore Kernel A (double):"), 
                "info_label":_get_info_label(self.conditional_params_frame,tooltip_enabled=True),
                "widget": tk.Entry(self.conditional_params_frame, textvariable=self.kernel_a_term_var, width=20),
                "variable": self.kernel_a_term_var,
                "guide":_get_guide_label(self.conditional_params_frame,"Il valore kernel A è la distanza di base tra il pixel che deve essere migliorato e il lato del kernel.\nAumentare questo valore fa incrementare la grandezza del kernel vicino al nucleo, senza un impatto significativo sulla grandezza del kernel lontano dal nucleo\nVisita la guida completa per una spiegazione più dettagliata sul kernel.\nQuesto valore è di solito nell'ordine di 1")
            },
            "kernel_b_term": {
                "label": tk.Label(self.conditional_params_frame, text="valore Kernel B (double):"), 
                "info_label":_get_info_label(self.conditional_params_frame,tooltip_enabled=True),
                "widget": tk.Entry(self.conditional_params_frame, textvariable=self.kernel_b_term_var, width=20),
                "variable": self.kernel_b_term_var,
                "guide":_get_guide_label(self.conditional_params_frame,"Il valore kernel B scala linearmente la grandezza del kernel in base alla distanza dal nucleo.\nAumentare questo valore fa incrementare la grandezza del kernel lungo tutte le distanze dal nucleo.\nVisita la guida completa per una spiegazione più dettagliata sul kernel.\nQuesto valore è di solito nell'ordine di 1")
            },
            "kernel_n_term": {
                "label": tk.Label(self.conditional_params_frame, text="valore Kernel N (double):"), 
                "info_label":_get_info_label(self.conditional_params_frame,tooltip_enabled=True),
                "widget": tk.Entry(self.conditional_params_frame, textvariable=self.kernel_n_term_var, width=20),
                "variable": self.kernel_n_term_var,
                "guide":_get_guide_label(self.conditional_params_frame,"Il valore N del kernel scala esponenzialmente la grandezza del kernel in base alla distanza dal nucleo.\nAumentare questo valore fa aumentare la grandezza del kernel quando lontano dal nucleo, senza avere un impatto significativo sulla grandezza del kernel vicino al nucleo.\nVisita la guida completa per una spiegazione più dettagliata sul kernel.\nQuesto valore è di solito nell'ordine di 0.1")
            },
            "transform_log": {
                "label": tk.Label(self.conditional_params_frame, text="Trasformare l'immagine input in scala log-10 prima di migliorarla?"), 
                "info_label":_get_info_label(self.conditional_params_frame),
                "widget": tk.Checkbutton(self.conditional_params_frame, variable=self.transform_log_var, text="Sì"),
                "variable": self.transform_log_var,
                "guide":_get_guide_label(self.conditional_params_frame,"Questa opzione permette di riscalare l'immagine nello spazio log-10 prima che la procedura di miglioramento venga avviata.\nLa presenza di pixel negativi non permetterà la trasformazione logaritmica.\nPuoi controllare il valore 'logarithmic image' nell'header dell'immagine migliorata per avere la certezza che la trasformazione sia stata applicata.")
            }
        }
        

        # Mappa le opzioni ai widget necessari e alla loro riga di griglia
        self.option_to_widgets_config = {
            OPTIONS[0]: [ # Division by Azimuthal Average
                ("rho_pixels", 0), ("theta_pixels", 1), ("std_dev_theta", 2)
            ],
            OPTIONS[1]: [ # Division by Azimuthal Median
                ("rho_pixels", 0), ("theta_pixels", 1)
            ],
            OPTIONS[2]: [ # Azimuthal Renormalization
                ("rho_pixels", 0), ("theta_pixels", 1), ("std_dev_theta", 2), ("min_max_std_dev", 3)
            ],
            OPTIONS[3]: [ # Division by 1/rho profile  
                #("rho_pixels", 0), ("theta_pixels", 1)
            ],
            OPTIONS[4]: [ # Radially Variable Spatial Filtering
                ("kernel_a_term", 0), ("kernel_b_term", 1), ("kernel_n_term", 2), ("transform_log", 3)
            ]
        }
        
        # Pulsante Submit (inizialmente nascosto e disabilitato)
        grande_font = tkFont.Font(family="Arial", size=20, weight="bold")
        self.submit_button = tk.Button(root, text="Submit",font=grande_font, command=self._on_submit_button_click,cursor="hand2") # Modificato il command

        # Inizializzazione degli stati
        self._validate_all_inputs()
    def params_process(self,option):
        o = Params
        if self.image_path:
            o.input_path=self.image_path
        o.option = option

        o.y_upper_lim = int(self.ymax_var.get())
        o.x_upper_lim = int(self.xmax_var.get())
        o.y_lower_lim = int(self.ymin_var.get())-1
        o.x_lower_lim = int(self.xmin_var.get())-1

        
        o.ynuc=float(self.entry_center_y.get())-1.0
        o.xnuc=float(self.entry_center_x.get())-1.0
        #[hdul,imold,log_abeled]=comet_pack.get_input_data(self.image_path)
        o.hdul=self.hdul
        (o.NROW,o.NCOL)= self.imold.shape
        
        if option==OPTIONS[0]:
            o.nrad=int(self.rho_pixels_var.get())
            o.ntheta=int(self.theta_pixels_var.get())

            try:
                o.rejsig=1/float(self.std_dev_theta_var.get())
            except ZeroDivisionError as e:
                o.rejsig = float('inf')

        if selected_option == OPTIONS[1]:
            o.nrad=int(self.rho_pixels_var.get())
            o.ntheta=int(self.theta_pixels_var.get()) 

        if option==OPTIONS[2]:
            o.nrad=int(self.rho_pixels_var.get())
            o.ntheta=int(self.theta_pixels_var.get())
            o.rejsig=1/float(self.std_dev_theta_var.get())
            o.nsig=float(self.min_max_std_dev_var.get())
        
        if option==OPTIONS[4]:
            o.A = float(self.kernel_a_term_var.get())
            o.B = float(self.kernel_b_term_var.get())
            o.N = float(self.kernel_n_term_var.get())
            o.NUMLOG = bool(self.transform_log_var.get())

        
        return o

    def _select_image(self):
        """Apre la finestra di dialogo per scegliere un'immagine FITS.""" # Modificato il testo
        file_path = filedialog.askopenfilename(
            title="Seleziona un'immagine FITS",
            filetypes=[("Immagini FITS", "*.fit *.fits")] # Modificato il filetype
        )
        if file_path:
            print('QUI1')
            self.image_path = file_path
            print('QUI2')
            try:
                with fits.open(file_path) as hdul: #hdul,imold,_upper_,y_upper_lim
                    self.hdul = hdul.copy()
                    self.imold = hdul[0].data.copy()
                    self.image_width=hdul[0].header['NAXIS1']
                    self.image_height=hdul[0].header['NAXIS2']
                    self.imold[np.isnan(self.imold)] = np.min(self.imold[~np.isnan(self.imold)])
                    if np.min(self.imold)<=0:
                        self.transform_log_var.set(False)
                        self.log_abeled = False
                    
                self.image_info_label.config(text=f"Immagine: {self.image_path.split('/')[-1]} ({self.image_width}x{self.image_height})")
                
                # Abilita TUTTE le caselle principali
                self.entry_center_x.config(state=tk.NORMAL)
                self.entry_center_y.config(state=tk.NORMAL)
                self.entry_xmin.config(state=tk.NORMAL)
                self.entry_xmax.config(state=tk.NORMAL)
                self.entry_ymin.config(state=tk.NORMAL)
                self.entry_ymax.config(state=tk.NORMAL)
                
                # Resetta i valori delle caselle principali a vuoto
                self.center_x_var.set("")
                self.center_y_var.set("")
                self.xmin_var.set("1")
                self.xmax_var.set(str(self.image_width))
                self.ymin_var.set("1")
                self.ymax_var.set(str(self.image_height))

                # Resetta e aggiorna lo stato dei campi condizionali
                self._update_conditional_fields()
                for e in self.all_entries:
                    e.info_label.tooltip.enable()

                self._validate_all_inputs() # Forza una validazione completa
            except Exception as e:
                messagebox.showerror("Errore Caricamento Immagine", f"Impossibile caricare l'immagine: {e}")
                self._reset_image_state()
        else:
            self._reset_image_state()

    def _update_conditional_fields(self):
        """Nasconde tutti i widget condizionali e abilita/mostra solo quelli necessari."""
        # 1. Nasconde e disabilita tutti i widget condizionali e resetta le loro variabili
        for key, components in self.conditional_widgets_map.items():
            label = components["label"]
            info_label = components["info_label"]
            widget = components["widget"]
            variable = components["variable"] 
            guide = components["guide"]
            
            info_label.grid_remove()
            label.grid_remove() 
            widget.grid_remove()
            guide.grid_remove()
            
            if isinstance(widget, tk.Entry):
                widget.config(state=tk.DISABLED)
                #variable.set("") # Resetta la StringVar
            elif isinstance(widget, tk.Checkbutton):
                widget.config(state=tk.DISABLED)
                variable.set(False) # Resetta BooleanVar

        # 2. Mostra/Nasconde il FRAME e popola i widget in base alla selezione
        selected_option = self.combobox_choice.get()
        print(selected_option,"AAA")

        if selected_option == OPTIONS[3]:
            # Opzione [3]: Nascondi il frame
            self.conditional_params_frame.pack_forget()
        
        elif selected_option in self.option_to_widgets_config:
            # Altre opzioni valide: Mostra frame e popola
            self.conditional_params_frame.pack(pady=10, padx=10, fill=tk.X) # MOSTRA
            self.conditional_params_frame.config(text=f"Parametri per '{selected_option}'")
            
            for widget_key, row_num in self.option_to_widgets_config[selected_option]:
                components = self.conditional_widgets_map[widget_key]
                label = components["label"]
                info_label=components["info_label"]
                widget = components["widget"]  
                guide = components["guide"]
                
                info_label.grid(row=row_num,column=0,padx=0,pady=0,sticky="w")
                label.grid(row=row_num, column=1, padx=5, pady=5, sticky="w")
                widget.grid(row=row_num, column=2, padx=5, pady=5, sticky="ew")
                widget.config(state=tk.NORMAL)
                guide.grid(row=row_num,column=3,padx=5,pady=5,sticky="ew")
            if self.log_abeled==False:
                self.conditional_widgets_map["transform_log"]["widget"].config(state=tk.DISABLED)
        else:
            # Caso fallback (es. opzione non in config, ma non è op[3])
            # o se l'opzione [3] non fosse gestita esplicitamente sopra
            self.conditional_params_frame.pack(pady=10, padx=10, fill=tk.X) # MOSTRA
            self.conditional_params_frame.config(text="Nessun parametro specifico per questa opzione.")


        # Nasconde il pulsante Submit (lo stato sarà gestito da _validate_all_inputs)
        self.submit_button.pack_forget()
        #self.submit_button.config(state=tk.DISABLED)


    def _reset_image_state(self):
        """Resetta lo stato relativo all'immagine e disabilita tutti gli input."""
        self.image_path = None
        self.image_width = 0
        self.image_height = 0
        self.image_info_label.config(text="Nessuna immagine selezionata.")
        
        # Disabilita tutte le caselle principali
        self.entry_center_x.config(state=tk.DISABLED)
        self.entry_center_y.config(state=tk.DISABLED)
        self.entry_xmin.config(state=tk.DISABLED)
        self.entry_xmax.config(state=tk.DISABLED)
        self.entry_ymin.config(state=tk.DISABLED)
        self.entry_ymax.config(state=tk.DISABLED)
        
        # Svuota i contenuti delle caselle principali
        self.center_x_var.set("")
        self.center_y_var.set("")
        self.xmin_var.set("")
        self.xmax_var.set("")
        self.ymin_var.set("")
        self.ymax_var.set("")

        for e in self.all_entries:
            e.info_label.config(fg='gray',cursor="arrow")
            e.info_label.tooltip.disable()
        # Aggiorna lo stato dei campi condizionali
        self._update_conditional_fields()

        self.conditional_params_frame.pack_forget() # Nasconde anche il frame contenitore
        self._validate_all_inputs() # Forza una validazione completa


    def _validate_double_input(self, value_var, min_val=None, max_val=None, allow_zero=True):
        """
        Valida se una stringa è un double e rientra nel range specificato.
        Se allow_zero è False, il valore deve essere > 0.
        """
        value_str = value_var.get()
        if not value_str: 
            stringa = "è assente, deve essere maggiore o uguale a "+str(min_val)
            if max_val is not None:
                stringa = stringa + " , minore o uguale a "+str(max_val)
            if not allow_zero:
                stringa = stringa +" e diverso da 0"
            return (False,stringa)
        try:
            val = float(value_str)
            if not allow_zero and val <= 0:
                return (False,"è minore o uguale a zero")
            
            if min_val is not None and val < min_val:
                return (False,"è minore di "+str(min_val))
            if max_val is not None and val > max_val:
                return (False,"è maggiore di "+str(max_val))
            
            return (True,"correttamente inserito")
        except ValueError:
            return (False,"non rappresenta un numero a virgola mobile")
    def _validate_rho_input_and_approximate(self,all_center_ok:bool,all_limits_ok:bool,limits_nucleus_ok:bool):
        if (not all_center_ok) and (not all_limits_ok):
            return (False,"dati sul Centro Nucleo Cometa e sui limiti di elaborazione non correttamente inseriti o assenti")
        if not all_limits_ok:
            return (False,"dati sui limiti di elaborazione non correttamente inseriti o assenti")
        if not all_center_ok:
            return (False,"dati sul Centro Nucleo Cometa non correttamente inseriti o assenti")
        if not limits_nucleus_ok:
            return (False,"i limiti di elaborazione non sono coerenti con le coordinate del Nucleo Cometa")
        
        return self._validate_int_input_and_approximate(self.rho_pixels_var,1,self.rho_max)
    

    def _validate_int_input_and_approximate(self, value_var, min_val=1, max_val=None):
        """
        Valida se una stringa può essere convertita a int (con approssimazione)
        e se rientra nel range [min_val, max_val] (se max_val è fornito).
        Aggiorna la StringVar con l'intero approssimato.
        """
        value_str = value_var.get()
        if not value_str: 
            stringa = "è assente, deve essere maggiore o uguale a "+str(min_val)
            if max_val is not None:
                stringa = stringa + " e minore o uguale a "+str(max_val)
            return (False,stringa)

        try:
            temp_val = float(value_str)
            int_val = round(temp_val)
            
            # Se il valore convertito a float non è un intero esatto, non è valido.
            # Questo controllo è per assicurarsi che l'utente inserisca numeri interi per campi int
            # ma permette la digitazione temporanea di decimali (es. "12." o "12.0")
            if abs(temp_val - int_val) > 1e-9: # Piccola tolleranza per float comparison
                 return (False,"non rappresenta un numero intero")
            
            # Solo aggiorna se il valore effettivo della StringVar è diverso
            if value_var.get() != str(int_val):
                value_var.set(str(int_val)) 

            if min_val is not None and int_val < min_val:
                return (False,"è minore di "+str(min_val))
            if max_val is not None and int_val > max_val:
                return (False,"è maggiore di "+str(max_val))
                
            
            return (True,"correttamente inserito")
        except ValueError:
            return (False,"non rappresenta un numero intero")
    
    
    def _validate_all_inputs(self, *args):#TODO: da rivedere
        """Funzione principale di validazione che gestisce lo stato di tutti i widget."""
        # args è presente qui a causa dei trace callback, ma non viene usato
        self.validations = Params
        all_ok = []
        if not self.image_path:
            self.conditional_params_frame.pack_forget() 
            self.submit_button.pack_forget() 
            #self.submit_button.config(state=tk.DISABLED)
            return

        # Validazione dei campi principali (Centro e Limiti)
        center_x_ok = self._validate_double_input(self.center_x_var, min_val=1, max_val=self.image_width)
        all_ok.append(center_x_ok)
        center_y_ok = self._validate_double_input(self.center_y_var, min_val=1, max_val=self.image_height)
        all_center_ok = center_x_ok[0] and center_y_ok[0]
        all_ok.append(center_y_ok)
        self.validations.all_center_ok = all_center_ok
        



        
        

        xmin_ok = self._validate_int_input_and_approximate(self.xmin_var, min_val=1, max_val=self.image_width)
        all_ok.append(xmin_ok)
        xmax_ok = self._validate_int_input_and_approximate(self.xmax_var, min_val=1, max_val=self.image_width)
        all_ok.append(xmax_ok)
        ymin_ok = self._validate_int_input_and_approximate(self.ymin_var, min_val=1, max_val=self.image_height)
        all_ok.append(ymin_ok)
        ymax_ok = self._validate_int_input_and_approximate(self.ymax_var, min_val=1, max_val=self.image_height)
        all_ok.append(ymax_ok)
        limits_x_ok = xmin_ok[0] and xmax_ok[0]
        limits_y_ok = ymin_ok[0] and ymax_ok[0]
        all_limits_ok = limits_x_ok and limits_y_ok

        
        xnuc = -1.0
        ynuc = -1.0
        xfin=-1.0
        xint = -1.0
        yfin=-1.0
        yint = -1.0
        
        self.limits_nucleus_ok = True
        if all_center_ok:
            xnuc = float(self.center_x_var.get())
            ynuc = float(self.center_y_var.get())


            limits_nucleus_x_ok=True
            if limits_x_ok:
                xfin=float(self.xmax_var.get())
                xint = float(self.xmin_var.get())
                limits_nucleus_x_ok = xnuc<=xfin and xnuc >= xint 
            limits_nucleus_y_ok=True
            if limits_y_ok:
                yfin=float(self.ymax_var.get())
                yint = float(self.ymin_var.get())
                limits_nucleus_y_ok = ynuc <= yfin and ynuc >= yint
            self.limits_nucleus_ok = limits_nucleus_x_ok and limits_nucleus_y_ok
            if not self.limits_nucleus_ok:
                self.limits_nucleus_x_ok=limits_nucleus_x_ok
                self.limits_nucleus_y_ok=limits_nucleus_y_ok

            if all_limits_ok:
                xmax = max((xfin-xnuc),(xnuc-xint))
                ymax =max((yfin-ynuc),(ynuc-yint))
                self.rho_max=int(np.floor((xmax**2+ymax**2)**0.5))
                self.rho_default= int(min((xfin-xnuc),(xnuc-xint),(yfin-ynuc),(ynuc-yint)))

        

        for i, entry_widget in enumerate(self.all_entries):
            conditional_config(entry_widget.info_label,all_ok[i])


        min_max_ok = False
        self.validations.all_limits_ok = all_limits_ok
        
        if all_limits_ok:
            x_min_max_ok = (int(self.xmin_var.get()) < int(self.xmax_var.get()))
            if not x_min_max_ok:
                conditional_config(self.entry_xmax.info_label,(False,"Limite superiore X minore o uguale al limite inferiore"))
                conditional_config(self.entry_xmin.info_label,(False,"Limite inferiore X maggiore o uguale al limite superiore"))
            y_min_max_ok=(int(self.ymin_var.get()) < int(self.ymax_var.get()))
            if not y_min_max_ok:
                conditional_config(self.entry_ymax.info_label,(False,"Limite superiore Y minore o uguale al limite inferiore"))
                conditional_config(self.entry_ymin.info_label,(False,"Limite inferiore X maggiore o uguale al limite superiore"))
            min_max_ok = x_min_max_ok and y_min_max_ok
            self.validations.min_max_ok = min_max_ok
        if not self.limits_nucleus_ok:
            if not self.limits_nucleus_x_ok:
                conditional_config(self.entry_xmax.info_label,(False,"Il range di valori sull'asse X non contniene il Nucleo Cometa"))
                conditional_config(self.entry_xmin.info_label,(False,"Il range di valori sull'asse X non contniene il Nucleo Cometa"))
            if not self.limits_nucleus_y_ok:
                conditional_config(self.entry_ymax.info_label,(False,"Il range di valori sull'asse Y non contniene il Nucleo Cometa"))
                conditional_config(self.entry_ymin.info_label,(False,"Il range di valori sull'asse Y non contniene il Nucleo Cometa"))

        main_inputs_valid = all_center_ok and all_limits_ok and min_max_ok and self.limits_nucleus_ok

        self.validations.main_inputs_valid = main_inputs_valid
        if not main_inputs_valid:
            self.validations.all_center_ok = all_center_ok
            if not all_center_ok:
                self.validations.center_x_ok = center_x_ok
                self.validations.center_y_ok=center_y_ok

            self.validations.all_limits_ok = all_limits_ok
            if not all_limits_ok:
                self.validations.xmin_ok=xmin_ok
                self.validations.xmax_ok = xmax_ok
                self.validations.ymin_ok=ymin_ok
                self.validations.ymax_ok = ymax_ok
            if all_limits_ok:
                if not min_max_ok:
                    self.validations.x_min_max_ok = x_min_max_ok
                    self.validations.y_min_max_ok = y_min_max_ok
            


        # Mostra sempre il frame contenitore dei parametri condizionali se un'immagine è caricata
        #self.conditional_params_frame.pack(pady=10, padx=10, fill=tk.X)
        
        
        selected_option = self.combobox_choice.get()
        #TODO: VEDERE SE SI PUO OTTIMIZZARE QUI IL CODICE
        # Validazione dei campi condizionali attualmente visibili
        
        
        if selected_option == OPTIONS[0]: # Division by Azimuthal Average
            rho_pixels_ok = self._validate_rho_input_and_approximate(all_center_ok,all_limits_ok,self.limits_nucleus_ok)
            conditional_config(self.conditional_widgets_map["rho_pixels"]["info_label"],rho_pixels_ok)
            theta_pixels_ok = self._validate_int_input_and_approximate(self.theta_pixels_var, min_val=1)
            conditional_config(self.conditional_widgets_map["theta_pixels"]["info_label"],theta_pixels_ok)
             
            std_dev_theta_ok= self._validate_double_input(self.std_dev_theta_var,min_val=0)
            conditional_config(self.conditional_widgets_map["std_dev_theta"]["info_label"],std_dev_theta_ok)            
            all_conditional_ok = rho_pixels_ok[0] and theta_pixels_ok[0] and std_dev_theta_ok[0]
            self.validations.option = 0
            self.validations.all_conditional_ok = all_conditional_ok
            if not all_conditional_ok:
                self.validations.rho_pixels_ok = rho_pixels_ok
                self.validations.theta_pixels_ok = theta_pixels_ok
                self.validations.std_dev_theta_ok = std_dev_theta_ok

        elif selected_option == OPTIONS[1]: # Azimuthal Median / Division by 1/rho profile
            rho_pixels_ok = self._validate_rho_input_and_approximate(all_center_ok,all_limits_ok,self.limits_nucleus_ok)
            conditional_config(self.conditional_widgets_map["rho_pixels"]["info_label"],rho_pixels_ok)
            theta_pixels_ok = self._validate_int_input_and_approximate(self.theta_pixels_var, min_val=1)
            conditional_config(self.conditional_widgets_map["theta_pixels"]["info_label"],theta_pixels_ok)
            all_conditional_ok = rho_pixels_ok[0] and theta_pixels_ok[0]
            self.validations.option = 1
            self.validations.all_conditional_ok = all_conditional_ok
            if not all_conditional_ok:
                self.validations.rho_pixels_ok = rho_pixels_ok
                self.validations.theta_pixels_ok = theta_pixels_ok
        elif selected_option==OPTIONS[3]:
            self.validations.option = 3
            all_conditional_ok = True
                
            
        elif selected_option == OPTIONS[2]: # Azimuthal Renormalization
            rho_pixels_ok = self._validate_rho_input_and_approximate(all_center_ok,all_limits_ok,self.limits_nucleus_ok)
            conditional_config(self.conditional_widgets_map["rho_pixels"]["info_label"],rho_pixels_ok)
            theta_pixels_ok = self._validate_int_input_and_approximate(self.theta_pixels_var, min_val=1)
            conditional_config(self.conditional_widgets_map["theta_pixels"]["info_label"],theta_pixels_ok)
            std_dev_theta_ok = self._validate_double_input(self.std_dev_theta_var, min_val=-1e-15, allow_zero=True)
            conditional_config(self.conditional_widgets_map["std_dev_theta"]["info_label"],std_dev_theta_ok)
            min_max_std_dev_ok = self._validate_double_input(self.min_max_std_dev_var, min_val=-1e-15, allow_zero=True)
            conditional_config(self.conditional_widgets_map["min_max_std_dev"]["info_label"],min_max_std_dev_ok)
            all_conditional_ok = rho_pixels_ok[0] and theta_pixels_ok[0] and std_dev_theta_ok[0] and min_max_std_dev_ok[0]

            self.validations.option = 2
            self.validations.all_conditional_ok = all_conditional_ok
            if not all_conditional_ok:
                self.validations.rho_pixels_ok = rho_pixels_ok
                self.validations.theta_pixels_ok = theta_pixels_ok
                self.validations.std_dev_theta_ok = std_dev_theta_ok
                self.validations.min_max_std_dev_ok = min_max_std_dev_ok


        elif selected_option == OPTIONS[4]: # Radially Variable Spatial Filtering
            kernel_a_ok = self._validate_double_input(self.kernel_a_term_var,min_val=0, allow_zero=False)
            conditional_config(self.conditional_widgets_map["kernel_a_term"]["info_label"],kernel_a_ok)

            kernel_b_ok = self._validate_double_input(self.kernel_b_term_var,min_val=0, allow_zero=False)
            conditional_config(self.conditional_widgets_map["kernel_b_term"]["info_label"],kernel_b_ok)

            kernel_n_ok = self._validate_double_input(self.kernel_n_term_var,min_val=0, allow_zero=False)
            conditional_config(self.conditional_widgets_map["kernel_n_term"]["info_label"],kernel_n_ok)

            all_conditional_ok = kernel_a_ok[0] and kernel_b_ok[0] and kernel_n_ok[0]

            self.validations.option = 4
            self.validations.all_conditional_ok = all_conditional_ok
            if not all_conditional_ok:
                self.validations.kernel_a_ok = kernel_a_ok
                self.validations.kernel_b_ok = kernel_b_ok
                self.validations.kernel_n_ok = kernel_n_ok
                

        else: # Nessun parametro specifico per altre opzioni (se ce ne fossero)
            self.validations.option = -1
            all_conditional_ok = True
        self.validations.all_data_ok = main_inputs_valid and all_conditional_ok
        
        # Gestione del pulsante Submit
        self.submit_button.pack(pady=15)
        # if main_inputs_valid and all_conditional_ok:
        #     self.submit_button.config(state=tk.NORMAL) 
        # else:
        #     #self.submit_button.pack_forget() 
        #     self.submit_button.config(state=tk.DISABLED) 


    def _on_combobox_selected_event_for_bind(self, event):
        """Funzione richiamata quando cambia la selezione del Combobox (da bind)."""
        self._update_conditional_fields() # Aggiorna la visibilità e lo stato dei campi
        self._validate_all_inputs() # Forza una validazione completa


    def _data_warning(self):
        detail_string=""
        if not self.validations.main_inputs_valid:

            if not self.validations.all_center_ok:
                if not self.validations.center_x_ok[0]:
                    detail_string = detail_string + "-Centro X "+self.validations.center_x_ok[1]+"\n"
                if not self.validations.center_y_ok[0]:
                    detail_string = detail_string + "-Centro Y "+self.validations.center_y_ok[1]+"\n"
            if not self.validations.all_limits_ok:
                detail_string = detail_string+"\n"
                if not self.validations.xmin_ok[0]:
                    detail_string = detail_string + "-X Min "+self.validations.xmin_ok[1]+"\n"
                if not self.validations.xmax_ok[0]:
                    detail_string = detail_string + "-X Max "+self.validations.xmax_ok[1]+"\n"
                if not self.validations.ymin_ok[0]:
                    detail_string = detail_string + "-Y min "+self.validations.ymin_ok[1]+"\n"
                if not self.validations.ymax_ok[0]:
                    detail_string = detail_string + "-Y max "+self.validations.ymax_ok[1]+"\n"
            else:
                if not self.validations.min_max_ok:
                    detail_string = detail_string+"\n"
                    if not self.validations.x_min_max_ok:
                        detail_string = detail_string + "-X Min >= X Max\n"
                    if not self.validations.y_min_max_ok:
                        detail_string = detail_string + "-Y Min >= Y Max\n"
            if not self.limits_nucleus_ok:
                detail_string = detail_string+"\n"
                if not self.limits_nucleus_x_ok:
                    detail_string = detail_string + "-il range di valori sull'asse X non contniene il Nucleo Cometa\n"
                if not self.limits_nucleus_y_ok:
                    detail_string = detail_string + "-il range di valori sull'asse Y non contniene il Nucleo Cometa\n"
        #TODO: VEDERE SE SI PUO OTTIMIZZARE IL CODICE ANCHE QUI

        if not self.validations.all_conditional_ok:
            detail_string = detail_string+"\n"
            match self.validations.option:
                case 0:
                    if not self.validations.rho_pixels_ok[0]:
                        #TODO: riprodurre su tutto rho
                        if (not self.validations.all_center_ok) or not(self.validations.all_limits_ok) or (not self.limits_nucleus_ok):
                            detail_string=detail_string + "-Limiti per numero di pixel asse ρ (rho) non calcolabili\n"
                        else:
                            detail_string = detail_string + "-Numero di pixel asse ρ (rho) "+self.validations.rho_pixels_ok[1]+", sarà impostato "+str(self.rho_default)+" che corrisponde alla più grande circonfernza contenuta nell'immagine a partire dal nucelo\n"
                            self.rho_pixels_var.set(str(self.rho_default))
                    
                    if not self.validations.theta_pixels_ok[0]:
                        detail_string = detail_string + "-Numero di pixel asse θ "+self.validations.theta_pixels_ok[1]+"\n"
                    if not self.validations.std_dev_theta_ok[0]:
                        detail_string = detail_string + "-Deviazioni standard  pixel asse θ "+self.validations.std_dev_theta_ok[1]+"\n"
                case 1:
                    if not self.validations.rho_pixels_ok[0]:
                        detail_string = detail_string + "-Numero di pixel asse ρ (rho) "+self.validations.rho_pixels_ok[1]+"\n"
                    if not self.validations.theta_pixels_ok[0]:
                        detail_string = detail_string + "-Numero di pixel asse θ "+self.validations.theta_pixels_ok[1]+"\n"
                case 2:
                    if not self.validations.rho_pixels_ok[0]:
                        detail_string = detail_string + "-Numero di pixel asse ρ (rho) "+self.validations.rho_pixels_ok[1]+"\n"
                    if not self.validations.theta_pixels_ok[0]:
                        detail_string = detail_string + "-Numero di pixel asse θ "+self.validations.theta_pixels_ok[1]+"\n"
                    if not self.validations.std_dev_theta_ok[0]:
                        detail_string = detail_string + "-Deviazioni standard pixel asse θ "+self.validations.std_dev_theta_ok[1]+"\n"
                    if not self.validations.min_max_std_dev_ok[0]:
                        detail_string = detail_string + "-Deviazioni standard pixel asse θ "+self.validations.min_max_std_dev_ok[1]+"\n"
                case 4:
                    if not self.validations.kernel_a_ok[0]:
                        detail_string = detail_string + "valore Kernel A "+self.validations.kernel_a_ok[1]+"\n"
                    if not self.validations.kernel_b_ok[0]:
                        detail_string = detail_string + "valore Kernel B "+self.validations.kernel_b_ok[1]+"\n"
                    if not self.validations.kernel_n_ok[0]:
                        detail_string = detail_string + "valore Kernel N "+self.validations.kernel_n_ok[1]+"\n"

        messagebox.showwarning("ATTENZIONE","Dati errati o mancanti\nElaborazione non possibile",detail=detail_string)
    
       

    def _on_submit_button_click(self):
        """Gestisce il click del pulsante Submit."""
        if self.validations.all_data_ok:
            self.submitted_successfully = True # Imposta il flag a True
            self.root.quit() # Chiude la finestra Tkinter per permettere al programma principale di continuare
        else:
            self._data_warning()
            


    def run(self):
        self.root.mainloop() 




def interactive_image_viewer(p:Params, gamma_step=0.05):
    """
    Mostra un'immagine astronomica in modo interattivo con controlli per lo stretching 
    dei livelli (RangeSlider verticale) e la correzione Gamma (scroll del mouse).

    Args:
        image_data (np.ndarray): Array NumPy (2D) contenente i dati dell'immagine.
        gamma_step (float): Il passo di incremento/decremento della correzione Gamma.
    """
    
    # Rimuovi i valori NaN e assicurati che i dati siano float32
    data = p.imn[0:p.x_upper_lim-p.x_lower_lim,0:p.y_upper_lim-p.y_lower_lim]
    data = cv2.flip(data.astype(np.float32),0)
    data[np.isnan(data)] = np.min(data[~np.isnan(data)])

    # --- Calcolo Range Iniziale ---
    data_min = np.min(data)
    data_max = np.max(data)

    print('data_min',data_min,'data_max',data_max)

    if not(data_max>data_min):
        messagebox.showwarning("Attenzione","L'elaborazione ha prodotto un'immagine monocromatica")
        plt.imshow(data,'gray')
        plt.show()
        return
    
    # --- Variabili di Stato Locali ---
    # Usiamo un wrapper per tenere traccia dello stato di Gamma
    class State:
        gamma_val = 1.0
        img_data = data
    
    state = State()
    

    
    

    # --- Funzioni di Elaborazione (Nested) ---    
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
        out_path_recommended=name+"_enhanced"+"Division by 1_rho profile"
    if p.option == OPTIONS[2]:
        out_path_recommended=out_path_recommended+"_nrad_"+str(p.nrad)+"_ntheta_"+str(p.ntheta)+"_stdtheta_"+str(1/p.rejsig)+"_nsig_"+str(p.nsig)+"_"
    if p.option == OPTIONS[4]:
        out_path_recommended=out_path_recommended+"_A_"+str(p.A)+"_B_"+str(p.B)+"_N_"+str(p.N)+"_NUMLOG_"+str(p.NUMLOG)+"_"

    

    #root = tk.Tk()
    #root.withdraw() # Nasconde la finestra principale vuota

    # 2. Definisci i tipi di file filtrabili
    filetypes = [("Immagini FITS", "*.fit *.fits")] 
    
    # 3. Apri la finestra di dialogo "Salva con nome"
    filepath = filedialog.asksaveasfilename(
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

    p.hdul[0].header['INP_IM']=(p.input_path,'input image name')
    p.hdul[0].header['NUC-X']=(p.xnuc,'optocenter X pixel value')
    p.hdul[0].header['NUC-Y']=(p.ynuc,'optocenter Y pixel value')


          
        
    if p.option ==OPTIONS[4]:
        if p.NUMLOG:
            p.hdul[0].header['NUMLOG']=(1,'numlog=1 => logarithmic image')
        else:
            p.hdul[0].header['NUMLOG']=(0,'numlog ne 1 => no logarithmic image')
        p.hdul[0].header['A']=(p.A,'coeffcicent a')
        p.hdul[0].header['B']=(p.B,'coefficent b')
        
        p.hdul[0].header['A']=(p.A,'coeffcicent a')
    elif p.option != OPTIONS[3]:
        p.hdul[0].header['ENH-RAD']=(p.nrad,'diameter=1+(2*radius)')




    p.hdul.writeto(filepath,overwrite=True)




if __name__ == "__main__":
    try:
        root = tk.Tk()
        root.resizable(False, False)
        app = ImageProcessingGUI(root)
        while True:
            app.submitted_successfully=False
            app.run()
            
            # --- Recupero dei valori dopo che la GUI è stata chiusa (dal submit) ---
            if app.submitted_successfully:
                
                print("\n--- Dati recuperati dalla GUI ---")
                print(f"Percorso immagine: {app.image_path}")
                print(f"Larghezza immagine: {app.image_width}")
                print(f"Altezza immagine: {app.image_height}")
                print(f"Centro X: {float(app.center_x_var.get())}")
                print(f"Centro Y: {float(app.center_y_var.get())}")
                print(f"X Min: {int(app.xmin_var.get())}")
                print(f"X Max: {int(app.xmax_var.get())}")
                print(f"Y Min: {int(app.ymin_var.get())}")
                print(f"Y Max: {int(app.ymax_var.get())}")
                print(f"Opzione selezionata: {app.combobox_choice.get()}")

                # Recupera i parametri condizionali in base all'opzione scelta
                selected_option = app.combobox_choice.get()
                
                o=app.params_process(selected_option)
                
                (imold_preprocessed, xnuc_rel, ynuc_rel) = comet_pack.preprocess_and_normalize_crop(app.imold,o.xnuc,o.ynuc,o.x_lower_lim,o.x_upper_lim,o.y_lower_lim,o.y_upper_lim)
                print(imold_preprocessed.shape,"GAUD")
                if selected_option == OPTIONS[0]:
                    print(f"  Rho Pixels: {int(app.rho_pixels_var.get())}")
                    print(f"  Theta Pixels: {int(app.theta_pixels_var.get())}")
                    print(f"  Std Dev Theta: {float(app.std_dev_theta_var.get())}")
                    
                    #TODO: AGGIUNGERE CHE AL CROP IL NUCLEO DEVE STARE DENTRO
                    imun = comet_pack.polarize(imold_preprocessed,o.nrad,o.ntheta,xnuc_rel,xnuc_rel)          
                    imien=comet_pack.azimuthal_average_division_vectorized(imun,o.rejsig)
                    
                    o.imn=comet_pack.reconstruct_from_polar(imien,o.NCOL,o.NROW,xnuc_rel,ynuc_rel)
                    
                elif selected_option == OPTIONS[1]:
                    print(f"  Rho Pixels: {int(app.rho_pixels_var.get())}")
                    print(f"  Theta Pixels: {int(app.theta_pixels_var.get())}") 
                    
                    imun = comet_pack.polarize(imold_preprocessed,o.nrad,o.ntheta,xnuc_rel,ynuc_rel)          
                    imien = comet_pack.azimuthal_median_division_vectorized(imun) 
                    o.imn=comet_pack.reconstruct_from_polar(imien,o.NCOL,o.NROW,xnuc_rel,ynuc_rel)

                        
                
                elif selected_option == OPTIONS[2]:#RENORMALIZATION
                    print(f"  Rho Pixels: {int(app.rho_pixels_var.get())}")
                    print(f"  Theta Pixels: {int(app.theta_pixels_var.get())}")
                    print(f"  Std Dev Theta: {float(app.std_dev_theta_var.get())}")
                    print(f"  Min/Max Std Dev: {float(app.min_max_std_dev_var.get())}")
                    
                    imiun = comet_pack.polarize(imold_preprocessed,o.nrad,o.ntheta,xnuc_rel,ynuc_rel)            
                    imien=comet_pack.azimuthal_renormalization_vectorized(imiun,o.rejsig,o.nsig)
                    o.imn=comet_pack.reconstruct_from_polar(imien,o.NCOL,o.NROW,xnuc_rel,ynuc_rel)

                elif selected_option == OPTIONS[3]:
                    o.imn=comet_pack.inverserho_vectorized(imold_preprocessed,xnuc_rel,ynuc_rel,o.NCOL,o.NROW)    

                elif selected_option == OPTIONS[4]:
                    print(f"  Kernel A Term: {float(app.kernel_a_term_var.get())}")
                    print(f"  Kernel B Term: {float(app.kernel_b_term_var.get())}")
                    print(f"  Kernel N Term: {float(app.kernel_n_term_var.get())}")
                    print(f"  Transform Log: {app.transform_log_var.get()}")
                    o.imn = comet_pack.radially_variable_spatial_filtering_vectorized(imold_preprocessed,o.A,o.B,o.N,o.NUMLOG,xnuc_rel,ynuc_rel,o.NCOL,o.NROW)
                app.submit_button.config(state=tk.DISABLED,cursor="arrow")
                #print("SHAPE USCITA",o.imn.shape) 
                interactive_image_viewer(o)
                #TODO: GESTIRE ERRORE DI CHIUSURA APP PRIMA DEL VISUALIZZATORE
                #ho notato che se si chiama interactive_image_viewer senza cambiare stato al bottone
                #l'applicazione continua a lavorare
                #magari si può gestire la chiamata all'elaborazione direttamente dalla funzione del bottone
                #così da rendere il flusso dell'applicazione indipendente dal resto


                try:
                    app.submit_button.config(state=tk.NORMAL,cursor="hand2")
                    
                except tk.TclError as e:
                    print("\nL'utente ha chiuso la finestra")
                    sys.exit()
                
                
                
                
                
            else:
                print("\nL'utente ha chiuso la finestra")
                sys.exit()
    except MemoryError:
        messagebox.showerror("Errore","La memoria del sistema non è sufficiente a completare l'elaborazione")
        print("Errore: Memoria esaurita! (MemoryError)")
    except Exception as e:
        messagebox.showerror("Errore",f"Si è verificato un altro errore: {e}")
        print(f"Si è verificato un altro errore: {e}")
    


#TODO: rho deve essere limitato anche dal crop
