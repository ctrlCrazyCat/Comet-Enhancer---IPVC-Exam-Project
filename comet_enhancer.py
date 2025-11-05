import tkinter as tk
from tkinter import filedialog, messagebox, ttk,font as tkFont

from astropy.io import fits
import matplotlib.pyplot as plt
from comet_pack import Params
import comet_pack
import numpy as np
from matplotlib.widgets import RangeSlider, Slider
import os
from matplotlib.cm import get_cmap
import matplotlib
import sys
matplotlib.use('TkAgg')

from tooltip import Tooltip

def conditional_config(entry_widget_info,configuration):
        (condition,text) = configuration
        entry_widget_info.config(fg='green')
        if not condition:
            entry_widget_info.config(fg='red')
        entry_widget_info.tooltip.text = text
def _get_info_label(frame,tooltip_enabled = False):
    out = tk.Label(frame, text="l",font=tkFont.Font(family="Wingdings", size=12, weight="bold"),fg="gray")
    out.tooltip= Tooltip(out,"",tooltip_enabled)
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
        self.theta_pixels_var = tk.StringVar()
        self.std_dev_theta_var = tk.StringVar()
        self.min_max_std_dev_var = tk.StringVar()
        self.kernel_a_term_var = tk.StringVar()
        self.kernel_b_term_var = tk.StringVar()
        self.kernel_n_term_var = tk.StringVar()
        self.transform_log_var = tk.BooleanVar()

        # --- Variabile per il Combobox ---
        self.combobox_choice = tk.StringVar()
        self.options = comet_pack.get_options()
        self.combobox_choice.set(self.options[0])

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
        self.select_image_button = tk.Button(root, text="Scegli Immagine FIT/FITS", command=self._select_image) # Modificato il testo
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

        # Frame per i Limiti (Int)
        self.limits_frame = ttk.LabelFrame(root, text="Limiti di Elaborazione (X Min/Max, Y Min/Max - int)")
        self.limits_frame.pack(pady=10, padx=10, fill=tk.X)

        # Riga X Min/Max
        self.xminmax_frame = tk.Frame(self.limits_frame)
        self.xminmax_frame.pack(fill=tk.X, pady=2)

             
        
        self.entry_xmin = _get_custom_entry(self.xminmax_frame,self.xmin_var,"X Min:")
        self.all_entries.append(self.entry_xmin)

        tk.Label(self.xminmax_frame).pack(side=tk.LEFT,padx=30)

        self.entry_xmax = _get_custom_entry(self.xminmax_frame,self.xmax_var,"X Max:")
        self.all_entries.append(self.entry_xmax)
        

        # Riga Y Min/Max
        self.yminmax_frame = tk.Frame(self.limits_frame)
        self.yminmax_frame.pack(fill=tk.X, pady=2)
        
        self.entry_ymin = _get_custom_entry(self.yminmax_frame,self.ymin_var,"Y Min:")
        self.all_entries.append(self.entry_ymin)

        tk.Label(self.yminmax_frame).pack(side=tk.LEFT,padx=30)

        self.entry_ymax = _get_custom_entry(self.yminmax_frame,self.ymax_var,"Y Max:")  
        self.all_entries.append(self.entry_ymax)

        # Combobox
        self.style = ttk.Style()
        self.style.configure('TCombobox', fieldbackground='white', background='lightgrey')
        self.combobox = ttk.Combobox(root, textvariable=self.combobox_choice, values=self.options, state='readonly', width=40)
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
                "variable": self.rho_pixels_var
            },
            "theta_pixels": {
                "label": tk.Label(self.conditional_params_frame, text="Numero di pixel asse θ:"), 
                "info_label":_get_info_label(self.conditional_params_frame,tooltip_enabled=True),
                "widget": tk.Entry(self.conditional_params_frame, textvariable=self.theta_pixels_var, width=20),
                "variable": self.theta_pixels_var
            },
            "std_dev_theta": {
                "label": tk.Label(self.conditional_params_frame, text="Quante deviazioni standard dovrebbero essere accettate per i pixel asse θ?"), 
                "info_label":_get_info_label(self.conditional_params_frame,tooltip_enabled=True),
                "widget": tk.Entry(self.conditional_params_frame, textvariable=self.std_dev_theta_var, width=20),
                "variable": self.std_dev_theta_var
            },
            "min_max_std_dev": {
                "label": tk.Label(self.conditional_params_frame, text="How many standard deviations from the mean should\nthe minimum and maximum pixel values be?"), 
                "info_label":_get_info_label(self.conditional_params_frame,tooltip_enabled=True),
                "widget": tk.Entry(self.conditional_params_frame, textvariable=self.min_max_std_dev_var, width=20),
                "variable": self.min_max_std_dev_var
            },
            "kernel_a_term": {
                "label": tk.Label(self.conditional_params_frame, text="valore Kernel A (double):"), 
                "info_label":_get_info_label(self.conditional_params_frame,tooltip_enabled=True),
                "widget": tk.Entry(self.conditional_params_frame, textvariable=self.kernel_a_term_var, width=20),
                "variable": self.kernel_a_term_var
            },
            "kernel_b_term": {
                "label": tk.Label(self.conditional_params_frame, text="valore Kernel B (double):"), 
                "info_label":_get_info_label(self.conditional_params_frame,tooltip_enabled=True),
                "widget": tk.Entry(self.conditional_params_frame, textvariable=self.kernel_b_term_var, width=20),
                "variable": self.kernel_b_term_var
            },
            "kernel_n_term": {
                "label": tk.Label(self.conditional_params_frame, text="valore Kernel N (double):"), 
                "info_label":_get_info_label(self.conditional_params_frame,tooltip_enabled=True),
                "widget": tk.Entry(self.conditional_params_frame, textvariable=self.kernel_n_term_var, width=20),
                "variable": self.kernel_n_term_var
            },
            "transform_log": {
                "label": tk.Label(self.conditional_params_frame, text="Trasformare l'immagine input in scala log-10 prima di migliorarla?"), 
                "info_label":_get_info_label(self.conditional_params_frame),
                "widget": tk.Checkbutton(self.conditional_params_frame, variable=self.transform_log_var, text="Sì"),
                "variable": self.transform_log_var
            }
        }
        #TODO: CREARE UN CONDITIONAL CONFIG SOLO PER RHO PIXEL

        # Mappa le opzioni ai widget necessari e alla loro riga di griglia
        self.option_to_widgets_config = {
            self.options[0]: [ # Division by Azimuthal Average
                ("rho_pixels", 0), ("theta_pixels", 1), ("std_dev_theta", 2)
            ],
            self.options[1]: [ # Division by Azimuthal Median
                ("rho_pixels", 0), ("theta_pixels", 1)
            ],
            self.options[2]: [ # Azimuthal Renormalization
                ("rho_pixels", 0), ("theta_pixels", 1), ("std_dev_theta", 2), ("min_max_std_dev", 3)
            ],
            self.options[3]: [ # Division by 1/rho profile #TODO: Modificare questa cosa 
                #("rho_pixels", 0), ("theta_pixels", 1)
            ],
            self.options[4]: [ # Radially Variable Spatial Filtering
                ("kernel_a_term", 0), ("kernel_b_term", 1), ("kernel_n_term", 2), ("transform_log", 3)
            ]
        }
        
        # Pulsante Submit (inizialmente nascosto e disabilitato)
        grande_font = tkFont.Font(family="Arial", size=20, weight="bold")
        self.submit_button = tk.Button(root, text="Submit",font=grande_font, command=self._on_submit_button_click) # Modificato il command

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
        [hdul,imold]=comet_pack.get_input_data(self.image_path)
        o.hdul = hdul
        o.imold = imold
        (o.NROW,o.NCOL)= imold.shape
        if option==self.options[0]:
            o.nrad=int(self.rho_pixels_var.get())
            o.ntheta=int(self.theta_pixels_var.get())

            try:
                o.rejsig=1/float(self.std_dev_theta_var.get())
            except ZeroDivisionError as e:
                o.rejsig = float('inf')

        if selected_option == app.options[1]:
            o.nrad=int(self.rho_pixels_var.get())
            o.ntheta=int(self.theta_pixels_var.get()) 

        if option==self.options[2]:
            o.nrad=int(self.rho_pixels_var.get())
            o.ntheta=int(self.theta_pixels_var.get())
            o.rejsig=1/float(self.std_dev_theta_var.get())
            o.nsig=float(self.min_max_std_dev_var.get())
        
        if option==self.options[4]:
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
            self.image_path = file_path
            try:
                with fits.open(file_path) as hdul: #hdul,imold,_upper_,y_upper_lim
                    self.hdul = hdul.copy()
                    self.imold = hdul[0].data.copy()
                    self.image_width=hdul[0].header['NAXIS1']
                    self.image_height=hdul[0].header['NAXIS2']
                    
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
            
            info_label.grid_remove()
            label.grid_remove() 
            widget.grid_remove()
            
            if isinstance(widget, tk.Entry):
                widget.config(state=tk.DISABLED)
                #variable.set("") # Resetta la StringVar
            elif isinstance(widget, tk.Checkbutton):
                widget.config(state=tk.DISABLED)
                variable.set(False) # Resetta BooleanVar

        # 2. Mostra/Nasconde il FRAME e popola i widget in base alla selezione
        selected_option = self.combobox_choice.get()
        print(selected_option,"AAA")

        if selected_option == self.options[3]:
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
                
                info_label.grid(row=row_num,column=0,padx=0,pady=0,sticky="w")
                label.grid(row=row_num, column=1, padx=5, pady=5, sticky="w")
                widget.grid(row=row_num, column=2, padx=5, pady=5, sticky="ew")
                widget.config(state=tk.NORMAL)
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
            e.info_label.config(fg='gray')
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
        if not value_str: return (False,"è assente")
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
    def _validate_rho_input_and_approximate(self,all_center_ok:bool):
        if not all_center_ok:
            return (False,"dati sul Centro Nucleo Cometa non correttamente inseriti o assenti")
        return self._validate_int_input_and_approximate(self.rho_pixels_var,1,self.rho_max)
    

    def _validate_int_input_and_approximate(self, value_var, min_val=1, max_val=None):
        """
        Valida se una stringa può essere convertita a int (con approssimazione)
        e se rientra nel range [min_val, max_val] (se max_val è fornito).
        Aggiorna la StringVar con l'intero approssimato.
        """
        value_str = value_var.get()
        if not value_str: return (False,"è assente")

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
        if all_center_ok:
            xnuc = float(self.center_x_var.get())
            ynuc = float(self.center_y_var.get())
            xmax = max((self.image_width-xnuc),xnuc)
            ymax =max((self.image_height-ynuc),ynuc)
            self.rho_max=int(np.floor((xmax**2+ymax**2)**0.5))
        
        

        xmin_ok = self._validate_int_input_and_approximate(self.xmin_var, min_val=1, max_val=self.image_width)
        all_ok.append(xmin_ok)
        xmax_ok = self._validate_int_input_and_approximate(self.xmax_var, min_val=1, max_val=self.image_width)
        all_ok.append(xmax_ok)
        ymin_ok = self._validate_int_input_and_approximate(self.ymin_var, min_val=1, max_val=self.image_height)
        all_ok.append(ymin_ok)
        ymax_ok = self._validate_int_input_and_approximate(self.ymax_var, min_val=1, max_val=self.image_height)
        all_ok.append(ymax_ok)
        all_limits_ok = xmin_ok[0] and xmax_ok[0] and ymin_ok[0] and ymax_ok[0]
        

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

        main_inputs_valid = all_center_ok and all_limits_ok and min_max_ok

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
        #TODO: modificare condizione di rho in modo che veda anche i paramtri del centro 
        
        if selected_option == self.options[0]: # Division by Azimuthal Average
            rho_pixels_ok = self._validate_rho_input_and_approximate(all_center_ok)
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

        elif selected_option == self.options[1]: # Azimuthal Median / Division by 1/rho profile
            rho_pixels_ok = self._validate_rho_input_and_approximate(all_center_ok)
            conditional_config(self.conditional_widgets_map["rho_pixels"]["info_label"],rho_pixels_ok)
            theta_pixels_ok = self._validate_int_input_and_approximate(self.theta_pixels_var, min_val=1)
            conditional_config(self.conditional_widgets_map["theta_pixels"]["info_label"],theta_pixels_ok)
            all_conditional_ok = rho_pixels_ok[0] and theta_pixels_ok[0]
            self.validations.option = 1
            self.validations.all_conditional_ok = all_conditional_ok
            if not all_conditional_ok:
                self.validations.rho_pixels_ok = rho_pixels_ok
                self.validations.theta_pixels_ok = theta_pixels_ok
        elif selected_option==self.options[3]:
            self.validations.option = 3
            all_conditional_ok = True
                
            
        elif selected_option == self.options[2]: # Azimuthal Renormalization
            rho_pixels_ok = self._validate_rho_input_and_approximate(all_center_ok)
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


        elif selected_option == self.options[4]: # Radially Variable Spatial Filtering
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
        #TODO: VEDERE SE SI PUO OTTIMIZZARE IL CODICE ANCHE QUI

        if not self.validations.all_conditional_ok:
            detail_string = detail_string+"\n"
            match self.validations.option:
                case 0:
                    if not self.validations.rho_pixels_ok[0]:
                        detail_string = detail_string + "-Numero di pixel asse ρ (rho) "+self.validations.rho_pixels_ok[1]+"\n"
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

if __name__ == "__main__":
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
            if selected_option == app.options[0]:
                print(f"  Rho Pixels: {int(app.rho_pixels_var.get())}")
                print(f"  Theta Pixels: {int(app.theta_pixels_var.get())}")
                print(f"  Std Dev Theta: {float(app.std_dev_theta_var.get())}")
                print("SIZE",o.imold.shape)
                imun = comet_pack.polarize(o.imold,o.nrad,o.ntheta,o.xnuc,o.ynuc)            
                print("SIZE",imun.shape)
                imien=comet_pack.azimuthal_average_division(imun,o.rejsig)
                print("SIZE",imien.shape)
                o.imn=comet_pack.reconstruct_from_polar(imien,o.NCOL,o.NROW,o.xnuc,o.ynuc,o.x_lower_lim,o.x_upper_lim,o.y_lower_lim,o.y_upper_lim)
                print("SIZE",o.imn.shape)

            elif selected_option == app.options[1]:
                print(f"  Rho Pixels: {int(app.rho_pixels_var.get())}")
                print(f"  Theta Pixels: {int(app.theta_pixels_var.get())}") 
                imiun = comet_pack.polarize(o.imold,o.nrad,o.ntheta,o.xnuc,o.ynuc)     
                imien = comet_pack.azimuthal_median_division(imiun)
                o.imn=comet_pack.reconstruct_from_polar(imien,o.NCOL,o.NROW,o.xnuc,o.ynuc,o.x_lower_lim,o.x_upper_lim,o.y_lower_lim,o.y_upper_lim)

                    
            
            elif selected_option == app.options[2]:#RENORMALIZATION
                print(f"  Rho Pixels: {int(app.rho_pixels_var.get())}")
                print(f"  Theta Pixels: {int(app.theta_pixels_var.get())}")
                print(f"  Std Dev Theta: {float(app.std_dev_theta_var.get())}")
                print(f"  Min/Max Std Dev: {float(app.min_max_std_dev_var.get())}")
                
                imiun = comet_pack.polarize(o.imold,o.nrad,o.ntheta,o.xnuc,o.ynuc)            
                imien=comet_pack.azimuthal_renormalization(imiun,o.rejsig,o.nsig)
                o.imn=comet_pack.reconstruct_from_polar(imien,o.NCOL,o.NROW,o.xnuc,o.ynuc,o.x_lower_lim,o.x_upper_lim,o.y_lower_lim,o.y_upper_lim)

            elif selected_option == app.options[3]:
                o.imn=comet_pack.enhance_inverserho_vectorized(o.imold,o.xnuc,o.ynuc,o.x_lower_lim,o.x_upper_lim,o.y_lower_lim,o.y_upper_lim)    

            elif selected_option == app.options[4]:
                print(f"  Kernel A Term: {float(app.kernel_a_term_var.get())}")
                print(f"  Kernel B Term: {float(app.kernel_b_term_var.get())}")
                print(f"  Kernel N Term: {float(app.kernel_n_term_var.get())}")
                print(f"  Transform Log: {app.transform_log_var.get()}")
                o.imn = comet_pack.radially_variable_spatial_filtering(o.imold,o.A,o.B,o.N,o.NUMLOG,o.xnuc,o.ynuc,o.x_lower_lim,o.x_upper_lim,o.y_lower_lim,o.y_upper_lim)
            app.submit_button.config(state=tk.DISABLED)
            print("SHAPE USCITA",o.imn.shape) 
            comet_pack.interactive_image_viewer(o)
            #TODO: GESTIRE ERRORE DI CHIUSURA APP PRIMA DEL VISUALIZZATORE
            #ho notato che se si chiama interactive_image_viewer senza cambiare stato al bottone
            #l'applicazione continua a lavorare
            #magari si può gestire la chiamata all'elaborazione direttamente dalla funzione del bottone
            #così da rendere il flusso dell'applicazione indipendente dal resto


            try:
                app.submit_button.config(state=tk.NORMAL)
            except tk.TclError as e:
                print("\nL'utente ha chiuso la finestra")
                sys.exit()
            
            
            
            
            
        else:
            print("\nL'utente ha chiuso la finestra")
            sys.exit()
    