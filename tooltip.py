import tkinter as tk

class Tooltip:
    """
    Classe Tooltip rivista con logica di attivazione/disattivazione 
    basata sullo stato 'is_enabled' per evitare problemi di unbind.
    """
    def __init__(self, widget, text,tooltip_enabled = True):
        self.widget = widget
        self.text = text
        self.tw = None
        self.id = None
        self.delay = 1000
        self.is_enabled = True # Stato di attivazione/disattivazione

        # E' fondamentale che i metodi 'schedule' e 'hide' siano definiti
        # PRIMA di queste righe. Nell'implementazione Python standard, questo è
        # garantito se i metodi sono definiti nell'ordine che segue l'__init__.

        # I binding restano ATTIVI. Il comportamento è controllato da self.is_enabled
        self.widget.bind("<Enter>", self.schedule)
        self.widget.bind("<Leave>", self.hide)
        self.widget.bind("<ButtonPress>", self.hide)
        if not tooltip_enabled:
            self.hide() # Assicura che la finestra sia chiusa
            self.is_enabled = False
        
    def unschedule(self, event=None):
        """Cancella il timer pendente (se presente)."""
        if self.id:
            self.widget.after_cancel(self.id)
            self.id = None

    def hide(self, event=None):
        """Nasconde e distrugge la finestra del tooltip."""
        # Deve essere definito prima di essere chiamato in __init__
        self.unschedule()
        
        if self.tw:
            self.tw.destroy()
        self.tw = None
        
    def schedule(self, event=None):
        """Pianifica la comparsa del tooltip solo se attivo."""
        if not self.is_enabled: 
            return
            
        self.unschedule()
        self.id = self.widget.after(self.delay, self.show)


    def show(self, event=None):
        """Crea e mostra la finestra pop-up del tooltip."""
        if self.tw:
            return

        x = self.widget.winfo_rootx() + 10
        #y = self.widget.winfo_rooty() + self.widget.winfo_height() + 1
        y = self.widget.winfo_rooty() - 10

        self.tw = tk.Toplevel(self.widget)
        self.tw.wm_overrideredirect(True)
        self.tw.wm_geometry(f"+{x}+{y}")
        
        label = tk.Label(self.tw, 
                          text=self.text, 
                          justify='left',
                          background="#ffffe0",
                          relief='solid', 
                          borderwidth=1,
                          font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    # --- Metodi di Attivazione/Disattivazione ---

    def disable(self):
        """Disattiva temporaneamente il tooltip."""
        self.hide() # Assicura che la finestra sia chiusa
        self.is_enabled = False
        print(f"Tooltip per {self.widget} disattivato.")

    def enable(self):
        """Riattiva il tooltip."""
        self.is_enabled = True
        print(f"Tooltip per {self.widget} riattivato.")
        
    def destroy_tooltip(self):
        """Rimuove permanentemente il tooltip."""
        self.disable()


# ====================================================================
# Esempio di Utilizzo per Attivare/Disattivare
# ====================================================================

if __name__ == '__main__':
    root = tk.Tk()
    root.title("Tkinter Tooltip Toggle Example")

    # Widget con Tooltip da controllare
    control_widget = tk.Label(root, text="Controlla il Tooltip Qui!", padx=10, pady=10, bg="lightblue")
    control_widget.pack(padx=50, pady=20)
    
    tooltip_text = "Questo tooltip può essere attivato e disattivato a piacere."
    # Salviamo l'istanza del Tooltip per chiamare i metodi enable/disable
    controlled_tooltip = Tooltip(control_widget, tooltip_text) 

    # --- Funzioni di Controllo ---
    def toggle_tooltip():
        """Alterna lo stato di attivazione/disattivazione del tooltip."""
        if controlled_tooltip.is_enabled:
            controlled_tooltip.disable()
            toggle_btn.config(text="Riattiva Tooltip")
            #control_widget.config(bg="gray")
        else:
            controlled_tooltip.enable()
            toggle_btn.config(text="Disattiva Tooltip")
            #control_widget.config(bg="lightblue")


    # Bottone per attivare/disattivare
    toggle_btn = tk.Button(root, text="Disattiva Tooltip", command=toggle_tooltip)
    toggle_btn.pack(pady=10)
    
    root.mainloop()