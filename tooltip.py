import tkinter as tk

class Tooltip:
    """
    Tooltip semplice (immediato) con stato enable/disable.
    
    Appare immediatamente all'evento <Enter> e scompare
    all'evento <Leave>, ma solo se lo stato 'is_enabled' è True.
    """
    def __init__(self, widget, text:str, tooltip_enabled=True):
        self.widget = widget
        self.text = text
        self.n_lines = 1 + text.count('\n')
        self.tw = None # Riferimento alla finestra Toplevel del tooltip

        self.estimated_line_height = 15 
        self.estimated_height = (self.n_lines * self.estimated_line_height) + 5 # + padding
        
        # Stato di attivazione
        self.is_enabled = tooltip_enabled 

        # I binding restano sempre attivi; il controllo
        # viene fatto all'interno dei metodi.
        self.widget.bind("<Enter>", self.show)
        self.widget.bind("<Leave>", self.hide)
        self.widget.bind("<ButtonPress>", self.hide)

    def show(self, event=None):
        """Mostra il tooltip solo se è abilitato."""
        
        # --- CONTROLLO CHIAVE ---
        # Se disabilitato, o se già mostrato, esci.
        if not self.is_enabled or self.tw:
            return

        # Calcola la posizione
        x = self.widget.winfo_rootx() + 10
        
        # Posiziona SOPRA:
        # Prendi la Y della cima del widget (winfo_rooty)
        # e sottrai l'altezza stimata del tooltip.
        y = self.widget.winfo_rooty() - self.estimated_height 

        # Controllo di sicurezza: se esce dallo schermo (in alto),
        # mettilo sotto.
        if y < 0:
            y = self.widget.winfo_rooty() + self.widget.winfo_height() + 1

        # Crea la finestra Toplevel
        self.tw = tk.Toplevel(self.widget)
        self.tw.wm_overrideredirect(True) # Rimuove bordi, titolo, ecc.
        self.tw.wm_geometry(f"+{x}+{y}")
        
        # Crea l'etichetta interna
        label = tk.Label(self.tw, 
                         text=self.text, 
                         justify='left',
                         background="#ffffe0", # Sfondo giallo pallido
                         relief='solid', 
                         borderwidth=1,
                         font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hide(self, event=None):
        """Nasconde e distrugge il tooltip."""
        # Questo metodo funziona indipendentemente dallo stato 'is_enabled'
        if self.tw:
            self.tw.destroy()
        self.tw = None

    # --- Metodi di Attivazione/Disattivazione ---

    def enable(self):
        """Attiva il tooltip."""
        self.is_enabled = True
        #print("Tooltip abilitato")

    def disable(self):
        """Disattiva il tooltip e lo nasconde se attualmente visibile."""
        self.is_enabled = False
        self.hide() # Nasconde se era già attivo
        #print("Tooltip disabilitato")


# ====================================================================
# Esempio di Utilizzo per Attivare/Disattivare
# ====================================================================

if __name__ == '__main__':
    root = tk.Tk()
    root.title("Tooltip Semplice (Abilita/Disabilita)")
    root.geometry("300x200")

    # Widget con Tooltip da controllare
    control_widget = tk.Label(root, text="Passa il mouse qui!", padx=20, pady=20, bg="lightblue")
    control_widget.pack(padx=50, pady=20)
    
    tooltip_text = "Questo è un tooltip immediato!"
    
    # Salviamo l'istanza del Tooltip
    controlled_tooltip = Tooltip(control_widget, tooltip_text) 

    # --- Funzioni di Controllo ---
    def toggle_tooltip():
        """Alterna lo stato di attivazione/disattivazione del tooltip."""
        if controlled_tooltip.is_enabled:
            controlled_tooltip.disable()
            toggle_btn.config(text="Attiva Tooltip")
            control_widget.config(bg="gray")
        else:
            controlled_tooltip.enable()
            toggle_btn.config(text="Disattiva Tooltip")
            control_widget.config(bg="lightblue")

    # Bottone per attivare/disattivare
    toggle_btn = tk.Button(root, text="Disattiva Tooltip", command=toggle_tooltip)
    toggle_btn.pack(pady=10)
    
    root.mainloop()