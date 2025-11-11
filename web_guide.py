import tkinter as tk
from tkinter import font
import webbrowser

# --- Funzione per il browser ---
def _apri_link(url):
    """Apre l'URL specificato nel browser."""
    webbrowser.open_new_tab(url)

# --- Funzione che CREA la finestra dei link (MODIFICATA) ---
def apri_finestra_link(root):
    """Crea e visualizza la finestra secondaria con i link."""
    
    # --- MODIFICA 1: Disabilita e sbiadisci root ---
    root.attributes('-disabled', True) 
    root.attributes('-alpha', 0.95)
    
    
    window_link = tk.Toplevel(root)
    window_link.title("Guida")
    window_link.geometry("400x180")
    window_link.resizable(False, False)
    window_link.transient(root) 
    window_link.grab_set()
    
    # --- Creazione del Widget Text ---
    bg_color = window_link.cget('bg')
    text_widget = tk.Text(
        window_link, 
        wrap="word",
        borderwidth=0,
        highlightthickness=0,
        bg=bg_color,
        font=("Arial", 11)
    )
    text_widget.pack(padx=15, pady=15, fill="both", expand=True)

    # --- Definiamo gli stili (Tag) ---
    link_font_sottolineato = font.Font(family="Arial", size=11, underline=True)

    # Tag 1
    text_widget.tag_configure("link1_tag", foreground="blue")
    text_widget.tag_bind("link1_tag", "<Button-1>", lambda e: _apri_link("https://cie.psi.edu/techniques.html"))
    text_widget.tag_bind("link1_tag", "<Enter>", lambda e: text_widget.config(cursor="hand2"))
    text_widget.tag_bind("link1_tag", "<Leave>", lambda e: text_widget.config(cursor=""))

    # Tag 2
    text_widget.tag_configure("link2_tag", foreground="blue", font=link_font_sottolineato)
    text_widget.tag_bind("link2_tag", "<Button-1>", lambda e: _apri_link("https://www.sciencedirect.com/science/article/abs/pii/S003206331500197X"))
    text_widget.tag_bind("link2_tag", "<Enter>", lambda e: text_widget.config(cursor="hand2"))
    text_widget.tag_bind("link2_tag", "<Leave>", lambda e: text_widget.config(cursor=""))

    # Tag 3
    text_widget.tag_configure("link3_tag", foreground="blue", font=link_font_sottolineato)
    text_widget.tag_bind("link3_tag", "<Button-1>", lambda e: _apri_link("https://cie.psi.edu/Description%20of%20the%20Enhancement%20Routines%20used%20in%20the%20Online%20Facility.pdf"))
    text_widget.tag_bind("link3_tag", "<Enter>", lambda e: text_widget.config(cursor="hand2"))
    text_widget.tag_bind("link3_tag", "<Leave>", lambda e: text_widget.config(cursor=""))

    # --- Inserimento del testo e applicazione dei Tag ---
    text_widget.insert(tk.END, "Questa applicazione è una reimplementazione in Python degli algoritmi presenti sul sito ")
    text_widget.insert(tk.END, "https://cie.psi.edu/", "link1_tag")
    text_widget.insert(tk.END, " e pubblicati in ")
    text_widget.insert(tk.END, "questo articolo", "link2_tag")
    text_widget.insert(tk.END, ".\n\nPer un approfondimento sugli algoritmi presentati cliclla il ")
    text_widget.insert(tk.END, "seguente link", "link3_tag")
    text_widget.insert(tk.END, ".")
    
    text_widget.config(state="disabled")
    
    
    # --- MODIFICA 2: Attendi che la finestra venga chiusa ---
    root.wait_window(window_link)
    
    
    # --- MODIFICA 3: Riabilita e ripristina root ---
    root.attributes('-disabled', False) 
    root.attributes('-alpha', 1.0)

    # --- MODIFICA 4: Forza la finestra root in primo piano ---
    root.lift()
    root.focus_force()
