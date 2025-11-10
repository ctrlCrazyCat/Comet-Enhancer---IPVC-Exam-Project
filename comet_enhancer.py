
from  gui_app import ImageProcessingGUI,OPTIONS,interactive_image_viewer
from tkinter import  messagebox
import tkinter as tk
import comet_pack
import sys
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
                
                (imold_preprocessed, xnuc_rel, ynuc_rel) = comet_pack.preprocess_normalize_and_crop(app.imold,o.xnuc,o.ynuc,o.x_lower_lim,o.x_upper_lim,o.y_lower_lim,o.y_upper_lim)
                
                if selected_option == OPTIONS[0]:
                    print(f"  Rho Pixels: {int(app.rho_pixels_var.get())}")
                    print(f"  Theta Pixels: {int(app.theta_pixels_var.get())}")
                    print(f"  Std Dev Theta: {float(app.std_dev_theta_var.get())}")
                    

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

                interactive_image_viewer(o)



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