import tkinter as tk
from tkinter import messagebox

from components.input import *
from components.image_generator import *
from components.browser import *
from components.output import *

from components.shared_globals import *

class MainWindow:
    def __init__(self, master):
        self.master = master
        self.master.title("SMILES to Image Generator")

        shared_globals.app = self

        # Frame for Input class
        self.input_frame = tk.Frame(master)
        self.input_frame.pack(pady=10)
        self.input_app = Input(self.input_frame)

        # Button to proceed to ImageGenerator
        self.image_generator_app = None
        self.generate_button = tk.Button(master, text="Generate Images", command=self.open_image_generator)
        self.generate_button.pack(pady=10)

        # Storing the browser window
        self.browser_app = None

        # Button to re-open a closed browser window
        self.browser_button = tk.Button(master, text = "Open image browser", command = self.open_browser_popup)
        self.browser_button.pack(pady = 10)

        # Frame for the Output class
        self.output_frame = tk.Frame(master)
        self.output_frame.pack(pady = 10)
        self.output_app = Output(self.output_frame)

    def open_image_generator(self):
        if self.input_app.csv_file:
            # Create the ImageGenerator class
            self.image_generator_app = ImageGenerator(self.master, self.input_app.csv_file)
            self.image_generator_app.create_molecules()

            # Open the Browser in a new pop-out window
            self.open_browser_popup()
        else:
            messagebox.showerror("Error", "Please select a CSV file first.")

    def open_browser_popup(self):
        if self.image_generator_app:
            # Create a new pop-up window for the browser
            browser_window = tk.Toplevel(self.master)
            browser_window.title("Molecule Browser")

            # Create the browser application in the new window
            self.browser_app = Browser(browser_window, self.image_generator_app.molecules)

# Main Tkinter loop
if __name__ == "__main__":
    root = tk.Tk()
    app = MainWindow(root)
    root.mainloop()