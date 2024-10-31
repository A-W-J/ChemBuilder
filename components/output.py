import tkinter as tk
from tkinter import filedialog
import os

from rdkit import Chem
from rdkit.Chem import Draw

class Output:
    def __init__(self, master):
        self.master = master
        self.molecules = None
        self.output_directory = None
        self.label = tk.Label(master, text = "path to output directory")
        self.label.pack(pady = 5)
        self.directory_entry = tk.Entry(master, width = 50)
        self.directory_entry.pack(pady = 10)
        #button to open the filepath browser dialog
        self.browse_button = tk.Button(master, text = "browse", command = self.browse_directory)
        self.browse_button.pack(pady = 10)
        #button to download all molecule images to the given directory
        self.download_all_button = tk.Button(master, text = "download all molecules to given directory", command = self.download_all)
        self.download_all_button.pack(pady = 10)
    def browse_directory(self):
        self.output_directory = filedialog.askdirectory()
        if self.output_directory:
            self.directory_entry.delete(0, tk.END)
            self.directory_entry.insert(0, self.output_directory)
    @staticmethod
    def make_filename(mol_name, counter):
        id = str(counter)
        translation_table = str.maketrans('<>:"/\\|?*', '_' * len('<>:"/\\|?*'))
        cleaned_string = mol_name.translate(translation_table)
        string_with_id = cleaned_string + '_' + id
        filename = '%s.png' % string_with_id
        return filename
    
    def download_all(self):
        if self.molecules and self.output_directory:
            os.chdir(self.output_directory)
            counter = 0
            for mol in self.molecules:
                filename = Output.make_filename(mol.name, counter)
                Draw.MolToFile(mol.structure, filename, size = (1000, 1000))
                counter = counter + 1





