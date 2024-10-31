import pandas as pd
import tkinter as tk
from tkinter import messagebox
from rdkit import Chem
from rdkit.Chem import Draw

from components.shared_globals import *

class Molecule:
    def __init__ (self, name, smile, structure, image):
        self.name = name
        self.smile = smile
        self.structure = structure
        self.image = image

class ImageGenerator:
    def __init__(self, master, csv_file):
        self.master = master
        self.csv_file = csv_file #filepath (from Input class)
        self.data = None #dataframe
        self.smiles_column = None #name of the column containing the SMILES strucutres
        self.name_column = None #name of the column containing the molecular names
        self.molecules = None #list of molecule objects assembled from the CSV file

        # Load the CSV file
        self.load_csv()
        # Try to automatically detect columns
        self.detect_columns()

        # If auto-detection fails, show the GUI for manual column input
        if not self.smiles_column or not self.name_column:
            self.create_manual_input_gui()
    def load_csv(self):
        """Load the CSV file and display its columns."""
        self.data = pd.read_csv(self.csv_file)
        print("Columns detected:", self.data.columns.tolist())
    def detect_columns(self):
        """Try to automatically detect SMILES and Name columns."""
        possible_smiles = ['smiles', 'SMILES', 'Smiles']
        possible_names = ['name', 'Name', 'Molecule', 'molecule']

        # Check for potential SMILES column
        for col in self.data.columns:
            if col in possible_smiles:
                self.smiles_column = col
                break

        # Check for potential Name column
        for col in self.data.columns:
            if col in possible_names:
                self.name_column = col
                break

        if self.smiles_column and self.name_column:
            print(f"Detected SMILES column: {self.smiles_column}")
            print(f"Detected Name column: {self.name_column}")
        else:
            print("Automatic detection failed.")
    def create_manual_input_gui(self):
        """Create GUI for manual input of SMILES and Name columns."""
        self.label_smiles = tk.Label(self.master, text="Enter the SMILES column name:")
        self.label_smiles.pack(pady=5)

        self.entry_smiles = tk.Entry(self.master)
        self.entry_smiles.pack(pady=5)

        self.label_name = tk.Label(self.master, text="Enter the Name column name:")
        self.label_name.pack(pady=5)

        self.entry_name = tk.Entry(self.master)
        self.entry_name.pack(pady=5)

        self.submit_button = tk.Button(self.master, text="Submit", command=self.submit_columns)
        self.submit_button.pack(pady=5)
    def submit_columns(self):
        """Store manually entered column names and check validity."""
        self.smiles_column = self.entry_smiles.get()
        self.name_column = self.entry_name.get()

        if self.smiles_column not in self.data.columns or self.name_column not in self.data.columns:
            messagebox.showerror("Error", "Invalid column names. Please check the CSV file.")
        else:
            messagebox.showinfo("Success", "Columns successfully set!")
            print(f"Manually set SMILES column: {self.smiles_column}")
            print(f"Manually set Name column: {self.name_column}")
    def create_molecules(self):
        global app
        """Extract SMILES and names from the specified columns."""
        if not self.smiles_column or not self.name_column:
            raise ValueError("SMILES or Name columns are not set.")

        smiles = self.data[self.smiles_column]
        names = self.data[self.name_column]

        molecules = []

        for smile, name in zip(smiles, names):
            structure = Chem.MolFromSmiles(str(smile))
            image = Draw.MolToImage(structure)
            molecules.append(Molecule(name, smile, structure, image))
        self.molecules = molecules

        #we can then attribute the same set of molecules to the output portion of the app
        shared_globals.app.output_app.molecules = molecules
        
