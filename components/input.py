import tkinter as tk
from tkinter import filedialog

class Input:
    def __init__(self, master):
        self.master = master

        # Create label and entry for CSV file path
        self.label = tk.Label(master, text="CSV File Path:")
        self.label.pack(pady=5)

        self.file_entry = tk.Entry(master, width=50)
        self.file_entry.pack(pady=5)

        # Create button to browse for CSV file
        self.browse_button = tk.Button(master, text="Browse", command=self.browse_file)
        self.browse_button.pack(pady=5)

        # Slot to store the selected CSV file
        self.csv_file = None

    def browse_file(self):
        # Open file dialog to select a CSV file
        self.csv_file = filedialog.askopenfilename(
            filetypes=[("CSV files", "*.csv")],
            title="Select a CSV file"
        )
        self.file_entry.delete(0, tk.END)
        self.file_entry.insert(0, self.csv_file)

# Main Tkinter loop
if __name__ == "__main__":
    root = tk.Tk()
    app = Input(root)
    root.mainloop()