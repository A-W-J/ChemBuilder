import tkinter as tk
from tkinter import Scrollbar, Canvas, Frame
from PIL import Image, ImageTk

class Browser:
    def __init__(self, master, molecules):
        self.master = master
        self.molecules = molecules
        self.image_labels = []

        # Create a frame to hold the canvas and scrollbars
        self.browser_frame = Frame(master)
        self.browser_frame.pack(fill="both", expand=True)

        # Create the canvas for image display
        self.canvas = Canvas(self.browser_frame)

        # Create the scrollbars
        self.scrollbar_y = Scrollbar(self.browser_frame, orient="vertical", command=self.canvas.yview)
        self.scrollbar_x = Scrollbar(self.browser_frame, orient="horizontal", command=self.canvas.xview)

        # Create a frame inside the canvas for scrolling content
        self.scrollable_frame = Frame(self.canvas)

        # Bind the frame's configuration to update the scroll region
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        )

        # Create the window inside the canvas to hold the scrollable frame
        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")

        # Configure canvas scroll commands
        self.canvas.configure(yscrollcommand=self.scrollbar_y.set, xscrollcommand=self.scrollbar_x.set)

        # Pack the canvas and scrollbars
        self.canvas.pack(side="left", fill="both", expand=True)
        self.scrollbar_y.pack(side="right", fill="y")
        self.scrollbar_x.pack(side="bottom", fill="x")

        # Display the images in a grid
        self.display_images()

    def display_images(self):
        """Display images in rows of three, each 300x300 pixels with a name caption."""
        max_image_size = (300, 300)  # Set the desired image size

        # Iterate through the molecules and place them in a grid
        row, col = 0, 0  # Initialize row and column for grid placement
        for i, molecule in enumerate(self.molecules):
            # Resize the image to 300x300 pixels
            resized_image = molecule.image.resize(max_image_size, Image.Resampling.LANCZOS)
            photo = ImageTk.PhotoImage(resized_image)

            # Create label for the molecule name (caption)
            label_name = tk.Label(self.scrollable_frame, text=molecule.name)
            label_name.grid(row=row*2, column=col, pady=5)  # Place the name above the image

            # Create label for the image
            image_label = tk.Label(self.scrollable_frame, image=photo)
            image_label.image = photo  # Keep a reference to avoid garbage collection
            image_label.grid(row=row*2 + 1, column=col, pady=5)  # Place image below name

            # Move to the next column (next "cell" in the row)
            col += 1

            # After placing 3 items in a row, move to the next row
            if col >= 3:
                col = 0
                row += 1

        # Update the scroll region to encompass all images and labels
        self.scrollable_frame.update_idletasks()
        self.canvas.config(scrollregion=self.canvas.bbox("all"))

        print(f"Displayed {len(self.molecules)} molecules in grid format.")