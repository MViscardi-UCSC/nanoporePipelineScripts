"""
settingsFileBuilder.py
Marcus Viscardi,    March 27, 2023

Simple goal of this script is to provide a fast and foolproof way
to build the settings files that I use for my nanopore pipeline.

Currently, even with my knowledge of what goes into this, it is an
annoying process of copying and pasting data paths, output paths,
and other details.

This script will provide a GUI based method to build the file(s).
"""
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
from functools import partial


# ChatGPT provided the framework I used to write the below class! Crazy.
class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.checkbox_dict = None
        self.submit_button = None
        self.cancel_button = None
        self.files_and_dirs_dict = None
        self.steps_dict = None
        self.checkbox_header = None
        self.lbl_header = None
        self.master = master
        self.pack(side="top", fill="both", expand=True)
        self.create_widgets()

    def create_widgets(self):
        current_row = 0
        
        # Add header
        self.lbl_header = tk.Label(self, text="Nanopore Pipeline Settings File Builder", font=("Arial", 16, "bold underline"))
        self.lbl_header.grid(row=0, column=0, columnspan=4, pady=10)
        current_row += 1

        # Add checkbox section
        self.checkbox_header = tk.Label(self, text="Select Steps to Run:", font=("Arial", 12, "bold"))
        self.checkbox_header.grid(row=current_row, column=0, columnspan=4)
        current_row += 1

        self.steps_dict = {"G": ["Guppy Basecalling", True],
                           "A": ["Filtering Alt. Genomes (currently pretty rough and slow)", False],
                           "T": ["Trimming TERA-Seq Adapters", False],
                           "M": ["Minimap2 and SamTools", True],
                           "C": ["Concatenate Files", True],
                           "N": ["Nanopolish Index and polyA Calling", True],
                           "F": ["FeatureCounts", True],
                           "P": ["Merging Results w/ Pandas", True],
                           "S": ["Mapping Old Standards (still experimental!!)", False],
                           "L": ["Calling Transcripts w/ Flair", True],
                           "X": ["Running random eXtra steps", False],
                           }
        self.checkbox_dict = {}

        for i, (step, (description, default_state)) in enumerate(self.steps_dict.items()):
            var = tk.BooleanVar(value=default_state)
            checkbox = tk.Checkbutton(self, text=f"{step}: {description}", variable=var)
            self.checkbox_dict[step] = var
            checkbox.grid(row=i + current_row, column=0, columnspan=4, sticky="w", padx=(20, 0), pady=(0, 5))
        current_row += len(self.checkbox_dict)
        
        sep = ttk.Separator(self, orient='horizontal')
        sep.grid(row=current_row, column=0, columnspan=4, sticky="ew", padx=(20, 0), pady=(0, 5))
        current_row += 1

        # Per path selection:
        files_and_dirs_to_select_dict = {'settingsFile': '',
                                         'dataDir': '',
                                         'genomeDir': '/data16/marcus/genomes/plus_cerENO2_elegansRelease100',
                                         'outputDir': ''}
        self.files_and_dirs_dict = {}  # each of the above items will get a label, entry box, and browse button
        for file_or_dir_target, default_path in files_and_dirs_to_select_dict.items():
            lbl = tk.Label(self, text=f"Enter path for {file_or_dir_target}:")
            lbl.grid(row=current_row, column=0, sticky="w", padx=(20, 0), pady=(10, 5))
            ent = tk.Entry(self)
            ent.grid(row=current_row, column=1, columnspan=2, sticky="ew", pady=(10, 5))
            if default_path:
                ent.insert(0, default_path)
            file_or_dir_specific_func = partial(self.browse_file, file_or_dir_target)
            btn = tk.Button(self, text="Browse", command=file_or_dir_specific_func)
            btn.grid(row=current_row, column=3, padx=(0, 20), pady=(10, 5))
            self.files_and_dirs_dict[file_or_dir_target] = [lbl, ent, btn, file_or_dir_specific_func]
            current_row += 1

        self.cancel_button = tk.Button(self, text="Cancel", command=self.cancel)
        self.cancel_button.grid(row=current_row, column=1, pady=(0, 20))
        self.submit_button = tk.Button(self, text="Submit", command=self.submit)
        self.submit_button.grid(row=current_row, column=2, pady=(0, 20))
        current_row += 1

        # Configure the root window to expand
        self.master.columnconfigure(0, weight=1)
        self.master.rowconfigure(0, weight=1)

        # Set weight for rows and columns to make the frame scale with the window
        for i in range(current_row):
            self.grid_rowconfigure(i, weight=1)
        for i in range(4):
            self.grid_columnconfigure(i, weight=1)

    def browse_file(self, file_or_dir_trigger: str):
        if file_or_dir_trigger.endswith("Dir"):
            file_path = filedialog.askdirectory()
        elif file_or_dir_trigger.endswith("File"):
            file_path = filedialog.askopenfilename()
        else:
            file_path = None
        self.files_and_dirs_dict[file_or_dir_trigger][1].delete(0, tk.END)
        self.files_and_dirs_dict[file_or_dir_trigger][1].insert(tk.END, file_path)

    def submit(self):
        # print(f"Steps To Run:")
        stepsToRun = ''
        for stepToRun, stepVar in self.checkbox_dict.items():
            if stepVar.get():
                stepsToRun += stepToRun
        print(f"stepsToRun|{stepsToRun}")
        
        for file_or_dir, (_, ent, _, _) in self.files_and_dirs_dict.items():
            enter_box_content = ent.get()
            if enter_box_content:
                print(f"{file_or_dir}|{enter_box_content}")
        self.master.destroy()

    def cancel(self):
        self.master.destroy()


if __name__ == '__main__':
    root = tk.Tk()
    app = Application(master=root)
    app.mainloop()
