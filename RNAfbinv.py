#!/usr/bin/env python3
'''
GUI wrapper to the RNAfbinv2.0 package
'''

from tkinter.ttk import Progressbar
from tkinter import filedialog
from tkinter import messagebox
from PIL import ImageTk, Image
import tkinter as tk
import configparser
import threading
import logging
import time
import os

from rnafbinv import shapiro_generator, sfb_designer, RNAfbinvCL, vienna, IUPAC
import varna_generator


# A tk Text widget with an additional custom event for text modification
# Currently not in use
class CustomText(tk.Text):
    def __init__(self, name, *args, **kwargs):
        """A text widget that report on internal widget commands"""
        tk.Text.__init__(self, *args, **kwargs)

        # create a proxy for the underlying widget
        self._orig = self._w + "_orig"
        self.tk.call("rename", self._w, self._orig)
        self.tk.createcommand(self._w, self._proxy)
        self.name = name

    def _proxy(self, command, *args):
        cmd = (self._orig, command) + args
        result = self.tk.call(cmd)

        if command in ("insert", "delete", "replace"):
            self.event_generate("<<TextModified>>")

        return result


# Image generation
BASE_IMG = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'img')
NO_STRUCTURE_IMG = os.path.join(BASE_IMG, 'NoImage.png')
ERROR_IMG = os.path.join(BASE_IMG, 'NoImage.jpg')
IMAGE_SIZE = (350, 400)


def get_no_img():
    img = Image.open(NO_STRUCTURE_IMG)
    resized = img.resize(IMAGE_SIZE, Image.ANTIALIAS)
    img = ImageTk.PhotoImage(resized)
    return img


def generate_image(structure: str, sequence: str, index_list=None):
    img_path = None
    try:
        if not index_list:
            index_list = None
        img_path = varna_generator.generate_image(structure, sequence, index_list)
        img = Image.open(img_path)
    except:
        img = Image.open(ERROR_IMG)
    finally:
        resized = img.resize(IMAGE_SIZE, Image.ANTIALIAS)
        img = ImageTk.PhotoImage(resized)
        if img_path is not None:
            try:
                os.remove(img_path)
            except OSError:
                pass
    return img


# image popup
def popup_image(res_no: int, img):
    win = tk.Toplevel()
    win.config(width=IMAGE_SIZE[0], height=IMAGE_SIZE[1])
    res_text = "Result {}".format(res_no)
    win.wm_title(res_text)
    win.resizable(0, 0)
    # ensure a consistent GUI size
    win.grid_propagate(False)

    img_frame = tk.Frame(win, borderwidth=2, relief=tk.GROOVE)
    img_frame.pack_propagate(0)
    img_frame.grid(row=0, rowspan=4, column=2, padx=10, pady=10)
    image_canvas = tk.Canvas(img_frame, width=IMAGE_SIZE[0], height=IMAGE_SIZE[1], bg='white')
    image_canvas.grid(row=0, column=0)
    image_canvas.create_image(0, 0, image=img, anchor='nw', tag=res_text)


# Frame with results
class RunFrame:
    all_frames = []

    def __init__(self, master, total_runs):
        self.run_index = len(RunFrame.all_frames)
        self.main_grid = tk.Frame(master, bd=1, relief=tk.GROOVE)
        self.master_canvas = master
        self.total_runs = total_runs
        self.label = tk.Label(self.main_grid, text='Run {}:'.format(self.run_index + 1))
        self.label.grid(row=0, column=0, sticky='nw')
        self.progress = Progressbar(self.main_grid, orient="horizontal", length=600, mode="determinate")
        self.progress.grid(row=1, column=0, sticky='nw')
        self.progress['value'] = 0
        self.progress["maximum"] = total_runs
        self.design_result = None
        # self.main_grid.grid(row=run_index, column=0, sticky='news')
        root.update()
        # calculate item box above
        height_above = self.master_canvas.bbox('all')
        logging.debug("creating windows bbox: {} run id: {}".format(height_above, self.run_index))
        if height_above is None:
            height_above = 0
        else:
            height_above = height_above[3]
        logging.debug("creating windows height: {} run id: {}".format(height_above, self.run_index))
        self.item_index = self.master_canvas.create_window(0, height_above, width=self.main_grid.winfo_reqwidth(),
                                                           height=self.main_grid.winfo_reqheight(),
                                                           window=self.main_grid, anchor='nw')
        RunFrame.all_frames.append(self.item_index)
        self.images = [None] * total_runs

    def update(self, current_runs):
        self.progress['value'] = current_runs

    def update_canvas(self):
        height_above = self.master_canvas.bbox(RunFrame.all_frames[self.run_index])[3]
        items_after = RunFrame.all_frames[self.run_index + 1:]
        for item in items_after:
            self.master_canvas.coords(item, [0, height_above])
            height_above = self.master_canvas.bbox(item)[3]
        self.master_canvas.config(scrollregion=self.master_canvas.bbox("all"))
        return

    def update_res(self, design_result):
        # save old height
        item_frame = RunFrame.all_frames[self.run_index]
        old_box = self.master_canvas.bbox(item_frame)
        # generate new data frame
        self.master_canvas.delete(item_frame)
        self.main_grid.destroy()
        self.main_grid = tk.Frame(self.master_canvas, bd=3, relief=tk.GROOVE)
        left_frame = tk.Frame(self.main_grid)
        left_frame.grid(row=0, column=0, sticky='wnes')
        self.text = tk.Text(left_frame, height=6, width=125, background='white smoke', bd=0)
        self.text.insert(tk.END, str(design_result))
        self.text['state'] = tk.DISABLED
        self.text.grid(row=0, column=0, sticky='wnes')
        scrollb = tk.Scrollbar(left_frame, command=self.text.yview)
        scrollb.grid(row=0, column=1, sticky='ns')
        self.text['yscrollcommand'] = scrollb.set
        right_frame = tk.Frame(self.main_grid)
        right_frame.grid(row=0, column=1, sticky='wnes')
        self.design_result = design_result
        right_frame.columnconfigure(1, weight=1)
        generate_image_button = tk.Button(right_frame, text="Generate image",
                                          command=lambda: self.show_design(self.run_index,
                                                                           design_result.sequence,
                                                                           design_result.structure))
        generate_image_button.grid(row=0, column=0, sticky='e')
        # self.main_grid.grid(row=run_index, column=0, sticky='wnes')
        # update window
        root.update()
        RunFrame.all_frames[self.run_index] = self.master_canvas.create_window(0, old_box[1], anchor='nw',
                                                                               width=self.main_grid.winfo_reqwidth(),
                                                                               height=self.main_grid.winfo_reqheight(),
                                                                               window=self.main_grid)

        # fix windows after
        self.update_canvas()

    def update_fail(self):
        # save old height
        item_frame = RunFrame.all_frames[self.run_index]
        old_box = self.master_canvas.bbox(item_frame)
        # generate new data frame
        self.master_canvas.delete(item_frame)
        self.main_grid.destroy()
        self.main_grid = tk.Frame(self.master_canvas, bd=1, relief=tk.GROOVE)
        self.label = tk.Label(self.main_grid, text='Run {} Failed!'.format(self.run_index + 1), justify=tk.LEFT)
        self.label.grid(row=0, column=0, sticky='wnes')
        # self.main_grid.grid(row=self.run_index, column=0, sticky='news')
        # update window
        root.update()
        RunFrame.all_frames[self.run_index] = self.master_canvas.create_window(0, old_box[1], anchor='nw',
                                                                               width=self.main_grid.winfo_reqwidth(),
                                                                               height=self.main_grid.winfo_reqheight(),
                                                                               window=self.main_grid)
        self.update_canvas()

    def show_design(self, res_no: int, sequence: str, structure: str):
        if self.images[res_no] is None:
            img = generate_image(structure, sequence)
            self.images[res_no] = img
        else:
            img = self.images[res_no]
        popup_image(res_no + 1, img)


MOTIF_NAME_MAP = {'H': 'Hairpin-loop', 'S': 'Stem', 'E': 'External',
                  'I': 'Internal-loop', 'B': 'Bulge', 'M': 'Multi-loop'}


# Main GUI class
def update_motif_list(motif_list, shapiro_list):
    shapiro_list = shapiro_list[::-1]
    motif_list.delete(0, tk.END)
    for i in range(0, len(shapiro_list)):
        item = shapiro_list[i]
        motif_str = '{} - {}'.format(MOTIF_NAME_MAP.get(item[0]), item)
        motif_list.insert(0, motif_str)


class RNAfbinvGUI(tk.Frame):
    MAIN_SCREEN_TITLE = "rnafbinv - Fragment based RNA design"
    FRAME_MIN_SIZE = (1150, 700)
    MODE_QUERY = 0
    MODE_RESULTS = 1

    def __init__(self, master=None):
        self.mode = None
        self.run_thread = None
        self.result_frame = None
        self.query_frame = None
        super().__init__(master)
        self.info_componenets = {}
        master.title(RNAfbinvGUI.MAIN_SCREEN_TITLE)
        master.minsize(*self.FRAME_MIN_SIZE)
        self.pack(fill='both', expand=True)
        # ensure a consistent GUI size
        self.grid_propagate(False)
        # implement stretchability
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        # create content
        self.create_query_widgets()
        self.keep_running = False
        self.arguments = None

    def create_query_widgets(self):
        def create_label_tb(parent, label_text):
            iframe = tk.Frame(parent)
            label = tk.Label(iframe, text=label_text, justify=tk.LEFT)
            label.grid(row=0, column=0, sticky='w', padx=2, pady=2)

            tb = tk.Text(iframe, height=5, width=80, font=('consolas', 12))
            tb.bind("<FocusOut>", self.on_change_custom_text)
            tb.grid(row=1, column=0, sticky='nsew', padx=2, pady=2)

            scrollb = tk.Scrollbar(iframe, command=tb.yview)
            scrollb.grid(row=1, column=1, sticky='nsew')
            tb['yscrollcommand'] = scrollb.set
            return iframe, tb

        def entry_check_combo(parent, name, text, validation, change_func, default_value=None):
            iframe = tk.Frame(parent)
            selval = tk.IntVar(0)
            self.info_componenets['is_{}'.format(name)] = selval
            cb = tk.Checkbutton(iframe, text=text, command=change_func, variable=selval)
            cb.grid(row=0, column=0, sticky='nw', padx=2, pady=5)
            vcmd = iframe.register(validation)
            var = tk.StringVar()
            entry = tk.Entry(iframe, validate='key', validatecommand=(vcmd, '%P', '%S'), width=6,
                             state='disabled', textvariable=var)
            if default_value is not None:
                var.set(default_value)
            self.info_componenets[name] = entry
            entry.grid(row=0, column=1, sticky='nw', padx=2, pady=5)
            return iframe

        self.mode = self.MODE_QUERY
        # check for prior installation
        if self.result_frame is not None:
            self.result_frame.grid_forget()
        if self.query_frame is not None:
            self.query_frame.grid(row=0, column=0, sticky='wens')
            return

        self.query_frame = tk.Frame(self)
        # ensure a consistent GUI size
        self.query_frame.grid_propagate(False)
        # implement stretchability
        self.query_frame.grid_rowconfigure(0, weight=1)
        self.query_frame.grid_columnconfigure(0, weight=1)
        self.query_frame.grid(row=0, column=0, sticky='wens')
        # This frame will hold parameter components
        top_frame = tk.Frame(self.query_frame)
        top_frame.grid_rowconfigure(0, weight=0)
        top_frame.grid_columnconfigure(0, weight=1)
        top_frame.grid(row=0, column=0, sticky='n', padx=10)

        # Left side, Parameters
        structure_frame, structure_tb = create_label_tb(top_frame, 'Target Structure:')
        self.info_componenets['target_structure'] = structure_tb
        structure_frame.grid(row=0, column=0, columnspan=2, sticky='nesw')
        sequence_frame, sequence_tb = create_label_tb(top_frame, 'Target sequence:')
        self.info_componenets['target_sequence'] = sequence_tb
        sequence_frame.grid(row=1, column=0, columnspan=2, sticky='nesw')

        output_no_frame = tk.Frame(top_frame)
        output_no_frame.grid(row=2, column=0, sticky='wn')
        output_no_label = tk.Label(output_no_frame, text='No. output sequences (1-1000): ', justify=tk.LEFT)
        output_no_label.grid(row=0, column=0, sticky='w')
        vcmd = output_no_frame.register(self.validate_out_no)
        output_no_val = tk.IntVar(output_no_frame)
        output_no_spinbox = tk.Spinbox(output_no_frame, from_=1, to=1000, validate='key',
                                       validatecommand=(vcmd, '%P'), width=5, textvariable=output_no_val)
        output_no_val.set(10)
        self.info_componenets['num_output'] = output_no_val
        output_no_spinbox.grid(row=0, column=1, sticky='wn')

        show_advanced_var = tk.IntVar(0)
        self.info_componenets['show_advanced'] = show_advanced_var
        show_advanced_check = tk.Checkbutton(top_frame, text='Show advanced options', command=self.show_advanced_change,
                                             variable=show_advanced_var)
        show_advanced_check.grid(row=3, column=0, sticky='wn')

        # advanced options frame
        advanced_frame = tk.Frame(top_frame, borderwidth=3, relief=tk.SUNKEN)
        self.info_componenets['advanced_frame'] = advanced_frame

        frame = entry_check_combo(advanced_frame, 'target_energy', 'Target energy:', self.validate_energy,
                                  self.target_energy_change)
        frame.grid(row=0, column=0, sticky='nsew')
        frame = entry_check_combo(advanced_frame, 'mr', 'Target mutational robustness:', self.validate_mr,
                                  self.mr_change)
        frame.grid(row=0, column=1, sticky='nsew')
        frame = entry_check_combo(advanced_frame, 'iter_no', 'No. of iterations:',
                                  self.validate_iteration, self.iteration_change, default_value=100)
        frame.grid(row=1, column=0, sticky='nsew')
        frame = entry_check_combo(advanced_frame, 'variable_length', 'Maximum length change:',
                                  self.validate_lc, self.lc_change, default_value=0)
        frame.grid(row=1, column=1, sticky='nsew')
        frame = entry_check_combo(advanced_frame, 'max_bi', 'Bulge and interior loop no penalty(max size):',
                                  self.validate_maxbi, self.maxbi_change, default_value=0)
        frame.grid(row=2, column=0, sticky='nsew')
        seq_motif_var = tk.IntVar(0)
        self.info_componenets['seq_motif'] = seq_motif_var
        show_advanced_check = tk.Checkbutton(advanced_frame, text='Allow sequence motifs (lower case)',
                                             variable=seq_motif_var)
        show_advanced_check.grid(row=2, column=1, sticky='nw')

        starting_seq_frame = tk.Frame(advanced_frame)
        starting_seq_var = tk.IntVar(0)
        self.info_componenets['is_start_seq'] = starting_seq_var
        starting_seq_label = tk.Checkbutton(starting_seq_frame, text='Starting sequence:',
                                            command=self.starting_seq_change, variable=starting_seq_var)
        starting_seq_label.grid(row=0, column=0, sticky='w', padx=0, pady=2)

        starting_seq_tb = tk.Text(starting_seq_frame, height=4, width=77, font=('consolas', 12),
                                  background='white smoke')
        starting_seq_tb.bind("<1>", lambda event: starting_seq_tb.focus_set())
        starting_seq_tb.configure(state='disabled')
        self.info_componenets['start_seq'] = starting_seq_tb
        starting_seq_tb.grid(row=1, column=0, sticky='nsew', padx=2, pady=2)
        starting_seq_frame.grid(row=3, column=0, columnspan=2, sticky='nsew')

        scrollb = tk.Scrollbar(starting_seq_frame, command=starting_seq_tb.yview)
        scrollb.grid(row=1, column=1, sticky='nsew')
        starting_seq_tb['yscrollcommand'] = scrollb.set

        # Right side, Image and motif selection
        img_frame = tk.Frame(top_frame, borderwidth=2, relief=tk.GROOVE)
        img_frame.pack_propagate(0)
        img_frame.grid(row=0, rowspan=4, column=2, padx=10, pady=10)
        no_img = get_no_img()
        self.info_componenets['image'] = no_img
        image_canvas = tk.Canvas(img_frame, width=IMAGE_SIZE[0], height=IMAGE_SIZE[1], bg='white')
        self.info_componenets['image_canvas'] = image_canvas
        image_canvas.grid(row=0, column=0)
        image_canvas.create_image(0, 0, image=no_img, anchor='nw', tag='FRONT_IMG')

        motif_frame = tk.Frame(top_frame)
        motif_frame.grid(row=4, column=2, pady=5, padx=10)
        motif_label = tk.Label(motif_frame, text='Motif selection (optional):')
        motif_label.grid(row=0, column=0, sticky='w')
        motif_list = tk.Listbox(motif_frame, selectmode=tk.MULTIPLE)
        motif_list.grid(row=1, column=0, sticky='wnes')
        self.info_componenets['motif_list'] = motif_list
        motif_list.bind('<<ListboxSelect>>', self.motif_selected)
        motif_scrollb = tk.Scrollbar(motif_frame, command=motif_list.yview)
        motif_scrollb.grid(row=1, column=1, sticky='nsew')
        motif_list['yscrollcommand'] = scrollb.set

        # This frame will hold the send / run buttons
        bottom_frame = tk.Frame(self.query_frame)
        bottom_frame.grid_rowconfigure(0, weight=1)
        bottom_frame.grid_rowconfigure(0, weight=1)
        bottom_frame.grid(row=1, column=0, columnspan=2, sticky='n', padx=10)

        run_button = tk.Button(bottom_frame, text="Run", command=self.run_design)
        run_button.grid(row=0, column=0, pady=10)

    def on_change_custom_text(self, event):
        if self.mode != self.MODE_QUERY:
            return
        logging.debug(
            'change event: name: "{}" value: "{}"'.format(event.widget, event.widget.get('1.0', tk.END).strip()))
        fold_image = None
        real_stuct = RNAfbinvCL.bracket_changer(self.info_componenets['target_structure'].get('1.0', tk.END).strip())
        target_seq = self.info_componenets['target_sequence'].get('1.0', tk.END).strip()
        motif_list = self.info_componenets.get('motif_list')
        if real_stuct != '' and RNAfbinvCL.is_valid_structure(real_stuct):
            if len(target_seq) != len(real_stuct):
                target_seq = 'N' * len(real_stuct)
            image_data = self.info_componenets.get('image_data')
            shapiro_list = shapiro_generator.get_shapiro(real_stuct).shapiro.replace('(', '').split(')')[:-2][::-1]
            if image_data is None or image_data != (real_stuct, target_seq, shapiro_list):
                self.info_componenets['image_data'] = (real_stuct, target_seq, shapiro_list)
                fold_image = generate_image(real_stuct, target_seq)
                self.info_componenets['image'] = fold_image
                # generate shapiro list from shapiro without 'R'
                update_motif_list(motif_list, shapiro_list)
        else:
            update_motif_list(motif_list, [])
            fold_image = get_no_img()
            self.info_componenets['image'] = fold_image
        if fold_image is not None:
            self.info_componenets.get('image_canvas').delete("FRONT_IMG")
            self.info_componenets.get('image_canvas').create_image(0, 0, image=fold_image, anchor='nw', tag='FRONT_IMG')

    def show_advanced_change(self):
        value = self.info_componenets.get('show_advanced').get()
        advanced_frame = self.info_componenets.get('advanced_frame')
        if value == 0:
            advanced_frame.grid_forget()
        else:
            advanced_frame.grid(row=4, column=0, sticky='wn')

    def starting_seq_change(self):
        value = self.info_componenets.get('is_start_seq').get()
        starting_seq_text = self.info_componenets.get('start_seq')
        if value == 0:
            starting_seq_text.configure(state='disabled', background='white smoke')
        else:
            starting_seq_text.configure(state='normal', background='white')

    def target_energy_change(self):
        value = self.info_componenets.get('is_target_energy').get()
        target_energy_entry = self.info_componenets.get('target_energy')
        if value == 0:
            target_energy_entry.configure(state='disabled')
        else:
            target_energy_entry.configure(state='normal')

    def validate_energy(self, value, added):
        try:
            if value != '-' or value != '':
                float(value)
            return True
        except ValueError:
            return False

    def validate_maxbi(self, value, added):
        try:
            ival = int(value)
            if ival < 0:
                return False
            return True
        except ValueError:
            return False

    def maxbi_change(self):
        value = self.info_componenets.get('is_max_bi').get()
        maxbi_entry = self.info_componenets.get('max_bi')
        if value == 0:
            maxbi_entry.configure(state='disabled')
        else:
            maxbi_entry.configure(state='normal')

    def mr_change(self):
        value = self.info_componenets.get('is_mr').get()
        mr_entry = self.info_componenets.get('mr')
        if value == 0:
            mr_entry.configure(state='disabled')
        else:
            mr_entry.configure(state='normal')

    def validate_mr(self, value, added):
        if added == '-':
            return False
        try:
            if value != '':
                ival = float(value)
                if ival < 0 or ival > 1:
                    return False
            return True
        except ValueError:
            return False

    def iteration_change(self):
        value = self.info_componenets.get('is_iter_no').get()
        iteration_entry = self.info_componenets.get('iter_no')
        if value == 0:
            iteration_entry.configure(state='disabled')
        else:
            iteration_entry.configure(state='normal')

    def validate_iteration(self, value, added):
        if added == '-':
            return False
        try:
            if value != '':
                ival = int(value)
                if ival < 0:
                    return False
            return True
        except ValueError:
            return False

    def lc_change(self):
        value = self.info_componenets.get('is_variable_length').get()
        lc_entry = self.info_componenets.get('variable_length')
        if value == 0:
            lc_entry.configure(state='disabled')
        else:
            lc_entry.configure(state='normal')

    def validate_lc(self, value, added):
        if added == '-':
            return False
        try:
            if value != '':
                ival = int(value)
                if ival < 0:
                    return False
            return True
        except ValueError:
            return False

    def validate_out_no(self, value):
        try:
            if value != '':
                ival = int(value)
                if ival < 1 or ival > 1000:
                    return False
            return True
        except ValueError:
            return False

    def create_run_screen(self, arguments):
        if self.query_frame is not None:
            self.query_frame.grid_forget()

        # check for prior installation
        if self.result_frame is not None:
            self.result_frame.destroy()
        self.result_frame = tk.Frame(self)
        # ensure a consistent GUI size
        self.result_frame.grid_propagate(False)
        # implement stretchability
        self.result_frame.grid_rowconfigure(0, weight=1)
        self.result_frame.grid_columnconfigure(0, weight=1)
        self.result_frame.grid(row=0, column=0, sticky='wnes')

        top_frame = tk.Frame(self.result_frame, bd=1, relief=tk.SUNKEN, width=1050, height=600)
        top_frame.rowconfigure(0, weight=1)
        top_frame.columnconfigure(0, weight=1)
        top_frame.grid(row=0, column=0, sticky='wens', padx=2, pady=5)
        top_canvas = tk.Canvas(top_frame, width=1050, height=600)
        top_canvas.grid(row=0, column=0, sticky='wens')
        RunFrame.all_frames = []
        progression_list = [RunFrame(top_canvas, arguments.get('iter')) for i in
                            range(0, arguments.get('output_num'))]
        vbar = tk.Scrollbar(top_frame, orient=tk.VERTICAL, command=top_canvas.yview)
        vbar.grid(row=0, column=1, sticky='ns')
        top_canvas.config(width=1050, height=600)
        top_canvas.config(yscrollcommand=vbar.set)
        top_canvas.config(scrollregion=top_canvas.bbox("all"))
        self.canvas = top_canvas
        self.canvas.bind("<Configure>", self.on_canvas_configure)
        bottom_frame = tk.Frame(self.result_frame)
        bottom_frame.grid(row=1, column=0, padx=2, pady=5)
        back_button = tk.Button(bottom_frame, text="Done", command=self.done_results)
        back_button.grid(row=0, column=0, padx=2, pady=5)
        back_button = tk.Button(bottom_frame, text="Export", command=self.export_results, state=tk.DISABLED)
        self.info_componenets['export_button'] = back_button
        back_button.grid(row=0, column=1, padx=2, pady=5)
        return progression_list

    def on_canvas_configure(self, event):
        self.canvas.config(scrollregion=self.canvas.bbox("all"))

    def done_results(self):
        answer = messagebox.askokcancel('Return to Query screen',
                                        'Returning to Query screen will result in removal of current query '
                                        '(cancelled if still running).\nIf the design process is still running '
                                        'this may take a few seconds')
        if answer is not None and answer:
            try:
                self.keep_running = False
                sfb_designer.stop(self.arguments)
                # sleep 100 milliseconds to let main algorithm reach it's last update call
                time.sleep(0.1)
                if self.run_thread is not None:
                    # Force UI update to release other thread so it can end
                    # Deadlock might occur if this runs after a call to update progress bar
                    root.update()
                    self.run_thread.join()
            except RuntimeError:
                pass
            self.create_query_widgets()

    def run_all(self, arguments, progression_list):
        rna_folder = vienna.LiveRNAfold(logging)
        rna_folder.start(False)
        arguments['RNAfold'] = rna_folder
        arguments['logger'] = logging
        self.info_componenets['result_list'] = []
        for item in progression_list:
            if not self.keep_running:
                break
            # run simulated annealing
            arguments['updater'] = item
            self.arguments = arguments
            logging.debug("Starting simulated_annealing\nArguments: {}".format(arguments))
            designed_sequence = sfb_designer.simulated_annealing(arguments)
            logging.debug("Finished simulated_annealing\nSequence: {}".format(designed_sequence))
            if self.keep_running:
                if designed_sequence is not None:
                    logging.info("Finished simulated annealing, resulting sequence: {}".format(designed_sequence))
                    result_object = sfb_designer.generate_res_object(designed_sequence, arguments)
                    item.update_res(result_object)
                    self.info_componenets['result_list'].append(result_object)
                else:
                    item.update_fail()
        self.info_componenets['export_button']['state'] = tk.NORMAL
        rna_folder.close()

    def export_results(self):
        files = [('All Files', '*.*'),
                 ('Text Document', '*.txt')]
        with filedialog.asksaveasfile(filetypes=files, defaultextension=files) as out_file:
            i = 0
            out_file.write("\tscore\tbase pair distance\tshapiro distance\tfree energy\tmutational robustness\t"
                           "aligned tree\n")
            for result in self.info_componenets['result_list']:
                i += 1
                out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(i, result.score, result.bp_dist,
                                                                     result.tree_edit_distance, result.energy,
                                                                     result.mutational_robustness, result.align_tree))

    def run_design(self):
        arguments = self.test_components()
        if arguments.get('error') is not None:
            messagebox.showerror("Error", arguments.get('error'))
            return
        else:
            self.mode = self.MODE_RESULTS
            arguments['stop'] = None
            progression_list = self.create_run_screen(arguments)
            self.run_thread = threading.Thread(target=self.run_all, args=(arguments, progression_list))
            self.keep_running = True
            self.run_thread.start()

    def test_components(self):
        arguments = {}
        # Target structure
        real_stuct = RNAfbinvCL.bracket_changer(self.info_componenets['target_structure'].get('1.0', tk.END).strip())
        if RNAfbinvCL.is_valid_structure(real_stuct) and len(real_stuct) > 0:
            arguments['target_structure'] = real_stuct
        else:
            arguments['error'] = 'Target structure must be a balanced dot bracket structure'
            return arguments
        # Target sequence
        target_seq = self.info_componenets['target_sequence'].get('1.0', tk.END).strip()
        if len(target_seq) > 0 and len(target_seq) == len(real_stuct) and IUPAC.is_valid_sequence(target_seq):
            arguments['target_sequence'] = target_seq if self.info_componenets['seq_motif'].get() != 0 else \
                target_seq.upper()
        elif target_seq == '':
            arguments['target_sequence'] = 'N' * len(real_stuct)
        else:
            arguments['error'] = 'Target sequence must be legal FASTA IUPAC sequence and be the same length as' \
                                 ' target structure'
            return arguments
        # no of outputs
        num_of_output = self.info_componenets.get('num_output').get()
        if num_of_output is not None and num_of_output <= 0:
            arguments['error'] = arguments['error'] = 'Number of output sequences should be a positive number'
            return arguments
        arguments['output_num'] = num_of_output
        # Target energy
        if self.info_componenets['is_target_energy'].get() != 0:
            try:
                target_energy = float(self.info_componenets['target_energy'].get())
                arguments['target_energy'] = target_energy
            except ValueError:
                arguments['error'] = 'Failed to read target energy (floating number)'
                return arguments
        # Target neutrality
        if self.info_componenets['is_mr'].get() != 0:
            try:
                mr = float(self.info_componenets['mr'].get())
                arguments['target_neutrality'] = mr
            except ValueError:
                arguments['error'] = 'Failed to read mutational robustness (0-1)'
                return arguments
        # number of iterations
        if self.info_componenets['is_iter_no'].get() != 0:
            try:
                no_iterations = int(self.info_componenets['iter_no'].get())
                arguments['iter'] = no_iterations
            except ValueError:
                arguments['error'] = 'Failed to read no of iterations (integer)'
                return arguments
        else:
            arguments['iter'] = RNAfbinvCL.DEF_NO_ITER
        # allowed variable length
        if self.info_componenets['is_variable_length'].get() != 0:
            try:
                no_iterations = int(self.info_componenets['variable_length'].get())
                arguments['vlength'] = no_iterations
            except ValueError:
                arguments['error'] = 'Failed to read maximum length change (integer)'
                return arguments
        else:
            arguments['vlength'] = 0
        # Starting sequence
        if self.info_componenets['is_start_seq'].get() != 0:
            starting_sequence = self.info_componenets['start_seq'].get(1.0, tk.END).strip()
            if IUPAC.is_valid_sequence(starting_sequence):
                arguments['starting_sequence'] = starting_sequence
            else:
                arguments['error'] = 'Starting sequence must be legal FASTA IUPAC sequence and be the same length as' \
                                     ' target structure'
                return arguments
        # Reduced penalty BI
        if self.info_componenets['is_max_bi'].get() != 0:
            arguments['reduced_bi'] = int(self.info_componenets['max_bi'].get())
        else:
            arguments['reduced_bi'] = 0
        # Motifs
        motif_list = []
        motif_widget = self.info_componenets['motif_list']
        selected = motif_widget.curselection()
        if selected:
            for i in selected:
                motif_str = motif_widget.get(i).rsplit(' ', 1)[1].strip()
                motif_list.append({'index': 1 + i, 'name': motif_str[0],
                                   'length': int(motif_str[1:])})
        arguments['motifs'] = motif_list
        # fold - currently not an option
        arguments['fold'] = 'MFE'
        # look_ahead - currently not an option
        arguments['look_ahead'] = 4
        return arguments

    def motif_selected(self, event):
        def str_list_to_int_list(lst):
            return [int(lst_item.strip()) for lst_item in lst]

        # gather seq + struct
        image_data = self.info_componenets.get('image_data')
        if image_data is not None and event.widget.size() > 0:
            real_stuct = image_data[0]
            target_seq = image_data[1]
            indexes = []
            # gather motif index list
            shapiro_indexes = shapiro_generator.get_shapiro(real_stuct).shapiro_indexes. \
                                  replace('(', '').split(')')[:-2][::-1]
            shapiro_indexes = [str_list_to_int_list(item[1:-1].split(',')) for item in shapiro_indexes]
            widget = event.widget
            selected = widget.curselection()
            for i in selected:
                indexes += shapiro_indexes[i]
            # generate image and replace
            img = generate_image(real_stuct, target_seq, [item + 1 for item in indexes])
            self.info_componenets['image'] = img
            self.info_componenets.get('image_canvas').delete("FRONT_IMG")
            self.info_componenets.get('image_canvas').create_image(0, 0, image=img, anchor='nw', tag='FRONT_IMG')


def read_config():
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.ini')
    config = configparser.ConfigParser()
    config.read(config_path)
    path_section = config['PATH']
    vienna.set_vienna_path(path_section.get('VIENNA', ''))
    java = path_section.get('JAVA')
    if java is not None and java != '':
        varna_generator.set_java_path(java)
    varna = path_section.get('VARNA')
    if varna is not None and varna != '':
        varna_generator.set_varna_path(varna)


if __name__ == '__main__':
    # logging.basicConfig(level=logging.DEBUG)
    read_config()
    root = tk.Tk()
    root.grid_propagate(False)
    my_gui = RNAfbinvGUI(root)
    root.mainloop()
