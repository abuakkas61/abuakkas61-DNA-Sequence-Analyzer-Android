import os
from collections import Counter

from kivy.app import App
from kivy.metrics import dp
from kivy.core.window import Window
from kivy.graphics import Color, Line, Rectangle
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.gridlayout import GridLayout
from kivy.uix.scrollview import ScrollView
from kivy.uix.label import Label
from kivy.uix.button import Button
from kivy.uix.textinput import TextInput
from kivy.uix.spinner import Spinner
from kivy.uix.popup import Popup
from kivy.uix.filechooser import FileChooserListView
from kivy.uix.widget import Widget

from sequence_engine import SequenceEngine, ORFEngine, CODON_TABLE

APP_TITLE = "DNA Sequence Analyzer"
BG = (0.965, 0.975, 0.985, 1)
NAVY = (0.055, 0.12, 0.22, 1)
BLUE = (0.08, 0.34, 0.62, 1)
TEXT = (0.10, 0.15, 0.21, 1)
MUTED = (0.35, 0.41, 0.48, 1)


class PlotWidget(Widget):
    """Lightweight Kivy-native plot; avoids desktop TkAgg/Matplotlib backends."""
    def __init__(self, title="", x_values=None, y_series=None, **kwargs):
        super().__init__(**kwargs)
        self.title = title
        self.x_values = x_values or []
        self.y_series = y_series or []
        self.bind(size=self.redraw, pos=self.redraw)
        self.redraw()

    def redraw(self, *args):
        self.canvas.clear()
        with self.canvas:
            Color(1, 1, 1, 1)
            Rectangle(pos=self.pos, size=self.size)
            left, right = self.x + dp(42), self.right - dp(12)
            bottom, top = self.y + dp(28), self.top - dp(30)
            Color(.78, .82, .87, 1)
            Line(points=[left, bottom, right, bottom, right, top], width=1)
            for _, values, rgba in self.y_series:
                if not values:
                    continue
                lo, hi = min(values), max(values)
                if hi == lo:
                    hi = lo + 1
                points = []
                n = len(values)
                for i, val in enumerate(values):
                    px = left + (right-left) * (i / max(1, n-1))
                    py = bottom + (top-bottom) * ((val-lo)/(hi-lo))
                    points.extend([px, py])
                if len(points) >= 4:
                    Color(*rgba)
                    Line(points=points, width=1.6)
        self.canvas.ask_update()


class DNAAnalyzerAndroid(App):
    def build(self):
        self.title = APP_TITLE
        Window.clearcolor = BG
        self.sequence = ""
        self.last_result = ""
        root = BoxLayout(orientation="vertical", spacing=dp(6), padding=dp(8))
        header = BoxLayout(size_hint_y=None, height=dp(58), padding=(dp(10), 0))
        with header.canvas.before:
            Color(*NAVY)
            self.header_rect = Rectangle(pos=header.pos, size=header.size)
        header.bind(pos=self._sync_header, size=self._sync_header)
        title_box = BoxLayout(orientation="vertical", spacing=0)
        title_box.add_widget(Label(text="[b]DNA Sequence Analyzer[/b]", markup=True,
                                   color=(1,1,1,1), font_size="19sp",
                                   halign="left", valign="middle"))
        title_box.add_widget(Label(text="Android edition  |  V12 analysis engine",
                                   color=(.82,.89,.96,1), font_size="11sp",
                                   halign="left", valign="middle"))
        header.add_widget(title_box)
        root.add_widget(header)

        self.mode = Spinner(text="Overview", values=(
            "Overview", "Nucleotide", "Sequence Tools", "Translation",
            "ORFs", "Visualization", "Comparison"),
            size_hint_y=None, height=dp(44), background_color=BLUE,
            color=(1,1,1,1))
        self.mode.bind(text=self.change_mode)
        root.add_widget(self.mode)

        root.add_widget(self.make_label("Input sequence (FASTA or raw DNA/RNA)", 14, True))
        self.input = TextInput(hint_text="Paste your DNA/RNA sequence here...",
                               multiline=True, size_hint_y=None, height=dp(150),
                               font_size="14sp", background_color=(1,1,1,1),
                               foreground_color=TEXT, padding=dp(10))
        root.add_widget(self.input)

        actions = GridLayout(cols=3, spacing=dp(5), size_hint_y=None, height=dp(88))
        for label, callback in [
            ("Analyze", self.analyze), ("Open FASTA", self.open_file),
            ("Save FASTA", self.save_fasta), ("Reverse", lambda *_: self.run_tool("reverse")),
            ("Complement", lambda *_: self.run_tool("complement")),
            ("Rev. complement", lambda *_: self.run_tool("reverse_complement")),
        ]:
            actions.add_widget(self.make_button(label, callback))
        root.add_widget(actions)

        self.content_scroll = ScrollView(do_scroll_x=False)
        self.content = BoxLayout(orientation="vertical", spacing=dp(8),
                                 padding=dp(5), size_hint_y=None)
        self.content.bind(minimum_height=self.content.setter("height"))
        self.content_scroll.add_widget(self.content)
        root.add_widget(self.content_scroll, 1)

        footer = Label(text="IUPAC ambiguity codes supported  •  Local analysis",
                       size_hint_y=None, height=dp(24), color=MUTED, font_size="10sp")
        root.add_widget(footer)
        self.show_overview()
        return root

    def _sync_header(self, widget, value):
        self.header_rect.pos = widget.pos
        self.header_rect.size = widget.size

    def make_label(self, text, size=13, bold=False):
        return Label(text=("[b]"+text+"[/b]") if bold else text,
                     markup=bold, color=TEXT, font_size=f"{size}sp",
                     size_hint_y=None, height=dp(30), halign="left",
                     valign="middle")

    def make_button(self, text, callback):
        b = Button(text=text, size_hint_y=None, height=dp(40),
                   background_normal="", background_color=BLUE,
                   color=(1,1,1,1), font_size="12sp")
        b.bind(on_release=callback)
        return b

    def clear_content(self):
        self.content.clear_widgets()

    def add_text(self, text, height=None, font=13, color=TEXT):
        lines = max(1, text.count("\n")+1)
        h = dp(height if height else min(520, 24 + lines*20))
        label = Label(text=str(text), color=color, font_size=f"{font}sp",
                      size_hint_y=None, height=h, halign="left",
                      valign="top", markup=True)
        label.bind(texture_size=lambda w, s: setattr(w, "text_size", (self.content.width-dp(20), None)))
        label.text_size = (self.content.width-dp(20), None)
        self.content.add_widget(label)
        return label

    def add_card(self, title, body):
        box = BoxLayout(orientation="vertical", padding=dp(10), spacing=dp(4),
                        size_hint_y=None)
        box.bind(minimum_height=box.setter("height"))
        box.add_widget(self.make_label(title, 14, True))
        self.add_text(body)
        # Body is added to content directly; card is intentionally avoided
        # to keep dynamic text sizing reliable on narrow screens.

    def get_valid_sequence(self):
        seq = SequenceEngine.clean(self.input.text)
        ok, message = SequenceEngine.validate(seq)
        if not ok:
            self.notice("Sequence error", message)
            return None
        self.sequence = seq
        return seq

    def notice(self, title, message):
        Popup(title=title, content=Label(text=message, color=TEXT),
              size_hint=(.88, .35)).open()

    def change_mode(self, spinner, value):
        if value == "Overview": self.show_overview()
        elif value == "Nucleotide": self.show_nucleotide()
        elif value == "Sequence Tools": self.show_tools()
        elif value == "Translation": self.show_translation()
        elif value == "ORFs": self.show_orfs()
        elif value == "Visualization": self.show_visualization()
        elif value == "Comparison": self.show_comparison()

    def show_overview(self):
        self.clear_content()
        self.add_text("[b]Welcome[/b]\nPaste a sequence and tap Analyze. The analyzer accepts DNA/RNA IUPAC symbols and FASTA input.", 78)
        self.add_text("[b]Available modules[/b]\n• Nucleotide composition and sequence statistics\n• Reverse, complement, reverse complement and transcription\n• Standard genetic-code translation and six-frame translation\n• ORF detection\n• Sliding-window GC, skew and entropy plots\n• FASTA open/save and text report export", 170)
        if self.sequence:
            self.add_text(f"[b]Loaded sequence:[/b] {len(self.sequence):,} nucleotides")

    def analyze(self, *_):
        seq = self.get_valid_sequence()
        if seq is None: return
        s = SequenceEngine.statistics(seq)
        self.last_result = (
            f"DNA/RNA type: {'RNA' if SequenceEngine.is_rna(seq) else 'DNA'}\n"
            f"Length: {s['length']}\nA: {s['A']}   T: {s['T']}   G: {s['G']}   C: {s['C']}   N: {s['N']}\n"
            f"Ambiguous bases: {s['ambiguous']}\nGC content: {s['GC']:.2f}%\n"
            f"AT content: {s['AT']:.2f}%\nPurines: {s['purines']}   Pyrimidines: {s['pyrimidines']}\n"
            f"GC skew: {s['gc_skew']:.4f}   AT skew: {s['at_skew']:.4f}\n"
            f"Shannon entropy: {s['entropy']:.4f} bits"
        )
        self.clear_content()
        self.add_text("[b]Sequence analysis report[/b]", 32)
        self.add_text(self.last_result, 220)
        self.add_text("Use the module selector for detailed nucleotide, translation, ORF and visualization results.", 55)
        self.mode.text = "Overview"

    def show_nucleotide(self):
        seq = self.get_valid_sequence()
        self.clear_content()
        if seq is None:
            self.add_text("Enter a valid sequence, then open Nucleotide again.")
            return
        s = SequenceEngine.statistics(seq)
        self.add_text("[b]Nucleotide composition[/b]", 32)
        self.add_text("\n".join(f"{k}: {s[k]}" for k in ("A","T","G","C","N","ambiguous","canonical")), 160)
        self.add_text(f"GC: {s['GC']:.2f}%   AT: {s['AT']:.2f}%\nPurines: {s['purines']}   Pyrimidines: {s['pyrimidines']}\nGC skew: {s['gc_skew']:.4f}\nAT skew: {s['at_skew']:.4f}\nEntropy: {s['entropy']:.4f} bits", 150)

    def show_tools(self):
        self.clear_content()
        self.add_text("[b]Sequence transformations[/b]\nChoose a transformation. Results appear below.", 60)
        grid = GridLayout(cols=2, spacing=dp(5), size_hint_y=None, height=dp(180))
        for name, key in [("Reverse","reverse"),("Complement","complement"),
                          ("Reverse complement","reverse_complement"),
                          ("Transcribe DNA → RNA","transcribe"),
                          ("Uppercase","upper"),("Clear input","clear")]:
            grid.add_widget(self.make_button(name, lambda inst, k=key: self.run_tool(k)))
        self.content.add_widget(grid)
        self.result = TextInput(text="", readonly=True, multiline=True,
                                size_hint_y=None, height=dp(180), font_size="13sp")
        self.content.add_widget(self.result)

    def run_tool(self, key):
        seq = self.get_valid_sequence()
        if key == "clear":
            self.input.text = ""
            if hasattr(self, "result"): self.result.text = ""
            return
        if seq is None: return
        funcs = {
            "reverse": SequenceEngine.reverse,
            "complement": SequenceEngine.complement,
            "reverse_complement": SequenceEngine.reverse_complement,
            "transcribe": SequenceEngine.transcribe,
            "upper": str.upper
        }
        out = funcs[key](seq)
        self.last_result = out
        self.mode.text = "Sequence Tools"
        if hasattr(self, "result"):
            self.result.text = out
        else:
            self.clear_content()
            self.add_text(f"[b]{key.replace('_',' ').title()}[/b]", 32)
            self.add_text(out, min(500, 30 + (len(out)//60+1)*20))

    def show_translation(self):
        self.clear_content()
        seq = self.get_valid_sequence()
        if seq is None:
            self.add_text("Enter a valid DNA/RNA sequence, then select Translation.")
            return
        dna = seq.replace("U","T")
        self.add_text("[b]Translation[/b]", 32)
        frames = SequenceEngine.six_frames(dna)
        for frame, protein in frames.items():
            self.add_text(f"[b]Frame {frame}[/b]  ({len(protein)} aa)\n{protein}", min(160, 45+(len(protein)//55+1)*18))
        self.add_text("[b]Codon usage[/b]", 32)
        counts = SequenceEngine.codon_counts(dna)
        self.add_text("\n".join(f"{c}: {n}" for c,n in counts.most_common(30)) or "No complete canonical codons.", 250)
        self.add_text("Select a frame's protein sequence above for interpretation. X denotes an unrecognized codon.", 48)

    def show_orfs(self):
        self.clear_content()
        seq = self.get_valid_sequence()
        if seq is None:
            self.add_text("Enter a valid sequence, then select ORFs.")
            return
        self.add_text("[b]Open reading frame detection[/b]", 32)
        grid = GridLayout(cols=2, size_hint_y=None, height=dp(48), spacing=dp(5))
        self.min_orf = TextInput(text="10", multiline=False, hint_text="Minimum amino acids",
                                 input_filter="int", size_hint_y=None, height=dp(42))
        grid.add_widget(self.min_orf)
        grid.add_widget(self.make_button("Find ORFs", self._find_orfs))
        self.content.add_widget(grid)
        self.orf_results = TextInput(readonly=True, multiline=True, size_hint_y=None,
                                     height=dp(260), font_size="12sp")
        self.content.add_widget(self.orf_results)
        self._find_orfs()

    def _find_orfs(self, *_):
        seq = self.get_valid_sequence()
        if seq is None: return
        try: minimum = max(1, int(self.min_orf.text or "10"))
        except ValueError: minimum = 10
        orfs = ORFEngine.find_orfs(seq, minimum)
        if not orfs:
            text = "No ORFs met the selected minimum length."
        else:
            rows = [f"ORFs found: {len(orfs)}\n"]
            for i, item in enumerate(orfs, 1):
                rows.append(f"{i}. Strand {item['strand']} | Frame {item['frame']} | "
                            f"Start {item['start']} | End {item['end']} | "
                            f"{item['length']} aa\n{item['protein']}\n")
            text = "\n".join(rows)
        self.orf_results.text = text

    def show_visualization(self):
        self.clear_content()
        seq = self.get_valid_sequence()
        if seq is None:
            self.add_text("Enter a valid sequence, then select Visualization.")
            return
        self.add_text("[b]Sequence visualization[/b]\nChoose window and step sizes. Plots are drawn with Kivy graphics.", 65)
        row = GridLayout(cols=2, spacing=dp(5), size_hint_y=None, height=dp(45))
        self.window_input = TextInput(text=str(min(100, max(10, len(seq)//5))), multiline=False,
                                      input_filter="int", hint_text="Window")
        self.step_input = TextInput(text="10", multiline=False, input_filter="int", hint_text="Step")
        row.add_widget(self.window_input); row.add_widget(self.step_input)
        self.content.add_widget(row)
        self.content.add_widget(self.make_button("Generate plots", self._make_plots))
        self.plot_area = BoxLayout(orientation="vertical", spacing=dp(8), size_hint_y=None)
        self.plot_area.bind(minimum_height=self.plot_area.setter("height"))
        self.content.add_widget(self.plot_area)
        self._make_plots()

    def _make_plots(self, *_):
        seq = self.get_valid_sequence()
        if seq is None: return
        try:
            window = max(2, int(self.window_input.text))
            step = max(1, int(self.step_input.text))
        except ValueError:
            self.notice("Input error", "Window and step must be positive integers.")
            return
        if window > len(seq):
            self.notice("Input error", "Window size cannot exceed sequence length.")
            return
        self.plot_area.clear_widgets()
        x, gc = SequenceEngine.sliding_gc(seq, window, step)
        self.add_plot("Sliding-window GC content (%)", x, [("GC%",gc,(.08,.34,.72,1))])
        x, gcs, ats = SequenceEngine.sliding_skew(seq, window, step)
        self.add_plot("Sliding-window skew", x, [("GC skew",gcs,(.08,.34,.72,1)),("AT skew",ats,(.8,.25,.18,1))])
        x, ent = SequenceEngine.sliding_entropy(seq, window, step)
        self.add_plot("Sliding-window Shannon entropy", x, [("Entropy",ent,(.12,.55,.35,1))])

    def add_plot(self, title, x, series):
        self.plot_area.add_widget(self.make_label(title, 13, True))
        chart = PlotWidget(title, x, series, size_hint_y=None, height=dp(210))
        self.plot_area.add_widget(chart)

    def show_comparison(self):
        self.clear_content()
        self.add_text("[b]Positional sequence comparison[/b]\nPaste a second sequence below. This compares symbols at matching positions; it is not an alignment algorithm.", 80)
        self.compare_input = TextInput(hint_text="Second DNA/RNA sequence...", multiline=True,
                                       size_hint_y=None, height=dp(120), font_size="13sp")
        self.content.add_widget(self.compare_input)
        self.content.add_widget(self.make_button("Compare sequences", self.compare_sequences))
        self.compare_result = TextInput(readonly=True, multiline=True,
                                        size_hint_y=None, height=dp(180), font_size="12sp")
        self.content.add_widget(self.compare_result)

    def compare_sequences(self, *_):
        first = self.get_valid_sequence()
        if first is None: return
        second = SequenceEngine.clean(self.compare_input.text)
        ok, msg = SequenceEngine.validate(second)
        if not ok:
            self.notice("Second sequence error", msg); return
        n = min(len(first), len(second))
        matches = sum(a == b for a,b in zip(first[:n], second[:n]))
        identity = matches/n*100 if n else 0
        self.compare_result.text = (f"Length 1: {len(first)}\nLength 2: {len(second)}\n"
                                    f"Compared positions: {n}\nMatches: {matches}\n"
                                    f"Mismatches: {n-matches}\nPositional identity: {identity:.2f}%\n"
                                    f"Length difference: {abs(len(first)-len(second))}")

    def open_file(self, *_):
        self.file_popup("Open FASTA / sequence", self._load_selected, ["*.fasta","*.fa","*.fna","*.txt","*"])

    def save_fasta(self, *_):
        seq = self.get_valid_sequence()
        if seq is None: return
        self.file_popup("Save FASTA", self._save_selected, ["*"])

    def file_popup(self, title, callback, filters):
        chooser = FileChooserListView(path=App.get_running_app().user_data_dir,
                                      filters=filters)
        layout = BoxLayout(orientation="vertical", spacing=dp(5), padding=dp(6))
        layout.add_widget(chooser)
        buttons = BoxLayout(size_hint_y=None, height=dp(44), spacing=dp(5))
        buttons.add_widget(self.make_button("Cancel", lambda *_: popup.dismiss()))
        buttons.add_widget(self.make_button("Select", lambda *_: callback(chooser.selection, popup)))
        layout.add_widget(buttons)
        popup = Popup(title=title, content=layout, size_hint=(.96,.92))
        popup.open()

    def _load_selected(self, selection, popup):
        if not selection:
            self.notice("File", "Select a file first."); return
        try:
            with open(selection[0], "r", encoding="utf-8") as f:
                self.input.text = f.read()
            popup.dismiss()
            self.analyze()
        except Exception as exc:
            self.notice("Open failed", str(exc))

    def _save_selected(self, selection, popup):
        if not selection:
            self.notice("File", "Choose a folder or file location first."); return
        path = selection[0]
        if os.path.isdir(path):
            path = os.path.join(path, "sequence.fasta")
        elif not path.lower().endswith((".fa",".fasta",".fna")):
            path += ".fasta"
        try:
            with open(path, "w", encoding="utf-8") as f:
                f.write(">DNA_Sequence_Analyzer_V12\n")
                seq = self.sequence
                for i in range(0, len(seq), 70): f.write(seq[i:i+70]+"\n")
            popup.dismiss()
            self.notice("Saved", f"FASTA saved to:\n{path}")
        except Exception as exc:
            self.notice("Save failed", str(exc))

    def export_report(self, *_):
        seq = self.get_valid_sequence()
        if seq is None: return
        stats = SequenceEngine.statistics(seq)
        report = "DNA SEQUENCE ANALYZER V12 REPORT\n\n" + "\n".join(
            f"{k}: {v}" for k,v in stats.items())
        path = os.path.join(self.user_data_dir, "DNA_Sequence_Analyzer_Report.txt")
        with open(path, "w", encoding="utf-8") as f: f.write(report)
        self.notice("Report exported", f"Report saved in app storage:\n{path}")


if __name__ == "__main__":
    DNAAnalyzerAndroid().run()
