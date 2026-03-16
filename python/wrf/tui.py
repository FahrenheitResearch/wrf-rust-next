"""
wrf-rust Terminal UI.

Browse WRF files, select variables, then explicitly compute/export/plot.
Nothing computes until you press a button.

Launch:
    python -m wrf tui [directory_or_file]
"""

from __future__ import annotations

import glob
import os
import sys
import time

import numpy as np

from textual import on, work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.widgets import (
    Footer,
    Header,
    Input,
    Label,
    OptionList,
    Static,
    Button,
    ProgressBar,
)
from textual.widgets.option_list import Option

from rich.panel import Panel
from rich.table import Table
from rich import box


# ── Helpers ──────────────────────────────────────────────────────────────────

def _find_wrf_files(path: str) -> list[str]:
    if os.path.isfile(path):
        return [path]
    files = []
    for pattern in ("wrfout*", "wrfout_*", "*.nc", "*.nc4"):
        files.extend(glob.glob(os.path.join(path, pattern)))
    return sorted(set(files))


def _load_wrf(path: str):
    from wrf import WrfFile
    return WrfFile(path)


def _get_var_list() -> list[dict]:
    from wrf import list_variables
    return list_variables()


def _parse_timesteps(text: str, nt: int) -> list[int]:
    """Parse timestep input like '0', '0-5', '0,2,4', 'all'."""
    text = text.strip().lower()
    if text in ("all", "*", ""):
        return list(range(nt))
    if "-" in text and "," not in text:
        parts = text.split("-")
        start = int(parts[0])
        end = int(parts[1])
        return list(range(start, min(end + 1, nt)))
    if "," in text:
        return [int(x.strip()) for x in text.split(",") if x.strip().isdigit()]
    return [int(text)]


# ── App ──────────────────────────────────────────────────────────────────────

class WrfTui(App):
    """WRF file browser and variable selector."""

    CSS = """
    Screen {
        layout: grid;
        grid-size: 3 1;
        grid-columns: 1fr 1fr 1fr;
        grid-gutter: 1;
    }

    #files-panel { height: 100%; padding: 0 1; }
    #vars-panel  { height: 100%; padding: 0 1; }
    #action-panel { height: 100%; padding: 0 1; }

    #file-list { height: 1fr; }
    #var-list  { height: 1fr; }
    #file-info { height: auto; margin-bottom: 1; }
    #var-detail { height: auto; margin-bottom: 1; }

    #selected-list { height: 1fr; margin-bottom: 1; }

    .panel-title {
        text-style: bold;
        margin-bottom: 1;
        color: $accent;
    }

    #progress-bar {
        height: auto;
        margin: 1 0;
    }

    #progress-label {
        height: auto;
        color: $text-muted;
    }

    #output-log {
        height: auto;
        max-height: 12;
        overflow-y: auto;
        margin-top: 1;
    }

    .action-btn {
        margin-bottom: 1;
        width: 100%;
    }
    """

    BINDINGS = [
        Binding("q", "quit", "Quit"),
        Binding("a", "select_all", "Select all"),
        Binding("c", "clear_selected", "Clear"),
    ]

    TITLE = "wrf-rust"

    def __init__(self, start_path: str | None = None):
        super().__init__()
        self.start_path = start_path or os.getcwd()
        self.wf = None
        self.wf_path: str | None = None
        self.all_vars: list[dict] = []
        self.selected_vars: list[str] = []

    def compose(self) -> ComposeResult:
        yield Header()

        # Left: file browser
        with Vertical(id="files-panel"):
            yield Label("[bold]Files[/bold]", classes="panel-title")
            yield Input(placeholder="Filter files...", id="file-filter")
            yield OptionList(id="file-list")
            yield Static(id="file-info")

        # Center: variable picker
        with Vertical(id="vars-panel"):
            yield Label("[bold]Variables[/bold]  [dim]Enter = toggle select[/dim]", classes="panel-title")
            yield Input(placeholder="Filter variables...", id="var-filter")
            yield OptionList(id="var-list")
            yield Static(id="var-detail")

        # Right: selected + actions
        with Vertical(id="action-panel"):
            yield Label("[bold]Selected[/bold]", classes="panel-title")
            yield OptionList(id="selected-list")
            yield Label("Timesteps [dim](e.g. 0  or  0-5  or  all)[/dim]:", classes="panel-title")
            yield Input(value="0", id="timestep-input")
            yield Button("Export to .npy", id="btn-export", variant="primary", classes="action-btn")
            yield Button("Plot to .png", id="btn-plot", variant="default", classes="action-btn")
            yield Button("Plot to .gif", id="btn-gif", variant="default", classes="action-btn")
            yield Button("Compute stats", id="btn-stats", variant="default", classes="action-btn")
            yield Label("", id="progress-label")
            yield ProgressBar(id="progress-bar", total=100, show_eta=False)
            yield Static(id="output-log")

        yield Footer()

    def on_mount(self) -> None:
        self.all_vars = _get_var_list()
        self._populate_var_list(self.all_vars)
        self._scan_files(self.start_path)
        self._refresh_selected_list()
        # Hide progress bar initially
        self.query_one("#progress-bar", ProgressBar).update(total=100, progress=0)

    # ── File list ──

    def _scan_files(self, path: str) -> None:
        file_list = self.query_one("#file-list", OptionList)
        file_list.clear_options()
        files = _find_wrf_files(path)
        if not files:
            file_list.add_option(Option("[dim]No WRF files found[/dim]", id="__none__"))
            return
        for f in files:
            file_list.add_option(Option(os.path.basename(f), id=f))

    @on(Input.Changed, "#file-filter")
    def _on_file_filter(self, event: Input.Changed) -> None:
        q = event.value.lower().strip()
        file_list = self.query_one("#file-list", OptionList)
        file_list.clear_options()
        for f in _find_wrf_files(self.start_path):
            name = os.path.basename(f)
            if q and q not in name.lower():
                continue
            file_list.add_option(Option(name, id=f))

    @on(OptionList.OptionSelected, "#file-list")
    def _on_file_select(self, event: OptionList.OptionSelected) -> None:
        if not event.option or event.option.id == "__none__":
            return
        self._load_file(str(event.option.id))

    @work(thread=True)
    def _load_file(self, path: str) -> None:
        self.call_from_thread(self._set_progress_label, f"Loading {os.path.basename(path)}...")
        try:
            wf = _load_wrf(path)
            self.wf = wf
            self.wf_path = path

            try:
                times = wf.times()
            except Exception:
                times = []

            tbl = Table(box=box.SIMPLE, show_header=False, padding=(0, 1))
            tbl.add_column("", style="bold cyan", width=8)
            tbl.add_column("")
            tbl.add_row("Grid", f"{wf.nx} x {wf.ny} x {wf.nz}")
            tbl.add_row("Times", str(wf.nt))
            tbl.add_row("dx", f"{wf.dx:g} m")
            if times:
                tbl.add_row("Start", times[0])
                if len(times) > 1:
                    tbl.add_row("End", times[-1])

            self.call_from_thread(
                self.query_one("#file-info", Static).update,
                Panel(tbl, title=os.path.basename(path), border_style="green"),
            )
            self.call_from_thread(self._set_subtitle, os.path.basename(path))
            self.call_from_thread(self._set_progress_label, f"Loaded {os.path.basename(path)}")

        except Exception as e:
            self.call_from_thread(
                self.query_one("#file-info", Static).update,
                Panel(f"[red]{e}[/red]", title="Error", border_style="red"),
            )
            self.call_from_thread(self._set_progress_label, "")

    def _set_subtitle(self, text: str) -> None:
        self.sub_title = text

    def _set_progress_label(self, text: str) -> None:
        self.query_one("#progress-label", Label).update(text)

    # ── Variable list (just browse + toggle, NO computation) ──

    def _populate_var_list(self, vars_list: list[dict]) -> None:
        var_list = self.query_one("#var-list", OptionList)
        var_list.clear_options()
        for v in vars_list:
            marker = "[green]\u2713[/green] " if v["name"] in self.selected_vars else "  "
            label = f"{marker}[bold]{v['name']}[/bold]  [dim]{v['units']}[/dim]"
            var_list.add_option(Option(label, id=v["name"]))

    @on(Input.Changed, "#var-filter")
    def _on_var_filter(self, event: Input.Changed) -> None:
        q = event.value.lower().strip()
        if q:
            filtered = [
                v for v in self.all_vars
                if q in v["name"].lower() or q in v["description"].lower()
                or q in v["units"].lower()
            ]
        else:
            filtered = list(self.all_vars)
        self._populate_var_list(filtered)

    @on(OptionList.OptionHighlighted, "#var-list")
    def _on_var_highlight(self, event: OptionList.OptionHighlighted) -> None:
        """Just show metadata. Zero computation."""
        if not event.option or not event.option.id:
            return
        varname = str(event.option.id)
        info = next((v for v in self.all_vars if v["name"] == varname), None)
        if not info:
            return
        detail = self.query_one("#var-detail", Static)
        detail.update(Panel(
            f"[bold]{info['name']}[/bold]\n{info['description']}\nUnits: {info['units']}",
            border_style="blue",
        ))

    @on(OptionList.OptionSelected, "#var-list")
    def _on_var_toggle(self, event: OptionList.OptionSelected) -> None:
        """Toggle checkmark. Zero computation."""
        if not event.option or not event.option.id:
            return
        varname = str(event.option.id)
        if varname in self.selected_vars:
            self.selected_vars.remove(varname)
        else:
            self.selected_vars.append(varname)
        self._refresh_var_marks()
        self._refresh_selected_list()

    @on(OptionList.OptionSelected, "#selected-list")
    def _on_deselect(self, event: OptionList.OptionSelected) -> None:
        """Remove from selected list on click."""
        if not event.option or event.option.id == "__none__":
            return
        varname = str(event.option.id)
        if varname in self.selected_vars:
            self.selected_vars.remove(varname)
            self._refresh_var_marks()
            self._refresh_selected_list()

    def _refresh_var_marks(self) -> None:
        q = self.query_one("#var-filter", Input).value.lower().strip()
        if q:
            filtered = [
                v for v in self.all_vars
                if q in v["name"].lower() or q in v["description"].lower()
            ]
        else:
            filtered = list(self.all_vars)
        self._populate_var_list(filtered)

    def _refresh_selected_list(self) -> None:
        sel_list = self.query_one("#selected-list", OptionList)
        sel_list.clear_options()
        if not self.selected_vars:
            sel_list.add_option(Option("[dim]None selected[/dim]", id="__none__"))
            return
        for name in self.selected_vars:
            info = next((v for v in self.all_vars if v["name"] == name), None)
            units = info["units"] if info else ""
            sel_list.add_option(Option(f"{name}  [dim]{units}[/dim]", id=name))

    def action_select_all(self) -> None:
        self.selected_vars = [v["name"] for v in self.all_vars]
        self._refresh_var_marks()
        self._refresh_selected_list()
        self.notify(f"Selected {len(self.selected_vars)} variables")

    def action_clear_selected(self) -> None:
        self.selected_vars.clear()
        self._refresh_var_marks()
        self._refresh_selected_list()
        self.notify("Cleared")

    # ── Actions (only compute when button is pressed) ──

    def _pre_action_check(self) -> bool:
        if not self.wf:
            self.notify("Open a file first", severity="warning")
            return False
        if not self.selected_vars:
            self.notify("Select variables first", severity="warning")
            return False
        return True

    def _get_timesteps(self) -> list[int]:
        text = self.query_one("#timestep-input", Input).value
        nt = self.wf.nt if self.wf else 1
        try:
            return _parse_timesteps(text, nt)
        except Exception:
            return [0]

    @on(Button.Pressed, "#btn-export")
    def _on_export(self) -> None:
        if self._pre_action_check():
            self._run_export()

    @on(Button.Pressed, "#btn-plot")
    def _on_plot(self) -> None:
        if self._pre_action_check():
            self._run_plot()

    @on(Button.Pressed, "#btn-gif")
    def _on_gif(self) -> None:
        if self._pre_action_check():
            ts = self._get_timesteps()
            if len(ts) < 2:
                self.notify("Need 2+ timesteps for GIF (e.g. '0-5' or 'all')", severity="warning")
                return
            self._run_gif()

    @on(Button.Pressed, "#btn-stats")
    def _on_stats(self) -> None:
        if self._pre_action_check():
            self._run_stats()

    @work(thread=True)
    def _run_export(self) -> None:
        from wrf import getvar
        outdir = os.path.dirname(self.wf_path) if self.wf_path else "."
        timesteps = self._get_timesteps()
        jobs = [(v, t) for v in self.selected_vars for t in timesteps]
        total = len(jobs)
        log_lines = []

        self.call_from_thread(self._reset_progress, total)

        for i, (varname, t) in enumerate(jobs):
            self.call_from_thread(self._set_progress_label,
                                  f"Exporting {varname} t={t}  ({i+1}/{total})")
            suffix = f"_t{t:04d}" if len(timesteps) > 1 else ""
            outpath = os.path.join(outdir, f"{varname}{suffix}.npy")
            try:
                data = getvar(self.wf, varname, timeidx=t)
                np.save(outpath, data)
                log_lines.append(f"[green]\u2713[/green] {varname} t={t}  {data.shape}")
            except Exception as e:
                log_lines.append(f"[red]\u2717[/red] {varname} t={t}: {e}")
            self.call_from_thread(self._advance_progress, i + 1, total)

        self.call_from_thread(self._set_progress_label, f"Done - {total} exports to {outdir}")
        self.call_from_thread(self._set_log, "\n".join(log_lines))

    @work(thread=True)
    def _run_plot(self) -> None:
        try:
            from wrf.plot import plot_field
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            self.call_from_thread(self.notify, "matplotlib not installed", severity="error")
            return

        outdir = os.path.dirname(self.wf_path) if self.wf_path else "."
        timesteps = self._get_timesteps()
        jobs = [(v, t) for v in self.selected_vars for t in timesteps]
        total = len(jobs)
        log_lines = []

        self.call_from_thread(self._reset_progress, total)

        for i, (varname, t) in enumerate(jobs):
            self.call_from_thread(self._set_progress_label,
                                  f"Plotting {varname} t={t}  ({i+1}/{total})")
            suffix = f"_t{t:04d}" if len(timesteps) > 1 else ""
            outpath = os.path.join(outdir, f"{varname}{suffix}.png")
            try:
                fig, _ = plot_field(self.wf, varname, timeidx=t)
                fig.savefig(outpath, dpi=150, bbox_inches="tight")
                plt.close(fig)
                log_lines.append(f"[green]\u2713[/green] {varname} t={t}")
            except Exception as e:
                log_lines.append(f"[red]\u2717[/red] {varname} t={t}: {e}")
            self.call_from_thread(self._advance_progress, i + 1, total)

        self.call_from_thread(self._set_progress_label, f"Done - {total} plots to {outdir}")
        self.call_from_thread(self._set_log, "\n".join(log_lines))

    @work(thread=True)
    def _run_gif(self) -> None:
        try:
            from wrf.plot import render_timesteps
            import matplotlib
            matplotlib.use("Agg")
        except ImportError:
            self.call_from_thread(self.notify, "matplotlib not installed", severity="error")
            return

        outdir = os.path.dirname(self.wf_path) if self.wf_path else "."
        timesteps = self._get_timesteps()
        total = len(self.selected_vars)
        log_lines = []

        self.call_from_thread(self._reset_progress, total)

        for i, varname in enumerate(self.selected_vars):
            self.call_from_thread(self._set_progress_label,
                                  f"Rendering {varname} ({len(timesteps)} frames)  ({i+1}/{total})")
            gif_path = os.path.join(outdir, f"{varname}.gif")
            try:
                def _progress(cur, tot, vn):
                    self.call_from_thread(self._set_progress_label,
                                          f"Rendering {vn} frame {cur}/{tot}  (var {i+1}/{total})")

                render_timesteps(
                    self.wf, varname,
                    timesteps=timesteps,
                    outdir=outdir,
                    gif=True,
                    gif_path=gif_path,
                    fixed_scale=True,
                    progress_callback=_progress,
                )
                log_lines.append(f"[green]\u2713[/green] {varname}  -> {os.path.basename(gif_path)}")
            except Exception as e:
                log_lines.append(f"[red]\u2717[/red] {varname}: {e}")

            self.call_from_thread(self._advance_progress, i + 1, total)

        self.call_from_thread(self._set_progress_label, f"Done - {total} GIFs to {outdir}")
        self.call_from_thread(self._set_log, "\n".join(log_lines))

    @work(thread=True)
    def _run_stats(self) -> None:
        from wrf import getvar
        timesteps = self._get_timesteps()
        jobs = [(v, t) for v in self.selected_vars for t in timesteps]
        total = len(jobs)

        self.call_from_thread(self._reset_progress, total)

        tbl = Table(box=box.ROUNDED, padding=(0, 1))
        tbl.add_column("Variable", style="bold")
        if len(timesteps) > 1:
            tbl.add_column("t", justify="right")
        tbl.add_column("Shape")
        tbl.add_column("Min", justify="right")
        tbl.add_column("Max", justify="right")
        tbl.add_column("Mean", justify="right")
        tbl.add_column("Std", justify="right")
        tbl.add_column("Units", style="dim")

        for i, (varname, t) in enumerate(jobs):
            self.call_from_thread(self._set_progress_label,
                                  f"Computing {varname} t={t}  ({i+1}/{total})")
            info = next((v for v in self.all_vars if v["name"] == varname), None)
            units = info["units"] if info else ""
            try:
                data = getvar(self.wf, varname, timeidx=t)
                valid = data[np.isfinite(data)]
                row = [varname]
                if len(timesteps) > 1:
                    row.append(str(t))
                if len(valid) > 0:
                    row.extend([
                        "x".join(str(d) for d in data.shape),
                        f"{valid.min():.4g}",
                        f"{valid.max():.4g}",
                        f"{valid.mean():.4g}",
                        f"{valid.std():.4g}",
                        units,
                    ])
                else:
                    row.extend(["x".join(str(d) for d in data.shape), "", "", "", "", "[dim]no data[/dim]"])
                tbl.add_row(*row)
            except Exception as e:
                row = [varname]
                if len(timesteps) > 1:
                    row.append(str(t))
                row.extend(["", "", "", "", "", f"[red]{e}[/red]"])
                tbl.add_row(*row)

            self.call_from_thread(self._advance_progress, i + 1, total)

        self.call_from_thread(self._set_progress_label, f"Done - {total} computations")
        self.call_from_thread(self._set_log, tbl)

    # ── Progress helpers ──

    def _reset_progress(self, total: int) -> None:
        pb = self.query_one("#progress-bar", ProgressBar)
        pb.update(total=total, progress=0)
        self.query_one("#output-log", Static).update("")

    def _advance_progress(self, current: int, total: int) -> None:
        pb = self.query_one("#progress-bar", ProgressBar)
        pb.update(total=total, progress=current)

    def _set_log(self, content) -> None:
        self.query_one("#output-log", Static).update(content)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
    app = WrfTui(path)
    app.run()


if __name__ == "__main__":
    main()
