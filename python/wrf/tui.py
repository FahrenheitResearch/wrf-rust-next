"""
wrf-rust Terminal UI.

Launch with:
    python -m wrf tui [wrfout_file]
    python -m wrf.tui [wrfout_file]
"""

from __future__ import annotations

import os
import sys
from typing import Optional

import numpy as np

from textual import on, work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.reactive import reactive
from textual.widgets import (
    Footer,
    Header,
    Input,
    Label,
    OptionList,
    Static,
)
from textual.widgets.option_list import Option
from textual.widget import Widget

from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich import box


# ── Helpers ──────────────────────────────────────────────────────────────────

def _load_wrf(path: str):
    """Open a WRF file, returning (WrfFile, error_string)."""
    try:
        from wrf import WrfFile
        return WrfFile(path), None
    except Exception as e:
        return None, str(e)


def _get_var_list():
    """Get list of variable dicts from wrf-core."""
    try:
        from wrf import list_variables
        return list_variables()
    except Exception:
        return []


def _compute_var(wf, varname: str, timeidx: int, level: int | None, units: str | None):
    """Compute a variable and return (data_2d, stats_dict, error)."""
    try:
        from wrf import getvar
        kwargs = {}
        if units:
            kwargs["units"] = units
        data = getvar(wf, varname, timeidx=timeidx, **kwargs)
        is_3d = data.ndim == 3
        if is_3d and level is not None:
            k = min(level, data.shape[0] - 1)
            data_2d = data[k]
        elif is_3d:
            data_2d = data[0]
        else:
            data_2d = data

        valid = data_2d[np.isfinite(data_2d)]
        stats = {}
        if len(valid) > 0:
            stats = {
                "min": float(np.min(valid)),
                "max": float(np.max(valid)),
                "mean": float(np.mean(valid)),
                "std": float(np.std(valid)),
                "p5": float(np.percentile(valid, 5)),
                "p25": float(np.percentile(valid, 25)),
                "p50": float(np.percentile(valid, 50)),
                "p75": float(np.percentile(valid, 75)),
                "p95": float(np.percentile(valid, 95)),
            }
        return data_2d, stats, is_3d, None
    except Exception as e:
        return None, {}, False, str(e)


def _spark_histogram(data: np.ndarray, width: int = 40) -> str:
    """Make a sparkline-style histogram from 2D data."""
    valid = data[np.isfinite(data)].ravel()
    if len(valid) == 0:
        return ""
    bars = " ▁▂▃▄▅▆▇█"
    counts, _ = np.histogram(valid, bins=width)
    if counts.max() == 0:
        return bars[0] * width
    scaled = (counts / counts.max() * (len(bars) - 1)).astype(int)
    return "".join(bars[s] for s in scaled)


def _ascii_mini_map(data: np.ndarray, width: int = 60, height: int = 20) -> str:
    """Downsample 2D data to a tiny ASCII representation using block chars."""
    if data is None:
        return ""
    ny, nx = data.shape

    # Downsample
    step_y = max(1, ny // height)
    step_x = max(1, nx // width)
    small = data[::step_y, ::step_x]

    valid = small[np.isfinite(small)]
    if len(valid) == 0:
        return "No valid data"

    vmin = np.percentile(valid, 2)
    vmax = np.percentile(valid, 98)
    if vmax == vmin:
        vmax = vmin + 1

    # Normalize to 0-1
    norm = np.clip((small - vmin) / (vmax - vmin), 0, 1)

    # Shade characters (light to dark)
    shades = " ░▒▓█"
    lines = []
    sh, sw = norm.shape
    for j in range(sh):
        row = ""
        for i in range(sw):
            v = norm[j, i]
            if not np.isfinite(v):
                row += " "
            else:
                idx = int(v * (len(shades) - 1))
                row += shades[idx]
        lines.append(row)

    return "\n".join(lines)


# ── Widgets ──────────────────────────────────────────────────────────────────

class FileInfo(Static):
    """Displays file metadata."""

    def render(self) -> Text:
        return Text("No file loaded", style="dim")


class VarFilter(Input):
    """Filter input for variable list."""
    pass


class StatsPanel(Static):
    """Displays statistics for the current variable."""

    def set_stats(self, varname: str, stats: dict, units: str, is_3d: bool, level: int):
        if not stats:
            self.update(Panel("No data", title="Stats", border_style="dim"))
            return

        tbl = Table(box=box.SIMPLE, show_header=False, padding=(0, 1))
        tbl.add_column("Stat", style="bold cyan", width=6)
        tbl.add_column("Value", justify="right")

        tbl.add_row("Min", f"{stats['min']:.4g}")
        tbl.add_row("P5", f"{stats['p5']:.4g}")
        tbl.add_row("P25", f"{stats['p25']:.4g}")
        tbl.add_row("Med", f"{stats['p50']:.4g}")
        tbl.add_row("P75", f"{stats['p75']:.4g}")
        tbl.add_row("P95", f"{stats['p95']:.4g}")
        tbl.add_row("Max", f"{stats['max']:.4g}")
        tbl.add_row("Mean", f"{stats['mean']:.4g}")
        tbl.add_row("Std", f"{stats['std']:.4g}")

        shape_str = "3D" if is_3d else "2D"
        title = f"{varname} ({units}) [{shape_str}, k={level}]" if is_3d else f"{varname} ({units}) [{shape_str}]"
        self.update(Panel(tbl, title=title, border_style="green"))

    def set_error(self, varname: str, err: str):
        self.update(Panel(f"[red]{err}[/red]", title=varname, border_style="red"))


class MapPreview(Static):
    """ASCII map preview of the data field."""

    def set_data(self, data: np.ndarray | None, varname: str = ""):
        if data is None:
            self.update(Panel("No data", title="Preview", border_style="dim"))
            return

        ascii_map = _ascii_mini_map(data)
        hist = _spark_histogram(data)

        content = f"{ascii_map}\n\n[cyan]Distribution:[/cyan] {hist}"
        self.update(Panel(content, title=f"Preview: {varname}", border_style="blue"))


# ── Main App ─────────────────────────────────────────────────────────────────

class WrfTui(App):
    """WRF Terminal Explorer."""

    CSS = """
    Screen {
        layout: grid;
        grid-size: 3 1;
        grid-columns: 28 1fr 36;
        grid-gutter: 1;
    }

    #left-panel {
        height: 100%;
        border: solid $primary-background;
        padding: 0 1;
    }

    #center-panel {
        height: 100%;
    }

    #right-panel {
        height: 100%;
        border: solid $primary-background;
        padding: 0 1;
    }

    #var-list {
        height: 1fr;
        margin-top: 1;
    }

    #filter-input {
        margin-bottom: 0;
    }

    #file-info {
        height: auto;
        margin-bottom: 1;
        padding: 1;
        background: $surface;
    }

    #map-preview {
        height: 1fr;
    }

    #stats-panel {
        height: auto;
        margin-bottom: 1;
    }

    #controls {
        height: auto;
        padding: 1;
        background: $surface;
        margin-bottom: 1;
    }

    .control-label {
        margin-top: 1;
        color: $text-muted;
    }
    """

    BINDINGS = [
        Binding("q", "quit", "Quit"),
        Binding("r", "refresh", "Refresh"),
        Binding("up", "level_up", "Level +", show=False),
        Binding("down", "level_down", "Level -", show=False),
        Binding("left", "time_prev", "Time -", show=False),
        Binding("right", "time_next", "Time +", show=False),
        Binding("j", "level_down", "Level -", show=False),
        Binding("k", "level_up", "Level +", show=False),
        Binding("h", "time_prev", "Time -", show=False),
        Binding("l", "time_next", "Time +", show=False),
    ]

    TITLE = "wrf-rust Explorer"

    # Reactive state
    current_var: reactive[str] = reactive("")
    current_time: reactive[int] = reactive(0)
    current_level: reactive[int] = reactive(0)
    current_units: reactive[str] = reactive("")

    def __init__(self, filepath: str | None = None):
        super().__init__()
        self.filepath = filepath
        self.wf = None
        self.all_vars = []
        self.filtered_vars = []

    def compose(self) -> ComposeResult:
        yield Header()

        # Left: variable list
        with Vertical(id="left-panel"):
            yield Label("[bold]Variables[/bold]")
            yield Input(placeholder="Filter...", id="filter-input")
            yield OptionList(id="var-list")

        # Center: map preview
        with Vertical(id="center-panel"):
            yield Static(id="file-info")
            yield MapPreview(id="map-preview")

        # Right: stats + controls
        with Vertical(id="right-panel"):
            yield StatsPanel(id="stats-panel")
            with Vertical(id="controls"):
                yield Label("[bold]Controls[/bold]")
                yield Label("Time step:", classes="control-label")
                yield Input(value="0", id="time-input", type="integer")
                yield Label("Level (3D vars):", classes="control-label")
                yield Input(value="0", id="level-input", type="integer")
                yield Label("Units:", classes="control-label")
                yield Input(placeholder="e.g. degC, knots, hPa", id="units-input")

        yield Footer()

    def on_mount(self) -> None:
        self.all_vars = _get_var_list()
        self.filtered_vars = list(self.all_vars)
        self._populate_var_list(self.all_vars)

        if self.filepath:
            self._open_file(self.filepath)

    def _open_file(self, path: str) -> None:
        wf, err = _load_wrf(path)
        if err:
            self.notify(f"Error: {err}", severity="error", timeout=5)
            return
        self.wf = wf

        # Update file info
        info_widget = self.query_one("#file-info", Static)
        tbl = Table(box=box.SIMPLE, show_header=False, padding=(0, 1))
        tbl.add_column("Key", style="bold")
        tbl.add_column("Value")
        tbl.add_row("File", os.path.basename(path))
        tbl.add_row("Grid", f"{wf.nx} x {wf.ny} x {wf.nz}")
        tbl.add_row("Times", str(wf.nt))
        tbl.add_row("Spacing", f"{wf.dx:g} m")
        info_widget.update(Panel(tbl, title="WRF File", border_style="cyan"))

        self.sub_title = f"{os.path.basename(path)}  [{wf.nx}x{wf.ny}x{wf.nz}]"

        # Auto-select first variable
        if self.all_vars:
            var_list = self.query_one("#var-list", OptionList)
            if var_list.option_count > 0:
                var_list.highlighted = 0
                self.current_var = self.all_vars[0]["name"]
                self._do_compute()

    def _populate_var_list(self, vars_list: list[dict]) -> None:
        var_list = self.query_one("#var-list", OptionList)
        var_list.clear_options()
        for v in vars_list:
            label = f"[bold]{v['name']}[/bold]  [dim]{v['units']}[/dim]\n  {v['description']}"
            var_list.add_option(Option(label, id=v["name"]))

    @on(Input.Changed, "#filter-input")
    def _on_filter(self, event: Input.Changed) -> None:
        query = event.value.lower().strip()
        if query:
            self.filtered_vars = [
                v for v in self.all_vars
                if query in v["name"].lower() or query in v["description"].lower()
            ]
        else:
            self.filtered_vars = list(self.all_vars)
        self._populate_var_list(self.filtered_vars)

    @on(OptionList.OptionHighlighted, "#var-list")
    def _on_var_select(self, event: OptionList.OptionHighlighted) -> None:
        if event.option and event.option.id:
            self.current_var = str(event.option.id)
            self._do_compute()

    @on(Input.Submitted, "#time-input")
    def _on_time_submit(self, event: Input.Submitted) -> None:
        try:
            self.current_time = int(event.value)
        except ValueError:
            pass
        self._do_compute()

    @on(Input.Submitted, "#level-input")
    def _on_level_submit(self, event: Input.Submitted) -> None:
        try:
            self.current_level = int(event.value)
        except ValueError:
            pass
        self._do_compute()

    @on(Input.Submitted, "#units-input")
    def _on_units_submit(self, event: Input.Submitted) -> None:
        self.current_units = event.value.strip()
        self._do_compute()

    def action_refresh(self) -> None:
        self._do_compute()

    def action_level_up(self) -> None:
        lvl_input = self.query_one("#level-input", Input)
        self.current_level = min(self.current_level + 1, (self.wf.nz - 1) if self.wf else 0)
        lvl_input.value = str(self.current_level)
        self._do_compute()

    def action_level_down(self) -> None:
        lvl_input = self.query_one("#level-input", Input)
        self.current_level = max(0, self.current_level - 1)
        lvl_input.value = str(self.current_level)
        self._do_compute()

    def action_time_next(self) -> None:
        t_input = self.query_one("#time-input", Input)
        self.current_time = min(self.current_time + 1, (self.wf.nt - 1) if self.wf else 0)
        t_input.value = str(self.current_time)
        self._do_compute()

    def action_time_prev(self) -> None:
        t_input = self.query_one("#time-input", Input)
        self.current_time = max(0, self.current_time - 1)
        t_input.value = str(self.current_time)
        self._do_compute()

    @work(thread=True)
    def _do_compute(self) -> None:
        """Compute the selected variable in a background thread."""
        if not self.wf or not self.current_var:
            return

        varname = self.current_var
        timeidx = self.current_time
        level = self.current_level
        units = self.current_units or None

        data, stats, is_3d, err = _compute_var(self.wf, varname, timeidx, level, units)

        # Find default units
        var_info = next((v for v in self.all_vars if v["name"] == varname), None)
        display_units = units or (var_info["units"] if var_info else "")

        # Update UI on the app thread
        self.call_from_thread(self._update_display, varname, data, stats, display_units, is_3d, level, err)

    def _update_display(self, varname, data, stats, units, is_3d, level, err):
        stats_panel = self.query_one("#stats-panel", StatsPanel)
        map_preview = self.query_one("#map-preview", MapPreview)

        if err:
            stats_panel.set_error(varname, err)
            map_preview.set_data(None)
        else:
            stats_panel.set_stats(varname, stats, units, is_3d, level)
            map_preview.set_data(data, varname)


def main():
    filepath = sys.argv[1] if len(sys.argv) > 1 else None
    app = WrfTui(filepath)
    app.run()


if __name__ == "__main__":
    main()
