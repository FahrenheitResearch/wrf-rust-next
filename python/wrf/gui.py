"""
Standalone WRF Explorer GUI.

Launch with:
    python -m wrf.gui [wrfout_file]
    python -m wrf gui [wrfout_file]
"""

import os
import sys
import tkinter as tk
from tkinter import ttk, filedialog
import threading

import numpy as np


def _get_dll_dirs():
    """Find HDF5/NetCDF DLL directories."""
    dirs = []
    home = os.path.expanduser("~")
    for base in ("miniforge3", "miniconda3", "anaconda3"):
        for env in ("wrfplot", "base", "wxsection"):
            d = os.path.join(home, base, "envs", env, "Library", "bin")
            if os.path.isdir(d):
                dirs.append(d)
        d = os.path.join(home, base, "Library", "bin")
        if os.path.isdir(d):
            dirs.append(d)
    return dirs


# Variable -> (colormap, vmin, vmax) defaults
_VAR_STYLES = {
    # CAPE
    "sbcape": ("YlOrRd", 0, 4000), "mlcape": ("YlOrRd", 0, 4000),
    "mucape": ("YlOrRd", 0, 4000), "cape": ("YlOrRd", 0, 4000),
    "effective_cape": ("YlOrRd", 0, 4000),
    # CIN
    "sbcin": ("Blues_r", -300, 0), "mlcin": ("Blues_r", -300, 0),
    "mucin": ("Blues_r", -300, 0), "cin": ("Blues_r", -300, 0),
    # SRH
    "srh1": ("RdPu", 0, 500), "srh3": ("RdPu", 0, 500),
    "srh": ("RdPu", 0, 500), "effective_srh": ("RdPu", 0, 500),
    # Severe
    "stp": ("Reds", 0, 10), "scp": ("Reds", 0, 10),
    "ehi": ("Reds", 0, 5), "ship": ("Reds", 0, 3),
    "critical_angle": ("RdYlBu_r", 0, 180),
    # Temperature
    "temp": ("RdBu_r", None, None), "tc": ("RdBu_r", None, None),
    "tk": ("RdBu_r", None, None),
    "theta": ("Spectral_r", None, None), "theta_e": ("Spectral_r", None, None),
    "tv": ("RdBu_r", None, None), "twb": ("RdBu_r", None, None),
    "td": ("YlGnBu", None, None),
    # Pressure
    "slp": ("viridis", None, None), "pressure": ("viridis", None, None),
    # Wind
    "wspd": ("YlGnBu", 0, 40), "wspd10": ("YlGnBu", 0, 40),
    "wdir": ("hsv", 0, 360), "wdir10": ("hsv", 0, 360),
    "shear_0_1km": ("YlOrRd", 0, 30), "shear_0_6km": ("YlOrRd", 0, 40),
    "bulk_shear": ("YlOrRd", 0, 40),
    # Radar
    "dbz": ("gist_ncar", -10, 75), "maxdbz": ("gist_ncar", -10, 75),
    # Moisture
    "rh": ("YlGnBu", 0, 100), "rh2m": ("YlGnBu", 0, 100),
    "pw": ("YlGnBu", 0, 60),
    # Height/terrain
    "terrain": ("terrain", None, None), "height": ("terrain", None, None),
    "height_agl": ("terrain", None, None),
    "freezing_level": ("cool", None, None), "wet_bulb_0": ("cool", None, None),
    # Vorticity
    "avo": ("RdBu_r", None, None), "pvo": ("RdBu_r", None, None),
    # Lapse rates
    "lapse_rate_700_500": ("RdYlBu_r", 4, 10),
    "lapse_rate_0_3km": ("RdYlBu_r", 4, 10),
    "lapse_rate": ("RdYlBu_r", 4, 10),
    # Cloud
    "ctt": ("Greys_r", -80, 20), "cloudfrac": ("Greys", 0, 100),
    # Fire
    "fosberg": ("YlOrRd", 0, 100), "haines": ("YlOrRd", 1, 6),
    "hdw": ("YlOrRd", 0, 200),
    # Helicity
    "uhel": ("RdPu", 0, 200),
    # Levels
    "lcl": ("YlGnBu_r", 0, 3000), "lfc": ("YlGnBu_r", 0, 5000),
    "el": ("YlGnBu_r", 0, 15000),
    "geopt": ("terrain", None, None),
    "omega": ("RdBu_r", None, None),
    "mixing_ratio": ("YlGnBu", 0, 0.02),
    "specific_humidity": ("YlGnBu", 0, 0.02),
    "dp2m": ("YlGnBu", None, None),
}


def _style_for(varname):
    """Get (cmap, vmin, vmax) for a variable name."""
    return _VAR_STYLES.get(varname, ("viridis", None, None))


class WrfExplorerApp:
    """Tkinter-based WRF file explorer with live-updating plots."""

    def __init__(self, filepath=None):
        self.root = tk.Tk()
        self.root.title("wrf-rust Explorer")
        self.root.geometry("1200x800")
        self.root.minsize(900, 600)

        # State
        self.wf = None
        self.var_list = []
        self.current_data = None
        self.lat = None
        self.lon = None

        self._build_ui()

        if filepath:
            self._open_file(filepath)

    def _build_ui(self):
        # ── Top toolbar ──
        toolbar = ttk.Frame(self.root, padding=5)
        toolbar.pack(fill=tk.X, side=tk.TOP)

        ttk.Button(toolbar, text="Open File...", command=self._browse_file).pack(side=tk.LEFT, padx=2)
        self.file_label = ttk.Label(toolbar, text="No file loaded", font=("Segoe UI", 9))
        self.file_label.pack(side=tk.LEFT, padx=10)

        # ── Main split: controls on left, plot on right ──
        paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Left panel: controls
        ctrl_frame = ttk.Frame(paned, padding=5)
        paned.add(ctrl_frame, weight=0)

        # Variable selector
        ttk.Label(ctrl_frame, text="Variable:", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W)
        self.var_combo = ttk.Combobox(ctrl_frame, state="readonly", width=30)
        self.var_combo.pack(fill=tk.X, pady=(0, 8))
        self.var_combo.bind("<<ComboboxSelected>>", self._on_var_change)

        # Search filter
        ttk.Label(ctrl_frame, text="Filter:").pack(anchor=tk.W)
        self.filter_var = tk.StringVar()
        self.filter_var.trace_add("write", self._on_filter_change)
        ttk.Entry(ctrl_frame, textvariable=self.filter_var, width=30).pack(fill=tk.X, pady=(0, 8))

        # Variable info
        self.var_info = tk.Text(ctrl_frame, height=3, width=30, font=("Consolas", 9),
                                wrap=tk.WORD, state=tk.DISABLED, bg="#f0f0f0")
        self.var_info.pack(fill=tk.X, pady=(0, 8))

        # Time slider
        ttk.Label(ctrl_frame, text="Time step:", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W)
        self.time_var = tk.IntVar(value=0)
        self.time_slider = ttk.Scale(ctrl_frame, from_=0, to=0, variable=self.time_var,
                                      orient=tk.HORIZONTAL, command=self._on_slider_change)
        self.time_slider.pack(fill=tk.X, pady=(0, 2))
        self.time_label = ttk.Label(ctrl_frame, text="t=0")
        self.time_label.pack(anchor=tk.W, pady=(0, 8))

        # Level slider
        ttk.Label(ctrl_frame, text="Vertical level:", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W)
        self.level_var = tk.IntVar(value=0)
        self.level_slider = ttk.Scale(ctrl_frame, from_=0, to=0, variable=self.level_var,
                                       orient=tk.HORIZONTAL, command=self._on_slider_change)
        self.level_slider.pack(fill=tk.X, pady=(0, 2))
        self.level_label = ttk.Label(ctrl_frame, text="k=0 (lowest)")
        self.level_label.pack(anchor=tk.W, pady=(0, 8))

        # Units
        ttk.Label(ctrl_frame, text="Units (optional):", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W)
        self.units_var = tk.StringVar()
        units_entry = ttk.Entry(ctrl_frame, textvariable=self.units_var, width=20)
        units_entry.pack(fill=tk.X, pady=(0, 8))
        units_entry.bind("<Return>", self._on_slider_change)

        # Update button
        ttk.Button(ctrl_frame, text="Update Plot", command=self._update_plot).pack(fill=tk.X, pady=5)

        # Stats display
        ttk.Label(ctrl_frame, text="Stats:", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W, pady=(10, 0))
        self.stats_text = tk.Text(ctrl_frame, height=5, width=30, font=("Consolas", 9),
                                   wrap=tk.WORD, state=tk.DISABLED, bg="#f0f0f0")
        self.stats_text.pack(fill=tk.X)

        # Right panel: plot
        plot_frame = ttk.Frame(paned)
        paned.add(plot_frame, weight=1)

        # Matplotlib canvas
        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        from matplotlib.figure import Figure

        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_facecolor("#f8f8f8")
        self.ax.text(0.5, 0.5, "Open a WRF file to begin",
                     ha="center", va="center", fontsize=14, color="#888",
                     transform=self.ax.transAxes)

        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.draw()

        # Toolbar
        nav_frame = ttk.Frame(plot_frame)
        nav_frame.pack(side=tk.BOTTOM, fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.canvas, nav_frame)
        self.toolbar.update()

        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Status bar
        self.status = ttk.Label(self.root, text="Ready", relief=tk.SUNKEN, padding=2)
        self.status.pack(fill=tk.X, side=tk.BOTTOM)

    def _browse_file(self):
        path = filedialog.askopenfilename(
            title="Open WRF Output File",
            filetypes=[
                ("WRF Output", "wrfout*"),
                ("NetCDF Files", "*.nc *.nc4"),
                ("All Files", "*.*"),
            ],
        )
        if path:
            self._open_file(path)

    def _open_file(self, path):
        self._set_status(f"Opening {os.path.basename(path)}...")
        try:
            from wrf import WrfFile, list_variables, getvar
            self.wf = WrfFile(path)
            self.file_label.config(text=f"{os.path.basename(path)}  ({self.wf.nx}x{self.wf.ny}x{self.wf.nz}, {self.wf.nt}t)")

            # Get lat/lon for georeferenced plotting
            try:
                self.lat = getvar(self.wf, "lat", timeidx=0)
                self.lon = getvar(self.wf, "lon", timeidx=0)
            except Exception:
                self.lat = None
                self.lon = None

            # Populate variable list
            self.var_list = list_variables()
            var_names = [f"{v['name']:20s} {v['description']}" for v in self.var_list]
            self.var_combo["values"] = var_names
            if var_names:
                self.var_combo.current(0)

            # Set slider ranges
            self.time_slider.config(to=max(0, self.wf.nt - 1))
            self.level_slider.config(to=max(0, self.wf.nz - 1))

            self._set_status(f"Loaded: {self.wf.nx}x{self.wf.ny}x{self.wf.nz} grid, {self.wf.nt} time steps")
            self._on_var_change(None)

        except Exception as e:
            self._set_status(f"Error: {e}")
            import traceback
            traceback.print_exc()

    def _get_selected_varname(self):
        sel = self.var_combo.get()
        if not sel:
            return None
        return sel.split()[0].strip()

    def _on_filter_change(self, *_args):
        filt = self.filter_var.get().lower()
        if not self.var_list:
            return
        if filt:
            filtered = [f"{v['name']:20s} {v['description']}"
                        for v in self.var_list
                        if filt in v["name"].lower() or filt in v["description"].lower()]
        else:
            filtered = [f"{v['name']:20s} {v['description']}" for v in self.var_list]
        self.var_combo["values"] = filtered
        if filtered:
            self.var_combo.current(0)

    def _on_var_change(self, _event):
        varname = self._get_selected_varname()
        if not varname:
            return
        # Update info box
        info = next((v for v in self.var_list if v["name"] == varname), None)
        self.var_info.config(state=tk.NORMAL)
        self.var_info.delete("1.0", tk.END)
        if info:
            self.var_info.insert(tk.END, f"{info['description']}\nUnits: {info['units']}")
        self.var_info.config(state=tk.DISABLED)

        self._update_plot()

    def _on_slider_change(self, *_args):
        t = int(self.time_var.get())
        k = int(self.level_var.get())
        self.time_label.config(text=f"t={t}")
        self.level_label.config(text=f"k={k}" + (" (lowest)" if k == 0 else ""))

    def _update_plot(self):
        varname = self._get_selected_varname()
        if not varname or not self.wf:
            return

        self._set_status(f"Computing {varname}...")
        # Run computation in background to keep UI responsive
        threading.Thread(target=self._do_plot, args=(varname,), daemon=True).start()

    def _do_plot(self, varname):
        try:
            from wrf import getvar

            t = int(self.time_var.get())
            k = int(self.level_var.get())
            units = self.units_var.get().strip() or None

            kwargs = {}
            if units:
                kwargs["units"] = units

            data = getvar(self.wf, varname, timeidx=t, **kwargs)

            # Handle 3D data: take a single level
            is_3d = data.ndim == 3
            if is_3d:
                level = min(k, data.shape[0] - 1)
                plot_data = data[level]
                level_str = f" [k={level}]"
            else:
                plot_data = data
                level_str = ""

            # Update level slider visibility info
            self.root.after(0, lambda: self.level_label.config(
                text=f"k={k}" + (" (3D)" if is_3d else " (2D, ignored)")
            ))

            # Get style
            cmap, vmin, vmax = _style_for(varname)
            if vmin is None:
                valid = plot_data[np.isfinite(plot_data)]
                if len(valid) > 0:
                    vmin = np.percentile(valid, 2)
                    vmax = np.percentile(valid, 98)
                else:
                    vmin, vmax = 0, 1

            # Compute stats
            valid = plot_data[np.isfinite(plot_data)]
            if len(valid) > 0:
                stats_str = (
                    f"Min:  {valid.min():.2f}\n"
                    f"Max:  {valid.max():.2f}\n"
                    f"Mean: {valid.mean():.2f}\n"
                    f"Std:  {valid.std():.2f}"
                )
            else:
                stats_str = "No valid data"

            # Schedule UI update on main thread
            self.root.after(0, lambda: self._render_plot(
                plot_data, varname, level_str, cmap, vmin, vmax, stats_str, units, t
            ))

        except Exception as e:
            self.root.after(0, lambda: self._set_status(f"Error: {e}"))
            import traceback
            traceback.print_exc()

    def _render_plot(self, data, varname, level_str, cmap, vmin, vmax, stats_str, units, t):
        self.fig.clear()
        ax = self.fig.add_subplot(111)

        if self.lat is not None and self.lon is not None:
            im = ax.pcolormesh(self.lon, self.lat, data, cmap=cmap,
                               vmin=vmin, vmax=vmax, shading="auto")
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
        else:
            im = ax.pcolormesh(data, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
            ax.set_xlabel("West-East")
            ax.set_ylabel("South-North")

        # Colorbar
        unit_str = units or next(
            (v["units"] for v in self.var_list if v["name"] == varname), ""
        )
        cbar = self.fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
        if unit_str:
            cbar.set_label(unit_str, fontsize=10)

        # Title
        desc = next((v["description"] for v in self.var_list if v["name"] == varname), varname)
        ax.set_title(f"{desc}{level_str}\nt={t}", fontsize=12)
        ax.set_aspect("equal" if self.lat is not None else "auto")

        self.fig.tight_layout()
        self.canvas.draw()

        # Update stats
        self.stats_text.config(state=tk.NORMAL)
        self.stats_text.delete("1.0", tk.END)
        self.stats_text.insert(tk.END, stats_str)
        self.stats_text.config(state=tk.DISABLED)

        self._set_status(f"{varname}: {data.shape[0]}x{data.shape[1]}")

    def _set_status(self, text):
        self.status.config(text=text)

    def run(self):
        self.root.mainloop()


def main():
    filepath = sys.argv[1] if len(sys.argv) > 1 else None
    app = WrfExplorerApp(filepath)
    app.run()


if __name__ == "__main__":
    main()
