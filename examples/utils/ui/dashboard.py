"""A generic matplotlib dashboard shell.

`Dashboard` is the entry point: you declare the knobs, say how the figure is
divided into panels, and fill each panel.  The shell owns the boring parts --
building the figure, laying out one widget per parameter, wiring callbacks,
and re-rendering on change -- so a script only writes the parts that are
specific to it.

    class MyView(Dashboard):
        title = "X_1(N) over F_p"

        def params(self):
            return [Param("N", 4, 30, 11, step=1),
                    Param("p", 2, 97, 13, step=1)]

        def panels(self):
            return ["locus", "points"]

        def draw(self, name, ax):
            ...                       # draw panel `name` onto `ax`

    MyView().show()                   # interactive window
    MyView().save("out.png")          # headless render

Parameters are `ui.base.Param`, the same objects `Drawable` uses, so a
Dashboard and a Drawable can share a parameter set.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Callable, Sequence

from ui.base import Param

__all__ = ["Dashboard", "Param", "TablePanel"]


class Dashboard(ABC):
    """A live matplotlib dashboard: parameter widgets plus one or more panels.

    Subclasses override `params()`, `panels()` and `draw()`.  Everything to do
    with figures, widgets and callbacks is handled here.
    """

    title: str = "dashboard"
    figsize: tuple[float, float] = (13.0, 7.0)
    control_width: float = 0.22

    def __init__(self, **overrides) -> None:
        self._params: dict[str, Param] = {p.name: p for p in self.params()}
        self._actions: dict[str, Callable[[], None]] = dict(self.actions())
        self._axes: dict[str, object] = {}
        self._widgets: dict[str, object] = {}
        self._fig = None
        self._live = False
        self._status: str = ""
        for key, value in overrides.items():
            self.set(key, value, redraw=False)

    def params(self) -> list[Param]:
        return []

    def panels(self) -> Sequence[str]:
        return ["main"]

    def actions(self) -> dict[str, Callable[[], None]]:
        return {}

    @abstractmethod
    def draw(self, name: str, ax) -> None:
        """Draw panel `name` onto `ax`.  Called for each panel on every update."""

    def header(self) -> str:
        return self.title

    def on_change(self, name: str, value) -> None:
        """Hook called after a parameter changes, before re-rendering."""

    @property
    def values(self) -> dict:
        return {name: p.value for name, p in self._params.items()}

    def __getitem__(self, key: str):
        return self._params[key].value

    def get(self, key: str, default=None):
        p = self._params.get(key)
        return default if p is None else p.value

    def set(self, key: str, value, redraw: bool = True):
        if key not in self._params:
            raise KeyError(f"unknown parameter {key!r}; have {sorted(self._params)}")
        new = self._params[key].set(value)
        self.on_change(key, new)
        if redraw:
            self.refresh()
        return new

    def update(self, **kwargs):
        for key, value in kwargs.items():
            self.set(key, value, redraw=False)
        self.refresh()
        return self

    def reset(self):
        for p in self._params.values():
            p.reset()
        self._sync_widgets()
        self.refresh()
        return self

    def status(self, message: str) -> None:
        self._status = message
        widget = self._widgets.get("__status__")
        if widget is not None:
            widget.set_text(message)
            self._draw_idle()

    def refresh(self) -> None:
        if self._fig is None:
            return
        for name, ax in self._axes.items():
            ax.cla()
            self.draw(name, ax)
        self._fig.suptitle(self.header())
        self._draw_idle()

    def _draw_idle(self) -> None:
        if self._fig is not None and self._fig.canvas is not None:
            self._fig.canvas.draw_idle()

    def build(self, interactive: bool = True):
        import matplotlib.pyplot as plt

        self._fig = plt.figure(figsize=self.figsize)
        names = list(self.panels())

        reserve = self.control_width if (interactive and self._has_widgets()) else 0.02
        grid = self._fig.add_gridspec(
            1, len(names), left=0.06, right=1.0 - reserve - 0.04, wspace=0.28
        )
        for i, name in enumerate(names):
            self._axes[name] = self._fig.add_subplot(grid[0, i])

        if interactive and self._has_widgets():
            self._build_widgets()

        self.refresh()
        return self._fig

    def show(self):
        import matplotlib.pyplot as plt

        self.build(interactive=True)
        self._live = True
        try:
            plt.show()
        finally:
            self._live = False
        return self._fig

    def save(self, path: str, dpi: int = 150):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        self.build(interactive=False)
        self._fig.savefig(path, dpi=dpi)
        plt.close(self._fig)
        self._fig = None
        self._axes.clear()
        return path

    def _slider_params(self) -> list[Param]:
        return [
            p for p in self._params.values() if p.enable_ui and p.type in (int, float)
        ]

    def _has_widgets(self) -> bool:
        return bool(self._slider_params() or self._actions)

    def _build_widgets(self) -> None:
        from matplotlib.widgets import Button, Slider

        left = 1.0 - self.control_width
        width = self.control_width - 0.04
        height, gap = 0.035, 0.018
        top = 0.88

        for p in self._slider_params():
            ax = self._fig.add_axes([left, top, width, height])
            slider = Slider(
                ax,
                p.label,
                p.slider_min,
                p.slider_max,
                valinit=p.slider_value,
                valstep=p.slider_step,
            )
            slider.valtext.set_text(p.format_value())
            slider.on_changed(self._slider_callback(p.name))
            self._widgets[p.name] = slider
            top -= height + gap

        top -= gap
        for label, callback in self._actions.items():
            ax = self._fig.add_axes([left, top, width, height])
            button = Button(ax, label)
            button.on_clicked(self._button_callback(callback))
            self._widgets[f"__button__{label}"] = button
            top -= height + gap

        self._widgets["__status__"] = self._fig.text(
            left, top - 0.01, self._status, fontsize=8, va="top", wrap=True
        )

    def _slider_callback(self, name: str):
        def _update(value):
            self._params[name].set_from_ui(value)
            self._widgets[name].valtext.set_text(self._params[name].format_value())
            self.on_change(name, self._params[name].value)
            self.refresh()

        return _update

    def _button_callback(self, callback: Callable[[], None]):
        def _clicked(_event):
            callback()
            self.refresh()

        return _clicked

    def _sync_widgets(self) -> None:
        for name, p in self._params.items():
            widget = self._widgets.get(name)
            if widget is not None and widget.val != p.slider_value:
                widget.set_val(p.slider_value)
            if widget is not None:
                widget.valtext.set_text(p.format_value())


class TablePanel:
    """Helper for rendering rows of text into a panel instead of a plot."""

    def __init__(self, ax, fontsize: int = 9) -> None:
        self.ax = ax
        self.fontsize = fontsize
        ax.set_axis_off()

    def rows(self, headers: Sequence[str], rows: Sequence[Sequence], max_rows: int = 25):
        body = [[str(cell) for cell in row] for row in list(rows)[:max_rows]]
        if not body:
            self.ax.text(
                0.5,
                0.5,
                "no rows",
                ha="center",
                va="center",
                fontsize=self.fontsize,
            )
            return self
        table = self.ax.table(
            cellText=body,
            colLabels=list(headers),
            loc="upper center",
            cellLoc="right",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(self.fontsize)
        return self

    def text(self, lines: Sequence[str]):
        self.ax.text(
            0.01,
            0.99,
            "\n".join(lines),
            family="monospace",
            fontsize=self.fontsize,
            va="top",
            ha="left",
            transform=self.ax.transAxes,
        )
        return self
