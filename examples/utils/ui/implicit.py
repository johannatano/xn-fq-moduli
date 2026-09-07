from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
import inspect

from ui.base import Drawable, Param

__all__ = [
    "ImplicitFunctionView",
    "PointSet2DView",
    "plot",
]


class PointSet2DView(Drawable):
    """A generic finite point set drawn in the (x, y) plane."""

    def __init__(
        self,
        point_source,
        *,
        x_coord: Callable | None = None,
        y_coord: Callable | None = None,
        title_text: str = "point set",
        title_builder: Callable[[dict], str] | None = None,
        xlabel: str = "x",
        ylabel: str = "y",
        params: Sequence[Param] | None = None,
        group_key: Callable | None = None,
        group_label: Callable | None = None,
        xlim=None,
        ylim=None,
        aspect: str | float = "equal",
        point_size: float = 45,
    ) -> None:
        self._point_source = point_source
        self._x_coord = x_coord or (lambda point: point[0])
        self._y_coord = y_coord or (lambda point: point[1])
        self._title_text = title_text
        self._title_builder = title_builder
        self.xlabel = xlabel
        self.ylabel = ylabel
        self._params_def = list(params or [])
        self._group_key = group_key
        self._group_label = group_label
        self._xlim = xlim
        self._ylim = ylim
        self._aspect = aspect
        self._point_size = point_size
        super().__init__()

    @classmethod
    def from_moduli_problem(cls, problem, F=None, p: int = 11) -> "PointSet2DView":
        from nt.moduli import Fp

        if F is not None:
            p = getattr(F, "p", p)

        params = [Param("p", 2, 97, p, step=1, label="prime p")]

        def point_source(values: dict):
            field = Fp(int(values["p"]))
            return list(problem.objects(field))

        def group_key(point, values: dict):
            field = Fp(int(values["p"]))
            return point.j(field)

        def group_label(group, points, values: dict):
            field = Fp(int(values["p"]))
            return f"j={group}, a={points[0].trace(field)}, {len(points)} pts"

        def title_builder(values: dict) -> str:
            return f"{problem!s} over F_{values['p']}"

        def xlim(values: dict):
            p_value = int(values["p"])
            return (-1, p_value)

        def ylim(values: dict):
            p_value = int(values["p"])
            return (-1, p_value)

        return cls(
            point_source,
            x_coord=lambda point: point.b,
            y_coord=lambda point: point.c,
            title_builder=title_builder,
            xlabel="x",
            ylabel="y",
            params=params,
            group_key=group_key,
            group_label=group_label,
            xlim=xlim,
            ylim=ylim,
        )

    def params(self) -> list[Param]:
        return list(self._params_def)

    def title(self) -> str:
        if self._title_builder is not None:
            return self._title_builder(self.param_values)
        return self._title_text

    def render(self, ax) -> None:
        import matplotlib.pyplot as plt
        from collections import defaultdict

        points = list(self._resolve_point_source())
        grouped = self._group_points(points)

        ordered = list(grouped.items())
        for (group, pts), color in zip(ordered, plt.cm.tab10.colors * 8):
            kwargs = {"s": self._point_size, "color": color}
            label = None
            if self._group_label is not None:
                label = self._group_label(group, pts, self.param_values)
            elif group is not None:
                label = str(group)
            if label:
                kwargs["label"] = label
            ax.scatter(
                [self._x_coord(point) for point in pts],
                [self._y_coord(point) for point in pts],
                **kwargs,
            )

        ax.set_aspect(self._aspect)
        xlim = self._resolve_limit(self._xlim)
        ylim = self._resolve_limit(self._ylim)
        if xlim is not None:
            ax.set_xlim(*xlim)
        if ylim is not None:
            ax.set_ylim(*ylim)
        if grouped and any(label for label in ax.get_legend_handles_labels()[1]):
            ax.legend(fontsize=7, loc="upper right")

    def _resolve_point_source(self):
        if callable(self._point_source):
            return self._call_with_optional_values(self._point_source)
        return self._point_source

    def _group_points(self, points: Iterable) -> dict:
        if self._group_key is None:
            return {None: list(points)}
        grouped = defaultdict(list)
        for point in points:
            grouped[self._group_key(point, self.param_values)].append(point)
        return dict(grouped)

    def _resolve_limit(self, limit):
        if limit is None:
            return None
        if callable(limit):
            return self._call_with_optional_values(limit)
        return limit

    def _call_with_optional_values(self, func, *args):
        try:
            arity = len(inspect.signature(func).parameters)
        except (TypeError, ValueError):
            arity = len(args) + 1
        if arity <= len(args):
            return func(*args)
        return func(*args, self.param_values)


class ImplicitFunctionView(Drawable):
    """The real zero-locus of an implicit equation poly(x, y) = 0."""

    def __init__(
        self,
        func,
        *,
        box=(-3.0, 3.0),
        grid: int = 800,
        title_text: str = "implicit locus",
        title_builder: Callable[[dict], str] | None = None,
        xlabel: str = "x",
        ylabel: str = "y",
        params: Sequence[Param] | None = None,
        level=0,
    ) -> None:
        self._func = func
        self.box = box
        self.grid = grid
        self._title_text = title_text
        self._title_builder = title_builder
        self.xlabel = xlabel
        self.ylabel = ylabel
        self._params_def = list(params or [])
        self.level = level
        super().__init__()

    @classmethod
    def from_equation(
        cls,
        poly,
        symbols,
        *,
        box=(-3.0, 3.0),
        grid: int = 800,
        title_text: str = "implicit locus",
        xlabel: str = "x",
        ylabel: str = "y",
        params: Sequence[Param] | None = None,
        title_builder: Callable[[dict], str] | None = None,
    ) -> "ImplicitFunctionView":
        import sympy as sp

        raw = sp.lambdify(symbols, poly, "numpy")

        def func(X, Y, _values):
            return raw(X, Y)

        return cls(
            func,
            box=box,
            grid=grid,
            title_text=title_text,
            title_builder=title_builder,
            xlabel=xlabel,
            ylabel=ylabel,
            params=params,
        )

    @classmethod
    def of(
        cls, problem, box=(-3.0, 3.0), grid: int = 800
    ) -> "ImplicitFunctionView":
        poly, symbols = problem.relation()
        return cls.from_equation(
            poly,
            symbols,
            box=box,
            grid=grid,
            title_builder=lambda _values: f"{problem!s} over R",
        )

    def params(self) -> list[Param]:
        return list(self._params_def)

    def title(self) -> str:
        if self._title_builder is not None:
            return self._title_builder(self.param_values)
        return self._title_text

    def render(self, ax) -> None:
        import numpy as np

        g = np.linspace(*self.box, self.grid)
        X, Y = np.meshgrid(g, g)
        values = self._call_with_optional_values(self._func, X, Y)
        ax.contour(X, Y, values, levels=[self.level], colors="black", linewidths=1.2)
        ax.axhline(0, lw=0.4, c="gray")
        ax.axvline(0, lw=0.4, c="gray")
        ax.set_aspect("equal")

    def _call_with_optional_values(self, func, *args):
        try:
            arity = len(inspect.signature(func).parameters)
        except (TypeError, ValueError):
            arity = len(args) + 1
        if arity <= len(args):
            return func(*args)
        return func(*args, self.param_values)


def plot(
    problem=None,
    F=None,
    poly=None,
    symbols=None,
    box=(-3.0, 3.0),
    path: str = "x1.png",
    grid: int = 800,
):
    """Draw a planar implicit locus and/or finite point set in the (x, y) plane."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels: list[Drawable] = []
    if poly is not None:
        panels.append(
            ImplicitFunctionView.from_equation(poly, symbols, box=box, grid=grid)
        )
    if F is not None:
        panels.append(PointSet2DView.from_moduli_problem(problem, F=F))

    if not panels:
        raise ValueError("nothing to draw: pass poly+symbols and/or problem+F")

    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(5.5 * n, 5), squeeze=False)
    for panel, ax in zip(panels, axes[0]):
        panel.draw(ax)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
