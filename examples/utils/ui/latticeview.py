from __future__ import annotations

from dataclasses import dataclass
from math import hypot

from .base import Drawable

__all__ = [
	"DisplayBasis2D",
	"LatticePoint",
	"LatticeView",
	"Sublattice",
]


@dataclass(frozen=True)
class DisplayBasis2D:
	"""Two display vectors that realize a rank-2 lattice in the plane."""

	e1: tuple[float, float] = (1.0, 0.0)
	e2: tuple[float, float] = (0.0, 1.0)
	labels: tuple[str, str] = ("e1", "e2")

	def map(self, coords: tuple[int | float, int | float]) -> tuple[float, float]:
		a, b = coords
		return (
			a * self.e1[0] + b * self.e2[0],
			a * self.e1[1] + b * self.e2[1],
		)


@dataclass(frozen=True)
class Sublattice:
	"""A rank-2 sublattice given by two generators in ambient coordinates."""

	generators: tuple[tuple[int, int], tuple[int, int]]
	label: str | None = None
	color: str = "C0"
	linewidth: float = 1.2
	alpha: float = 0.9
	linestyle: str = "-"
	zorder: int = 2


@dataclass(frozen=True)
class LatticePoint:
	"""A marked lattice point in ambient coordinates."""

	coords: tuple[int | float, int | float]
	label: str | None = None
	color: str = "crimson"
	size: float = 90.0
	marker: str = "o"
	zorder: int = 9
	annotate_offset: tuple[float, float] = (10.0, 8.0)


class LatticeView(Drawable):
	"""A generic planar view of a rank-2 lattice with sublattice overlays."""

	xlabel = "x"
	ylabel = "y"

	def __init__(
		self,
		display_basis: DisplayBasis2D | None = None,
		radius: float = 8.0,
		title_text: str = "lattice view",
		metadata: list[str] | None = None,
		show_ambient: bool = True,
		show_basis: bool = True,
		show_axes: bool = True,
		max_ambient_lines: int = 120,
	) -> None:
		self.display_basis = display_basis or DisplayBasis2D()
		self.radius = radius
		self.title_text = title_text
		self.metadata = list(metadata or [])
		self.show_ambient = show_ambient
		self.show_basis = show_basis
		self.show_axes = show_axes
		self.max_ambient_lines = max_ambient_lines
		self.sublattices: list[Sublattice] = []
		self.points: list[LatticePoint] = []
		super().__init__()

	def title(self) -> str:
		if not self.metadata:
			return self.title_text
		return self.title_text + "\n" + "\n".join(self.metadata)

	def add_sublattice(
		self,
		generators: tuple[tuple[int, int], tuple[int, int]],
		label: str | None = None,
		color: str = "C0",
		linewidth: float = 1.2,
		alpha: float = 0.9,
		linestyle: str = "-",
		zorder: int = 2,
	) -> "LatticeView":
		self.sublattices.append(
			Sublattice(generators, label, color, linewidth, alpha, linestyle, zorder)
		)
		return self

	def add_point(
		self,
		coords: tuple[int | float, int | float],
		label: str | None = None,
		color: str = "crimson",
		size: float = 90.0,
		marker: str = "o",
		zorder: int = 9,
		annotate_offset: tuple[float, float] = (10.0, 8.0),
	) -> "LatticeView":
		self.points.append(
			LatticePoint(coords, label, color, size, marker, zorder, annotate_offset)
		)
		return self

	def set_metadata(self, *lines: str) -> "LatticeView":
		self.metadata = list(lines)
		return self

	def render(self, ax) -> None:
		if self.show_ambient and self._line_count(((1, 0), (0, 1))) <= self.max_ambient_lines:
			self._draw_grid(
				ax,
				((1, 0), (0, 1)),
				color="0.88",
				linewidth=0.5,
				alpha=1.0,
				linestyle="-",
				zorder=0,
			)

		for sublattice in self.sublattices:
			self._draw_grid(
				ax,
				sublattice.generators,
				color=sublattice.color,
				linewidth=sublattice.linewidth,
				alpha=sublattice.alpha,
				linestyle=sublattice.linestyle,
				zorder=sublattice.zorder,
				label=sublattice.label,
			)

		for point in self.points:
			x, y = self.display_basis.map(point.coords)
			ax.scatter(
				[x],
				[y],
				s=point.size,
				marker=point.marker,
				color=point.color,
				zorder=point.zorder,
				edgecolors="k",
				linewidths=0.7,
				label=point.label,
			)
			if point.label:
				ax.annotate(
					point.label,
					(x, y),
					textcoords="offset points",
					xytext=point.annotate_offset,
					color=point.color,
					fontsize=11,
				)

		if self.show_axes:
			ax.axhline(0, lw=0.7, c="steelblue", zorder=1)
			ax.axvline(0, lw=0.7, c="steelblue", zorder=1)

		if self.show_basis:
			self._draw_basis(ax)

		ax.set_aspect("equal")
		ax.set_xlim(-self.radius, self.radius)
		ax.set_ylim(-self.radius, self.radius)
		handles, labels = ax.get_legend_handles_labels()
		if labels:
			ax.legend(fontsize=8, loc="upper left", framealpha=0.92)

	def _draw_basis(self, ax) -> None:
		for coords, label in zip(((1, 0), (0, 1)), self.display_basis.labels):
			x, y = self.display_basis.map(coords)
			ax.annotate(
				"",
				xy=(x, y),
				xytext=(0, 0),
				arrowprops={"arrowstyle": "->", "lw": 1.0, "color": "black"},
				zorder=6,
			)
			ax.text(x, y, label, fontsize=10, ha="left", va="bottom")

	def _draw_grid(
		self,
		ax,
		generators: tuple[tuple[int, int], tuple[int, int]],
		*,
		color: str,
		linewidth: float,
		alpha: float,
		linestyle: str,
		zorder: int,
		label: str | None = None,
	) -> None:
		u = self.display_basis.map(generators[0])
		v = self.display_basis.map(generators[1])
		span = 2 * self.radius + 2
		ku = int(span / max(self._norm(v), 1e-9)) + 2
		kv = int(span / max(self._norm(u), 1e-9)) + 2

		first = True
		for k in range(-ku, ku + 1):
			base = (k * v[0], k * v[1])
			start = (base[0] - span * u[0], base[1] - span * u[1])
			end = (base[0] + span * u[0], base[1] + span * u[1])
			ax.plot(
				[start[0], end[0]],
				[start[1], end[1]],
				color=color,
				lw=linewidth,
				alpha=alpha,
				ls=linestyle,
				zorder=zorder,
				label=label if first else None,
			)
			first = False

		for k in range(-kv, kv + 1):
			base = (k * u[0], k * u[1])
			start = (base[0] - span * v[0], base[1] - span * v[1])
			end = (base[0] + span * v[0], base[1] + span * v[1])
			ax.plot(
				[start[0], end[0]],
				[start[1], end[1]],
				color=color,
				lw=linewidth,
				alpha=alpha,
				ls=linestyle,
				zorder=zorder,
			)

	def _line_count(
		self, generators: tuple[tuple[int, int], tuple[int, int]]
	) -> int:
		u = self.display_basis.map(generators[0])
		v = self.display_basis.map(generators[1])
		span = 2 * self.radius + 2
		ku = int(span / max(self._norm(v), 1e-9)) + 2
		kv = int(span / max(self._norm(u), 1e-9)) + 2
		return (2 * ku + 1) + (2 * kv + 1)

	@staticmethod
	def _norm(vec: tuple[float, float]) -> float:
		return hypot(vec[0], vec[1])
