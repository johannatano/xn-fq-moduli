from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence

__all__ = [
	"Drawable",
	"Param",
]


class Param:
	"""A tunable parameter: its bounds, current value, and UI hints."""

	def __init__(
		self,
		name: str,
		min_val: float,
		max_val: float,
		default_val,
		step: float = 0.01,
		label: str | None = None,
		enable_ui: bool = True,
		choices: Sequence | None = None,
	) -> None:
		self.name = name
		self.choices = tuple(choices) if choices is not None else None
		if self.choices is not None and not self.choices:
			raise ValueError(f"Param {name!r} choices must not be empty")
		if self.choices is not None and default_val not in self.choices:
			raise ValueError(
				f"Param {name!r} default {default_val!r} must be in choices"
			)
		self.min_val = 0 if self.choices is not None else min_val
		self.max_val = (
			len(self.choices) - 1 if self.choices is not None else max_val
		)
		self.default_val = default_val
		self.value = default_val
		self.step = step
		self.label = label or name
		self.type = type(default_val)
		self.enable_ui = enable_ui

	def set(self, value):
		"""Coerce to the declared type and clamp into range."""
		if self.choices is not None:
			if value not in self.choices:
				raise ValueError(
					f"invalid value {value!r} for {self.name!r}; "
					f"expected one of {list(self.choices)!r}"
				)
			self.value = value
			return self.value
		value = self.type(value)
		if isinstance(value, (int, float)) and not isinstance(value, bool):
			value = max(self.min_val, min(self.max_val, value))
		self.value = value
		return self.value

	def set_from_ui(self, value):
		if self.choices is None:
			return self.set(value)
		index = int(round(value))
		index = max(0, min(len(self.choices) - 1, index))
		self.value = self.choices[index]
		return self.value

	@property
	def slider_min(self):
		return 0 if self.choices is not None else self.min_val

	@property
	def slider_max(self):
		return len(self.choices) - 1 if self.choices is not None else self.max_val

	@property
	def slider_step(self):
		if self.choices is not None:
			return 1
		return 1 if self.type is int else self.step

	@property
	def slider_value(self):
		if self.choices is None:
			return self.value
		return self.choices.index(self.value)

	def format_value(self, value=None) -> str:
		value = self.value if value is None else value
		return str(value)

	def reset(self):
		self.value = self.default_val
		return self.value

	def __repr__(self) -> str:
		return (
			f"Param({self.name}={self.value!r}, range=[{self.min_val}, {self.max_val}])"
		)


class Drawable(ABC):
	"""Anything that can render itself onto a matplotlib axes."""

	xlabel: str = "x"
	ylabel: str = "y"

	def __init__(self) -> None:
		self._params = {p.name: p for p in self.params()}
		self._on_change = None

	def params(self) -> list[object]:
		return []

	@abstractmethod
	def render(self, ax) -> None:
		"""Draw the current state onto `ax`."""

	def title(self) -> str:
		return type(self).__name__

	@property
	def param_values(self) -> dict:
		return {name: p.value for name, p in self._params.items()}

	def get_param(self, key: str):
		return self._params[key].value

	def set_param(self, key: str, value, redraw: bool = True):
		if key not in self._params:
			raise KeyError(f"unknown parameter {key!r}; have {sorted(self._params)}")
		self._params[key].set(value)
		if redraw and self._on_change is not None:
			self._on_change()
		return self._params[key].value

	def set_params(self, **kwargs):
		for key, value in kwargs.items():
			self.set_param(key, value, redraw=False)
		if self._on_change is not None:
			self._on_change()
		return self

	def reset_params(self):
		for p in self._params.values():
			p.reset()
		if self._on_change is not None:
			self._on_change()
		return self

	def figure(self, figsize=(6, 6), with_sliders: bool = False):
		import matplotlib.pyplot as plt

		fig = plt.figure(figsize=figsize)
		if with_sliders and self._slider_params():
			ax = fig.add_axes([0.10, 0.10, 0.60, 0.80])
		else:
			ax = fig.add_subplot(111)
		return fig, ax

	def draw(self, ax) -> None:
		ax.cla()
		ax.set_title(self.title())
		ax.set_xlabel(self.xlabel)
		ax.set_ylabel(self.ylabel)
		self.render(ax)

	def plot(
		self,
		path: str | None = None,
		figsize=(6, 6),
		dpi: int = 150,
		show: bool = False,
	):
		import matplotlib

		if path is not None and not show:
			matplotlib.use("Agg")
		import matplotlib.pyplot as plt

		fig, ax = self.figure(figsize=figsize)
		self.draw(ax)
		fig.tight_layout()

		if path is not None:
			fig.savefig(path, dpi=dpi)
		if show:
			plt.show()
		elif path is not None:
			plt.close(fig)
			return path
		return fig

	def interactive(self, figsize=(11, 6)):
		import matplotlib.pyplot as plt
		from matplotlib.widgets import Slider

		fig, ax = self.figure(figsize=figsize, with_sliders=True)

		sliders = {}
		tunable = self._slider_params()
		for i, p in enumerate(tunable):
			height, spacing = 0.04, 0.02
			bottom = 0.85 - i * (height + spacing)
			ax_slider = fig.add_axes([0.78, bottom, 0.18, height])
			sliders[p.name] = Slider(
				ax_slider,
				p.label,
				p.slider_min,
				p.slider_max,
				valinit=p.slider_value,
				valstep=p.slider_step,
			)
			sliders[p.name].valtext.set_text(p.format_value())

		def rerender():
			self.draw(ax)
			fig.canvas.draw_idle()

		def make_callback(name):
			def _update(value):
				self._params[name].set_from_ui(value)
				sliders[name].valtext.set_text(self._params[name].format_value())
				rerender()

			return _update

		for p in tunable:
			sliders[p.name].on_changed(make_callback(p.name))

		fig._nt_sliders = sliders  # type: ignore[attr-defined]

		self._on_change = rerender
		try:
			rerender()
			plt.show()
		finally:
			self._on_change = None
		return fig

	def _slider_params(self) -> list[object]:
		return [
			p for p in self._params.values() if p.enable_ui and p.type in (int, float)
		]
