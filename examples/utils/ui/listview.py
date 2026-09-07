from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

from matplotlib.patches import Rectangle
from matplotlib.widgets import Button

__all__ = ["ListColumn", "SelectListView"]


@dataclass(frozen=True)
class ListColumn:
    key: str
    label: str
    getter: Callable[[object], object]
    width: float = 0.2


class SelectListView:
    """A scrollable, clickable list renderer for matplotlib axes."""

    def __init__(
        self,
        *,
        items: Sequence[object] | None = None,
        columns: Sequence[ListColumn] | None = None,
        sorters: dict[str, Callable[[object], object]] | None = None,
        title_text: str = "items",
        empty_text: str = "no rows",
        page_size: int = 8,
        on_select: Callable[[object], None] | None = None,
    ) -> None:
        self.items = list(items or [])
        self.columns = tuple(columns or [])
        self.sorters = dict(sorters or {})
        for column in self.columns:
            self.sorters.setdefault(column.key, column.getter)
        self.title_text = title_text
        self.empty_text = empty_text
        self.page_size = page_size
        self.on_select = on_select
        self.sort_key: str | None = None
        self.sort_reverse = False
        self.filter_predicate: Callable[[object], bool] | None = None
        self.filter_label = "all"
        self.offset = 0
        self._selected_item = None
        self._ax = None
        self._canvas = None
        self._button_cid = None
        self._scroll_cid = None
        self._pager_prev_bounds: tuple[float, float, float, float] | None = None
        self._pager_next_bounds: tuple[float, float, float, float] | None = None
        self._pager_prev_ax = None
        self._pager_next_ax = None
        self._pager_prev_button: Button | None = None
        self._pager_next_button: Button | None = None

    def set_items(self, items: Sequence[object]) -> "SelectListView":
        self.items = list(items)
        visible = self._visible_items()
        if self._selected_item not in visible:
            self._selected_item = visible[0] if visible else None
        self._clamp_offset(len(visible))
        return self

    def set_sort(self, key: str | None, reverse: bool = False) -> "SelectListView":
        self.sort_key = key
        self.sort_reverse = reverse
        visible = self._visible_items()
        if self._selected_item not in visible:
            self._selected_item = visible[0] if visible else None
        self._clamp_offset(len(visible))
        return self

    def set_filter(
        self,
        predicate: Callable[[object], bool] | None,
        label: str = "all",
    ) -> "SelectListView":
        self.filter_predicate = predicate
        self.filter_label = label
        visible = self._visible_items()
        if self._selected_item not in visible:
            self._selected_item = visible[0] if visible else None
        self._clamp_offset(len(visible))
        return self

    def selected_item(self):
        return self._selected_item

    def render(self, ax) -> "SelectListView":
        self._attach(ax)
        ax.set_axis_off()

        visible = self._visible_items()
        self._clamp_offset(len(visible))
        if self._selected_item not in visible:
            self._selected_item = visible[0] if visible else None

        ax.text(
            0.01,
            0.99,
            self.title_text,
            fontsize=10,
            fontweight="bold",
            ha="left",
            va="top",
            transform=ax.transAxes,
        )

        footer = (
            f"{len(visible)} rows | sort={self.sort_key or 'none'} | "
            f"filter={self.filter_label} | wheel/buttons=scroll, click=select"
        )
        ax.text(
            0.01,
            0.02,
            footer,
            fontsize=8,
            family="monospace",
            color="0.35",
            ha="left",
            va="bottom",
            transform=ax.transAxes,
        )

        if not visible:
            ax.text(
                0.5,
                0.5,
                self.empty_text,
                fontsize=11,
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            return self

        header_y = 0.87
        row_top = 0.82
        row_bottom = 0.10
        row_height = (row_top - row_bottom) / max(1, self.page_size)
        self._draw_pager(ax, len(visible))
        self._ensure_pager_widgets(ax, len(visible))
        x = 0.02
        for column in self.columns:
            ax.text(
                x,
                header_y,
                column.label,
                fontsize=9,
                fontweight="bold",
                ha="left",
                va="center",
                transform=ax.transAxes,
            )
            x += column.width

        rows = visible[self.offset : self.offset + self.page_size]
        for row_index, item in enumerate(rows):
            y = row_top - row_index * row_height
            is_selected = item == self._selected_item
            bg = "#dbeafe" if is_selected else ("#f8fafc" if row_index % 2 == 0 else "white")
            ax.add_patch(
                Rectangle(
                    (0.01, y - row_height + 0.003),
                    0.98,
                    row_height - 0.006,
                    facecolor=bg,
                    edgecolor="0.88",
                    linewidth=0.6,
                    transform=ax.transAxes,
                    zorder=0,
                )
            )
            x = 0.02
            for column in self.columns:
                value = column.getter(item)
                ax.text(
                    x,
                    y - 0.5 * row_height,
                    str(value),
                    fontsize=9,
                    family="monospace",
                    ha="left",
                    va="center",
                    transform=ax.transAxes,
                    zorder=1,
                )
                x += column.width

        return self

    def _draw_pager(self, ax, total_rows: int) -> None:
        prev_bounds = (0.82, 0.905, 0.07, 0.055)
        next_bounds = (0.90, 0.905, 0.07, 0.055)
        self._pager_prev_bounds = prev_bounds
        self._pager_next_bounds = next_bounds

        prev_active = self.offset > 0
        next_active = self.offset + self.page_size < total_rows
        for bounds, label, active in (
            (prev_bounds, "Prev", prev_active),
            (next_bounds, "Next", next_active),
        ):
            x, y, w, h = bounds
            ax.add_patch(
                Rectangle(
                    (x, y),
                    w,
                    h,
                    facecolor="#e2e8f0" if active else "#f8fafc",
                    edgecolor="0.75",
                    linewidth=0.8,
                    transform=ax.transAxes,
                    zorder=2,
                )
            )
            ax.text(
                x + 0.5 * w,
                y + 0.5 * h,
                label,
                fontsize=8,
                ha="center",
                va="center",
                color="black" if active else "0.6",
                transform=ax.transAxes,
                zorder=3,
            )

    def _visible_items(self) -> list[object]:
        visible = list(self.items)
        if self.filter_predicate is not None:
            visible = [item for item in visible if self.filter_predicate(item)]
        if self.sort_key is not None:
            sorter = self.sorters[self.sort_key]
            visible = sorted(visible, key=sorter, reverse=self.sort_reverse)
        return visible

    def _attach(self, ax) -> None:
        self._ax = ax
        canvas = ax.figure.canvas
        if self._canvas is canvas:
            return
        if self._canvas is not None:
            if self._button_cid is not None:
                self._canvas.mpl_disconnect(self._button_cid)
            if self._scroll_cid is not None:
                self._canvas.mpl_disconnect(self._scroll_cid)
        self._canvas = canvas
        self._button_cid = canvas.mpl_connect("button_press_event", self._on_click)
        self._scroll_cid = canvas.mpl_connect("scroll_event", self._on_scroll)

    def _clamp_offset(self, size: int) -> None:
        max_offset = max(0, size - self.page_size)
        self.offset = max(0, min(self.offset, max_offset))

    def _on_scroll(self, event) -> None:
        if not self._event_targets_list(event):
            return
        visible = self._visible_items()
        self._clamp_offset(len(visible))
        step = getattr(event, "step", 0)
        if event.button == "up" or step > 0:
            self.offset = max(0, self.offset - 1)
        elif event.button == "down" or step < 0:
            self.offset = min(max(0, len(visible) - self.page_size), self.offset + 1)
        self._rerender()

    def _on_click(self, event) -> None:
        if not self._event_targets_list(event):
            return
        x_ax, y_ax = self._ax.transAxes.inverted().transform((event.x, event.y))
        if self._in_bounds(x_ax, y_ax, self._pager_prev_bounds):
            self.offset = max(0, self.offset - 1)
            self._rerender()
            return
        if self._in_bounds(x_ax, y_ax, self._pager_next_bounds):
            visible = self._visible_items()
            self.offset = min(max(0, len(visible) - self.page_size), self.offset + 1)
            self._rerender()
            return
        if not (0.01 <= x_ax <= 0.99 and 0.10 <= y_ax <= 0.82):
            return

        visible = self._visible_items()
        rows = visible[self.offset : self.offset + self.page_size]
        if not rows:
            return

        row_top = 0.82
        row_bottom = 0.10
        row_height = (row_top - row_bottom) / max(1, self.page_size)
        row_index = int((row_top - y_ax) / row_height)
        if row_index < 0 or row_index >= len(rows):
            return

        self._selected_item = rows[row_index]
        if self.on_select is not None:
            self.on_select(self._selected_item)
        else:
            self._rerender()

    @staticmethod
    def _in_bounds(
        x_ax: float,
        y_ax: float,
        bounds: tuple[float, float, float, float] | None,
    ) -> bool:
        if bounds is None:
            return False
        x, y, w, h = bounds
        return x <= x_ax <= x + w and y <= y_ax <= y + h

    def _ensure_pager_widgets(self, ax, total_rows: int) -> None:
        if self._pager_prev_bounds is None or self._pager_next_bounds is None:
            return

        fig = ax.figure
        if self._pager_prev_ax is None or self._pager_prev_ax.figure is not fig:
            self._remove_pager_widgets()
            self._pager_prev_ax = fig.add_axes([0, 0, 0.01, 0.01], zorder=50)
            self._pager_next_ax = fig.add_axes([0, 0, 0.01, 0.01], zorder=50)
            self._pager_prev_button = Button(self._pager_prev_ax, "Prev")
            self._pager_next_button = Button(self._pager_next_ax, "Next")
            self._pager_prev_button.on_clicked(lambda _event: self._page(-1))
            self._pager_next_button.on_clicked(lambda _event: self._page(1))

        self._position_pager_axis(self._pager_prev_ax, ax, self._pager_prev_bounds)
        self._position_pager_axis(self._pager_next_ax, ax, self._pager_next_bounds)

        prev_active = self.offset > 0
        next_active = self.offset + self.page_size < total_rows
        self._style_pager_button(self._pager_prev_ax, self._pager_prev_button, prev_active)
        self._style_pager_button(self._pager_next_ax, self._pager_next_button, next_active)

    def _page(self, delta: int) -> None:
        visible = self._visible_items()
        self._clamp_offset(len(visible))
        self.offset = max(0, min(max(0, len(visible) - self.page_size), self.offset + delta))
        if self._canvas is not None:
            self._canvas.draw_idle()

    def _position_pager_axis(self, pager_ax, list_ax, bounds) -> None:
        x, y, w, h = bounds
        left_bottom = list_ax.transAxes.transform((x, y))
        right_top = list_ax.transAxes.transform((x + w, y + h))
        inv = list_ax.figure.transFigure.inverted()
        fig_left_bottom = inv.transform(left_bottom)
        fig_right_top = inv.transform(right_top)
        pager_ax.set_position(
            [
                fig_left_bottom[0],
                fig_left_bottom[1],
                fig_right_top[0] - fig_left_bottom[0],
                fig_right_top[1] - fig_left_bottom[1],
            ]
        )

    @staticmethod
    def _style_pager_button(pager_ax, button: Button | None, active: bool) -> None:
        if pager_ax is None or button is None:
            return
        pager_ax.set_facecolor("#e2e8f0" if active else "#f8fafc")
        button.label.set_color("black" if active else "0.6")

    def _remove_pager_widgets(self) -> None:
        for pager_ax in (self._pager_prev_ax, self._pager_next_ax):
            if pager_ax is not None:
                pager_ax.remove()
        self._pager_prev_ax = None
        self._pager_next_ax = None
        self._pager_prev_button = None
        self._pager_next_button = None

    def _event_targets_list(self, event) -> bool:
        if self._ax is None:
            return False
        if getattr(event, "inaxes", None) == self._ax:
            return True
        x = getattr(event, "x", None)
        y = getattr(event, "y", None)
        if x is None or y is None:
            return False
        return self._ax.bbox.contains(x, y)

    def _rerender(self) -> None:
        if self._ax is None:
            return
        self._ax.cla()
        self.render(self._ax)
        if self._canvas is not None:
            self._canvas.draw_idle()