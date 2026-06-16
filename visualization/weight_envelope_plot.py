from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from PIL import Image, ImageDraw, ImageFont


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = PROJECT_ROOT / "outputs" / "cases" / "ch3_weight_envelope_output.json"
OUTPUT_PATH = PROJECT_ROOT / "outputs" / "figures" / "ch3_weight_envelope.png"


def read_json(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8-sig") as f:
        return json.load(f)


def get_font(size: int, bold: bool = False) -> ImageFont.ImageFont:
    candidates = [
        "consolab.ttf" if bold else "consola.ttf",
        "arialbd.ttf" if bold else "arial.ttf",
        "DejaVuSansMono-Bold.ttf" if bold else "DejaVuSansMono.ttf",
        "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def as_xy(point: dict[str, Any]) -> tuple[float, float]:
    if "x" in point and "y" in point:
        return float(point["x"]), float(point["y"])
    return float(point["XBAR"]), float(point["WEIGHT"])


def collect_points(plot_data: dict[str, Any]) -> list[tuple[float, float]]:
    points: list[tuple[float, float]] = []
    for key in ["basic_points", "forward_edge", "aft_edge", "envelope_polygon"]:
        points.extend(as_xy(point) for point in plot_data.get(key, []))

    structural = plot_data.get("structural_limits") or {}
    for line in structural.get("plot_lines", []):
        points.extend(as_xy(point) for point in line.get("points", []))
    for point in structural.get("points", {}).values():
        points.append(as_xy(point))
    return points


def format_axis(value: float) -> str:
    if abs(value) >= 1000:
        return f"{value / 1000:.1f}E+3"
    if abs(value) >= 10:
        return f"{value / 10:.1f}E+1"
    return f"{value:.1f}"


def draw_polyline(
    draw: ImageDraw.ImageDraw,
    points: list[dict[str, Any]],
    mapper,
    fill: tuple[int, int, int],
    width: int = 3,
) -> None:
    if len(points) < 2:
        return
    draw.line([mapper(*as_xy(point)) for point in points], fill=fill, width=width, joint="curve")


def draw_labeled_points(
    draw: ImageDraw.ImageDraw,
    points: list[dict[str, Any]],
    mapper,
    font: ImageFont.ImageFont,
    fill: tuple[int, int, int],
    label_offsets: dict[str, tuple[int, int]] | None = None,
) -> None:
    label_offsets = label_offsets or {}
    for point in points:
        x, y = mapper(*as_xy(point))
        label = str(point.get("label", point.get("ADDED", "")))
        if "point_id" in point:
            label = f"[{point['point_id']}] {label}"
        elif "POINT_ID" in point:
            label = f"[{point['POINT_ID']}] {label}"
        draw.ellipse((x - 4, y - 4, x + 4, y + 4), fill=fill)
        offset_key = str(point.get("label", point.get("ADDED", ""))).upper()
        dx, dy = label_offsets.get(offset_key, (8, -15))
        draw.text((x + dx, y + dy), label, fill=fill, font=font)


def split_shared_endpoint_labels(
    forward_labels: list[dict[str, Any]],
    aft_labels: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    if not forward_labels or not aft_labels:
        return forward_labels, aft_labels, []

    forward_last = forward_labels[-1]
    aft_last = aft_labels[-1]
    fx, fy = as_xy(forward_last)
    ax, ay = as_xy(aft_last)
    same_coordinate = abs(fx - ax) < 1.0e-6 and abs(fy - ay) < 1.0e-6
    same_point_id = forward_last.get("point_id") == aft_last.get("point_id")
    if not same_coordinate or not same_point_id:
        return forward_labels, aft_labels, []

    forward_label = str(forward_last.get("label", ""))
    aft_label = str(aft_last.get("label", ""))
    if forward_label == aft_label:
        merged_label = forward_label
    else:
        merged_label = f"{forward_label} / {aft_label}"

    merged = dict(forward_last)
    merged["label"] = merged_label
    return forward_labels[:-1], aft_labels[:-1], [merged]


def plot_weight_envelope(
    input_path: str | Path = INPUT_PATH,
    output_path: str | Path = OUTPUT_PATH,
) -> Path:
    data = read_json(input_path)
    plot_data = data["envelope"]["plot_data"]
    structural = plot_data.get("structural_limits") or {}

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    width, height = 960, 710
    left, right, top, bottom = 110, 48, 58, 92
    plot_width = width - left - right
    plot_height = height - top - bottom
    axis_bottom = top + plot_height

    all_points = collect_points(plot_data)
    min_x = math.floor(min(x for x, _ in all_points) / 5.0) * 5.0
    max_x = math.ceil(max(x for x, _ in all_points) / 5.0) * 5.0
    min_y = math.floor(min(y for _, y in all_points) / 200.0) * 200.0
    max_y = math.ceil(max(y for _, y in all_points) / 200.0) * 200.0 + 200.0

    def sx(x_value: float) -> float:
        return left + (x_value - min_x) / (max_x - min_x) * plot_width

    def sy(y_value: float) -> float:
        return top + (max_y - y_value) / (max_y - min_y) * plot_height

    def mapper(x_value: float, y_value: float) -> tuple[float, float]:
        return sx(x_value), sy(y_value)

    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    font = get_font(15, bold=True)
    small_font = get_font(13, bold=True)
    axis_font = get_font(14, bold=False)
    title_font = get_font(18, bold=True)

    for x in range(0, width, 10):
        for y in range(0, height, 10):
            draw.point((x, y), fill=(210, 210, 210))

    draw.line((left, top, left, axis_bottom), fill=(0, 0, 0), width=2)
    draw.line((left, axis_bottom, left + plot_width, axis_bottom), fill=(0, 0, 0), width=2)

    for index in range(6):
        x_value = min_x + (max_x - min_x) * index / 5
        x = sx(x_value)
        draw.line((x, axis_bottom - 7, x, axis_bottom + 7), fill=(0, 0, 0), width=2)
        draw.text((x - 26, axis_bottom + 16), format_axis(x_value), fill=(0, 0, 0), font=axis_font)

    y_step = 200.0
    tick = min_y
    while tick <= max_y + 0.1:
        y = sy(tick)
        draw.line((left - 10, y, left, y), fill=(0, 0, 0), width=2)
        draw.text((left - 86, y - 9), format_axis(tick), fill=(0, 0, 0), font=axis_font)
        tick += y_step

    title = plot_data.get("title", "USEFUL LOAD ENVELOPE AND STRUCTURAL LIMITS")
    title_width = draw.textlength(title, font=title_font)
    draw.text(((width - title_width) / 2, 18), title, fill=(0, 0, 0), font=title_font)

    x_axis = plot_data.get("x_axis", "Fuselage Station")
    x_axis_width = draw.textlength(x_axis, font=axis_font)
    draw.text((left + (plot_width - x_axis_width) / 2, height - 36), x_axis, fill=(0, 0, 0), font=axis_font)
    y_axis = plot_data.get("y_axis", "Weight")
    y_label = Image.new("RGBA", (120, 28), (255, 255, 255, 0))
    y_label_draw = ImageDraw.Draw(y_label)
    y_label_draw.text((0, 0), y_axis, fill=(0, 0, 0), font=axis_font)
    y_label = y_label.rotate(90, expand=True)
    image.paste(y_label, (8, int(top + plot_height / 2 - y_label.height / 2)), y_label)

    envelope = plot_data.get("envelope_polygon", [])
    if envelope:
        closed_envelope = envelope + [envelope[0]]
        draw_polyline(draw, closed_envelope, mapper, fill=(0, 0, 0), width=3)

    draw_polyline(draw, plot_data.get("forward_edge", []), mapper, fill=(0, 0, 0), width=3)
    draw_polyline(draw, plot_data.get("aft_edge", []), mapper, fill=(0, 0, 0), width=3)
    draw_polyline(draw, plot_data.get("basic_points", []), mapper, fill=(0, 0, 0), width=3)

    if structural:
        points = structural.get("points", {})
        fwd_reduced = points.get("FWDRED")
        fwd_gross = points.get("FWDGROSS")
        aft_gross = points.get("AFTGROSS")
        boundary: list[dict[str, Any]] = []
        if fwd_reduced and fwd_gross and aft_gross:
            boundary = [
                {"x": fwd_reduced["XBAR"], "y": min_y, "label": "Fwd reduced lower"},
                {"x": fwd_reduced["XBAR"], "y": fwd_reduced["WEIGHT"], "label": "Fwd Reduced"},
                {"x": fwd_gross["XBAR"], "y": fwd_gross["WEIGHT"], "label": "Fwd Gross"},
                {"x": aft_gross["XBAR"], "y": aft_gross["WEIGHT"], "label": "Aft Gross"},
                {"x": aft_gross["XBAR"], "y": min_y, "label": "Aft gross lower"},
            ]
            draw_polyline(draw, boundary, mapper, fill=(0, 0, 0), width=3)

    offsets = {
        "EMPTY WEIGHT": (8, -2),
        "PILOT": (8, -18),
        "MINIMUM WEIGHT": (8, -14),
        "FUEL TO FULL": (8, -18),
        "BAGGAGE": (8, -16),
        "COPILOT": (8, -8),
        "3RD PERSON": (8, -16),
        "4TH PERSON": (8, 2),
        "5TH PERSON": (8, -16),
        "6TH PERSON": (8, -16),
    }
    forward_labels = plot_data.get("forward_edge", [])[1:]
    aft_labels = plot_data.get("aft_edge", [])[1:]
    forward_labels, aft_labels, shared_endpoint_labels = split_shared_endpoint_labels(
        forward_labels,
        aft_labels,
    )

    draw_labeled_points(draw, plot_data.get("basic_points", []), mapper, small_font, (0, 0, 0), offsets)
    draw_labeled_points(draw, forward_labels, mapper, small_font, (0, 0, 0), offsets)
    draw_labeled_points(draw, aft_labels, mapper, small_font, (0, 0, 0), offsets)
    draw_labeled_points(draw, shared_endpoint_labels, mapper, small_font, (0, 0, 0), offsets)

    if structural:
        points = structural.get("points", {})
        aft_gross = points.get("AFTGROSS")
        x0 = sx(float(aft_gross["XBAR"])) + 18 if aft_gross else left + plot_width * 0.62
        y0 = sy(2260.0)
        draw.text((x0, y0), "LIMITS", fill=(0, 0, 0), font=small_font)
        rows = [
            ("Fwd Gross", points.get("FWDGROSS")),
            ("Aft Gross", points.get("AFTGROSS")),
            ("Fwd Reduced", points.get("FWDRED")),
        ]
        for idx, (name, point) in enumerate(rows, start=1):
            if not point:
                continue
            draw.text(
                (x0, y0 + idx * 18),
                f"{name:<12} {point['WEIGHT']:.0f}  {point['XBAR']:.1f}",
                fill=(0, 0, 0),
                font=small_font,
            )

    image.save(output_path)
    return output_path


if __name__ == "__main__":
    saved_path = plot_weight_envelope()
    print(f"Saved figure: {saved_path}")
