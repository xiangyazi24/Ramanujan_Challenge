#!/usr/bin/env python3
"""Apply the audited reporting/ranking corrections before the final rerun."""

from pathlib import Path


MAIN = Path("problems/3.2/research/scripts/q7311_centered_crt_fourier.py")
SHAPE = Path("problems/3.2/research/scripts/q7311_shape_postprocess.py")


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count == 0 and new in text:
        return text
    if count != 1:
        raise RuntimeError(f"{label}: expected one source occurrence, found {count}")
    return text.replace(old, new, 1)


def patch_main() -> None:
    text = MAIN.read_text(encoding="utf-8")

    text = replace_once(
        text,
        '    abs_mass = sum(abs_values, Fraction())\n    signed_primitive = centered_degree[3]\n',
        '    abs_mass = sum(abs_values, Fraction())\n'
        '    ordered_abs_mass = 6 * abs_mass\n'
        '    signed_primitive = centered_degree[3]\n'
        '    if ordered_abs_mass < abs(signed_primitive):\n'
        '        raise AssertionError("ordered primitive absolute mass is below absolute signed mass")\n',
        "ordered primitive mass",
    )
    text = replace_once(
        text,
        '        "primitive_abs_mass": abs_mass,\n',
        '        "primitive_abs_mass": ordered_abs_mass,\n',
        "ordered primitive metric",
    )

    text = replace_once(
        text,
        '    keep = max(4 * top_modes, 40)\n'
        '    heap: List[Tuple[float, int, float, int, int, int]] = []\n',
        '    # Retain every scanned primitive pair; high-precision reranking then\n'
        '    # covers the complete scanned prefix rather than a floating heap.\n'
        '    scanned_candidates: List[Tuple[float, int, float, int, int, int]] = []\n',
        "frequency candidate storage",
    )
    text = replace_once(
        text,
        '            item = (abs(value), k, value, a, b, c)\n'
        '            if len(heap) < keep:\n'
        '                heapq.heappush(heap, item)\n'
        '            elif item[0] > heap[0][0]:\n'
        '                heapq.heapreplace(heap, item)\n',
        '            item = (abs(value), k, value, a, b, c)\n'
        '            scanned_candidates.append(item)\n',
        "frequency candidate append",
    )
    text = replace_once(
        text,
        '        if k % 64 == 0 and len(heap) >= top_modes:\n'
        '            ranked = sorted(heap, reverse=True)\n',
        '        if k % 64 == 0 and len(scanned_candidates) >= top_modes:\n'
        '            ranked = sorted(scanned_candidates, reverse=True)\n',
        "frequency checkpoint",
    )
    text = replace_once(
        text,
        '    candidates = sorted(heap, reverse=True)\n',
        '    candidates = sorted(scanned_candidates, reverse=True)\n',
        "frequency high-precision rerank",
    )

    # Make the ordered factor explicit in every raw-table field name.
    text = text.replace("primitive_signed_total", "ordered_primitive_signed_total")
    text = text.replace("primitive_abs_mass", "ordered_primitive_abs_mass")

    MAIN.write_text(text, encoding="utf-8")


def patch_shape() -> None:
    text = SHAPE.read_text(encoding="utf-8")
    text = replace_once(
        text,
        '(X, "primitive_abs_mass")',
        '(X, "ordered_primitive_abs_mass")',
        "shape random metric lookup",
    )
    SHAPE.write_text(text, encoding="utf-8")


def main() -> None:
    patch_main()
    patch_shape()


if __name__ == "__main__":
    main()
