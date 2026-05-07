"""Round-trip sweep across maps × (mode, layout_mode) combos.

Reports per-combo OK/FAIL counts plus per-map failures, prints
incrementally so progress is visible.
"""
import glob
import os
import sys
import tempfile
import warnings

import momapy.io.core
import pd2af


warnings.simplefilter("ignore")

_DEFAULT_MAPS = [
    "tests/maps/celldesigner/SNCA_expression.xml",
    "tests/maps/celldesigner/Apoptosis_pathway.xml",
    "tests/maps/celldesigner/LRRK2_activity.xml",
    "tests/maps/celldesigner/Iron_metabolism.xml",
]

if "--all" in sys.argv:
    MAPS = sorted(glob.glob("tests/maps/celldesigner/*.xml"))
else:
    MAPS = _DEFAULT_MAPS

COMBOS = [
    ("keep-species", "plain"),
    ("keep-species", "overlay"),
    ("keep-species", "auto"),
    ("keep-species-no-complex", "plain"),
    ("keep-species-no-complex", "overlay"),
    ("keep-species-no-complex", "auto"),
    ("no-complex", "auto"),
    ("normal", "auto"),
    ("casq", "plain"),
    ("casq", "overlay"),
    ("casq", "auto"),
]


def run_one(map_obj, mode, layout_mode):
    out = pd2af.transform(map_obj, mode=mode, layout_mode=layout_mode)
    with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as handle:
        path = handle.name
    try:
        momapy.io.core.write(out, path, writer="celldesigner")
        momapy.io.core.read(path)
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass


def main():
    counts = {}
    failures = {}
    for map_path in MAPS:
        name = os.path.basename(map_path)
        try:
            map_obj = momapy.io.core.read(map_path).obj
        except Exception as exc:
            print(f"SKIP read {name}: {type(exc).__name__} {exc}", flush=True)
            continue
        for mode, layout_mode in COMBOS:
            key = (mode, layout_mode)
            counts.setdefault(key, [0, 0])
            try:
                run_one(map_obj, mode, layout_mode)
                counts[key][0] += 1
            except Exception as exc:
                counts[key][1] += 1
                failures.setdefault(key, []).append(
                    (name, type(exc).__name__, str(exc)[:80])
                )
                print(
                    f"FAIL {name:50s} {mode:25s} {str(layout_mode):8s}"
                    f" {type(exc).__name__}: {str(exc)[:80]}",
                    flush=True,
                )
        print(f"done map {name}", flush=True)
    print("\n=== summary ===", flush=True)
    for (mode, layout_mode), (ok, fail) in sorted(counts.items()):
        total = ok + fail
        print(f"  {mode:25s} {str(layout_mode):8s} {ok:3d} / {total:3d} ok")


if __name__ == "__main__":
    main()
