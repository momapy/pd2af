import os.path
import glob
import shutil
import pathlib

import momapy.io.core
import momapy.rendering.skia
import momapy.rendering.core

import pd2af
import pd2af.casq


def remake_dir(path: str | pathlib.Path) -> pathlib.Path:
    path = pathlib.Path(path)
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True)
    return path


def _get_output_dir_path(
    input_file_path, input_main_dir_path, output_main_dir_path, mode
):
    output_dir_path = (
        output_main_dir_path
        / mode
        / input_file_path.parent.relative_to(input_main_dir_path)
    )
    return output_dir_path


def make_pd2af_results():
    INPUT_MAIN_DIR_PATH = pathlib.Path("data/celldesigner/")
    INPUT_FILE_PATHS = list(INPUT_MAIN_DIR_PATH.glob("**/*.xml"))
    OUTPUT_MAIN_DIR_PATH = pathlib.Path("results/pd2af/")
    output_dir_paths = set([])
    MODES = ["pd2af", "pd2af-no-complex"]
    LAYOUT_MODES = ("plain", "overlay", "auto")
    for mode in MODES:
        for input_file_path in INPUT_FILE_PATHS:
            output_dir_path = _get_output_dir_path(
                input_file_path, INPUT_MAIN_DIR_PATH, OUTPUT_MAIN_DIR_PATH, mode
            )
            output_dir_paths.add(output_dir_path)
    output_dir_paths = sorted(list(output_dir_paths))
    for output_dir_path in output_dir_paths:
        remake_dir(output_dir_path)
    for input_file_path in INPUT_FILE_PATHS:
        input_file_name = input_file_path.name
        output_file_name = f"{input_file_name.split('.')[0]}.pdf"
        print("--------------------------------------------------------------------")
        print(f"{input_file_name}")
        print("--------------------------------------------------------------------")
        print("* Reading map...")
        reader_result = momapy.io.core.read(input_file_path)
        map_ = reader_result.obj
        for mode in MODES:
            output_dir_path = _get_output_dir_path(
                input_file_path, INPUT_MAIN_DIR_PATH, OUTPUT_MAIN_DIR_PATH, mode
            )
            output_file_path = os.path.join(output_dir_path, output_file_name)
            print(f"* Transforming map with pd2af in {mode} mode...")
            new_maps = [
                pd2af.transform(map_, mode=mode, layout_mode=layout_mode)
                for layout_mode in LAYOUT_MODES
            ]
            print(f"* Rendering map(s) to {output_file_path}...")
            momapy.rendering.core.render_maps(
                new_maps,
                output_file_path,
                multi_pages=True,
            )


def make_casq_results():
    INPUT_FILE_PATHS = glob.glob("data/celldesigner/**/*.xml")
    OUTPUT_DIR_PATH = "results/casq/"
    LAYOUT_MODES = ("plain", "overlay", "auto")
    for input_file_path in INPUT_FILE_PATHS:
        input_file_name = input_file_path.split("/")[-1]
        output_file_name = f"{input_file_name.split('.')[0]}.pdf"
        output_file_path = os.path.join(OUTPUT_DIR_PATH, output_file_name)
        print("--------------------------------------------------------------------")
        print(f"{input_file_name}")
        print("--------------------------------------------------------------------")
        print("* Transforming map with casq...")
        new_maps = [
            pd2af.casq.transform(input_file_path, layout_mode=layout_mode)
            for layout_mode in LAYOUT_MODES
        ]
        print(f"* Rendering map(s) to {output_file_path}...")
        momapy.rendering.core.render_maps(
            new_maps,
            output_file_path,
            renderer="skia",
            format_="pdf",
            multi_pages=True,
        )


def main():
    make_pd2af_results()
    # make_casq_results()


if __name__ == "__main__":
    main()
