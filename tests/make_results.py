import os.path
import glob

import momapy.io.core
import momapy.rendering.skia
import momapy.rendering.core

import pd2af.core

if __name__ == "__main__":
    INPUT_FILE_PATHS = glob.glob("data/**/*.pickle")
    OUTPUT_DIR_PATH = "results/"
    OUTPUT_PD2AF_DIR_PATH = "results/pd2af/"
    LAYOUT_MODE = "all"
    for input_file_path in INPUT_FILE_PATHS:
        input_file_name = input_file_path.split("/")[-1]
        output_file_name = f"{input_file_name.split('.')[0]}.pdf"
        print("--------------------------------------------------------------------")
        print(f"{input_file_name}")
        print("--------------------------------------------------------------------")
        print("* Reading map...")
        reader_result = momapy.io.core.read(input_file_path)
        map_ = reader_result.obj
        for mode in ["pd2af", "casq"]:
            output_file_path = os.path.join(OUTPUT_DIR_PATH, mode, output_file_name)
            print(f"* Transforming map in {mode} mode...")
            new_maps = pd2af.core.transform_map(
                map_, mode=mode, layout_mode=LAYOUT_MODE
            )
            print(f"* Rendering map(s) to {output_file_path}...")
            momapy.rendering.core.render_maps(
                new_maps,
                output_file_path,
                renderer="skia",
                format_="pdf",
                multi_pages=True,
            )
