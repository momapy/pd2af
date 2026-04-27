import os.path

import momapy.io.core
import momapy.rendering.skia
import momapy.rendering.core

import pd2af

if __name__ == "__main__":
    # INPUT_FILE_PATH = "example.xml"
    INPUT_FILE_PATH = "./data/celldesigner/pd_dm/Histamine_signaling.xml"
    # MODE = "pd2af-no-complex"
    MODE = "pd2af"
    LAYOUT_MODE = "overlay"
    input_file_name = INPUT_FILE_PATH.split("/")[-1]
    output_file_name = f"{input_file_name.split('.')[0]}_{MODE}_mode.pdf"
    print("--------------------------------------------------------------------")
    print(f"{input_file_name}")
    print("--------------------------------------------------------------------")
    print("* Reading map...")
    reader_result = momapy.io.core.read(INPUT_FILE_PATH)
    map_ = reader_result.obj
    output_file_path = os.path.join(output_file_name)
    print(f"* Transforming map in {MODE} mode...")
    new_map = pd2af.transform(
        map_, mode=MODE, layout_mode=LAYOUT_MODE
    )
    print(f"* Rendering map to {output_file_path}...")
    momapy.rendering.core.render_maps(
        [new_map],
        output_file_path,
        renderer="skia",
        format_="pdf",
        multi_pages=True,
    )
