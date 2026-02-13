import os.path

import momapy.io.core
import momapy.rendering.skia
import momapy.rendering.core

import pd2af.core

if __name__ == "__main__":
    INPUT_FILE_PATH = "example.xml"
    MODE = "pd2af"
    ACTIVE = []
    LAYOUT_MODE = "auto"
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
    maps = pd2af.core.transform_map(
        map_, mode=MODE, active=ACTIVE, layout_mode=LAYOUT_MODE
    )
    print(f"* Rendering map to {output_file_path}...")
    momapy.rendering.core.render_maps(
        maps,
        output_file_path,
        renderer="skia",
        format_="pdf",
        multi_pages=True,
    )
