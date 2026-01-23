import os.path

import momapy.io.core
import momapy.rendering.skia
import momapy.rendering.core

import pd2af.core

if __name__ == "__main__":
    INPUT_FILE_PATH = "bug_cd.xml"
    input_file_name = INPUT_FILE_PATH.split("/")[-1]
    output_file_name = f"{input_file_name.split('.')[0]}.pdf"
    print("--------------------------------------------------------------------")
    print(f"{input_file_name}")
    print("--------------------------------------------------------------------")
    print("* Reading map...")
    reader_result = momapy.io.core.read(INPUT_FILE_PATH)
    map_ = reader_result.obj
    for mode in ["casq"]:
        output_file_path = os.path.join(output_file_name)
        print(f"* Transforming map in {mode} mode...")
        maps = pd2af.core.transform_map(map_, mode=mode)
        print(f"* Rendering map to {output_file_path}...")
        momapy.rendering.core.render_maps(
            maps,
            output_file_path,
            renderer="skia",
            format_="pdf",
            multi_pages=True,
        )
