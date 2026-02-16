import os.path

import momapy.io.core
import momapy.rendering.skia
import momapy.rendering.core

import pd2af.casq

if __name__ == "__main__":
    INPUT_FILE_PATH = "phenotype.xml"

    input_file_name = INPUT_FILE_PATH.split("/")[-1]
    output_file_name = f"{input_file_name.split('.')[0]}.pdf"
    print("--------------------------------------------------------------------")
    print(f"{input_file_name}")
    print("--------------------------------------------------------------------")
    print("* Reading map...")
    output_file_path = output_file_name
    print("* Transforming map...")
    maps = pd2af.casq.transform_map(INPUT_FILE_PATH)
    print(f"* Rendering map to {output_file_path}...")
    momapy.rendering.core.render_maps(
        maps,
        output_file_path,
        renderer="skia",
        format_="pdf",
        multi_pages=True,
    )
