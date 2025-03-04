import argparse
import os

import imageio
import napari
import pandas as pd
from tqdm import tqdm


def _to_tif(data_folder, name):
    pos = name[name.find("Pos"):].rstrip(".csv")
    return os.path.join(data_folder, f"MMStack_{pos}.ome.tif")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-e", "--use_enhancer", type=int, default=0)
    args = parser.parse_args()
    use_enhancer = bool(args.use_enhancer)

    annotation_folder = args.input
    assert os.path.exists(annotation_folder), annotation_folder
    timepoint, cycle = args.input.split("/")[-2:]

    data_folder = f"/g/kreshuk/data/marioni/shila/{timepoint}/{cycle}"
    assert os.path.exists(data_folder), data_folder

    seg_folder_nuclei = f"/g/kreshuk/data/marioni/shila/nucleus_segmentation/{timepoint}/{cycle}/watershed"
    assert os.path.exists(seg_folder_nuclei), seg_folder_nuclei

    seg_name = "enhancer" if use_enhancer else "vanilla"
    seg_folder_cells = f"/g/kreshuk/data/marioni/shila/cell_segmentation/{timepoint}/{cycle}/{seg_name}"
    have_cell_seg = os.path.exists(seg_folder_cells)

    names = os.listdir(annotation_folder)
    names.sort()
    for name in tqdm(names):
        annotation_path = os.path.join(annotation_folder, name)
        annotations = pd.read_csv(annotation_path).drop(columns=["index", "axis-1"])

        tif_path = _to_tif(data_folder, name)
        image = imageio.volread(tif_path)
        mem_image = image[:, 0]
        nuc_image = image[:, -1]

        nuc_path = _to_tif(seg_folder_nuclei, name)
        nuclei = imageio.volread(nuc_path)
        assert nuc_image.shape == nuclei.shape

        if have_cell_seg:
            cell_path = _to_tif(seg_folder_cells, name)
            assert os.path.exists(cell_path)
            cells = imageio.volread(cell_path)
            assert cells.shape == nuc_image.shape

        v = napari.Viewer()
        v.title = name
        v.add_image(mem_image)
        v.add_image(nuc_image)
        v.add_labels(nuclei)
        if have_cell_seg:
            v.add_labels(cells)
        v.add_points(annotations, size=30)
        napari.run()


if __name__ == "__main__":
    main()
