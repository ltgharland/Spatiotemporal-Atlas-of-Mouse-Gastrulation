import argparse
import os

import imageio
import napari
import pandas as pd
from tqdm import tqdm


def _to_tif(data_folder, name):
    pos = name[name.find("Pos"):].rstrip(".csv")
    return os.path.join(data_folder, f"MMStack_{pos}.ome.tif")


def show_annotations(annotation_folder, data_folder, names):
    for name in tqdm(names):
        csv_path = os.path.join(annotation_folder, name)
        annotations = pd.read_csv(csv_path).drop(columns=["index", "axis-1"])

        tif_path = _to_tif(data_folder, name)
        im = imageio.volread(tif_path)[:, -1]

        v = napari.Viewer()
        v.title = name
        v.add_image(im)
        v.add_points(annotations)
        napari.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", required=True)
    parser.add_argument("--names", type=int, nargs="+", default=None)
    args = parser.parse_args()

    input_ = args.input
    timepoint, cycle = input_.split("/")[-2:]
    names = args.names
    if names is None:
        names = os.listdir(input_)
        names.sort()

    data_folder = f"/g/kreshuk/data/marioni/shila/{timepoint}/{cycle}"
    assert os.path.exists(data_folder), data_folder
    show_annotations(input_, data_folder, names)


if __name__ == "__main__":
    main()
