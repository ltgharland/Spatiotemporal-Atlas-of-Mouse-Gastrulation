import os
import argparse

import imageio
import napari
import pandas as pd


# TODO size filter
def check_corrections(name, version, save, min_size=25):
    im_name = name.replace(".tif", ".ome.tif")
    im_path = os.path.join("./data/nuclei", im_name)
    assert os.path.exists(im_path), im_path
    image = imageio.imread(im_path)

    seg = imageio.volread(f"./data/corrected/v{version}/{name}")
    pos = name[name.find("Pos"):name.find("Pos")+4]
    points = pd.read_csv(f"./data/annotations/TimEmbryos-020420_points_{pos}.csv").values[:, 1:]
    cz = tuple(points[0, :2].astype("int"))
    seg = seg[cz].astype("int16")
    assert seg.shape == image.shape
    points = points[:, 2:]

    v = napari.Viewer()
    v.add_image(image)
    v.add_labels(seg)
    v.add_points(points)
    napari.run()

    if save:
        out_folder = f"./data/ground_truth/v{version}"
        os.makedirs(out_folder, exist_ok=True)
        out_path = os.path.join(out_folder, im_name)
        imageio.imsave(out_path, seg)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", required=True, type=int)
    parser.add_argument("-n", "--names", default=None, type=str, nargs="+")
    parser.add_argument("-s", "--save", type=int, default=0)

    args = parser.parse_args()
    if args.names is None:
        names = os.listdir(f"./data/corrected/v{args.version}")
    else:
        names = args.names
    for name in names:
        check_corrections(name, args.version, args.save)


if __name__ == "__main__":
    main()
