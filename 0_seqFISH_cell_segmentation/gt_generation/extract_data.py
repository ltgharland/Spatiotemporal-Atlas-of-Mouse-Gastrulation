import argparse
import os
from subprocess import run

import imageio
import pandas as pd


def _to_tif(name):
    pos = name[name.find("Pos"):name.find("Pos")+4]
    return f"MMStack_{pos}.ome.tif"


def copy_data(names):
    out_folder = "./data/raw"
    os.makedirs(out_folder, exist_ok=True)
    for name in names:
        tif_name = _to_tif(name)
        input_path = f"/g/kreshuk/data/marioni/shila/TimEmbryos-020420/HybCycle_29/{tif_name}"
        outp = os.path.join(out_folder, tif_name)
        inp = f"pape@gpu6.cluster.embl.de:{input_path}"
        cmd = ["scp", inp, outp]
        # print(cmd)
        run(cmd)


def extract_data(names):
    annotation_folder = "./data/annotations"
    raw_folder = "./data/raw"
    data_folder = "./data/nuclei"
    os.makedirs(data_folder, exist_ok=True)

    for name in names:
        annotations = pd.read_csv(os.path.join(annotation_folder, name))
        cz_coords = tuple(annotations.values[:, 1:3][0].astype("int"))
        tif_name = _to_tif(name)
        im_path = os.path.join(raw_folder, tif_name)
        im = imageio.volread(im_path)
        im = im[cz_coords]
        out_path = os.path.join(data_folder, tif_name)
        imageio.imsave(out_path, im)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--copy", type=int, default=0)
    names = os.listdir("./data/annotations")
    args = parser.parse_args()
    if args.copy:
        copy_data(names)
    extract_data(names)


if __name__ == "__main__":
    main()
