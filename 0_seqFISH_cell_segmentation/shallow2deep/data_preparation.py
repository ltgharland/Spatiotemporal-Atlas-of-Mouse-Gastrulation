import argparse
import os
from glob import glob

import imageio
import h5py
import napari
import numpy as np
from tifffile import TiffWriter
from tqdm import tqdm


def fix_data(vol):
    shape = vol.shape
    x = vol.reshape((shape[0] * shape[1], shape[2], shape[3]))
    y = np.zeros_like(vol)

    nc, nz = shape[1], shape[0]
    for c in range(nc):
        start = c * nz
        stop = start + nz
        y[:, c] = x[start:stop]

    return y


def load_seg(seg_files):
    seg = []
    for segf in seg_files:
        with h5py.File(segf, "r") as f:
            seg.append(f["exported_data"][:])
    seg = np.concatenate(seg)
    return seg


def view(image_path, cell_files, nuc_files):
    image = imageio.volread(image_path)
    cell = load_seg(cell_files)
    nuc = load_seg(nuc_files)
    name = os.path.basename(image_path)
    # the usual data screw-up
    if name != "MMStack_Pos0.ome.tif":
        print("Fix image:", name)
        image = fix_data(image)
    v = napari.Viewer()
    v.add_image(image)
    v.add_labels(cell[:, None], scale=(1, 1, 4, 4))
    v.add_labels(nuc[:, None], scale=(1, 1, 4, 4))
    v.title = name
    napari.run()


def save(image_path, image_out, cell_files, cell_out, nuc_files, nuc_out):
    image = imageio.volread(image_path)
    cell = load_seg(cell_files)
    nuc = load_seg(nuc_files)
    name = os.path.basename(image_path)
    # the usual data screw-up
    if name != "MMStack_Pos0.ome.tif":
        print("Fix image:", name)
        image = fix_data(image)
    with TiffWriter(image_out) as tif:
        tif.save(image)
    with TiffWriter(cell_out) as tif:
        tif.save(cell)
    with TiffWriter(nuc_out) as tif:
        tif.save(nuc)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("embryo")
    parser.add_argument("-v", "--view", type=int, default=0)
    parser.add_argument("-s", "--save", type=int, default=0)
    args = parser.parse_args()
    embryo = args.embryo
    assert embryo in ("embryo1_embryo2", "embryo3")

    data_folder = f"/g/kreshuk/data/marioni/shila/mouse-atlas-2020/{embryo}_raw"
    data_out_folder = f"/g/kreshuk/data/marioni/shila/mouse-atlas-2020/{embryo}"

    cell_in_folder = f"/g/kreshuk/data/marioni/shila/mouse-atlas-2020/segmentation_raw/{embryo}"
    cell_out_folder = f"/g/kreshuk/data/marioni/shila/mouse-atlas-2020/segmentation/{embryo}/cells"

    nuc_in_folder = f"/g/kreshuk/data/marioni/shila/mouse-atlas-2020/segmentation_raw/{embryo}"
    nuc_out_folder = f"/g/kreshuk/data/marioni/shila/mouse-atlas-2020/segmentation/{embryo}/nuclei"

    os.makedirs(cell_out_folder, exist_ok=True)
    os.makedirs(nuc_out_folder, exist_ok=True)
    os.makedirs(data_out_folder, exist_ok=True)

    images = glob(os.path.join(data_folder, "*.ome.tif"))
    images.sort()

    for image in tqdm(images):
        name = os.path.basename(image)
        pos = name[name.find("Pos"):].rstrip(".ome.tif")

        cell_pattern = os.path.join(cell_in_folder, f"boundaryseg_{pos}-*")
        cell_files = glob(cell_pattern)
        cell_files.sort()

        nuc_pattern = os.path.join(nuc_in_folder, f"dapi_{pos}-*")
        nuc_files = glob(nuc_pattern)
        nuc_files.sort()

        assert len(cell_files) == 6
        if len(nuc_files) != 6:
            print("Only", len(nuc_files), "nuclei slices for", name)
        if args.view:
            view(image, cell_files, nuc_files)
        if args.save:
            image_out = os.path.join(data_out_folder, name)
            cell_out = os.path.join(cell_out_folder, name)
            nuc_out = os.path.join(nuc_out_folder, name)
            save(image, image_out, cell_files, cell_out, nuc_files, nuc_out)


if __name__ == "__main__":
    main()
