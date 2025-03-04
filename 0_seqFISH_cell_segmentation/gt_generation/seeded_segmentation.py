import argparse
import os

import imageio
import h5py
import napari
import numpy as np
import pandas as pd
import tifffile

from scipy.ndimage.morphology import distance_transform_edt
from skimage.segmentation import watershed


def get_seeds(im, points, radius=6):
    coords = np.round(points).astype("int")
    coords = tuple(coords[:, i] for i in range(coords.shape[1]))

    # this is more correct
    seeds = np.zeros(im.shape, dtype="uint32")
    seed_ids = np.arange(1, len(points) + 1)
    seeds[coords] = seed_ids

    dist = distance_transform_edt(seeds == 0)
    seeds = watershed(dist.max() - dist, seeds, mask=dist < 4,
                      compactness=10)
    return seeds


def get_foreground_and_boundaries(pred_path, seeds, threshold=0.5, alpha=0.3):
    with h5py.File(pred_path, "r") as f:
        data = f["exported_data"][:]

    foreground = data[..., 2] < threshold  # 2 is the bg channel
    # enforce seeds in foreground
    foreground = np.logical_or(foreground, seeds > 0)

    boundaries = data[..., 1]  # 1 is the boundary channel
    dists = distance_transform_edt(seeds == 0)
    dists[~foreground] = 0
    max_dist = np.percentile(dists, 95)
    dists /= max_dist
    dists = np.clip(dists, 0, 1)

    boundaries = alpha * boundaries + (1. - alpha) * dists
    return foreground, boundaries


def get_nn_foreground_and_boundaries(pred_path, seeds, threshold=0.5, alpha=0.9):
    with h5py.File(pred_path, "r") as f:
        data = f["prediction"][:].squeeze()

    foreground = data[0] > threshold
    # enforce seeds in foreground
    foreground = np.logical_or(foreground, seeds > 0)

    boundaries = np.mean(data[1:5], axis=0)
    dists = distance_transform_edt(seeds == 0)
    dists[~foreground] = 0
    max_dist = np.percentile(dists, 95)
    dists /= max_dist
    dists = np.clip(dists, 0, 1)

    boundaries = alpha * boundaries + (1. - alpha) * dists
    return foreground, boundaries


def seeded_segmentation(im_path, pred_path, ann_path, view=False):
    im = imageio.imread(im_path)
    points = pd.read_csv(ann_path).values[:, 3:]

    seeds = get_seeds(im, points)
    if "Probabilities" in pred_path:
        foreground, boundaries = get_foreground_and_boundaries(pred_path, seeds)
    else:
        foreground, boundaries = get_nn_foreground_and_boundaries(pred_path, seeds)
    seg = watershed(boundaries, seeds, mask=foreground)

    if view:
        v = napari.Viewer()
        v.add_image(im)
        v.add_points(points)

        v.add_labels(seeds)
        v.add_labels(foreground, visible=False)
        v.add_image(boundaries)
        v.add_labels(seg)

        napari.run()

    return seg


def _to_tif(name):
    pos = name[name.find("Pos"):name.find("Pos")+4]
    return f"MMStack_{pos}.ome.tif"


def _to_pred(name):
    return name.replace(".ome.tif", ".ome_Probabilities Stage 2.h5")


def _to_nn_pred(name):
    return name.replace(".ome.tif", ".h5")


def write_ome_tiff(path, data):
    with tifffile.TiffWriter(path) as tif:
        tif.save(data)


def export_segmentation(im_name, ann_name, seg, version):
    if version is None:
        out_folder = "./data/segmentations"
    else:
        out_folder = f"./data/segmentations/v{version}"
    os.makedirs(out_folder, exist_ok=True)
    im = imageio.volread(os.path.join("./data/raw", im_name))
    full_seg = np.zeros(im.shape, dtype="int16")
    cz = pd.read_csv(os.path.join("./data/annotations", ann_name)).values[0, 1:3].astype("int")
    full_seg[cz[0], cz[1]] = seg
    # v = napari.Viewer()
    # v.add_image(im)
    # v.add_labels(full_seg)
    # napari.run()
    seg_path = os.path.join(out_folder, im_name)
    write_ome_tiff(seg_path, full_seg)


def segment_all(save, version):
    annotation_folder = "./data/annotations"
    names = os.listdir(annotation_folder)
    data_folder = "./data/nuclei"
    if version is None:
        pred_folder = "./data/predictions"
    else:
        pred_folder = f"./data/nn_predictions/v{version}"
    for name in names:
        tif_name = _to_tif(name)
        pred_name = _to_pred(tif_name) if version is None else _to_nn_pred(tif_name)
        im_path = os.path.join(data_folder, tif_name)
        pred_path = os.path.join(pred_folder, pred_name)
        ann_path = os.path.join(annotation_folder, name)
        seg = seeded_segmentation(im_path, pred_path, ann_path, view=not save)
        if save:
            export_segmentation(tif_name, name, seg, version)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--save", type=int, default=0)
    parser.add_argument("-v", "--version", default=None, type=int)
    args = parser.parse_args()
    segment_all(bool(args.save), args.version)


if __name__ == "__main__":
    main()
