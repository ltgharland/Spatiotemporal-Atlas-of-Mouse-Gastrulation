import argparse
import os

import bioimageio.core
import imageio
import numpy as np
import pandas as pd
import tifffile

from bioimageio.core.prediction import predict_with_tiling
from elf.segmentation.mutex_watershed import mutex_watershed as mws
from elf.segmentation.watershed import apply_size_filter
from scipy.ndimage.morphology import distance_transform_edt
from skimage.segmentation import watershed
from tqdm import tqdm


def predict_image(model, image):
    assert image.ndim == 2, f"{image.shape}"
    image = image[None, None]
    tiling = {"tile": {"x": 1024, "y": 1024}, "halo": {"x": 128, "y": 128}}
    pred = predict_with_tiling(model, image, tiling=tiling)[0]
    return np.array(pred[0, 0]), np.array(pred[0, 1:])


def _get_seeds(annotations, shape, radius=6):
    coords = np.round(annotations.values).astype("int")
    coords = tuple(coords[:, i] for i in range(coords.shape[1]))

    assert len(coords) == len(shape), f"{len(coords)}, {shape}"
    if any(coord.max() >= sh for coord, sh in zip(coords, shape)):
        print("Coords exceed shape:", [coord.max() for coord in coords], shape)
        print("Clipping to shape")
        coords = tuple(np.clip(coord, 0, sh - 1) for coord, sh in zip(coords, shape))

    seeds = np.zeros(shape, dtype="uint32")
    seed_ids = np.arange(1, len(annotations) + 1)
    seeds[coords] = seed_ids

    dist = distance_transform_edt(seeds == 0)
    seeds = watershed(dist.max() - dist, seeds, mask=dist < radius)
    return seeds


def _get_foreground(fg_pred, seeds, threshold):
    foreground = fg_pred > threshold
    # enforce seeds in foreground
    foreground = np.logical_or(foreground, seeds > 0)
    return foreground


def _get_hmap(affs, foreground, seeds, alpha=0.9):
    boundaries = np.mean(affs[:4], axis=0)
    dists = distance_transform_edt(seeds == 0)
    dists[~foreground] = 0
    max_dist = np.percentile(dists, 95)
    dists /= max_dist
    dists = np.clip(dists, 0, 1)
    boundaries = alpha * boundaries + (1. - alpha) * dists
    return boundaries


def segment_seeded_watershed(fg_pred, affs, annotations, threshold=0.5):
    seeds = _get_seeds(annotations, fg_pred.shape)
    foreground = _get_foreground(fg_pred, seeds, threshold)
    hmap = _get_hmap(affs, foreground, seeds)
    assert hmap.shape == foreground.shape

    # import napari
    # v = napari.Viewer()
    # v.add_image(foreground)
    # v.add_image(hmap)
    # napari.run()

    seg = watershed(hmap, seeds, mask=foreground)
    return seg


def segment_mws(fg, affs, threshold=0.5, strides=[3, 3], min_size=50):
    offsets = [[-1, 0], [0, -1], [-3, 0], [0, -3],
               [-9, 0], [0, -9], [-27, 0], [0, -27]]
    mask = fg > threshold
    seg = mws(affs, offsets, strides, randomize_strides=True, mask=mask).astype("uint32")

    if min_size > 0:
        seg = apply_size_filter(seg, affs[0], min_size, exclude=[0])[0]
        seg[~mask] = 0

    return seg


def _write_res(path, seg, full_shape):
    assert seg.shape == full_shape
    with tifffile.TiffWriter(path) as tif:
        tif.save(seg)


def segment_single_image(model, image, annotations, use_mws):
    fg, affs = predict_image(model, image)
    seg_ws = segment_seeded_watershed(fg, affs, annotations)
    seg_mws = segment_mws(fg, affs) if use_mws else None
    return seg_ws, seg_mws


def segment_image_stack(model, image, annotations, use_mws):
    seg_ws = np.zeros_like(image, dtype="uint16")
    seg_mws = np.zeros_like(image, dtype="uint16") if use_mws else None
    slices = np.unique(annotations["axis-0"].values).astype("int")
    for z in slices:
        zann = annotations[annotations["axis-0"] == z].drop(columns=["axis-0"])
        wsz, mwsz = segment_single_image(model, image[z], zann, use_mws)
        seg_ws[z] = wsz
        if use_mws:
            seg_mws[z] = mwsz
    return seg_ws, seg_mws


def segment_image(model, data_path, annotation_path, seg_root, name, use_mws=False):
    image = imageio.volread(data_path)[:, -1]
    annotations = pd.read_csv(annotation_path).drop(columns=["index", "axis-1"])
    seg_ws, seg_mws = segment_image_stack(model, image, annotations, use_mws)
    path_ws = _to_tif(os.path.join(seg_root, "watershed"), name)
    # return
    _write_res(path_ws, seg_ws, image.shape)
    if use_mws:
        path_mws = _to_tif(os.path.join(seg_root, "mutex_watershed"), name)
        _write_res(path_mws, seg_mws, image.shape)


def _load_model(version, device):
    ckpt = f"./checkpoints/nuclei_v{version}/bioimageio-model-affinities/custom-nucleus-segmentation.zip"
    assert os.path.exists(ckpt), ckpt
    resource = bioimageio.core.load_resource_description(ckpt)
    devices = None if device is None else [device]
    model = bioimageio.core.prediction_pipeline.create_prediction_pipeline(bioimageio_model=resource, devices=devices)
    return model


def _to_tif(data_folder, name):
    pos = name[name.find("Pos"):].rstrip(".csv")
    return os.path.join(data_folder, f"MMStack_{pos}.ome.tif")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-v", "--version", required=True, type=int)
    parser.add_argument("-d", "--device", default=None, type=str)
    args = parser.parse_args()

    annotation_folder = args.input
    assert os.path.exists(annotation_folder), annotation_folder
    timepoint, cycle = args.input.split("/")[-2:]
    names = os.listdir(annotation_folder)
    names.sort()

    data_folder = f"/g/kreshuk/data/marioni/shila/{timepoint}/{cycle}"
    assert os.path.exists(data_folder), data_folder

    seg_root = f"/g/kreshuk/data/marioni/shila/nucleus_segmentation/{timepoint}/{cycle}"
    os.makedirs(seg_root, exist_ok=True)
    os.makedirs(os.path.join(seg_root, "watershed"), exist_ok=True)
    # os.makedirs(os.path.join(seg_root, "mutex_watershed"), exist_ok=True)

    model = _load_model(args.version, args.device)

    for name in tqdm(names):
        annotation_path = os.path.join(annotation_folder, name)
        data_path = _to_tif(data_folder, name)
        segment_image(model, data_path, annotation_path, seg_root, name, use_mws=False)


if __name__ == "__main__":
    main()
