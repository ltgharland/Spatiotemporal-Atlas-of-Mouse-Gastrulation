import argparse
import os
import subprocess

import imageio
import h5py
import numpy as np
import tifffile
from scipy.ndimage.morphology import distance_transform_edt
from skimage.segmentation import watershed
from tqdm import tqdm

# from ilastik.experimental.api import from_project_file
# from xarray import DataArray


def predict_enhancer(model, input_):
    import torch
    from torch_em.util import get_trainer

    boundaries = np.zeros(input_.shape, dtype="float32")
    model = get_trainer(model, device="cuda").model
    model.eval()
    with torch.no_grad():
        for z in range(boundaries.shape[0]):
            input_z = input_[z][None, None]
            input_z = torch.from_numpy(input_z).to("cuda")
            bd = model(input_z)
            boundaries[z] = bd.detach().cpu().numpy().squeeze()
    return boundaries


def predict_with_ilastik(ilp, image):
    tmp_input = "./ilastik-projects/tmp-raw.h5"
    tmp_output = "./ilastik-projects/tmp-pred.h5"
    with h5py.File(tmp_input, "w") as f:
        f.create_dataset("raw", data=image, chunks=True)

    if os.path.exists(tmp_output):
        os.remove(tmp_output)

    ilastik_folder = "/g/kreshuk/pape/Work/software/src/ilastik/ilastik-1.4.0b8-Linux"
    ilastik_exe = os.path.join(ilastik_folder, "run_ilastik.sh")
    assert os.path.exists(ilastik_exe), ilastik_exe
    input_str = f"{tmp_input}/raw"
    cmd = [ilastik_exe, "--headless",
           "--project=%s" % ilp,
           "--output_format=hdf5",
           "--raw_data=%s" % input_str,
           "--output_filename_format=%s" % tmp_output,
           "--readonly=1"]
    print("Run ilastik prediction ...")
    subprocess.run(cmd)
    with h5py.File(tmp_output, "r") as f:
        pred = f["exported_data"][:]
    return pred


def segment_image(input_path, seed_path, ilp, model, out_path, view):
    image = imageio.volread(input_path)

    # print("Run prediction with ilastik ...")
    # FIXME bug in ilastik, for now we use the hard-coded prediction
    # ilp = from_project_file(ilp)
    # image = DataArray(image, dims=tuple("zcyx"))
    # pred = ilp.predict(image)
    pred = predict_with_ilastik(ilp, image)
    # note: Shila has labeled 1 as boundary and 0 as non-boundary
    boundaries = (pred[:, 0] - pred[:, 1]).astype("float32")

    seeds = imageio.volread(seed_path)
    assert boundaries.shape == seeds.shape, f"{boundaries.shape}, {seeds.shape}"

    if model is not None:
        boundaries = predict_enhancer(model, boundaries)

    # print("Run foreground segmentation ...")
    # segment foreground / background
    fg_seed = (seeds > 0).astype("uint8")
    bg_seed = np.ones(fg_seed.shape, dtype="bool")
    bg_seed[:, 5:-5, 5:-5] = 0
    bg_seed_exclusion = distance_transform_edt(fg_seed == 0) <= 5
    bg_seed[bg_seed_exclusion] = 0
    fg_seed[bg_seed] = 2
    fg_mask = np.zeros_like(seeds, dtype="bool")
    for z in range(fg_mask.shape[0]):
        fg_mask[z] = (watershed(boundaries[z], markers=fg_seed[z]) == 1)

    # print("Run segmentation ...")
    # segment cells
    cell_seg = np.zeros_like(seeds, dtype="uint16")
    for z in range(cell_seg.shape[0]):
        if (seeds[z] != 0).sum() == 0:
            continue
        cell_seg[z] = watershed(boundaries[z], markers=seeds[z], mask=fg_mask[z])

    with tifffile.TiffWriter(out_path) as tif:
        tif.save(cell_seg)

    if view:
        import napari
        v = napari.Viewer()
        v.add_image(image[:, -1])
        v.add_image(pred)
        v.add_image(boundaries)
        v.add_labels(seeds)
        v.add_labels(fg_mask.astype("uint8"))
        v.add_labels(cell_seg)
        napari.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-e", "--use_enhancer", type=int, default=0)
    parser.add_argument("-v", "--view", type=int, default=0)
    args = parser.parse_args()
    use_enhancer = bool(args.use_enhancer)

    ilp = "./ilastik-projects/shilaV2.ilp"

    model = "./shallow2deep/checkpoints/embryo-cell-boundaries" if use_enhancer else None
    seg_name = "enhancer" if use_enhancer else "vanilla"

    input_folder = args.input
    split = input_folder.find("shila") + len("shila")
    nucleus_seg_folder = os.path.join(
        input_folder[:split], "nucleus_segmentation", input_folder[split:].lstrip("/"), "watershed"
    )
    assert os.path.exists(nucleus_seg_folder), nucleus_seg_folder

    out_folder = nucleus_seg_folder.replace("nucleus_segmentation", "cell_segmentation")
    out_folder = out_folder.replace("watershed", seg_name)
    os.makedirs(out_folder, exist_ok=True)

    names = os.listdir(nucleus_seg_folder)
    for name in tqdm(names):
        input_path = os.path.join(input_folder, name)
        assert os.path.exists(input_path), input_path
        seed_path = os.path.join(nucleus_seg_folder, name)
        out_path = os.path.join(out_folder, name)
        segment_image(input_path, seed_path, ilp, model, out_path, bool(args.view))


if __name__ == "__main__":
    main()
