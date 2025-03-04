import os
from glob import glob
from shutil import copyfile

import imageio
import numpy as np
from tifffile import TiffWriter


def fix_data(path, out_path):
    vol = imageio.volread(path)
    shape = vol.shape
    x = vol.reshape((shape[0] * shape[1], shape[2], shape[3]))
    y = np.zeros_like(vol)

    nc, nz = shape[1], shape[0]
    for c in range(nc):
        start = c * nz
        stop = start + nz
        y[:, c] = x[start:stop]

    with TiffWriter(out_path) as tif:
        tif.save(y)


def fix_all_tiffs(in_folder, out_folder, exclude=None):
    input_files = glob(os.path.join(in_folder, "*.ome.tif"))
    os.makedirs(out_folder, exist_ok=True)
    for ff in input_files:
        name = os.path.basename(ff)
        out_path = os.path.join(out_folder, name)
        if exclude is not None and name in exclude:
            copyfile(ff, out_path)
            continue
        fix_data(ff, out_path)


def fix_020420_cycle0():
    in_folder = "/g/kreshuk/data/marioni/shila/TimEmbryos-020420/HybCycle_0_original"
    out_folder = "/g/kreshuk/data/marioni/shila/TimEmbryos-020420/HybCycle_0"
    exclude = ["MMStack_Pos0.ome.tif"]
    fix_all_tiffs(in_folder, out_folder, exclude)


def fix_030320_cycle0():
    in_folder = "/g/kreshuk/data/marioni/shila/TimEmbryos-030320/HybCycle_0_original"
    out_folder = "/g/kreshuk/data/marioni/shila/TimEmbryos-030320/HybCycle_0"
    exclude = ["MMStack_Pos0.ome.tif"]
    fix_all_tiffs(in_folder, out_folder, exclude)


if __name__ == "__main__":
    # fix_020420_cycle0()
    fix_030320_cycle0()
