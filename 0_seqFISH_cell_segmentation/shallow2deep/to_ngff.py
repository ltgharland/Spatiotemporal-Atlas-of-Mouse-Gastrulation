import argparse
import os
from glob import glob

import imageio
import ome_zarr
import ome_zarr.scale
import ome_zarr.writer
import zarr
from tqdm import tqdm


def get_axes_and_trafos(mip, axis_names, units, resolution):
    axes = []
    for ax in axis_names:
        axis = {"name": ax, "type": "channel" if ax == "c" else "space"}
        unit = units.get(ax, None)
        if unit is not None:
            axis["unit"] = unit
        axes.append(axis)

    is_scaled = {"c": False, "z": False, "y": True, "x": True}
    trafos = [
        [{
            "scale": [resolution[ax] * 2**scale_level if is_scaled[ax] else resolution[ax] for ax in axis_names],
            "type": "scale"
        }]
        for scale_level in range(len(mip))
    ]

    return axes, trafos


def get_storage_opts(axis_names):
    # provide additional storage options for zarr
    chunks = (1, 1, 512, 512) if len(axis_names) == 4 else (1, 512, 512)
    assert len(chunks) == len(axis_names)
    return {"chunks": chunks}


def convert_image_data(in_path, group, resolution, units, name):
    # load the input data from ome.tif
    vol = imageio.volread(in_path)

    # the data is stored as 'zcyx'. This is currently not allowed by ome.zarr, so we reorder to 'czyx'
    vol = vol.transpose((1, 0, 2, 3))

    # create scale pyramid
    scaler = ome_zarr.scale.Scaler()
    mip = scaler.local_mean(vol)

    # specify the axis metadata
    axis_names = tuple("czyx")
    axes, trafos = get_axes_and_trafos(mip, axis_names, units, resolution)
    storage_opts = get_storage_opts(axis_names)

    ome_zarr.writer.write_multiscale(
        mip, group, name=name,
        axes=axes, coordinate_transformations=trafos,
        storage_options=storage_opts
    )


def convert_label_data(in_path, group, resolution, units, label_name, colors=None):
    # load the input data from ome.tif
    vol = imageio.volread(in_path)
    if vol.ndim != 3:
        print("Labels have unexpected shape", vol.shape, "for", in_path, label_name)
        print("Adding these labels will be skipped")
        return False

    # create scale pyramid
    scaler = ome_zarr.scale.Scaler()
    mip = scaler.nearest(vol)

    # specify the axis metadata
    axis_names = tuple("zyx")
    axes, trafos = get_axes_and_trafos(mip, axis_names, units, resolution)
    storage_opts = get_storage_opts(axis_names)

    ome_zarr.writer.write_multiscale_labels(
        mip, group, label_name,
        axes=axes, coordinate_transformations=trafos,
        storage_options=storage_opts, colors=colors
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("embryo")
    args = parser.parse_args()
    embryo = args.embryo
    assert embryo in ("embryo1_embryo2", "embryo3"), embryo

    input_folder = f"/g/kreshuk/data/marioni/shila/mouse-atlas-2020/{embryo}"
    assert os.path.exists(input_folder), input_folder
    output_folder = f"/g/kreshuk/data/marioni/shila/mouse-atlas-2020/ngff/{embryo}"
    os.makedirs(output_folder, exist_ok=True)
    images = glob(os.path.join(input_folder, "*.ome.tif"))

    cell_segmentation_folder = f"/g/kreshuk/data/marioni/shila/mouse-atlas-2020/segmentation/{embryo}/cells"
    nucleus_segmentation_folder = f"/g/kreshuk/data/marioni/shila/mouse-atlas-2020/segmentation/{embryo}/nuclei"

    # the xy-resolution is different for the two embryos
    if embryo == "embryo3":
        resolution = {"c": 1.0, "z": 4.0, "y": 0.17, "x": 0.17}
    else:
        resolution = {"c": 1.0, "z": 4.0, "y": 0.11, "x": 0.11}

    # the labels are downsampled by a factor of 4 in xy
    label_resolution = {ax: res * 4 if ax in "xy" else res for ax, res in resolution.items()}
    units = {"c": None, "z": "micrometer", "y": "micrometer", "x": "micrometer"}

    # do we store each position as a single ome.zarr or do we merge accoriding to the grid info?
    # discuss with shila!
    for image in tqdm(images, desc=f"Convert images from {input_folder} to ngff"):
        name = os.path.basename(image)
        out_name = name.replace(".ome.tif", ".ome.zarr")
        out_path = os.path.join(output_folder, out_name)

        loc = ome_zarr.io.parse_url(out_path, mode="w")
        group = zarr.group(loc.store)
        convert_image_data(image, group, resolution, units, name="image")

        cell_segmentation = os.path.join(cell_segmentation_folder, name)
        assert os.path.exists(cell_segmentation)
        convert_label_data(cell_segmentation, group, label_resolution, units, label_name="cells")

        nucleus_segmentation = os.path.join(nucleus_segmentation_folder, name)
        assert os.path.exists(nucleus_segmentation)
        colors = [{"label-value": 1, "rgba": [0, 0, 255, 255]}]
        convert_label_data(nucleus_segmentation, group, label_resolution, units, label_name="nuclei", colors=colors)


if __name__ == "__main__":
    main()
