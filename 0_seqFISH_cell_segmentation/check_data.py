import argparse
import os
from glob import glob
import imageio
import napari


def check_inputs(input_folder):
    images = glob(os.path.join(input_folder, "*.ome.tif"))
    names = [os.path.basename(im) for im in images]
    images = [imageio.volread(im) for im in images]

    v = napari.Viewer()
    for name, im in zip(names, images):
        v.add_image(im, name=name)
    napari.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    args = parser.parse_args()
    check_inputs(args.input)
