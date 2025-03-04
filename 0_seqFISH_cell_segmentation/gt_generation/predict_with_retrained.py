import argparse
import os

import bioimageio.core
import imageio
import h5py
from bioimageio.core.prediction import predict
from tqdm import tqdm


def _predict(pp, in_path, out_path):
    input_ = imageio.imread(in_path)[None, None]
    output = predict(pp, input_)[0]
    with h5py.File(out_path, "w") as f:
        f.create_dataset("prediction", data=output, compression="gzip")


def run_prediction(version, device):
    input_folder = "./data/nuclei"
    output_folder = f"./data/nn_predictions/v{version}"
    os.makedirs(output_folder, exist_ok=True)
    names = os.listdir(input_folder)
    model = bioimageio.core.load_resource_description(
        f"../checkpoints/nuclei_v{version}/bioimageio-model/custom-nucleus-segmentation.zip"
    )
    devices = None if device is None else [device]
    pp = bioimageio.core.prediction_pipeline.create_prediction_pipeline(bioimageio_model=model, devices=devices)
    for name in tqdm(names):
        in_path = os.path.join(input_folder, name)
        out_name = name.replace(".ome.tif", ".h5")
        out_path = os.path.join(output_folder, out_name)
        _predict(pp, in_path, out_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", required=True, type=int)
    parser.add_argument("-d", "--device", default=None, type=str)
    args = parser.parse_args()
    version = args.version
    run_prediction(version, args.device)


if __name__ == "__main__":
    main()
