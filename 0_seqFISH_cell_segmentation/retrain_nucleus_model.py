import argparse
import os
import sys
from glob import glob

import imageio
from torch_em.util import export_bioimageio_model


def to_modelzoo(version, to_boundaries):
    input_data = imageio.imread("./gt_generation/data/nuclei/MMStack_Pos0.ome.tif")[:256, :256]
    checkpoint = f"./checkpoints/nuclei_v{version}"
    if to_boundaries:
        postprocessing = "affinities_with_foreground_to_boundaries2d"
        output = f"./checkpoints/nuclei_v{version}/bioimageio-model"
    else:
        postprocessing = None
        output = f"./checkpoints/nuclei_v{version}/bioimageio-model-affinities"
    export_bioimageio_model(
        checkpoint, output,
        input_data=input_data,
        name="custom-nucleus-segmentation",
        authors=[{"name": "Constantin Pape; @constantinpape"}],
        tags=["nucleus-segmentation"],
        license="CC-BY-4.0",
        documentation="Doc that",
        git_repo='https://github.com/constantinpape/torch-em.git',
        model_postprocessing=postprocessing,
        input_optional_parameters=False
    )


def retrain_for_nuclei(version, pretrained, batch_size, n_iterations):
    sys.path.append("../..")
    from segmentation.training.default_models import get_nucleus_trainer

    label_paths = glob(f"./gt_generation/data/ground_truth/v{version}/*.ome.tif")
    fnames = [os.path.split(lp)[1] for lp in label_paths]
    raw_paths = [os.path.join("./gt_generation/data/nuclei", fname) for fname in fnames]

    name = f"nuclei_v{version}"
    trainer = get_nucleus_trainer(name,
                                  raw_paths, label_paths,
                                  raw_paths, label_paths,
                                  batch_size=batch_size, patch_shape=(512, 512))
    trainer.fit(n_iterations)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", type=int, required=True)
    parser.add_argument("-c", "--convert", type=int, default=0)
    parser.add_argument("-t", "--to_boundaries", type=int, default=1)
    parser.add_argument("-p", "--pretrained", type=int, default=1)
    parser.add_argument("-b", "--batch_size", type=int, default=4)
    parser.add_argument("-n", "--n_iterations", type=int, default=10000)
    args = parser.parse_args()
    if bool(args.convert):
        to_modelzoo(args.version, bool(args.to_boundaries))
    else:
        retrain_for_nuclei(args.version, bool(args.pretrained), args.batch_size, args.n_iterations)


if __name__ == "__main__":
    main()
