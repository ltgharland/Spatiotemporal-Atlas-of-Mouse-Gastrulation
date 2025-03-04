import os
from glob import glob

import imageio
import h5py
from skimage.transform import resize


def prepare_training_data(embryo):
    raw_folder = f"/g/kreshuk/data/marioni/shila/mouse-atlas-2020/{embryo}"
    assert os.path.exists(raw_folder), raw_folder
    images = glob(os.path.join(raw_folder, "*.ome.tif"))
    images.sort()
    label_folder = f"/g/kreshuk/data/marioni/shila/mouse-atlas-2020/segmentation/{embryo}/cells"
    assert os.path.exists(label_folder), label_folder
    labels = glob(os.path.join(label_folder, "*.ome.tif"))
    labels.sort()
    assert len(images) == len(labels), f"{len(images)}, {len(labels)}"

    out_folder = "/g/kreshuk/data/marioni/shila/mouse-atlas-2020/training_data"
    os.makedirs(out_folder, exist_ok=True)

    for ii, (im, lab) in enumerate(zip(images, labels)):
        out_path = os.path.join(out_folder, f"{embryo}_{ii}.h5")
        image = imageio.volread(im)[:, 0]
        label = imageio.volread(lab)
        assert image.shape[0] == label.shape[0]
        image = resize(image, label.shape)
        assert image.shape == label.shape
        with h5py.File(out_path, "w") as f:
            f.create_dataset("raw", data=image, compression="gzip")
            f.create_dataset("labels", data=label, compression="gzip")


def main():
    # prepare_training_data("embryo1_embryo2")
    prepare_training_data("embryo3")


# extract relevant channels and downsample raw to have training data for shallow2deep training
if __name__ == "__main__":
    main()
