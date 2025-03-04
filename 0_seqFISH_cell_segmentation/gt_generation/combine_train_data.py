import os
from glob import glob
from shutil import copyfile

import imageio
import pandas as pd


def _to_tif(name):
    pos = name[name.find("Pos"):].rstrip(".csv")
    return f"MMStack_{pos}.ome.tif"


def combine_train_data():
    out_folder = "./data/ground_truth/v4"
    os.makedirs(out_folder, exist_ok=True)

    gt_image_folder = "./data/nuclei"
    gt_label_folder = "./data/ground_truth/v3"

    names = []
    gt_labels = glob(os.path.join(gt_label_folder, "*.ome.tif"))
    for gt_label in gt_labels:
        name = os.path.basename(gt_label)
        assert os.path.exists(os.path.join(gt_image_folder, name))
        copyfile(gt_label, os.path.join(out_folder, name))
        names.append(name)

    data_folder = "/g/kreshuk/data/marioni/shila/TimEmbryos-020420/HybCycle_29"
    seg_folder = "/g/kreshuk/data/marioni/shila/segmentation/watershed"
    annotation_folder = "../point_annotations/TimEmbryos-020420/HybCycle_29"

    all_names = os.listdir(annotation_folder)
    all_names.sort()
    for name in all_names:
        tif_name = _to_tif(name)
        if tif_name in names:
            print("SKIPPPING!", tif_name)
            continue
        annotations = pd.read_csv(os.path.join(annotation_folder, name)).iloc[:, 1:]
        cz = annotations.iloc[0, :2]
        annotations = annotations.iloc[:, 2:]
        im = imageio.volread(os.path.join(data_folder, tif_name))[int(cz[0]), int(cz[1])]
        seg_name = name.replace(".csv", ".ome.tif")
        seg = imageio.volread(os.path.join(seg_folder, seg_name))[int(cz[0]), int(cz[1])]

        out_labels = os.path.join(out_folder, tif_name)
        imageio.imwrite(out_labels, seg)

        out_image = os.path.join(gt_image_folder, tif_name)
        print(out_image)
        imageio.imwrite(out_image, im)


combine_train_data()
