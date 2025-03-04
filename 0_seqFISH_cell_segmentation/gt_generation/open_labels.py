import argparse
import imageio
import napari


def open_labels(image_path, label_path):
    im = imageio.volread(image_path)
    labels = imageio.volread(label_path)
    assert im.shape == labels.shape

    v = napari.Viewer()
    v.add_image(im)
    v.add_labels(labels)
    napari.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--image", help="Filepath to the tif with image data", required=True)
    parser.add_argument("-l", "--labels", help="Filepath to the tif with the label data (segmentation)", required=True)
    args = parser.parse_args()
    open_labels(args.image, args.labels)


if __name__ == "__main__":
    main()
