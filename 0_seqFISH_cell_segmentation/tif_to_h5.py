import argparse
import os
import imageio
import h5py


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    args = parser.parse_args()
    path = args.path
    assert os.path.exists(path)

    vol = imageio.volread(path)
    assert vol.ndim == 4

    folder, name = os.path.split(path)
    out_name = name.replace(".ome.tif", ".h5")
    out_path = os.path.join(folder, out_name)

    print("Save tif file:", path)
    print("as h5 file to:", out_path)

    with h5py.File(out_path, "w") as f:
        f.create_dataset("data", data=vol, compression="gzip", chunks=(1, 1, 512, 512))


if __name__ == "__main__":
    main()
