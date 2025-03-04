import imageio
from ilastik.experimental.api import from_project_file
from xarray import DataArray
input_path = "/g/kreshuk/data/marioni/shila/TimEmbryos-020420/HybCycle_0/MMStack_Pos0.ome.tif"
image = imageio.volread(input_path)
ilp = from_project_file("/g/kreshuk/pape/Work/my_projects/deep-micro-val/custom_data/shila/shila-boundaries.ilp")
image = DataArray(image, dims=tuple("zcyx"))
pred = ilp.predict(image)
