from yt.mods import *
import Scientific.IO.NetCDF
from numpy import *
from Scientific.IO.NetCDF import *
from yt.utilities.physical_constants import cm_per_kpc
import h5py

file = NetCDFFile('sod-mpi1007.nc','r')
data = {}
for k in file.variables.keys():
	data[k] = file.variables[k].getValue()

bbox = array([[0.0,1.0],[0.0,1.0],[0.0,1.0]])
pf = load_uniform_grid(data,data["Density"].shape,100.*cm_per_kpc,bbox=bbox,nprocs=1,periodicity=(False,False,False))
slc = SlicePlot(pf,"x",["Density"])
slc.save()
