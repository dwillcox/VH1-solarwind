from yt.mods import *
import Scientific.IO.NetCDF
from numpy import *
from Scientific.IO.NetCDF import *
from yt.utilities.physical_constants import cm_per_kpc
import h5py

filePrefix='sod-mpi'
fileSuffix='.nc'

startNum = 1000
endNum = 1050

dim = 'z'
field = 'Density'

bbox = array([[0.0,1.0],[0.0,1.0],[0.0,1.0]])
scale = 100.*cm_per_kpc

for k in range(startNum,endNum+1,1):
	fileName = filePrefix + str(k) + fileSuffix
	outputName = ('Slice_' + field + '_' + dim + 
			'_' + filePrefix + '_' + str(k))
	file = NetCDFFile(fileName,'r')
	data = {}
	for k in file.variables.keys():
		data[k] = file.variables[k].getValue()
	pf = load_uniform_grid(data,data[field].shape,
			scale,bbox=bbox,nprocs=1,
			periodicity=(False,False,False))
	slc = SlicePlot(pf,dim,[field],origin=('native'))
	slc.save(name=outputName)
