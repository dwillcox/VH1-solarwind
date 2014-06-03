# Make a movie of each of the variables shown below in fnames

# All this script requires is the directory containing the png files
# as a command line option when this script is run, i.e.:
# $> python movieddtvars.py /path/to/awesome/pngs/

# Note this script expects the png naming convention to follow the glob
# pattern in inputfilePrefix.

import os
import sys
import glob

# List of field names
fnames = ['Density']

# sys.argv[0] = script name
# sys.argv[1] = directory containing the png files

# The following two lines aren't needed unless you have a static dist.
ffmpegDir = '/home/dwillcox/codes/ffmpeg-2.1.3-64bit-static/'
os.chdir(ffmpegDir)

fileDirect = sys.argv[1]
if(fileDirect[-1]!='/'):
	fileDirect+='/'

fileNameTemplate = '*_Slice_*.png'

infileTemplate = glob.glob(fileDirect + fileNameTemplate)[0]

inputfilePrefix = 'sod-mpi_%04d_Slice_z_'
outputfilePrefix = 'Slice_z_'

cmd1 = './ffmpeg -y -start_number 1000 -i ' + fileDirect + inputfilePrefix

cmd2 = '.png -c:v libx264 -pix_fmt yuv420p -r 10 ' + fileDirect + outputfilePrefix

cmd3 = '.mp4'

for fn in fnames:
    cmd = cmd1 + fn + cmd2 + fn + cmd3
    os.system(cmd)


