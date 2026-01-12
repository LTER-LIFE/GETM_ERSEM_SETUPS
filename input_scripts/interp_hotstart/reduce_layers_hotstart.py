#! /usr/bin/env python

# python script to read an .nc file, halving number of layers 
# and write
# a new file
# Works on a restart file

# import the relevant packages
import sys
sys.path.append('/share/apps/python2.6.6/lib/python2.6/site-packages/pynetcdf')
from netCDF4 import Dataset
#from NetCDF import *
from numpy import *
from pylab import *
#import Numeric

############## settings ##################################################

# from script nest_bdy
#infname=os.environ['infname']
#ofname=os.environ['ofname']

# set hard (comment out if using script)
infname='/export/lv1/user/jvandermolen/model_output/active_runs/boundaries/dws_200m_nwes/restart_201501_hydro.nc'

ofname='/export/lv9/user/qzhan/home/model_input_files/restart/restart_201501_hydro_reducedlayers.nc'

##################################################################################
# Main routine
#

# Open input files
print('halving layers in nc file.')
print('Input files:')
print(infname)
infile=Dataset(infname,'r',format='NETCDF4') #NetCDFFile(infname,'r')

print('Output file: ',ofname)

# get dimensions and variables from file
alldimnames=list(infile.dimensions.keys())
varnames=list(infile.variables.keys())
print('Variables: ',varnames)

# open and initialise output file
#outfile=NetCDFFile(ofname,'w')
outfile=Dataset(ofname,'w',format='NETCDF3_CLASSIC')
ndims=len(alldimnames)
for idim in range(ndims):
  dimname=alldimnames[idim]
  print(dimname)
  dimvalue=infile.dimensions[alldimnames[idim]]
  print('dimvalue', dimvalue)

  if alldimnames[idim] == 'zax':
    lendim=1+(len(dimvalue)-1)//2
  else:
    lendim=len(dimvalue)

  outfile.createDimension(alldimnames[idim],lendim)

# process

for varname in varnames:
  # read variable
  print(varname)
  var=infile.variables[varname]
  dimnames=var.dimensions
#  units=var.units
  datavals=var[:]
#  data_attlist=dir(var)
  datatype=datavals.dtype.kind
##  print datatype
#  if datatype=='f':
#    datatype='d'
#  if ('_FillValue' in data_attlist):  
#    mv=getattr(var,'_FillValue')
#  else:
#    if ('missing_value' in data_attlist):  
#       mv=getattr(var,'missing_value')
#    else:
#       mv=None

#  print(data_attlist)
  # save time variable
  if varname=='timestemp':
    time_2d=datavals
    time_2d_units=units
    
  if len(dimnames)==3:
    # adjust
    sv=shape(datavals)
    out=zeros(((sv[0]-1)//2+1,sv[1],sv[2]))
    if varname=='ho' or varname=='hn' :   # add 2 levels
      out[0,:,:]=datavals[0,:,:]
      for nl in range(1,(sv[0]-1)//2+1):
        out[nl,:,:]=datavals[2*nl-1,:,:]+datavals[2*nl,:,:]
    else:                                         # average 2 levels
      out[0,:,:]=datavals[0,:,:]
      for nl in range(1,(sv[0]-1)//2+1):
        out[nl,:,:]=(datavals[2*nl-1,:,:]+datavals[2*nl,:,:])/2
  else:
    out=datavals

  if varname=='zax':
    print('yes')
    sv=datavals.shape
    newz = 1 + (sv[0]-1)//2
    out = np.arange(newz, dtype=datavals.dtype)  

  # write variable
#  outvar=outfile.createVariable(varnames[j],datatype,dimnames,fill_value=mv,chunksizes=out.shape)

  print(out)
  outvar=outfile.createVariable(varname,datatype,dimnames,chunksizes=out.shape)
  sout=list(shape(out))
  print(sout)
  outvar[:]=out
  #if len(sout)>1:
  #  outvar[:]=out
   

  #else:
  #  print(outvar[:])
  #  print(sout)
  #  print(out)
  #  outvar[:]=out
  #  print(not(sout))
  #  if not(sout):
  #    outvar[:]=out
   # else:
   #   outvar[:]=out

  ##units and other attributes should be written as well!!
  # copy attributes
  for att in var.ncattrs():
    if (att!='_FillValue') & (att!='assignValue') & (att!='getValue') & (att!='typecode'):
      setattr(outvar,att,getattr(var,att))
 
  
# close files  
infile.close()
outfile.close()

print('Done')

