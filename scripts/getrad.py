from netCDF4 import Dataset
import numpy as np
from windspharm.standard import VectorWind
from windspharm.tools import prep_data, recover_data, order_latdim
from scipy import interpolate

dset = "erain"
scheme = "UOM"

for yr in range(2005,2006):

    fnc  = "/Users/irudeva/work/DATA/%st/erain.mslp.%d.nc"%(dset,yr)
    ftrk = "../data/trk%s/trk_%s/ntrkdat.%s.%d"%(scheme,dset,dset,yr)
    frad = "../output/rad%s.%s.%d.trk"%(scheme,dset,yr)
    print "fnc  =",fnc
    print "ftrk =",ftrk
    print "frad =",frad

    # read netcdf
    print 'fnc reading:'
    dimnam=('longitude','latitude','time')
    varnam=['longitude','latitude','time','msl']

    nc = Dataset(fnc, 'r')
    v=0
    for var in varnam:
        if nc.variables[varnam[v]].name != var:
            print "Variables don't agree", var, nc.variables[varnam[v]].name, v
            exit()
        v += 1

    lons = nc.variables[varnam[0]][:]
    lats = nc.variables[varnam[1]][:]
    time = nc.variables[varnam[2]][:]
    mslp = nc.variables[varnam[3]][:]/100.
    mslp0=mslp/100.


    print "\nSLP Interpolation"
    # data from south to north

    #for t in range(time.size) :
    slpint = interpolate.interp2d(lons, lats, mslp[0,:,:], kind='cubic')
    print slpint(0,30)

    print "end interpolation"


    #read trk

    print 'ftrk reading:'
    max_ntrk = 20000
    max_trklength = 200
    lon  = np.zeros((max_ntrk,max_trklength))
    lat  = np.zeros_like(lon)
    cslp = np.zeros_like(lon)

    f = open(ftrk, 'r')

    for line in range(1,max_ntrk):
        break
        emptyline = f.readline()
        header    = f.readline()
        emptyline = f.readline()
        emptyline = f.readline()
        if header == '':
         print ' ftrk is at the eof'
         break
        else :
         print ' header: ', header.strip()

        columns = header.split()
        tmp=columns[1]
        if(tmp[-1] == "(") :
            ntrk = int(tmp[:-1])
            nit=int(columns[14])
        elif (tmp[-1] == ":") :
            ntrk = int(tmp[:-8])
            nit=int(columns[13])

        if nit > max_trklength :
            print " ERROR!!! Track length is more than max_trklength: nit = %d, max_trklength = %d" %(nit,max_trklength)
            quit()

        #print 'ntrk=',ntrk, 'nit=',nit

        for n in range(0,nit):
            l = f.readline()
            columns = l.split()
            lon[ntrk-1,n]=float(columns[7])
            lat[ntrk-1,n]=float(columns[8])
            cslp[ntrk-1,n]=float(columns[9])
            #print lon[ntrk-1,n],lat[ntrk-1,n],cslp[ntrk-1,n]

    f.close()
    print 'ftrk closed'


    quit()
