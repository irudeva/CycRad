from netCDF4 import Dataset
import numpy as np
from windspharm.standard import VectorWind
from windspharm.tools import prep_data, recover_data, order_latdim
from scipy import interpolate

def polerot(plat,plon,ilat,ilon,n):
    rad=6371032 #in m
    dtr = pi/180.
    rtd = 180./pi
    plat = (plat-90.)*dtr
    plon = plon*dtr
    ilon = ilon*dtr
    ilat = ilat*dtr

    do i in range(0,n) :

        # convert to cartesian coordinates
        # from bronstein p.217

        x=rad*np.cos(ilat(i))*np.cos(ilon(i))
        y=rad*np.cos(ilat(i))*np.sin(ilon(i))
        z=rad*np.sin(ilat(i))


        # turn XY
        xn=x*np.cos(plon)+y*np.sin(plon)
        yn=x*(-np.sin(plon))+y*np.cos(plon)
        zn=z


        # tilt XZ
        xnn=xn*cos(plat)+zn*sin(plat)
        ynn=yn
        znn=xn*(-np.sin(plat))+zn*np.cos(plat)

        # convert back to polar coordinates
        rr=np.sqrt(xnn*xnn+ynn*ynn+znn*znn)
        nlat(i)=np.arcsin(znn/rr)*rtd


        if xnn==0. :
            if ynn > 0. :
                nlon(i)=0.
        	else
        	    nlon(i)=180.

        elif ynn==0. :
            nlon(i)=90.
        else
            nlon(i)=np.arctan(ynn,xnn)*rtd


    return nlat,nlon




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
    np = nnp.zeros(max_ntrk)  # number of points in a track

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

            np[ntrk-1]=nit
        for n in range(0,nit):
            l = f.readline()
            columns = l.split()
            lon[ntrk-1,n]=float(columns[7])
            lat[ntrk-1,n]=float(columns[8])
            cslp[ntrk-1,n]=float(columns[9])
            #print lon[ntrk-1,n],lat[ntrk-1,n],cslp[ntrk-1,n]

    f.close()
    print 'ftrk closed'

    #get radius



    #pole rotation
    for itrk in range(0,ntrk):
        for ip in range(0,np(ntrk-1)):



    quit()
