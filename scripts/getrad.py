from netCDF4 import Dataset
import numpy as np
from windspharm.standard import VectorWind
from windspharm.tools import prep_data, recover_data, order_latdim
from scipy import interpolate

print np.pi


def polerot(plat,plon,ilat,ilon):
# the location of the new pole is (plat,plon)
# the new coordinates are in a Cartesian c. system where:
# real North pole is (y=90,x=0) in the new system if plat > 0 and the East is (y=0,x=90)
# real South pole is (90,0) in the new system if plat < 0 and the East is (y=0,x=-90)


    import numpy as np
    rad  =6371032 #in m
    pi   = np.pi
    dtr  = pi/180.
    rtd  = 180./pi

    plat0 = plat*1
    plon0 = plon*1.

    plat = (plat-90.)*dtr
    plon = plon*dtr

    print "ilat", ilat

    print 'ilon=',ilon
    ilon[:] = [x*dtr for x in ilon]
    ilat[:] = [x*dtr for x in ilat]

    nlat = np.zeros_like(ilat)
    nlon = np.zeros_like(ilon)

    for index, (lon,lat) in enumerate(zip(ilon, ilat)):

        # print lon*rtd,lat*rtd
        print 'ind=',index
        print "lat, lon = ",lat*rtd,lon*rtd
        print "plat, plon = ",plat*rtd,plon*rtd

        # convert to cartesian coordinates
        # from bronstein p.217

        cos_plon = np.cos(plon)
        if plon == 90*dtr or plon == -90*dtr:
            cos_plon = 0
        cos_lon = np.cos(lon)
        if lon == 90*dtr or lon == -90*dtr:
            cos_lon = 0

        x=rad*np.cos(lat)*cos_lon
        y=rad*np.cos(lat)*np.sin(lon)
        print "!!!!!!!!!np.cos(lon)", lon*rtd, 90, x
        z=rad*np.sin(lat)
        # print "z, rad, sin(lat) =",z, rad, np.sin(lat)
        print "x, y, z =", x,y,z


        # turn XY
        xn=x*cos_plon+y*np.sin(plon)
        yn=x*(-np.sin(plon))+y*cos_plon
        print "x*(-np.sin(plon))+y*np.cos(plon)",x,(-np.sin(plon)),y,cos_plon
        print "____yn", yn
        zn=z
        # print "xn, yn, zn =", xn,yn,zn

        # tilt XZ
        xnn=xn*np.cos(plat)+zn*np.sin(plat)
        # print "xn*np.cos(plat),zn*np.sin(plat)",xn*np.cos(plat),zn*np.sin(plat)
        ynn=yn
        # print ynn/xnn
        # print xnn/ynn
        znn=xn*(-np.sin(plat))+zn*np.cos(plat)
        # print "xnn, ynn, znn =", xnn,ynn,znn
        # print "plat, lat=",plat, lat
        # print ilat
        # print "zn,np.cos(plat)=",zn,np.cos(plat)

        # convert back to polar coordinates
        rr=np.sqrt(xnn*xnn+ynn*ynn+znn*znn)
        nlat[index]=np.arcsin(znn/rr)*rtd


        # if xnn==0.:
        #     print 'xnn = 0'
        #     if ynn > 0.:
        #       nlon[index]=0.
        #     else:
        # 	  nlon[index]=180.

        # elif ynn==0.:
        #     print 'ynn == 0'
        #     print plon, x,y,yn
        #     print x*(-np.sin(plon)),y*np.cos(plon)
        #     nlon[index]=90.
        # else:
        nlon[index]=np.arctan2(xnn,ynn)*rtd
        print "xnn ynn =", xnn,ynn
        print "xnn/ynn =", xnn/ynn

        print "res----", nlat[index],nlon[index]

    # if plat0>0 :
    #      nlon = -nlon
    # else :
    #      nlat = -nlat


    return nlat,nlon




dset = "erain"
scheme = "UOM"

for yr in range(2005,2006):

    fnc  = "/Users/Irina/work/DATA/%st/erain.mslp.%d.nc"%(dset,yr)
    ftrk = "../data/trk%s/trk_%s/ntrkdat.%s.%d"%(scheme,dset,dset,yr)
    frad = "../output/rad%s.%s.%d.trk"%(scheme,dset,yr)
    print "fnc  =",fnc
    print "ftrk =",ftrk
    print "frad =",frad

    # read netcdf
    # print 'fnc reading:'
    # dimnam=('longitude','latitude','time')
    # varnam=['longitude','latitude','time','msl']
    #
    # nc = Dataset(fnc, 'r')
    # v=0
    # for var in varnam:
    #     if nc.variables[varnam[v]].name != var:
    #         print "Variables don't agree", var, nc.variables[varnam[v]].name, v
    #         exit()
    #     v += 1
    #
    # lons = nc.variables[varnam[0]][:]
    # lats = nc.variables[varnam[1]][:]
    # time = nc.variables[varnam[2]][:]
    # mslp = nc.variables[varnam[3]][:]/100.
    # mslp0=mslp/100.


    print "\nSLP Interpolation"
    # data from south to north

    #for t in range(time.size) :
    # slpint = interpolate.interp2d(lons, lats, mslp[0,:,:], kind='cubic')
    # print slpint(0,30)

    print "end interpolation"


    #read trk

    print 'ftrk reading:'
    max_ntrk = 20000
    max_trklength = 200
    lon  = np.zeros((max_ntrk,max_trklength))
    lat  = np.zeros_like(lon)
    cslp = np.zeros_like(lon)
    np = np.zeros(max_ntrk, dtype=np.int)  # number of points in a track

    f = open(ftrk, 'r')

    for line in range(1,max_ntrk):
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
        print header
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

        break  # for testing - first tack only


    f.close()
    print 'ftrk closed'

    #get radius




    #pole rotation
    for itrk in range(0,1): #  ntrk):
        for ip in range(0,1): #np[itrk]):

        #  print 'lon = ', lon[itrk,ip]
        #  print 'lat = ', lat[itrk,ip]
         lon[itrk,ip] = 165
         lat[itrk,ip] = -11
        #  tmplat = [  90, -90, 15]
        #  tmplon = [ 165, 165, 175]
         tmplat = [  90, -90, -89.9, -89.9999]
         tmplon = [ 165,165, 165, -15]
         tmplat1 = tmplat*1
         tmplon1 = tmplon*1
         nlat, nlon = polerot(lat[itrk,ip],lon[itrk,ip],tmplat,tmplon)
         print 'lat = ', tmplat1
         print 'lon = ', tmplon1
         print 'nlat = ', nlat
         print 'nlon = ', nlon
        #  print nlat[0],nlon[0],nlat[1],nlon[1]
         print "place pole to (lat %f,lon %d)"%(nlat[0],nlon[0])
         nlat1, nlon1 = polerot(nlat[0], nlon[0],nlat,nlon)
         print 'nlat = ', nlat1
        #  print 'nlon = ', nlon1
        #  print 'nlon = ', nlon1+lon[itrk,ip] +90 # works for NH
         print 'nlon = ', nlon1
         print 'nlon = ', nlon1+90+lon[itrk,ip]
# quit()
