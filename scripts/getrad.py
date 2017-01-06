from netCDF4 import Dataset
import numpy as np
from windspharm.standard import VectorWind
from windspharm.tools import prep_data, recover_data, order_latdim
from scipy import interpolate
import datetime as datetime  # Python standard library datetime  module



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

    ilon[:] = [x*dtr for x in ilon]
    ilat[:] = [x*dtr for x in ilat]

    nlat = np.zeros_like(ilat)
    nlon = np.zeros_like(ilon)

    for index, (lon,lat) in enumerate(zip(ilon, ilat)):

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
        z=rad*np.sin(lat)


        # turn XY
        xn=x*cos_plon+y*np.sin(plon)
        yn=x*(-np.sin(plon))+y*cos_plon
        zn=z

        # tilt XZ
        xnn=xn*np.cos(plat)+zn*np.sin(plat)
        ynn=yn
        znn=xn*(-np.sin(plat))+zn*np.cos(plat)

        # convert back to polar coordinates
        rr=np.sqrt(xnn*xnn+ynn*ynn+znn*znn)
        nlat[index]=np.arcsin(znn/rr)*rtd


        nlon[index]=np.arctan2(xnn,ynn)*rtd

        # print "res----", nlat[index],nlon[index]

        return nlat,nlon




dset = "erain"
scheme = "UOM"

dt_nc = []

for yr in range(2005,2006):

    ftrk = "../data/trk%s/trk_%s/ntrkdat.%s.%d"%(scheme,dset,dset,yr+1)
    frad = "../output/rad%s.%s.%d.trk"%(scheme,dset,yr)
    print "ftrk =",ftrk
    print "frad =",frad

    # for iyr,cyr in enumerate(np.arange(yr,yr+1)):
    for iy,cyr in enumerate(np.arange(yr,yr+1)):
     print "loop", iy, cyr
     fnc  = "/Users/irudeva/work/DATA/%st/erain.mslp.%d.nc"%(dset,cyr)
     print "fnc  =",fnc

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
     #mslp0=mslp/100.
     if iy ==0:
         slp=np.tile(mslp,(3,1,1,1))
     slp[iy,:,:,:] = mslp

# create an array with 3 time series
     dt_nc = [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
           for t in time]


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
    date = np.zeros_like(lon,dtype = np.int)
    trktime  = np.zeros_like(date)
    it   = np.zeros_like(date)
    iyr  = np.zeros_like(date)
    print iyr.shape
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
        import numpy as np
        for n in range(0,nit):
            l = f.readline()
            columns = l.split()
            lon[ntrk-1,n]=float(columns[7])
            lat[ntrk-1,n]=float(columns[8])
            cslp[ntrk-1,n]=float(columns[9])
            date[ntrk-1,n]=columns[1]
            trktime[ntrk-1,n]=columns[2]

            iyr[ntrk-1,n] = -1
            for ind,year in enumerate(np.arange(yr-1,yr+1)):
                print 'trk yr: ',str(date[ntrk-1,n])[0:4]
                print 'year=',year
                if str(year) == str(date[ntrk-1,n])[0:4]:
                    print 'ok', ind, ntrk-1,n
                    print iyr.shape
                    iyr[ntrk-1,n] = ind
            if iyr[ntrk-1,n] == -1:
                print "!!!! Check years in trk file"
                quit()
            for ind,t in enumerate(dt_nc) :
                year = int(str(date[ntrk-1,n])[0:4])
                mon  = int(str(date[ntrk-1,n])[4:6])
                dat  = int(str(date[ntrk-1,n])[6:8])
                hr   = int(str(trktime[ntrk-1,n]/100))
                if t == datetime.datetime(year,mon,dat,hr):
                    it[ntrk-1,n] =ind
            if it[ntrk-1,n] == -1:
                print "!!!! Check time in trk file"
                quit()
            # print "current time", iyr, it

            #print lon[ntrk-1,n],lat[ntrk-1,n],cslp[ntrk-1,n]

        break  # for testing - first tack only


    f.close()
    print 'ftrk closed'

    #get radius




    #pole rotation
    for itrk in range(0,1): #  ntrk):
        for ip in range(0,1): #np[itrk]):

         print date[itrk,ip]
         print str(date[itrk,ip])[4:6]

        #  print 'lon = ', lon[itrk,ip]
        #  print 'lat = ', lat[itrk,ip]
        #  plon = [lon[itrk,ip]]
        #  plat = [lat[itrk,ip]]
         plon = [165]
         plat = [-11]
         nplat0 = plat[0]
         nplon0 = plon[0]
         print "place pole to cyc center (lat %2.f,lon %d)"%(plat[0], plon[0])
         nplat, nplon = polerot(plat[0],plon[0],[90],[nplon0])
        #  print 'lat = ', tmplat1
        #  print 'lon = ', tmplon1
        #  print 'nlat = ', nlat
        #  print 'nlon = ', nlon

         import numpy as np
         lonrange = np.arange(0., 360., 10.)
         latrange = np.arange(90., 69.5, -0.5)

         for ilon in lonrange :
             print "ilon=",ilon

             xnplat = nplat[0]
             xnplon = nplon[0]

             gridlat = np.copy(latrange)
             gridlon = np.zeros_like(gridlat)+ilon

             print 'gridlat=',gridlat
             print 'gridlon=',gridlon

             print "place pole back to NP"
             #  print "place pole to (lat %2.f,lon %d)"%(xnplat, xnplon)
             nlat, nlon = polerot(xnplat, xnplon,gridlat,gridlon)

             quit()

# quit()
