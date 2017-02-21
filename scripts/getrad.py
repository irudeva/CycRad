from netCDF4 import Dataset
import numpy as np
from windspharm.standard import VectorWind
from windspharm.tools import prep_data, recover_data, order_latdim
from scipy import interpolate
import datetime as datetime  # Python standard library datetime  module
import time as ftime


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

    plat0 = np.copy(plat)
    plon0 = np.copy(plon)

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


        nlon[index]=np.arctan2(ynn,xnn)*rtd

        # print "res----", nlat[index],nlon[index]

    return nlat,nlon


def rotated_grid_transform(plat, plon, ilat, ilon, option):

    import numpy as np

    rad  =6371032 #in m
    pi   = np.pi
    dtr  = pi/180.
    rtd  = 180./pi
    # print 'rotated'
    # print plat, plon

    ilon[:] = [x*dtr for x in ilon]
    ilat[:] = [x*dtr for x in ilat]

    nlat = np.zeros_like(ilat)
    nlon = np.zeros_like(ilon)

    theta = plat-90 # Rotation around y-axis
    phi = plon #  Rotation around z-axis

    phi = phi*dtr
    theta = theta*dtr
    if option == 2: # Rotated -> Regular
     phi = -phi
     theta = -theta

    cos_phi = np.cos(phi)
    if phi == 90*dtr or phi == -90*dtr:
        cos_phi = 0
    sin_phi = np.sin(phi)
    if phi == 0*dtr or phi == 180*dtr or phi == -180*dtr:
        sin_phi = 0
    cos_theta = np.cos(theta)
    if theta == 90*dtr or theta == -90*dtr:
        cos_theta = 0
    sin_theta = np.sin(theta)
    if theta == 0*dtr :
        sin_theta = 0

    # print cos_theta,cos_phi, sin_phi,sin_theta

    for index, (lon,lat) in enumerate(zip(ilon, ilat)):

        cos_lon = np.cos(lon)
        if lon == 90*dtr or lon == -90*dtr:
            cos_lon = 0
        sin_lon = np.sin(lon)
        if lon == 0*dtr or lon == 180*dtr or lon == -180*dtr:
            sin_lon = 0
        cos_lat = np.cos(lat)
        if lat == 90*dtr or lat == -90*dtr:
            cos_lat = 0
        sin_lat = np.sin(lat)
        if lat == 0*dtr :
            sin_lat = 0

        x = cos_lat*cos_lon
        y = cos_lat*sin_lon
        z = sin_lat

        # print cos_lat,cos_lon, sin_lon


        if option == 1: # Regular -> Rotated

            x_new = cos_theta*cos_phi*x + cos_theta*sin_phi*y + sin_theta*z
            y_new = -sin_phi*x + cos_phi*y
            z_new = -sin_theta*cos_phi*x - sin_theta*sin_phi*y + cos_theta*z


        elif option == 2: # Rotated -> Regular

            # phi = -phi
            # theta = -theta

            x_new = cos_theta*cos_phi*x + sin_phi*y + sin_theta*cos_phi*z
            y_new = -cos_theta*sin_phi*x + cos_phi*y - sin_theta*sin_phi*z
            z_new = -sin_theta*x + cos_theta*z


        lon_new = np.arctan2(y_new,x_new) # % Convert cartesian back to spherical coordinates
        lat_new = np.arcsin(z_new)

        nlon[index] = lon_new*rtd  # % Convert radians back to degrees
        nlat[index] = lat_new*rtd

    # print 'end rotation'
    return nlat,nlon

# test!!!!

# Rotate Pole
# nplat, nplon = rotated_grid_transform(11,165,[88],[50],1)
# print nplat, nplon
# Reverse rotation
# nplat, nplon = rotated_grid_transform(11,165,[nplat],[nplon],2)
# print nplat, nplon
# quit()

# end test!!!

dset = "erain"
scheme = "UOM"

dt_nc = [ datetime.datetime(1900, 1, 1, 0) for i in range(1464)]

for yr in range(2005,2006):

    ftrk = "../data/trk%s/trk_%s/ntrkdat.%s.%d"%(scheme,dset,dset,yr)
    frad = "../output/rad%s.%s.%d.trk"%(scheme,dset,yr)
    print "ftrk =",ftrk
    print "frad =",frad
    fout = open(frad, 'wb')

    # for iyr,cyr in enumerate(np.arange(yr,yr+1)):
    for iy,cyr in enumerate(np.arange(yr-1,yr+2)):
    #  print "loop", iy, cyr
     fnc  = "/Users/Irina/work/DATA/%st/erain.mslp.%d.nc"%(dset,cyr)
     print "fnc  =",fnc

    # read netcdf
     print 'fnc reading...'
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

     # fix for leap years
     if iy == 0:
         slp=np.tile(mslp,(3,1,1,1))
     slp[iy,:time.size,:,:] = mslp

# create an array with 3 time series
     dt_nc[0:time.size] = [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
           for t in time]
     if iy == 0:
        dtall=np.tile(dt_nc,(3,1))

     dtall[iy,:] = dt_nc
    #  print len(dtall)
    #  print dtall[:,:10]

    print "\nSLP Interpolation"
    # data from south to north

    #for t in range(time.size) :
    # slpint = interpolate.interp2d(lons, lats, mslp[0,:,:], kind='cubic')
    # print slpint(0,30)

    print "end interpolation"


    #read trk

    print 'ftrk reading:'
    max_ntrk =  20000
    max_trklength = 200
    npnt  = np.zeros(max_ntrk,dtype = np.int)
    lon  = np.zeros((max_ntrk,max_trklength))
    lat  = np.zeros_like(lon)
    cslp = np.zeros_like(lon)
    date = np.zeros_like(lon,dtype = np.int)
    trktime  = np.zeros_like(date)
    it   = np.zeros_like(date)
    iyr  = np.zeros_like(date)
    print iyr.shape

    f = open(ftrk, 'r')

    for line in range(1,max_ntrk):
        emptyline = f.readline()
        header    = f.readline()
        emptyline = f.readline()
        smthline = f.readline()

        fout.write(emptyline)
        fout.write(header)
        fout.write(emptyline)
        fout.write(smthline)

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

        npnt[ntrk-1]=nit
        for n in range(0,nit):
            l = f.readline()
            fout.write(l)
            columns = l.split()
            lon[ntrk-1,n]=float(columns[7])
            lat[ntrk-1,n]=float(columns[8])
            cslp[ntrk-1,n]=float(columns[9])
            date[ntrk-1,n]=columns[1]
            trktime[ntrk-1,n]=columns[2]

            iyr[ntrk-1,n] = -1
            for ind,year in enumerate(np.arange(yr-1,yr+1)):
                if str(year) == str(date[ntrk-1,n])[0:4]:
                    iyr[ntrk-1,n] = ind
            if iyr[ntrk-1,n] == -1:
                print "!!!! Check years in trk file"
                quit()
            for ind,t in enumerate(dtall[iyr[ntrk-1,n],:]) :
                year = int(str(date[ntrk-1,n])[0:4])
                mon  = int(str(date[ntrk-1,n])[4:6])
                dat  = int(str(date[ntrk-1,n])[6:8])
                hr   = int(str(trktime[ntrk-1,n]/100))
                if t == datetime.datetime(year,mon,dat,hr):
                    it[ntrk-1,n] =ind
                    break
            if it[ntrk-1,n] == -1:
                print "!!!! Check time in trk file"
                quit()

            # start radius estimation
            plon = [lon[ntrk-1,n]]
            plat = [lat[ntrk-1,n]]
            #  plon = [162]
            #  plat = [33]
            # nplat0 = plat[0]
            # nplon0 = plon[0]

            # grid around cyclone center
            dlon = 10
            dlat = 0.5
            lonrange = np.arange(0., 360., dlon)
            latrange = np.arange(90., 64.5, -dlat)

            rslp     = np.zeros_like(latrange)
            dslp     = np.zeros_like(latrange[:-1])
            lslp     = np.zeros_like(lonrange)
            flslp    = np.copy(lslp)
            flat     = np.copy(lslp)
            flon     = np.copy(lslp)
            rad      = np.copy(lslp)

            nlon     = np.zeros((latrange.size,lonrange.size))
            nlat     = np.zeros_like(nlon)

            mr = 4  # min radius = mr*dlat (deg.lat)

            for i,ilon in enumerate(lonrange) :
                #  print "t=",n,"/",nit,"   ilon=",ilon

                 gridlat = np.copy(latrange)
                 gridlon = np.zeros_like(gridlat)+ilon

                 nlat[:,i], nlon[:,i] = rotated_grid_transform(plat[0],plon[0],gridlat,gridlon,2)

                 slpint = interpolate.interp2d(lons, lats, slp[iyr[ntrk-1,n],it[ntrk-1,n],:,:], kind='cubic')

                 for j,jlat in enumerate(latrange[:-1]) :
                     rslp[j] = slpint(nlon[j,i],nlat[j,i])[0]
                     dslp[j] = slpint(nlon[j+1,i],nlat[j+1,i]) - slpint(nlon[j,i],nlat[j,i])

                #  print rslp

                 if all(dslp[0:mr-1]) < 0. :
                     print "!!!! slp bug"
                     print dslp
                    #  ftime.sleep(20)
                     quit()

                 for j in range(mr,latrange.size-1):
                     if dslp[j] < 0 or j == latrange.size-2:
                         slp0 = slpint(nlon[0,i],nlat[0,i])
                         lslp[i] = slpint(nlon[j,i],nlat[j,i])[0]
                        #  fout.write(" {:>14}{:4.0f}{:>5}{:8.3f}{:>5} {:7.3f}{:>5} {:4.1f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format("angle=",ilon,"lon=",nlon[j+1,i],"lat=",nlat[j+1,i],"rad=",(j+1)*dlat,"cslp=",slp0[0],"fslp=",lslp[i]))
                        #  ftime.sleep(10)
                         break


            # find the last closed isobar (fslp)
            fslp = np.amin(lslp)
            cycdepth = fslp-slp0[0]

            # set rad = 0 for weak cyclones
            if cycdepth >= 1 :
                # find where last closed isobar fslp cross each radius
                gridlat = np.copy(latrange)
                for i,ilon in enumerate(lonrange) :
                    for j in range(mr,latrange.size-1):
                     slp2 = slpint(nlon[j+1,i],nlat[j+1,i])[0]
                     slp1 = slpint(nlon[j,i],nlat[j,i])[0]

                     if slp1 <= fslp and slp2 > fslp :
                        rad[i] = j*dlat + (fslp-slp1)/(slp2-slp1)*dlat
                        flat[i], flon[i] = rotated_grid_transform(plat[0],plon[0],[90-rad[i]],[ilon],2)
                        # print "flat=",flat, "flon=",flon
                        flslp[i] = slpint(flon[0],flat[0])[0]
                        #  fout.write(" {:>14}{:4.0f}{:>5}{:8.3f}{:>5} {:7.3f}{:>5} {:7.3f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format("angle=",ilon,"lon=",flon[0],"lat=",flat[0],"rad=",rad[i],"cslp=",slp0[0],"fslp=",flslp[i]))
                        break

                     if j == latrange.size-1:
                        if slpint(nlon[mr,i],nlat[mr,i])[0] > fslp :
                            rad[i] = mr*dlat
                            flat[i], flon[i] = rotated_grid_transform(plat[0],plon[0],rad[i],[ilon],2)
                            flslp[i] = slpint(flon[0],flat[0])[0]
                            #  fout.write(" {:>14}{:4.0f}{:>5}{:8.3f}{:>5} {:7.3f}{:>5} {:7.3f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format("angle=",ilon,"lon=",flon[0],"lat=",flat[0],"rad=",rad[i],"cslp=",slp0[0],"fslp=",flslp[i]))
                        elif slpint(nlon[j],nlat[j])[0] < fslp :
                            print "!!!! fslp bug"
                            print fslp, slpint(nlon[j],nlat[j])[0]
                            #  ftime.sleep(20)
                            quit()

                    # fout.write(" {:>14}{:4.0f}{:>5}{:8.3f}{:>5} {:7.3f}{:>5} {:7.3f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format
                    #           ("angle=",ilon,"lon=",flon[i],"lat=",flat[i],"rad=",rad[i],"cslp=",slp0[0],"fslp=",flslp[i]))


            else :
                for i,ilon in enumerate(lonrange) :
                    rad[i]   = 0
                    flat[i]  =  plat[0]
                    flon[i]  =  plon[0]
                    flslp[i]    = fslp
                    # fout.write(" {:>14}{:4.0f}{:>5}{:8.3f}{:>5} {:7.3f}{:>5} {:7.3f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format
                    #           ("angle=",ilon,"lon=",flon[i],"lat=",flat[i],"rad=",rad[i],"cslp=",slp0[0],"fslp=",flslp[i]))

            Area = 0
            for i,ilon in enumerate(lonrange) :
                Area = Area + rad[i]**2
            Area = Area/lonrange.size

            effrad = np.sqrt(Area)
            fout.write(" {:>14}{:7.3f}\n".format("effrad=",effrad))

            if Area > 0 :
                for i,ilon in enumerate(lonrange) :
                    fout.write(" {:>14}{:4.0f}{:>5}{:8.3f}{:>5} {:7.3f}{:>5} {:7.3f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format
                          ("angle=",ilon,"lon=",flon[i],"lat=",flat[i],"rad=",rad[i],"cslp=",slp0[0],"fslp=",flslp[i]))

            # end radius estimation
            # break  # for the first step of the track

        # break  # for testing - first tack only


    f.close()
    print 'ftrk closed'

    quit()
