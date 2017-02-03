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

    ilon[:] = [x*dtr for x in ilon]
    ilat[:] = [x*dtr for x in ilat]

    nlat = np.zeros_like(ilat)
    nlon = np.zeros_like(ilon)

    theta = -plat-90 # Rotation around y-axis
    phi = -plon #  Rotation around z-axis

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

    print cos_theta,cos_phi, sin_phi,sin_theta

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

        print cos_lat,cos_lon, sin_lon


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

    return nlat,nlon



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
     fnc  = "/Users/irudeva/work/DATA/%st/erain.mslp.%d.nc"%(dset,cyr)
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
    max_ntrk =  2 #20000
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
            nplat0 = plat[0]
            nplon0 = plon[0]
            # print "place pole to cyc center (lat %2.f,lon %d)"%(plat[0], plon[0])
            #nplat, nplon = polerot(plat[0],plon[0],[90],[nplon0])
            print "new pole", plat[0],plon[0]

            nplat, nplon = rotated_grid_transform(plat[0],plon[0],[90],[nplon0],1)
            print plat[0],plon[0]
            print "north pole bacomes", nplat, nplon
            nplat, nplon = rotated_grid_transform(plat[0],plon[0],[11.78],[180],2)
            print nplat, nplon
            nplat, nplon = rotated_grid_transform(plat[0],plon[0],[-90],[nplon0],1)
            print "south pole becomes",nplat, nplon
            nplat, nplon = rotated_grid_transform(plat[0],plon[0],[-11.78],[0],2)
            print nplat, nplon
            nplat, nplon = rotated_grid_transform(plat[0],plon[0],[-11.78],[90],2)
            print nplat, nplon
            # nplat, nplon = rotated_grid_transform(plat[0],plon[0],nplat, nplon,2)
            # print "back",nplat, nplon
            #
            #
            # nplat, nplon = rotated_grid_transform(plat[0],plon[0],[-90],[nplon0],1)
            # print nplat, nplon
            # nplat, nplon = rotated_grid_transform(plat[0],plon[0],nplat, nplon,2)
            # print "back",nplat, nplon
            # nplat, nplon = rotated_grid_transform(plat[0],plon[0],[-90],[nplon0],1)
            # nplat, nplon = rotated_grid_transform(nplat, nplon ,[90],[nplon0],1 )
            # print "back tweak",nplat, nplon
            # nplat, nplon = rotated_grid_transform(plat[0],plon[0],[11.78],[200],1)
            # print nplat, nplon
            # nplat, nplon = rotated_grid_transform(plat[0],plon[0],nplat, nplon,2)
            # print "back",nplat, nplon
            quit()

            # grid around cyclone center
            dlon = 10
            dlat = 0.5
            lonrange = np.arange(0., 360., dlon)
            latrange = np.arange(90., 69., -dlat)

            dslp     = np.zeros_like(latrange[:-1])
            lslp     = np.zeros_like(lonrange)
            flslp    = np.copy(lslp)

            mr = 4  # min radius = mr*dlat (deg.lat)

            for i,ilon in enumerate(lonrange) :
                 print "t=",n,"/",nit,"   ilon=",ilon

                 xnplat = nplat[0]
                 xnplon = nplon[0]

                 gridlat = np.copy(latrange)
                 gridlon = np.zeros_like(gridlat)+ilon

                #  print "place pole back to NP"
                 print gridlat
                 print gridlon
                 print plat[0],plon[0]
                 nlat, nlon = rotated_grid_transform(plat[0],plon[0],gridlat,gridlon,2)
                 print nlat
                 print nlon
                 print plat[0],plon[0]
                 quit()
                 #nlat, nlon = polerot(xnplat, xnplon,gridlat,gridlon)
                 #nlon = nlon+90+plon

                 slpint = interpolate.interp2d(lons, lats, slp[iyr[ntrk-1,n],it[ntrk-1,n],:,:], kind='cubic')

                 for j,jlat in enumerate(latrange[:-1]) :
                     dslp[j] = slpint(nlon[j+1],nlat[j+1]) - slpint(nlon[j],nlat[j])

                 if all(dslp[0:mr-1]) < 0. :
                     print "!!!! slp bug"
                     print dslp
                    #  ftime.sleep(20)
                     quit()

                 for j in range(mr,latrange.size-1):
                     if dslp[j] < 0 or j == latrange.size-2:
                         slp0 = slpint(nlon[0],nlat[0])
                         lslp[i] = slpint(nlon[j],nlat[j])[0]
                         fout.write(" {:>14}{:4.0f}{:>5}{:8.3f}{:>5} {:7.3f}{:>5} {:4.1f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format("angle=",ilon,"lon=",nlon[j+1],"lat=",nlat[j+1],"rad=",(j+1)*dlat,"cslp=",slp0[0],"fslp=",lslp[i]))
                        #  ftime.sleep(10)
                         break


            # find the last closed isobar (fslp)
            print "np.amin=", np.amin(lslp)
            print lslp
            fslp = np.amin(lslp)
            print fslp
            print fslp-slp0[0]  # add this critarium

         #  !!! skip if fslp-slp0 < 1.

            # find where last closed isobar fslp cross each radius
            gridlat = np.copy(latrange)
            for ilon in lonrange:
                for j in range(mr,latrange.size):
                 slp2 = slpint(nlon[j],nlat[j])[0]
                 slp1 = slpint(nlon[j-1],nlat[j-1])[0]

                 print ilon,j, fslp, slp1, slp2
                 if slp1 <= fslp and slp2 > fslp :
                     tlat = 90. - (j-1)*dlat - (fslp-slp1)/(slp2-slp1)*dlat
                     print ilon,j, 90. - (j-1)*dlat, gridlat[j-1]
                     flat, flon = polerot(xnplat, xnplon,[tlat,90],[ilon,ilon])
                    #  print "xnplat=",xnplat, "xnplon=",xnplon
                     print "flat=",flat, "flon=",flon
                     flslp[i] = slpint(flon[0],flat[0])[0]
                     fout.write(" {:>14}{:4.0f}{:>5}{:8.3f}{:>5} {:7.3f}{:>5} {:7.3f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format("angle=",ilon,"lon=",flon[0],"lat=",flat[0],"rad=",90.-tlat,"cslp=",slp0[0],"fslp=",flslp[i]))
                     break

                 if j == latrange.size-1 and slpint(nlon[mr-1],nlat[mr-1])[0] > fslp :
                     tlat = 90.-mr*dlat
                     flat, flon = polerot(xnplat, xnplon,[tlat,90],[ilon,ilon])
                     flslp[i] = slpint(flon[0],flat[0])[0]
                     fout.write(" {:>14}{:4.0f}{:>5}{:8.3f}{:>5} {:7.3f}{:>5} {:7.3f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format("angle=",ilon,"lon=",flon[0],"lat=",flat[0],"rad=",90.-tlat,"cslp=",slp0[0],"fslp=",flslp[i]))


            # end radius estimation

            break # for the first time step of a track
        break  # for testing - first tack only


    f.close()
    print 'ftrk closed'

    quit()

    #get radius




    #pole rotation
    for itrk in range(0,ntrk):
        fout.write("Track = %d\n"%itrk)
        for ip in range(npnt[itrk]):
             fout.write("t %d out of %d\n"%(ip+1,npnt[itrk]))

             print "start ip=",ip
             print date[itrk,ip]

            #  print 'lon = ', lon[itrk,ip]
            #  print 'lat = ', lat[itrk,ip]
             plon = [lon[itrk,ip]]
             plat = [lat[itrk,ip]]
            #  plon = [162]
            #  plat = [33]
             nplat0 = plat[0]
             nplon0 = plon[0]
             print "place pole to cyc center (lat %2.f,lon %d)"%(plat[0], plon[0])
             nplat, nplon = polerot(plat[0],plon[0],[90],[nplon0])
            #  print 'lat = ', tmplat1
            #  print 'lon = ', tmplon1
            #  print 'nlat = ', nlat
            #  print 'nlon = ', nlon

             lonrange = np.arange(0., 360., 10.)
             latrange = np.arange(90., 69., -0.5)
             dslp     = np.zeros_like(latrange[1:-1])

             for ilon in lonrange :
                 print "ilon=",ilon
                 fout.write("ilon=%d\n"%ilon)

                 xnplat = nplat[0]
                 xnplon = nplon[0]

                 gridlat = np.copy(latrange)
                 gridlon = np.zeros_like(gridlat)+ilon

                 print "place pole back to NP"
                 nlat, nlon = polerot(xnplat, xnplon,gridlat,gridlon)
                 nlon = nlon+90+plon

                #  print 'nlat=',nlat
                #  print 'nlon=',nlon

                 slpint = interpolate.interp2d(lons, lats, slp[iyr[itrk,ip],it[itrk,ip],:,:], kind='cubic')

                 for i,ilat in enumerate(latrange[:-2]) :
                     print i,nlat[0],nlon[0]
                     print slpint(nlon[i],nlat[i])
                     print i+2,nlat[2],nlon[2]
                     print slpint(nlon[i+2],nlat[i+2])
                     print 'cslp=',cslp[itrk,ip]
                    #  print slpint(nlon[i+2],nlat[i+2]) - slpint(nlon[i],nlat[i])
                     dslp[i] = slpint(nlon[i+2],nlat[i+2]) - slpint(nlon[i],nlat[i])
                    #  print  dslp[i]

                    #  print slpint(0,30)
                    #  print lats[60]
                    #  print slp[iyr[ntrk-1,n],it[ntrk-1,n],60,0]

                    #  print lats[57]
                    #  print lons[162]
                    #  print slp[iyr[itrk,ip],it[itrk,ip],57,162]
                    #  print slpint(162,33)
                     #
                     print iyr[itrk,ip],it[itrk,ip]
                    #  break

                 if all(dslp[0:4]) < 0. :
                     print "!!!! slp bug"
                     print dslp
                    #  ftime.sleep(20)
                     quit()

                 for i in range(5,latrange.size-2):
                     if dslp[i] < 0:
                         fout.write(" %d %d  %f %f\n"%(i,latrange[i+1],slpint(nlon[i+1],nlat[i+1]),slpint(nlon[0],nlat[0])))
                         fout.write(" %d %d\n"%(nlon[0],nlat[0]))
                         fout.write(" %d %d\n"%(nlon[i+1],nlat[i+1]))
                        #  print "rjad:",i,latrange[i+1],[slpint(nlon[j],nlat[j]) for j in range(latrange.size-2)]
                        #  print "dslp:",[dslp[j] for j in range(i+1)]
                        #  ftime.sleep(10)
                         break



                #  break  # loop over lonrange
            #  ftime.sleep(5)

        # break
        # ftime.sleep(5)


                #  quit()

# quit()
