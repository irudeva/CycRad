 	module radius
	real, parameter:: Rad=6371.032, sRad=(Rad**2)*3.141592654
        integer, parameter:: ne=36,me=36  
        real, parameter:: radstep=50. ! => R=ne*radstep=1800m
	real, parameter:: grad=111.195
	end module radius

      Program rad_calculating
!interpolates data from 181 grid on radii and calculates cyclone size

	use radius
      Implicit none

	real, parameter:: r=0.2,dr=0.1, rm=5.
c	real sin45
	
	character *70 file_trk, file_glo, file_out,file_area, maskdir
	integer lenmask

!coor
	real plat(ne*me),plon(ne*me),px(ne*me),py(ne*me)
		
      integer total, num
	character * 19 time
! total - number of cyclones per year
! num - number of points in one cyclone
	integer aix, aiy, ax,ay
	integer dat,press, widthx,widthy,x(181*181),y(181*181)
	real cslp
	real alon,alat,nlon(181*181),nlat(181*181),nx(181*181),ny(181*181)
	real glon(181*181),glat(181*181)

	real slp(181,181), slparea(181*181),slprad(ne*me)

!radius in directions
	real slpdir(me,ne), maxslp(me),minslp, rcycl(me)

	real s,rcount

	integer i, j, n,m, k,t 

c	sin45=sind(45.)
c----------making 181 coordinates--------------------------------------

c     create a lat/lon source grid
	k=0
	do 11 m=1,me
	do 11 n=1,ne  
	k=k+1  
      plat(k)=90.-(radstep/grad)*float(ne-n+1)
      plon(k)=(m-1)*10.
11	call coor181(plon(k),plat(k),px(k),py(k))
c--------reading shell and openning files------------------------------
      
	read(*,1111) file_trk
	read(*,1111) file_glo
	read(*,1111) file_out
	read(*,1111) maskdir
	lenmask=len_trim(maskdir)
	print*, file_trk,file_glo, file_out,maskdir

	open(10,file=file_trk, status='old')
	open(20, file=file_glo,status='old',
     &FORM='FORMATTED',ACCESS='DIRECT',RECL=181*7)
	open(30,file=file_out,status ='unknown')

c---------reading cyclone----------------------------------------------
c       read(10,*)	!for northen hemishere
c       read(10,'(10x,i10)') total
c	 write(*,*) 'total number of cyclones in NAA =', total
c	 write(30,'(a10,i10)') 'itog', total

c	do 22 t=1,total
c	write(*,*) t,'from',total
c	open (50, file='buf', status='unknown')

9991	 read (10,1000,end=9992) time
	 write (30,1000) time
	 read (10,1001) num
	 write (30,1001) num
	 write(*,*)time,num

 	do 221 i=1,num              ! every point of cyclone

   	read (10,1002) aix, aiy, dat, press
	cslp=float(press)/10.
	call coorint(aix,aiy,alat)
c---------reading data-------------------------------------------------
	 do 223 k=1,181
223	 read (20,1006,rec=181*(dat-1)+k ) slp(:,k)


c----open grid file----------------------------------------------------
	if(aiy<10)then
	 widthy=1
	elseif(aiy<100)then
	 widthy=2
	else
	 widthy=3
	endif


	if(aix<10)then
	 widthx=1
	elseif(aix<100)then
	 widthx=2
	else
	 widthx=3
	endif
	
	write(file_area,'(a<lenmask>,i<widthx>,a1,i<widthy>)')
     &maskdir,aix,'_',aiy
	open (40,file=file_area,status='old',form='unformatted',
     &action='read')

c	read(40,10031) ax,ay,alon,alat
c	if(ax/=aix.or.ay/=aiy) print*,'WRONG FILE TO OPEN!!!'
c	read(40,*) k;
	read(40) k
	read(40) ((x(j), y(j),glon(j),glat(j), nlon(j),nlat(j)),j=1,k)
	do 222 j=1,k;
c	read(40,10032) x(j), y(j),glon(j),glat(j), nlon(j),nlat(j)
	call coor181(nlon(j),nlat(j),nx(j),ny(j))
222	slparea(j)=slp(x(j),y(j))
c222	write(50, '(7f10.2)')
c     &glon(j),glat(j),nlon(j),nlat(j),nx(j),ny(j),slparea(j)
	close(40)


c****interpolation*****************************************************
	do 225 j=1,ne*me
225   call inter(k,slparea(1:k),nx(1:k),ny(1:k),dr,rm,
     &px(j),py(j),slprad(j),r)

c----radius in directions----------------------------------------------
	j=0
	do 224 m=1,me
	do 224 n=1,ne  
	j=j+1  
c224	write(50, *) plon(j),plat(j),slprad(j)
224	slpdir(m,n)=slprad(j)

!-----looking for maxslp in every direction----------------------------      
	do 331 m=1,me
	maxslp(m)=cslp
	do 331 n=ne,1,-1
331	 if(slpdir(m,n)>=maxslp(m))maxslp(m)=slpdir(m,n)

!-----looking for maxslp for closed isobars----------------------------      
	minslp=maxslp(1)
	do 332 m=2,me
332    if(maxslp(m).lt.minslp) minslp=maxslp(m)

!distance at every direction
	do 333 m=1,me
	rcycl(m)=float(ne)*radstep
	do 334 n=ne,1,-1
	if(n==ne) then
	 if(slpdir(m,n)>minslp)then
	 rcycl(m)=radstep*(minslp-cslp)/(slpdir(m,n)-cslp)
	 goto 333
	 endif
	endif
	 
	 if(slpdir(m,n)==minslp)then
	 rcycl(m)=radstep*(ne-n+1)
	 goto 333
	 elseif(slpdir(m,n)>minslp)then
c	 rcycl(m)=radstep*(ne-n+(minslp-cslp)/(slpdir(m,n)-cslp))
	 rcycl(m)=radstep*(ne-n+
     &(minslp-slpdir(m,n+1))/(slpdir(m,n)-slpdir(m,n+1)))
	 goto 333
334	 endif
333	continue

	s=0.
	do 335 m=1,me
c	rcycl(m)=rcycl(m)*sin45/sind(alat)
c335   s=s+2*sRad*(1-cos(rcycl(m)/Rad))/float(me)        !calculating area!!!!!!!
335   s=s+2*sRad*(1-cos(rcycl(m)/Rad))/float(me)        !calculating area!!!!!!!
	rcount = Rad*acos(1-(s/(2*sRad)))	

c335	write(50,*) (m-1)*10.,90-rcycl(m)/111.
c	write(50,*) minslp, rcount
c	close(50)
	if(any(rcycl==0.))then
	rcount=-999.99
	rcycl=-999.99
	endif
	write(30,1010) aix,aiy,dat,cslp,minslp,rcount,rcycl

c		go to 999
221	continue
c22	continue   
	goto 9991

9992    close(10)
        close(20)
        close(30)

1111  FORMAT (a70)	   
1000  FORMAT (a19)
1001  FORMAT (i4)
1002  FORMAT (i4,2i5,i6)
10031 FORMAT (2i5,2f7.2)
10032 FORMAT (2i4,4f10.2)
1006  FORMAT (181f7.1)

1010  FORMAT (3i5,<me+3>f8.2)

999   end


c----------------------------------------------------------------------
      subroutine coor181(lon, lat, nx, ny)
	
	IMPLICIT NONE
	real lon,lat, nx, ny, s
          
           s=90.*cosd(lat)
           nx=s*cosd(lon)
           ny=s*sind(lon)

      return
      end
c----------------------------------------------------------------------
      subroutine coorint(ax, ay, lat)
	
	IMPLICIT NONE
      real x, y, lon, lat
      real  ppp, gl
	integer ax,ay,err
	x=float(ax-91)
	y=float((-1)*ay+91)

 	  err=0

	 if (x.eq.0.) then
           if (y.gt.0.) then
		 lon=90.
	     elseif (y.eq.0.) then
		 lon=0.
           else
		 lon=270.
	     endif
		    lat=acosd(abs(y)/90.)
	 
	 else
           
		 gl = atand(y/x)

           if (x.lt.0.)   lon = 180.+gl

           if (x.ge.0..and.y.lt.0.) lon = 360.+gl

           if (x.ge.0..and.y.ge.0.) lon = gl


              ppp = x / 90. / cosd(lon)
              if (ppp.gt.999.) ppp = 1.

	  
	    if (abs(ppp).le.1.)then
           lat = acosd(ppp)
	    else
	     err = 1
	    endif

       endif

      return
      end

c**********************************************************************
      SUBROUTINE INTER(PASNUM,Q,LX,LY,DR,RM,XH,YH,QK,R)


c********************************************************

C     PASNUM - the number of points

C     Q      = function values in these points

C     LX, LY - coordinates of points

C     DR     = extention of radius

C     RM     = maximum of radius

C     XH, YH = coordinates of the point of interpolation

C     QK     = the result of interpolation

C     R       = initial radius

c********************************************************

      INTEGER PASNUM

      REAL Q(PASNUM),LX(PASNUM),LY(PASNUM)


      rw=r

3     continue

      NREAL=0

      A=0.

      B=0.
      DO 1 I=1,PASNUM


      IF(Q(I).LT.-999.0)GOTO 1

      RR=SQRT((XH-LX(I))*(XH-LX(i))+(YH-LY(I))*(YH-LY(I)))

      IF(RR.EQ.0.)THEN

      QK=Q(I)

      RETURN

      END IF

      IF(RR.LT.Rw)THEN

      NREAL=NREAL+1

      A=A+Q(I)/RR*RR

      B=B+1./RR*RR

      END IF

1     CONTINUE

      IF(NREAL.LT.4)THEN

      rw=rw+dr

      IF(RW.GE.RM)GOTO 2

      GOTO 3

      END IF

      QK=A/B

2     RETURN

      END

