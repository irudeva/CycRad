	Program Trk

!Create a file with selected tracks

      Implicit none

	character *100 file_trk,file_outnum,file_out

	character monch*3
	integer year, date,time,month,num,ncycl
	integer ax, ay,dat
	real cslp

	integer i,n

c--------reading shell and openning files------------------------------
      
	read(*,'(a100)') file_trk
	read(*,'(a100)') file_outnum
	read(*,'(a100)') file_out
	print*, file_trk,file_out

	open(10,file=file_trk, status='old',action='read')
c	read(10,*)
	open(20, file=file_outnum,action='write',
     &FORM='UNFORMATTED',ACCESS='DIRECT',RECL=7)
	open(21, file=file_out,action='write',
     &FORM='UNFORMATTED',ACCESS='DIRECT',RECL=5)

c---------reading cyclone----------------------------------------------
		ncycl=0;n=0

9991	 read (10,1000,end=999) date,monch,year,time
	 call monthnumsubr(monch,month)
	 ncycl=ncycl+1
	 read (10,1001) num
	 write(20,REC=ncycl) ncycl,n+1,year,month,date,time,num

 	do 221 i=1,num              ! every point of cyclone
	n=n+1
  	read (10,1002) ax,ay,dat,cslp
221	write(21,REC=n) ncycl,dat,ax,ay,cslp

	goto 9991
999	close (10);close(20);close(21)
	
c	open(20, file=file_outnum,action='read',
c     &FORM='UNFORMATTED',ACCESS='DIRECT',RECL=7)
c	read(20,rec=2567) ncycl,n,year,month,date,time,num
c	print*, ncycl,n,year,month,date,time,num

c	open(21, file=file_out,action='read',
c     &FORM='UNFORMATTED',ACCESS='DIRECT',RECL=5)
c	read(21,rec=n) ncycl,dat,ax,ay,cslp
c	print*, ncycl,dat,ax,ay,cslp

1000  FORMAT (1x,i2,1x,a3,1x,i4,1x,i2)
1001  FORMAT (i4)
1002  FORMAT (i4,2i5,f6.1)

	stop
	end

c---------------------------------------------------------
	subroutine monthnumsubr(time,month)
	character time*3
	integer month

	if(time.eq.'Jan')month=1
	if(time.eq.'Feb')month=2
	if(time.eq.'Mar')month=3
	if(time.eq.'Apr')month=4
	if(time.eq.'May')month=5
	if(time.eq.'Jun')month=6
	if(time.eq.'Jul')month=7
	if(time.eq.'Aug')month=8
	if(time.eq.'Sep')month=9
	if(time.eq.'Oct')month=10
	if(time.eq.'Nov')month=11
	if(time.eq.'Dec')month=12

	return
	end
c---------------------------------------------------------
	subroutine seasonsubr(time,sea)
	character time*3
	integer sea

	if(time.eq.'Jan'.or.time.eq.'Feb'.or.time.eq.'Mar')sea=1
	if(time.eq.'Apr'.or.time.eq.'May'.or.time.eq.'Jun')sea=2
	if(time.eq.'Jul'.or.time.eq.'Aug'.or.time.eq.'Sep')sea=3
	if(time.eq.'Oct'.or.time.eq.'Nov'.or.time.eq.'Dec')sea=4

	return
	end
