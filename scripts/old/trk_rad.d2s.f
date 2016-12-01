	Program Direct2sequencial
	
!to recover radii

      Implicit none

	character *100 file_dir,file_seq,fin,ftrk

	character mon*3
	integer ncycl,nline,year,month,date,time,num
	integer ax, ay,dat
	real cslp,fslp,rcount,rcycl(36)

	integer i,n

c--------reading shell and openning files------------------------------
      
	read(*,'(a100)') fin
	read(*,'(a100)') ftrk
	read(*,'(a100)') file_dir
	read(*,'(a100)') file_seq

	open(10,file=file_dir, action='read',				 
     &FORM='UNFORMATTED',ACCESS='DIRECT',RECL=40)
	open(20, file=file_seq,action='write',
     &status='unknown')
	open(11, file=fin,status='old',action='read',
     &FORM='UNFORMATTED',ACCESS='DIRECT',RECL=7)
	open(12, file=ftrk,action='read',
     &  FORM='UNFORMATTED',ACCESS='DIRECT',RECL=5)
     
       ncycl=1;print*,ncycl
9991   	read(11,rec=ncycl,err=999) n,nline,year,month,date,time,num
	print*, ncycl,nline,year,month,date,time,num
	 call monthtrans(month,mon)
	 if(time<12)then
	 write(20,'(i3,"-",a3,"-",i4," 0",i1,":00")') date,mon,year,time
	 write(*,'(i3,"-",a3,"-",i4," 0",i1,":00")') date,mon,year,time
	 else
	 write(20,'(i3,"-",a3,"-",i4,i3,":00")') date,mon,year,time	 
	 endif
	 write (20,1001) num

	do 1 i=1,num
	read(12,rec=nline-1+i)n,dat,ax,ay,cslp
	read(10,REC=nline-1+i) n,cslp,fslp,rcount,rcycl

1	write(20,1002) ax,ay,dat,cslp,fslp,rcount,rcycl

	ncycl=ncycl+1
	goto 9991
 
999	close(10);close(11);close(12);close(20)

1001  FORMAT (i4)
1002  FORMAT (3i5,39f8.2)

	stop
	end

c---------------------------------------------------------
	subroutine monthtrans(monnum,mon)
	character time*3
	integer month

	if(monnum==1)mon='Jan'
	if(monnum==2)mon='Feb'
	if(monnum==3)mon='Mar'
	if(monnum==4)mon='Apr'
	if(monnum==5)mon='May'
	if(monnum==6)mon='Jun'
	if(monnum==7)mon='Jul'
	if(monnum==8)mon='Aug'
	if(monnum==9)mon='Sep'
	if(monnum==10)mon='Oct'
	if(monnum==11)mon='Nov'
	if(monnum==12)mon='Dec'

	return
	end
