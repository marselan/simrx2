character*20 function int2str(k)
implicit none
integer :: k
write(int2str, *) k
int2str = adjustl(int2str)
end


subroutine state2det(proyNum)
implicit none
integer                         :: maxCount
integer                         :: proyNum


integer*4                       :: IR, IW
integer*4                       :: aux4
integer*8                       :: aux8
character*20                    :: int2str
integer, parameter :: imageHeight = 501
integer, parameter :: imageWidth = 501

integer*4, Dimension(imageHeight,imageWidth) :: Detector

IR = 10
IW = 11
Open(IR,File='state' // int2str(proyNum),Form='unformatted',Access='SEQUENTIAL')
Read(IR) aux8 ! skip first value (current photon count)
Read(IR) aux4 ! skip second value (ISEED1)
Read(IR) aux4 ! skip third value (ISEED2)
Read(IR) Detector

! flip up to down
Detector = Detector(imageHeight:1:-1, :)
Open(IW,File='detector' // int2str(ProyNum),Form='unformatted',Access='DIRECT',RECL=imageWidth*imageHeight * 4)
Write(IW, rec=1) Detector  

Close(Unit = IR)
Close(Unit = IW)

return
end


program main
	integer :: proyNum
	do proyNum = 1, 360
		call state2det(proyNum)
	end do	
end program

