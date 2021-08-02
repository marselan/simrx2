character*20 function int2str(k)
implicit none
integer :: k
write(int2str, *) k
int2str = adjustl(int2str)
end


subroutine readState(proyNum, currentPhotonCount, seed1, seed2, detector, height, width)
implicit none
integer           :: proyNum, width, height
integer*4, Dimension(height, width) :: detector
integer*4         :: seed1, seed2, fid
integer*8         :: currentPhotonCount
logical           :: fileExists
character*20 	  :: int2str

fid = 10
inquire(file='state' // int2str(ProyNum), exist=fileExists)
if( .NOT. fileExists ) then
	currentPhotonCount = 0
	seed1 = 8111117
	seed2 = 9111117
	detector = 0
	open(fid,File='state'// int2str(ProyNum),Form='unformatted',Access='SEQUENTIAL')
	write(fid) currentPhotonCount
	write(fid) seed1
	write(fid) seed2
	write(fid) detector 
	close(unit = fid)
else 
	open(fid,File='state'// int2str(ProyNum),Form='unformatted',Access='SEQUENTIAL')
	read(fid) currentPhotonCount
	read(fid) seed1
	read(fid) seed2
	read(fid) detector
	close(unit = fid)
end if
end

subroutine saveState(proyNum, currentPhotonCount, seed1, seed2, detector, height, width)
implicit none
integer           :: proyNum, width, height
integer*4, Dimension(height, width) :: detector
integer*4         :: seed1, seed2, fid
integer*8         :: currentPhotonCount
character*20 	  :: int2str

fid = 10
open(fid,File='state'// int2str(ProyNum),Form='unformatted',Access='SEQUENTIAL')
write(fid) currentPhotonCount
write(fid) seed1
write(fid) seed2
write(fid) detector 
close(unit = fid)
print *, 'State saved for proy:[', proyNum, '] currentPhotonCount:[', currentPhotonCount, '] seed1: [', seed1, '] seed2: [', seed2
end

subroutine simulate(ProyNum, NProy)
implicit none
integer                       :: ProyNum, NProy
character*20 :: int2str
Common/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
real*8            :: E,X,Y,Z,U,V,W,WGHT
integer*4         :: KPAR,IBODY,MAT,ILB

integer*4, Parameter :: MAXMAT = 10

integer*4, Parameter :: ELECTRON = 1
integer*4, Parameter :: PHOTON   = 2
integer*4, Parameter :: POSITRON = 3

real*8, Parameter :: CONIC_BEAM = 0d0
real*8, Parameter :: PARALLEL_BEAM = 1d0

Common/STATECOMMON/saveStateEvery
real*8            :: saveStateEvery

Common/PROCESSCOMMON/shutdownAfter,continueExecution
integer           :: shutdownAfter, continueExecution

Common/SCANNEDARCCOMMON/scannedAngle
Real*8            :: scannedAngle

Common/DETECTORCOMMON/detectorWidth,detectorHeight,detectorDistance
real*8            :: detectorWidth, detectorHeight, detectorDistance

Common/IMAGECOMMON/imageWidth,imageHeight,bitResolution
integer           :: imageWidth, imageHeight, bitResolution

Common/BEAMCOMMON/beamType,halfAngle,beamDistance,beamWidth,beamHeight,npmax,emax
real*8            :: beamType,halfAngle,beamDistance,beamWidth,beamHeight,npmax,emax

Common/GEOMETRYFILECOMMON/geometryFile
character*20      :: geometryFile

Common/MATERIALCOMMON/nmat,pmfile
integer                             :: nmat
character(len=100), dimension(10)   :: pmfile


Common/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),WCR(MAXMAT)
Real*8            :: EABS,C1,C2,WCC,WCR

Common/DSMAXCOMMON/dsmax
Real*8            :: dsmax(0:MAXMAT)

COMMON/STOKES/SP1,SP2,SP3,IPOL
double precision  :: SP1,SP2,SP3
integer*4         :: IPOL
integer*8           :: PhotonCount,NP,currentPhotonCount, isaveStateEvery

Common/RSEED/ISEED1,ISEED2
integer*4         :: ISEED1,ISEED2

integer*4         :: IRD, IWR, INFO, sfid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision            :: DS
double precision            :: DSEF,DE
integer*4                   :: NCROSS,ICOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer*4, parameter        :: NPINP = 0
double precision            :: PARINP(NPINP) 
integer*4                   :: NMATI4, NBOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8                      :: halfDetectorWidth
real*8                      :: halfDetectorHeight
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer*4, Dimension(imageHeight, imageWidth) :: Detector
integer*4                       :: MaxCount
double precision                :: DeltaAngle  
integer                         :: xi , yj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8, parameter :: Pi = 3.1415926535897932, DosPi = 2d0*Pi
real*8, parameter :: Deg2Rad = Pi/180d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8            :: Phi, CosTita, Suv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer*4         :: LEFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision                :: RAND
external RAND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PhotonCount = npmax
call readState(ProyNum, currentPhotonCount, ISEED1, ISEED2, Detector, imageHeight, imageWidth)
if (currentPhotonCount >= PhotonCount) then
	print *, 'Projection [', ProyNum, '] already simulated'
	return
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
isaveStateEvery = saveStateEvery

DSMAX(0) = 1d30 ! max free path in void region

WGHT = 1d0

IWR = 10
INFO = 0

NMATI4 = nmat

halfDetectorWidth = detectorWidth / 2d0
halfDetectorHeight = detectorHeight / 2d0

Open(Unit = IWR, File = 'report.txt')
call PEINIT(EMAX, NMATI4, IWR, INFO, PMFILE) 
Close(Unit = IWR)

IRD = 10
IWR = 11
Open(Unit = IRD, File = GeometryFile)
Open(Unit = IWR, File = 'geometry.rep')
call GEOMIN(PARINP,NPINP,NMATI4,NBOD,IRD,IWR)
Close(Unit = IRD)
Close(Unit = IWR)

deltaAngle = (scannedAngle / NProy) * (ProyNum - 1) * Deg2Rad  

! rotate sample	
call SampleRS(0d0, deltaAngle, 0d0, 0d0, 0d0, 0d0) 

print *, 'Simulating Projection: [', ProyNum, ']'	

currentPhotonCount = currentPhotonCount + 1

do NP = currentPhotonCount, PhotonCount
        Print *, 'tick', ProyNum, NP
	E = EMAX
	KPAR = 2
	if (beamtype == CONIC_BEAM) then

		! haz conico

		X = 0d0
		Y = 0d0
		Z = -beamDistance

		Phi = DosPi*Rand(1d0)

		CosTita = cos(halfAngle*Deg2Rad);
		W = rand(1d0) * (1d0 - CosTita) + CosTita !genera un nro uniforme entre [cos(tita), 1]
		Suv = Sqrt(1d0 - W*W)
		U = Suv*Cos(Phi)
		V = Suv*Sin(Phi)
	else
		! haz paralelo

		X = (rand(1d0)-0.5)*beamWidth
	      	Y = (rand(1d0)-0.5)*beamHeight
	      	Z = -beamDistance
		U = 0d0
		V = 0d0
		W = 1d0
	end if

	ILB = (/ 1, 0, 0, 0, 0 /)
	! IBODY and MAT are set by function LOCATE
	call LOCATE
	call CLEANS
	call START
	do
		call LOCATE
		if (MAT.EQ.0) then
			DS = DSMAX(0) ! void region
		else
			call JUMP(DSMAX(MAT), DS)
		endif

		call STEP(DS, DSEF, NCROSS)
		if (MAT.EQ.0) then
			! process secondary particles, if any
			call SECPAR(LEFT)
			if (LEFT.GT.0) then
				call START
				CYCLE
			else
				EXIT
			end if
		end if
		! NCROSS greater than 0 means the particle crossed an interface
		if (NCROSS.GT.0) then
			! IBODY equals 1 means that the particle hit the detector
			if (IBODY.EQ.1) then
                                xi= Nint((X+halfDetectorWidth)*(real(imageWidth)/detectorWidth))
                                if (xi.EQ.0) then
                                    xi = xi + 1
                                end if
                                if (xi.GT.imageWidth) then
                                    xi = imageWidth
                                end if
      				yj= Nint((Y+halfDetectorHeight)*(real(imageHeight)/detectorHeight)) + 1
                                if (yj.EQ.0) then
                                    yj = yj + 1
                                end if
                                if (yj.GT.imageHeight) then
                                    yj = imageHeight
                                end if
                                Detector(yj,xi) = Detector(yj,xi) + 1
				! process secondary particles, if any
	  			call SECPAR(LEFT)
				if (LEFT.GT.0) then
					call START
					CYCLE
				else
					EXIT
				end if
      			else
				call START
				CYCLE
			end if
		end if
		call KNOCK(DE, ICOL)
		if (E.LT.EABS(KPAR,MAT)) then
			! process seconday particles, if any		
			call SECPAR(LEFT)
			if (LEFT.GT.0) then
				call START
				CYCLE
			else
				EXIT
			end if
		else
			CYCLE
		end if
	end do
	! save state		
	if ( (mod(NP, iSaveStateEvery) == 0) .OR. continueExecution == 0) then
		call saveState(ProyNum, NP, ISEED1, ISEED2, Detector, imageHeight, imageWidth) 
		if (continueExecution == 0) then
			return
		end if
	end if
end do

MaxCount = maxval(Detector)
Open(IWR,File='maxvalue' // int2str(ProyNum),Form='unformatted',Access='DIRECT',RECL=4)
Write(IWR, rec=1) MaxCount  
Close(Unit = IWR)
print *, 'Writing MaxValue [', MaxCount, '] for Projection: [', ProyNum, ']'

return       
end

integer function getmax(Nproy)
implicit none
integer                         :: NProy, proyNum
integer*4                       :: IR
integer*4                       :: maxValue, value
character*20                    :: int2str

maxValue = 0
IR = 10
do proyNum = 1, Nproy
	Open(IR,File='maxvalue' // int2str(proyNum),Form='unformatted',Access='DIRECT',RECL=4)
	Read(IR, rec=1) value
	maxValue = max(maxValue, value)
	Close(Unit = IR, STATUS = 'DELETE') 
	print *, 'Read MaxValue [', value, '] for Projection: [', proyNum, ']'
end do
getmax = maxValue
end

subroutine det2img(maxCount, proyNum)
implicit none
integer                         :: maxCount
integer                         :: proyNum

Common/IMAGECOMMON/imageWidth,imageHeight,bitResolution
integer           :: imageWidth, imageHeight, bitResolution

integer*4                       :: IR, IW, ID
integer*4                       :: aux4
integer*8                       :: aux8
character*20                    :: int2str
integer*4, Dimension(imageHeight,imageWidth) :: Detector
integer*2, Dimension(imageHeight,imageWidth) :: Image
integer*2, Dimension(imageWidth,imageHeight) :: TImage

IR = 10
IW = 11
ID = 12
Open(IR,File='state' // int2str(proyNum),Form='unformatted',Access='SEQUENTIAL')
Read(IR) aux8 ! skip first value (current photon count)
Read(IR) aux4 ! skip second value (ISEED1)
Read(IR) aux4 ! skip third value (ISEED2)
Read(IR) Detector

! flip up to down
Detector = Detector(imageHeight:1:-1, :)
! correct contrast
Image = ( 1e0 - (maxCount - Detector) / real(maxCount) ) * (2**BitResolution - 1)
! transpose and save image
TImage = TRANSPOSE(Image)
Open(IW,File='image' // int2str(ProyNum),Form='unformatted',Access='DIRECT',RECL=imageWidth*imageHeight * 2)
Write(IW, rec=1) TImage  

Open(ID,File='detector' // int2str(ProyNum),Form='unformatted',Access='DIRECT',RECL=imageWidth*imageHeight * 4)
Write(ID, rec=1) Detector  

Close(Unit = IR)
Close(Unit = IW)
Close(Unit = ID)

return
end

subroutine printParams()
implicit none

integer              :: i
integer*4, Parameter :: MAXMAT = 10
integer*4, Parameter :: ELECTRON = 1
integer*4, Parameter :: PHOTON   = 2
integer*4, Parameter :: POSITRON = 3

real*8, Parameter :: CONIC_BEAM = 0d0
real*8, Parameter :: PARALLEL_BEAM = 1d0

Common/STATECOMMON/saveStateEvery
real*8           :: saveStateEvery

Common/SCANNEDARCCOMMON/scannedAngle
Real*8            :: scannedAngle

Common/DETECTORCOMMON/detectorWidth,detectorHeight,detectorDistance
real*8            :: detectorWidth, detectorHeight, detectorDistance

Common/IMAGECOMMON/imageWidth,imageHeight,bitResolution
integer           :: imageWidth, imageHeight, bitResolution

Common/BEAMCOMMON/beamType,halfAngle,beamDistance,beamWidth,beamHeight,npmax,emax
real*8            :: beamType,halfAngle,beamDistance,beamWidth,beamHeight,npmax,emax

Common/GEOMETRYFILECOMMON/geometryFile
character*20      :: geometryFile

Common/MATERIALCOMMON/nmat,pmfile
integer                             :: nmat
character(len=100), dimension(10)   :: pmfile
Common/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),WCR(MAXMAT)
Real*8            :: EABS,C1,C2,WCC,WCR

Common/DSMAXCOMMON/dsmax
Real*8            :: dsmax(0:MAXMAT)

Print *, 'Save state every: ', saveStateEvery, ' photons'
Print *, 'Scanned arc: ', scannedAngle
Print *, 'Detector width (cm): ', detectorWidth
Print *, 'Detector height (cm): ', detectorHeight
Print *, 'Detector distance from beam (cm): ', detectorDistance
Print *, 'Image Width: ', imageWidth
Print *, 'Image Height: ', imageHeight
Print *, 'Bit Resolution: ', bitResolution
if (beamType == CONIC_BEAM ) then
Print *, 'Beam geometry: [CONIC] - Half angle: [', halfAngle, '] - DFO: [', beamDistance, ']'
else
Print *, 'Beam geometry: [PARALLEL] - Width: [', beamWidth, '] - Height: [', beamHeight, '] - DFO: [', beamDistance, ']'
end if
Print *, 'Number of photons (npmax): ', npmax
Print *, 'Energy of (primary) photons (emax) :', emax
Print *, 'Geometry file: [', geometryFile, ']'
Print *, ''
Print *, 'MATERIALS (total: ', nmat, ')'
Print *, ''



do i = 1, nmat
	Print *, 'Material ', i
	Print '(A2A100A1)', ' [', pmfile(i),']'
	Print *, 'EABS Electrons: ', EABS(ELECTRON, i)
	Print *, 'EABS Photons: ', EABS(PHOTON, i)
	Print *, 'EABS Positrons: ', EABS(POSITRON, i)
	Print *, 'C1: ', C1(i)
	Print *, 'C2: ', C2(i)
	Print *, 'WCC: ', WCC(i)
	Print *, 'WCR: ', WCR(i)
	Print *, 'dSmax: ', dsmax(i)
	Print *, '' 
end do
return
end

