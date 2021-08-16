character*20 function int2str(k)
	implicit none
	integer :: k
	write(int2str, *) k
	int2str = adjustl(int2str)
end

subroutine stateFileName(workerId, name)
	implicit none
	integer           :: proyNum, workerId
	character*25 	  :: name
	character*20 	  :: int2str
    name = 'state' // int2str(workerId)
end


subroutine readState(workerId, currentPhotonCount, seed1, seed2, detector, height, width)
	implicit none
	integer           :: workerId, width, height
	integer*4, Dimension(height, width) :: detector
	integer*4         :: seed1, seed2, fid
	integer*8         :: currentPhotonCount
	logical           :: fileExists
	character*25 	  :: stateName

	call stateFileName(workerId, stateName)

	fid = 10
	inquire(file=stateName, exist=fileExists)
	if( .NOT. fileExists ) then
		currentPhotonCount = 0
		seed1 = 8111117
		seed2 = 9111117
		detector = 0
		open(fid,File=stateName, Form='unformatted', Access='SEQUENTIAL')
		write(fid) currentPhotonCount
		write(fid) seed1
		write(fid) seed2
		write(fid) detector 
		close(unit = fid)
	else 
		open(fid,File=stateName, Form='unformatted', Access='SEQUENTIAL')
		read(fid) currentPhotonCount
		read(fid) seed1
		read(fid) seed2
		read(fid) detector
		close(unit = fid)
	end if
end

subroutine saveState(workerId, currentPhotonCount, seed1, seed2, detector, height, width)
	implicit none
	integer           :: workerId, width, height
	integer*4, Dimension(height, width) :: detector
	integer*4         :: seed1, seed2, fid
	integer*8         :: currentPhotonCount
	character*25 	  :: stateName

	call stateFileName(workerId, stateName)

	fid = 10
	open(fid,File=stateName, Form='unformatted', Access='SEQUENTIAL')
	write(fid) currentPhotonCount
	write(fid) seed1
	write(fid) seed2
	write(fid) detector 
	close(unit = fid)
end

subroutine checkDetectorHit(Detector)
	implicit none
	Common/DETECTORCOMMON/detectorWidth,detectorHeight,detectorDistance
	real*8            :: detectorWidth, detectorHeight, detectorDistance

	Common/IMAGECOMMON/imageWidth,imageHeight,bitResolution
	integer           :: imageWidth, imageHeight, bitResolution

	Common/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
	real*8            :: E,X,Y,Z,U,V,W,WGHT
	integer*4         :: KPAR,IBODY,MAT,ILB

	real*8                      :: halfDetectorWidth
	real*8                      :: halfDetectorHeight
	integer                     :: xi , yj
	integer*4, Dimension(imageHeight, imageWidth) :: Detector

	halfDetectorWidth = detectorWidth / 2d0
	halfDetectorHeight = detectorHeight / 2d0

	xi = Nint((X+halfDetectorWidth)*(real(imageWidth)/detectorWidth))
	if (xi.EQ.0) then
		xi = xi + 1
	end if
	if (xi.GT.imageWidth) then
		xi = imageWidth
	end if
	yj = Nint((Y+halfDetectorHeight)*(real(imageHeight)/detectorHeight)) + 1
	if (yj.EQ.0) then
		yj = yj + 1
	end if
	if (yj.GT.imageHeight) then
		yj = imageHeight
	end if
	Detector(yj,xi) = Detector(yj,xi) + 1
end

subroutine simulate(NProy, ProyNum, workerId, photonCount)
	implicit none
	integer                       :: NProy, ProyNum, workerId, photonCount
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
	integer*8           :: NP,currentPhotonCount

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

	print *, 'Simulating Projection: [', ProyNum, '] workerID: [', workerId, ']'

	call readState(workerId, currentPhotonCount, ISEED1, ISEED2, Detector, imageHeight, imageWidth)

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

	currentPhotonCount = currentPhotonCount + 1

	do NP = currentPhotonCount, PhotonCount
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
		if (MAT.EQ.0) then
			call STEP(DSMAX(0), DSEF, NCROSS)
			if (MAT.EQ.0) then
				CYCLE
			endif
			if (IBODY.EQ.1) then
				call checkDetectorHit(Detector)
				CYCLE
			endif
		endif
		call CLEANS
		call START
		do
			call JUMP(DSMAX(MAT), DS)
			call STEP(DS, DSEF, NCROSS)
			! IBODY equals 1 means that the particle hit the detector
			if (IBODY.EQ.1) then
				call checkDetectorHit(Detector)
			else
				! NCROSS greater than 0 means the particle crossed an interface
				if (NCROSS.GT.0) then
					if (MAT.NE.0) then
						call START
						CYCLE
					endif
				else
					call KNOCK(DE, ICOL)
					if (E.GE.EABS(KPAR,MAT)) then
						CYCLE
					endif
				endif
			endif
			! process seconday particles, if any		
			call SECPAR(LEFT)
			if (LEFT.GT.0) then
				call START
				CYCLE
			else
				EXIT
			end if
		end do
		! exit when shutdown signal is received
		if (continueExecution == 0) then
			EXIT
		end if
	end do

	call saveState(workerId, NP, ISEED1, ISEED2, Detector, imageHeight, imageWidth) 
	print *, 'State saved for proy:[', proyNum, ']', ' workerId:[', workerId, '] currentPhotonCount:[', NP, ']'
	MaxCount = maxval(Detector)
	Open(IWR,File='maxvalue' // int2str(workerId),Form='unformatted',Access='DIRECT',RECL=4)
	Write(IWR, rec=1) MaxCount  
	Close(Unit = IWR)
	print *, 'Writing MaxValue [', MaxCount, '] for workerId: [', workerId, ']'

	return       
end

integer function getmax(Nproy)
implicit none
integer                         :: NProy, workerId
integer*4                       :: IR
integer*4                       :: maxValue, value
character*20                    :: int2str

maxValue = 0
IR = 10
do workerId = 1, Nproy
	Open(IR,File='maxvalue' // int2str(workerId),Form='unformatted',Access='DIRECT',RECL=4)
	Read(IR, rec=1) value
	maxValue = max(maxValue, value)
	Close(Unit = IR, STATUS = 'DELETE') 
	print *, 'Read MaxValue [', value, '] for workerId: [', workerId, ']'
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
Print *, 'DETECTOR [', proyNum, '] MAX VAL [', MAXVAL(Detector), ']'

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

