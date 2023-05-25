!Simple physics module. Models points along a 1d line from the centre to edge of a cloud. Assuming the cloud is spherical you can average
!over the points to get a 1d average and then assume the rest of sphere is the same.

MODULE physics
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    INTEGER :: dstep,points,ntime
    !Switches for processes are also here, 1 is on/0 is off.
    INTEGER :: collapse,switch,phase
    INTEGER :: h2desorb,crdesorb,crdesorb2,uvcr,desorb

    !evap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see initializeChemistry)
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    INTEGER :: evap,ion,solidflag,volcflag,coflag,tempindx

    double precision :: filetime,filedens,junk
    integer :: tstep,it,endcheck
    character(len=100) :: trajecfile
    double precision,allocatable :: timegrid(:),densgrid(:,:),tempgrid(:,:),avgrid(:,:),nh2grid(:,:),ncogrid(:,:)

    !variables either controlled by physics or that user may wish to change
    DOUBLE PRECISION :: initialDens,timeInYears,targetTime,currentTime,currentTimeold,finalDens,finalTime,grainRadius,initialTemp
    DOUBLE PRECISION :: cloudSize,rout,rin,baseAv,bc,olddens,maxTemp,vs
    
    !Arrays for phase 2 temp profiles. parameters for equation chosen by index
    !arrays go [1Msun,5, 10, 15, 25,60]
   
    DOUBLE PRECISION,PARAMETER :: tempa(6)=(/1.927d-1,4.8560d-2,7.8470d-3,9.6966d-4,1.706d-4,4.74d-7/)
    DOUBLE PRECISION,PARAMETER :: tempb(6)=(/0.5339,0.6255,0.8395,1.085,1.289,1.98/)
    DOUBLE PRECISION,PARAMETER :: solidtemp(6)=(/20.0,19.6,19.45,19.3,19.5,20.35/)
    DOUBLE PRECISION,PARAMETER :: volctemp(6)=(/84.0,86.3,88.2,89.5,90.4,92.2/)
    DOUBLE PRECISION,PARAMETER :: codestemp(6)=(/95.0,97.5,99.4,100.8,101.6,103.4/)

    DOUBLE PRECISION, allocatable :: av(:),coldens(:),temp(:),dens(:)
    !Everything should be in cgs units. Helpful constants and conversions below
    DOUBLE PRECISION,PARAMETER ::PI=3.141592654,mh=1.67e-24,kbolt=1.38d-16
    DOUBLE PRECISION, PARAMETER :: year=3.16455d-08,pc=3.086d18

    !variables for collapse modes
    DOUBLE PRECISION :: unitrho,unitr,unitt,c_s
    DOUBLE PRECISION :: dimrho,dimr,dimt,maxdimt
    DOUBLE PRECISION :: rho0,r0
    DOUBLE PRECISION,PARAMETER :: G_N = 6.674d-8

    integer :: ii,jj,kk

CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN.
!YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

    SUBROUTINE initializePhysics
        IF (ALLOCATED(av)) DEALLOCATE(av,coldens,temp,dens)
        ALLOCATE(av(points),coldens(points),temp(points),dens(points))
        cloudSize=(rout-rin)*pc
        dens=initialDens
        temp=initialTemp

      do ii=1,points
         read(66,*) timegrid(1),(junk,kk=1,3),densgrid(ii,1),(junk,kk=1,3),tempgrid(ii,1),(junk,kk=1,2),avgrid(ii,1),nh2grid(ii,1),ncogrid(ii,1)
!         dens(ii) = densgrid(ii,1)/(1.4*mh)
         do jj=2,ntime
            read(66,*) timegrid(jj),(junk,kk=1,3),densgrid(ii,jj),(junk,kk=1,3),tempgrid(ii,jj),(junk,kk=1,2),avgrid(ii,jj),nh2grid(ii,jj),ncogrid(ii,jj)
         end do
      end do

      coldens = 2e4*sqrt(dens/(1.4*mh*G_N))
      av = baseAv + coldens*6e-22
      
    END SUBROUTINE

    SUBROUTINE updateTargetTime
    !This is the time step for outputs from UCL_CHEM NOT the timestep for the integrator.
    !targetTime in seconds for DLSODE, timeInYears in years for output.

      targettime = timegrid(tstep) + 1./year
      
    END SUBROUTINE updateTargetTime

    SUBROUTINE updatePhysics

      dens(dstep) = densgrid(dstep,tstep)
      temp(dstep) = tempgrid(dstep,tstep)
      coldens(dstep) = avgrid(dstep,tstep)
      av(dstep) = 5.348e-22 * coldens(dstep)

    END SUBROUTINE updatePhysics

!This function works out the time derivative of the density, allowing DVODE to update density with the rest of our ODEs
!It get's called by F, the SUBROUTINE in chem.f90 that sets up the ODEs for DVODE
!Currently set to Rawlings 1992 freefall.
    pure FUNCTION densdot(density)
        DOUBLE PRECISION, INTENT(IN) :: density
        DOUBLE PRECISION :: densdot
        !Rawlings et al. 1992 freefall collapse. With factor bc for B-field etc
        IF (density .lt. finalDens) THEN
             densdot=bc*(density**4./initialDens)**0.33*&
             &(8.4d-30*initialDens*((density/initialDens)**0.33-1.))**0.5
        ELSE
            densdot=0.0
        ENDIF
    END FUNCTION densdot
END MODULE physics

!REQUIRED PHYSICS ENDS HERE, ANY ADDITIONAL PHYSICS CAN BE ADDED BELOW.
