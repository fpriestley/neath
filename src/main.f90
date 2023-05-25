! 2015 UCLCHEM by Serena Viti update by Jon Holdship
! Rewritten for FORTRAN 95 and modulised
PROGRAM uclchem
!Everything to do with physics should be stored in a physics module based on physics-template.f90
!UCLCHEM uses density and temperature in chemical calculations so these MUST be provided, everything else is
!user dependent
USE physics
USE chemistry
use OMP_LIB
IMPLICIT NONE
    INTEGER :: ios=0,line=0,fileNum=12,pos,readIndx
    CHARACTER (LEN=100):: paramFile, buffer,label
    CHARACTER (LEN=100):: abundFile,outputFile,columnFile
    LOGICAL :: columnFileRead=.False.
    integer :: maxthreads,threadno,iomp

    !All model parameters are given a default value in paramters.f90
    !Full explanations of those parameters in the comments of that file
    INCLUDE 'defaultparameters.f90'
    !Any subset of the parameters can be passed in a file on program start
    !see example.inp
    
    !read(*,*) trajecfile,points,ntime
    trajecfile = 'neath_example_data.out'
    points = 250
    ntime = 121
    outputFile='output/'//TRIM(trajecfile)
    
    INCLUDE 'readparameters.f90'

    dstep=1
    tstep = 1
    currentTime=0.0
    timeInYears=0.0

    write(*,*) 'Reading initial conditions'
    write(*,"(I5,' particles, ',I5,' times')") points,ntime

    open(unit=66,file='data/'//trajecfile,status='old')

    allocate(timegrid(ntime),densgrid(points,ntime),tempgrid(points,ntime),avgrid(points,ntime),nh2grid(points,ntime),ncogrid(points,ntime))

    !Set up with initial values. For chemistry this is setting initial abundances and assigning memory for ODE solver
    CALL initializePhysics
    CALL initializeChemistry
    
    !update physics
    DO dstep=1,points
        CALL updatePhysics
    END DO

    !loop until the end condition of the model is reached 
    DO WHILE (tstep .le. ntime)

        !store current time as starting point for each depth step
        IF (points .gt. 1) THEN
!            currentTimeold=targetTime
!            currentTime=currentTimeold
           currentTimeold = currentTime
        END IF
        !Each physics module has a subroutine to set the target time from the current time
        CALL updateTargetTime
 
        write(*,"('Timestep ',I4)") tstep
        !loop over parcels, counting from centre out to edge of cloud
        DO dstep=1,points

           !reset target for next depth point
           if (points .gt. 1)currentTime=currentTimeold
           
           !update chemistry
           CALL updateChemistry

            !set time to the final time of integrator rather than target     
           !targetTime=currentTime

           CALL updatePhysics
            
            !get time in years for output
            timeInYears= currentTime*year
            !write this depth step
            CALL output
         END DO

         tstep = tstep + 1
         
      END DO

      deallocate(timegrid,densgrid,tempgrid,avgrid,nh2grid,ncogrid)

      close(unit=66)
      
END PROGRAM uclchem
