Program MyNBody
  implicit none
  real :: Start, Finish, CPUTIME
  integer Numpart                                       ! Number of particles 
  integer SpaceDim                                      ! Number of Spatial Dimensions 
  double precision tStart, tEnd, tau, totalM, pi, rp    ! Time
  integer Numsteps                                      ! Number of time steps
  integer i,j                                           ! Do loops
  integer nt                                            ! Index for time step

!---------------------------------------------------------------------------------------------!
  !y(:,:) is an array of dimensions (2*Numpart, SpaceDim) speicifying the state of the system.
  !the first NumPart elements are positions and the second NumPart
  ! each element is considerd a vector of SpaceDim dimensions
!---------------------------------------------------------------------------------------------!
  double precision, allocatable :: y(:,:), rt1(:), rt2(:), MassArray(:), f2bt(:), yCM(:,:)
  double precision BoxL               !BoxL is the size of the intial volume of the box

!---------------------------------------------------------------------------------------------!
  call cpu_time(Start)
  pi=dacos(-1.d0)
 
  SpaceDim=3  ! Set the Spatial Dimension
 
  NumPart=3   ! Set the number of particles
  
  tStart = 0.d0
  tEnd = 2*pi

 
  NumSteps=5000  ! Set the number of RK4 steps:

  
  tau= (tEnd-tStart)/dble(NumSteps-1)  ! Set the time step:

  rp = 1.d0                            !(4.d0*pi**2)**(-1.d0/3.d0)
  allocate(y(2*NumPart,SpaceDim),rt1(SpaceDim), rt2(SpaceDim), MassArray(NumPart),f2bt(SpaceDim))
  allocate(yCM(2, SpaceDim))
  !--------------------------------------!
  ! Intial conditions for n-body problem !
  !--------------------------------------!
  MassArray(1)=3.0d-6
  MassArray(2)=3.694d-8
  MassArray(3)=1.d0
 

  y(1,:)=(/ -1.d0, 0.d0, 0.d0/)       ! Earth
  y(2,:)=(/-1.0025699d0, 0.d0, 0.d0/) ! Moon
  y(3,:)=(/0.d0, 0.d0, 0.d0/)         ! Sun

  y(1+NumPart,:)=(/0.d0, -1.d0, 0.d0 /)                    ! Velocity of Earth
  y(2+NumPart,:)=(/0.d0,-(12.0d0*0.0025699d0+1.0d0),0.d0/) ! Velocity of Moon
  y(3+NumPart,:)=(/0.d0,0.d0,0.d0/)                        ! Velocity of Sun
  
  !----------------------------------------------------------------------------------!
  !Calculates the Cm position and velocity, then subtracts them to form the CM frame !
  !----------------------------------------------------------------------------------!

   totalM=0.d0
   yCM(:,:)=0.d0
   do i =1, NumPart
      totalM=totalM+MassArray(i)
      yCM(1,:)=yCM(1,:)+MassArray(i)*y(i,:)          ! position of Center of Mass Transformation, actually totalM*yCM!

      yCM(2,:)=yCM(2,:)+MassArray(i)*y(i+NumPart,:)  ! Velocity of CM transformation, actually totalM*yCM !

   enddo
   yCM(:,:)=yCM(:,:)/totalM
   ! write(6,*) yCM(:,:)
   write(6,*) "rCM="
   write(6,*) yCM(1,:)
   write(6,*) "vCM="
   write(6,*) yCM(2,:)
   do i=1, NumPart
      y(i,:)=y(i,:)-yCM(1,:)                ! subtract off the center of mass position
      y(i+NumPart,:)= y(i+NumPart,:)-yCM(2,:) !subtract off the center of mass velocity
      !This is the transformation part of turning the positons of each particle and converting them into a simplified coordinate system.
   enddo

      totalM=0.d0
   yCM(:,:)=0.d0
   do i =1, NumPart
      totalM=totalM+MassArray(i)
      yCM(1,:)=yCM(1,:)+MassArray(i)*y(i,:)          ! position of Center of Mass Transformation, !actually totalM*yCM!

      yCM(2,:)=yCM(2,:)+MassArray(i)*y(i+NumPart,:)  ! Velocity of CM transformation,! actually totalM*yCM !

   enddo
   yCM(:,:)=yCM(:,:)/totalM
 ! write(6,*) yCM(:,:)
   write(6,*) "rCM="
   write(6,*) yCM(1,:)
   write(6,*) "vCM="
   write(6,*) yCM(2,:)
   
    do i=1,NumPart
      write(101,10) y(i,:)
   enddo
   
   nt=0
   do nt=1, NumSteps
     ! do i=1, NumPart
      !write(100,10) nt*tau,( y(i,:), y(i+NumPart,:), i=1,NumPart)
      write(100,10) nt*tau,( y(i,:), i=1,NumPart)
     ! enddo
      call RK4Step(y,tau,NumPart, SpaceDim, MassArray)
     ! call algo12(
   end do
   
   deallocate (y ,rt1,rt2,MassArray,f2bt,yCM)

10 format(1P,100e16.8)
   call cpu_time(Finish)
   CPUTIME= Finish-Start
   print *, "Run time =",CPUTIME,"seconds"
   print *, "Simulation time for entire Solar System is",CPUTIME*(3.1709792d-8)*(4.6d9)/(2*pi),"years"
 endprogram MyNBody
 
 !-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-!
 ! Subroutine f2b calcualte the force that particle j exerts on i   !
 ! in this case just particle 1 on particel 2 and vice versa        !
 ! The result is returned in a SpaceDim vector f2bi                 !
 !-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-!

 subroutine f2b(t,j,i,y,SpaceDim,NumPart,MassArray,f2bi) ! Force per mass-i that particle j exerts on i
   implicit none
   integer SpaceDim, NumPart, i, j
   double precision y(2*NumPart,SpaceDim),ri(SpaceDim), rj(SpaceDim), t, f2bi(SpaceDim), rvecji(SpaceDim),rji
   double precision MassArray(NumPart)
   double precision eps ! the Softening parameter in the force = -Gm1m2/(r**2+eps^2)set (G=1)
   eps=0.d0

   !rvecji points from j to i
   rvecji=y(i,:)-y(j,:)
   rji = dsqrt(dot_product(rvecji,rvecji))
   
   !Gravitational
   f2bi = -MassArray(j)*rvecji/(rji+eps)**3
   ! write(6,*) j, i, rji, f2bi*MassArray(i)
 end subroutine f2b

subroutine f2bfixedstar(t,j,i,y,SpaceDim,NumPart,MassArray,f2bi) ! Force per mass-i that particle j exerts on i
   implicit none
   integer SpaceDim, NumPart, i, j
   double precision y(2*NumPart,SpaceDim),ri(SpaceDim), rj(SpaceDim), t, f2bi(SpaceDim), rvecji(SpaceDim),rji
   double precision MassArray(NumPart)
   double precision eps ! the Softening parameter in the force = -Gm1m2/(r**2+eps^2)set (G=1)
   eps=0.d0
   !i,j =1,2, or 3
   !rvecji points from j to i
   rvecji=y(i,:)-y(j,:)
   rji = dsqrt(dot_product(rvecji,rvecji))
   if ((i.eq.2).or.(i.eq.3)) then
   f2bi=0.d0
   else if (i.eq.1) then
   f2bi = -MassArray(j)*rvecji/(rji+eps)**3
end if

 end subroutine f2bfixedstar

 
 !-------______________---------------________-----____---_--_-_--_-___--_-__-_--!
 !Subrountine Force_Oni calcuates the total force per mass-i acting on
 !particle i due to all other particles.
 !The result is returned in a SpaceDim Vector ffi
 !For MD(???) simulations involving short range interactions, this subroutine should be modified
 !to restrict the sum to only those particles within a definined radius of influence
 !------___________----_____________-----______----____---__-_-_-__--_--__--_--_-!

 subroutine ForceOn_i(t,i,y,SpaceDim,NumPart,MassArray,ffi)
   implicit none
   double precision y(2*NumPart,SpaceDim), t, ffi(SpaceDim), fji(SpaceDim) ! fji is the force j exerts on i!
   double precision MassArray(NumPart)
   integer NumPart, SpaceDim,i,j

   !________________-------------______________!
   !If particles get too close, then the G field
   !Between them explodes and gets too big.....!
   !----------------_____________--------------!

   ffi= 0.d0
   do j=1,i-1
      call f2b(t,j,i,y,SpaceDim,NumPart,MassArray,fji)
       !call f2bfixedstar(t,j,i,y,SpaceDim,NumPart,MassArray,fji)
      ffi= ffi+fji
   end do
   do j=i+1, NumPart
      call f2b(t,j,i,y,SpaceDim,NumPart,MassArray,fji)
      !call f2bfixedstar(t,j,i,y,SpaceDim,NumPart,MassArray,fji)
      ffi = ffi+fji
   end do
 ! print *, ffi
 end subroutine ForceOn_I
 !---------________________________--------------!
 ! ^^Looping over all the particels that arnt the one youre calculating the total force on
 ! ffi is a vector that is going to store the total force on particle i
 ! f2b then calcualtes the two body on particle i
 !__________----------------------_______________!

 subroutine dydt_f(t,y,SpaceDim,NumPart,MassArray,fresult)
   implicit none
   double precision y(2*NumPart, SpaceDim), t, ffi(SpaceDim), fji(SpaceDim)
   double precision MassArray(NumPart), fresult(2*NumPart,SpaceDim)
   integer NumPart, SpaceDim, i, j

   do i=1, NumPart
      fresult(i,:)=y(i+NumPart,:)
      call ForceOn_i(t,i,y,SpaceDim,NumPart,MassArray,ffi)
      fresult(i+NumPart,:)=ffi !Calculate the total force on particle i
    ! write(6,*) i, ffi
   end do

 end subroutine dydt_f
 !____________--------------_______________________-------------______________!
 subroutine RK4Step(y,tau,NumPart, SpaceDim, MassArray)
   implicit none
   integer n, i ,j, NumPart, SpaceDim
   double precision t
   double precision y(2*NumPart,SpaceDim), ytemp(2*NumPart,SpaceDim), tau, fres(2*NumPart,SpaceDim), MassArray(NumPart)
   double precision k1(2*NumPart, SpaceDim), k2(2*NumPart,SpaceDim),k3(2*NumPart,SpaceDim), k4(2*NumPart,SpaceDim)

   call dydt_f(t,y,SpaceDim, NumPart,MassArray,fres)
   k1 = tau*fres
   call dydt_f(t+0.5d0*tau, y+0.5d0*k1,SpaceDim, NumPart,MassArray,fres)
   k2=tau*fres
   call dydt_f(t+0.5d0*tau, y+0.5d0*k2,SpaceDim, NumPart,MassArray,fres)
   k3= tau*fres
   call dydt_f(t+tau, y+ k3 ,SpaceDim, NumPart,MassArray,fres)
   k4= tau*fres

   ytemp=y +k1/6.0d0 + k2/3.0d0 + k3/3.0d0 + k4/6.0d0
   y=ytemp
 end subroutine RK4Step
 
!  ************** testsymp.for **************************************************
  ! Program Suite testsymp (fixed form Fortran)
  !
  ! Simple driver and test program for symplectic integration routines
  ! of E. Teloy. See paper of A. Seiter and Ch. Schlier
  ! "High-Order Symplectic Integration: An Assessment"
  ! Computer Physics Communications, accepted
  !
  ! Written in standard Fortran 90 (fixed form). To change into free form
  ! module coefficients must be changed as follows: blancs in the numbers
  ! must be deleted, and the continuation codes (&) moved.
  !
  ! Select double or quadruple precision (mk = 8 or 16) in module
  ! coordinates. These kinds may be processor or compiler dependent.
  !
  ! Set dimension of problem (qdim) also in module coordinates.
  !
  ! One out of the documentd five integrators can be selected by
  ! uncommenting them in subroutine algo12.
  !
  ! The main program is placed below, since some compilers want to compile
  ! the modules first!
  !
  ! ---------------------------------------------------------------------
 
 ! xtime == tau
 ! xdt == NumSteps
 
  module coordinates

    !     Defines: mk, qdim, xq, xqd, xtime, xdt used by the other program units

    implicit none

    !     Select double or quadruple precision here:
    !     integer, parameter :: mk = 16
    integer, parameter :: mk = 8

    !     Dimension of coordinate space is defined here:
    integer, parameter :: qdim = 1000

    !     Coordinates (q:1-6, p:7-12) and forces ("right sides")
    real (kind = mk), dimension (1:2*qdim) :: xq, xqd
    ! Time and time step
    real (kind = mk) :: xtime, xdt

    ! Potential and other parameters could be added here

  end module coordinates

  !----------------------------------------------------------------------

  module coefficients

    ! Coefficients for five symplectic integration routines of E. Teloy
    ! See paper of A. Seiter and Ch. Schlier
    ! "High-Order Symplectic Integration: An Assessment"
    ! Computer Physics Communications, to be published
    !
    ! Listed is only the first half of the coefficients, ai(0) to ai(nc),
    ! where 2*nc is the number of substages of the integrator.
    ! The full set is computed on first start in subroutine algo12.
    !
    ! This module is correct only in fixed form Fortran 90!
    ! To change it to free form the blanks in the coefficients must be
    ! deleted, and the "&"s must be moved from before the continuation lines
    ! to the end of the lines which are to be continued.

    use coordinates, only: mk

    implicit none

    integer, parameter :: &
     typ6a = 1, nc6a = 9, &
     typ6b = 2, nc6b = 8, &
     typ8a = 1, nc8a = 17, &
     typ8b = 2, nc8b = 17, &
     typ8c = 2, nc8c = 17

    real (kind = mk), dimension(0:9), parameter :: &
    ai6a = (/ &
      0.09517625454177405267746114335519342_mk, &
      0.66629689399770780134207498907168068_mk, &
     -0.12795028552368677941219191621429411_mk, &
      0.02461890095210508713078430308713062_mk, &
      0.10597295345325113143793587608716998_mk, &
     -0.41072553361795113231992873918199025_mk, &
      0.44822227660082748416851634186561201_mk, &
      0.65772926205091317768935130009339042_mk, & 
     -0.02142119907216588887172144509368130_mk, &
     -0.87583904676554986768456370614042295_mk /)

    real (kind = mk), dimension(0:8), parameter :: &
     ai6b = (/ &
      0.06942944346252987735848865824703402_mk, &
      0.28487837717280084052745346456657828_mk, &
     -0.13315519831598209409961309951373512_mk, &
      0.32783975759612945412054678367325547_mk, &
      0.00129038917981078974230481746443284_mk, &
     -0.38122104271932629475622784374211274_mk, &
      0.42243536567364142699881962380226825_mk, &
      0.26850290795039600010822759550227899_mk, &
      0.28000000000000000000000000000000000_mk /)

    real (kind = mk), parameter :: &
     ai8a(0:17) = (/ &
      0.04020757626295627296653921454892367_mk, &
      0.10968252140081995880852111452131455_mk, &
      0.17023759564885894706257453906663563_mk, &
      0.36756158806337006433149757369026277_mk, &
      0.24370233998503432353195633486895307_mk, &
     -0.04544131419758065661437375963088864_mk, & 
      0.56601963795366046019899599701939548_mk, &
      0.00022167162169864039643822185570309_mk, &
     -0.58169695762497039518529999797620005_mk, &
      0.05519927098092328759679762829526377_mk, &
     -0.24138639830477987453171482029238617_mk, &
     -0.12513929981618023524050370745321727_mk, &
      0.36115097569793127373014000321599616_mk, &
     -0.04284389352937610255914308734324331_mk, &
     -0.53225450460377165284025446933453953_mk, &
     -0.00393367299329157410510456094858013_mk, &
      0.47401973498508064506706319888322175_mk, &
      0.36938625693923323477174115402677030_mk /)

    real (kind = mk), parameter :: &
     ai8b(0:17) = (/ &
      0.03676680389912337302666154929429291_mk, &
      0.11072655003739784175754797312279745_mk, &
      0.16040429374255560219395381214509780_mk, &
      0.61101267825171523627962718607785428_mk, & 
     -0.00472877643941287918639412436088645_mk, &
     -0.19202809069032535396838334049379558_mk, &
      0.02983098489335056954884440558763334_mk, &
     -0.25979073929811660257162833544861286_mk, &
      0.19135844311091097984885756175207225_mk, &
      0.38384564066882093754274499421236298_mk, &
     -0.03781968145745128677723635761417376_mk, &
      0.32661664886778120135972921761872954_mk, &
      0.00351845996378093605518443870229385_mk, & 
     -0.53463443374897025678663398242742174_mk, &
      0.13067013867271618676514580608303276_mk, &
     -0.39935632081078281354806842349635698_mk, &
     -0.01000066638557348147501709158936269_mk, &
      0.90721613344495961987012942166888585_mk /)

    real (kind = mk), parameter :: &
    ai8c(0:17) = (/  &
      0.04463795052359022755913999625733590_mk, &
      0.13593258071690959145543264213495574_mk, &
      0.21988440427147072254445535069606167_mk, &
      0.13024946780523828601621193778196846_mk, &
      0.10250365693975069608261241007779814_mk, &
      0.43234521869358547487983257884877035_mk, &
     -0.00477482916916881658022489063962934_mk, & 
     -0.58253476904040845493112837930861212_mk, &
     -0.03886264282111817697737420875189743_mk, &
      0.31548728537940479698273603797274199_mk, & 
      0.18681583743297155471526153503972746_mk, &
      0.26500275499062083398346002963079872_mk, &
     -0.02405084735747361993573587982407554_mk, & 
     -0.45040492499772251180922896712151891_mk, &
     -0.05897433015592386914575323926766330_mk, &
     -0.02168476171861335324934388684707580_mk, &
      0.07282080033590128173761892641234244_mk, &
      0.55121429634197067334405601381594315_mk /)

  end module coefficients
  !
  ! ---------------------------------------------------------------------
  ! qdim == SpaceDim
  ! xq == y(:,:)
  ! xqd == 
  !----------------------------------------------------------------------

  subroutine algo12

    ! This subroutine performs one full integration step.
    ! The increments xqdh are added to xqd only after the final substep to
    ! preserve as much accuracy as possisble.
    ! The use of the FSAL property has not been implemented.

    ! Uncomment one of the five following lines to determine which
    ! integrator you want to use:

    !      use coefficients, only : typ =>typ6a, nc => nc6a, ai => ai6a
    !      use coefficients, only : typ =>typ6b, nc => nc6b, ai => ai6b
    !      use coefficients, only : typ =>typ8a, nc => nc8a, ai => ai8a
    !      use coefficients, only : typ =>typ8b, nc => nc8b, ai => ai8b
    use coefficients, only : typ =>typ8c, nc => nc8c, ai => ai8c

    ! Exchange of parameters and coordinates with other program units:
    use coordinates,  only : mk, qdim, xq, xqd, xtime, xdt

    implicit none

    integer :: nn, ii

    real (kind = mk), dimension (0:2*nc) :: hhc
    real (kind = mk), dimension (1:2*qdim) :: xqh, xqdh

    external pots1, pots2

    entry algini

    write (*,*) "typ ", typ, " nc: ", nc
    write (9,*) "typ ", typ, " nc: ", nc
    ! Compute full set of coefficients:
    hhc(0:nc) = ai(0:nc)
    hhc(nc+1:2*nc) = ai(nc-1:0:-1)
    ! Check integrity of coefficients (sum must be 2.0):
    !      write (*,*)"sum of hh ", sum(hhc)
    !      write (*,*)

    hhc(0:2*nc) = xdt*hhc(0:2*nc)

    entry algrun

    xqh(1:2*qdim) = xq(1:2*qdim)
    xqdh(1:2*qdim) = 0.0e0_mk

    choice: if (typ == 1) then      ! Choose algorithm: typ = 1
       do ii = 0, 2*nc-2, 2
          call pots1      ! dT/dp    !replace with intial velocity vectors ?
          xqdh(1:qdim) = xqdh(1:qdim) + xqd(qdim+1:2*qdim)*hhc(ii)
          xq(1:qdim) = xqh(1:qdim) + xqdh(1:qdim)
          call pots2     ! -dV/dq    !replace with dy/dt_f
          xqdh((qdim+1):2*qdim)=xqdh((qdim+1):2*qdim)+xqd(1:qdim)*hhc(ii+1)
          xq( qdim+1:2*qdim) = xqh( qdim+1:2*qdim) + xqdh( qdim+1:2*qdim)
       enddo
       call pots1      ! dT/dp
       xqdh(1:qdim) = xqdh(1:qdim) + xqd(qdim+1:2*qdim)*hhc(2*nc)

    else  choice                    ! Choose algorithm: typ = 2

       do ii = 0, 2*nc-2, 2
          call pots2  ! -dV/dq     !replace with dy/dt_f
          xqdh( qdim+1:2*qdim) = xqdh( qdim+1:2*qdim) + xqd(1:qdim)*hhc(ii)
          xq( qdim+1:2*qdim) = xqh( qdim+1:2*qdim) + xqdh( qdim+1:2*qdim)
          call pots1   ! dT/dp     !replace with intial velocity vectors?
          xqdh(1:qdim) = xqdh(1:qdim) + xqd(qdim+1:2*qdim)*hhc(ii+1)
          xq(1:qdim) = xqh(1:qdim) + xqdh(1:qdim)
       enddo
       call pots2  ! -dV/dq replace with Dy/dt_f
       xqdh(qdim+1:2*qdim) = xqdh(qdim+1:2*qdim) + xqd(1:qdim)*hhc(2*nc)

    endif choice

    xq(1:2*qdim) = xqh(1:2*qdim) + xqdh(1:2*qdim)
    xtime = xtime + xdt

    return
  end subroutine algo12

  !----------------------------------------------------------------------

  ! Subroutines for the "right sides" of the harmonic oscillator with
  ! unit frequency:
  ! pots1 calculates the beginining of the calculation of the change in Kinetic Energy
  subroutine pots1        ! dT/dp
    use coordinates, only : mk, qdim, xq, xqd
    xqd(qdim+1:2*qdim) = xq(qdim+1:2*qdim)
    return
  end subroutine pots1
  ! pots2 is the potential
  subroutine pots2       ! -dV/dq
    use coordinates,  only : mk, qdim, xq, xqd
    xqd(1:qdim) = -xq(1:qdim)
    return
  end subroutine pots2
!  ****
! Want to replace pots1/ pots2 with our subroutines.
