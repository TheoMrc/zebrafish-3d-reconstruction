module TimeScheme
  use AdvectionProblem
  
  implicit none
  
contains  
  
!! Func to initialize the time scheme
subroutine SetInitialCondition(rho0, u, v, rho, rhou, rhov)
  real(pr), dimension(:,:), intent(in) :: rho0, u, v
  real(pr), dimension(:,:), intent(inout) :: rho, rhou, rhov
  alp(1) = 0._pr
  alp(2) = dt*0.5_pr
  alp(3) = dt*0.5_pr
  alp(4) = dt

  beta(1) = dt/6._pr
  beta(2) = dt/3._pr
  beta(3) = dt/3._pr
  beta(4) = dt/6._pr

  rho = rho0
  rhou = rho*u
  rhov = rho*v

end

!! Calculate the approximate solution at time t
subroutine AdvanceTime(nt, dt_, timeSchemeChoice, advectionSchemeChoice, interpChoice, BC, rho, u, v, rhou, rhov)
  integer, intent(in) :: nt, timeSchemeChoice, advectionSchemeChoice, interpChoice, BC
  real(pr), intent(in) :: dt_
  real(pr), dimension(:,:), intent(inout) :: rho, u, v, rhou, rhov
  real(pr),dimension(size(rho,1),size(rho,2)) :: u0, v0, rho00, rhou0, rhov0, rhoNext1, rhoNext2, un, zero
  integer :: i, j
  real(pr) :: maxLS, minLS, seuilLS
  u0 = u
  v0 = v
  un = 1._pr
  zero = 0._pr
  rho00 = rho
  rhou0 = rhou
  rhov0 = rhov

  if (timeSchemeChoice == 1) then
          call AdvanceEulerExp(nt, dt_, advectionSchemeChoice, BC, rho, u, v, rhou, rhov)
  else if (timeSchemeChoice == 3) then
          call AdvanceRK3(nt, dt_, advectionSchemeChoice, BC, rho, u, v, rhou, rhov)
  else if (timeSchemeChoice == 4) then
          call AdvanceRK4(nt, dt_, advectionSchemeChoice, BC, rho, u, v, rhou, rhov)
  else if (advectionSchemeChoice == 2) then

          !/!\ DIVERGENCE NULL /!\!
          call AdvanceLWperiodic(nt, dt_, rho, u, v)
          call AdvanceLWperiodic(nt, dt_, rhou, u, v)
          call AdvanceLWperiodic(nt, dt_, rhov, u, v)
  else if (advectionSchemeChoice == 3) then
          write(*,*) "no particle"
  else if (advectionSchemeChoice == 5) then
          write(*,*)"Lax-Friedrichs Flux Splitting NOT Implemented"
  endif

end

!/*************************
! * ExplicitEulerIterator *
! *************************/

!! Main function that advances the scheme in time
subroutine AdvanceEulerExp(n, tn, scheme, BC, rho, u, v, rhou, rhov)
  integer, intent(in) :: n
  real(pr), intent(in) :: tn
  integer, intent(in) :: scheme, BC
  real(pr), dimension(:,:), intent(inout) :: rho, u, v, rhou, rhov
  real(pr), dimension(size(rho,1),size(rho,2)) :: rhoNext, rhouNext, rhovNext, zero, un
  real(pr), dimension(size(rho,1),size(rho,2)) :: rhoNext1, rhoNext2
  real(pr), dimension(size(rho,1),size(rho,2)) :: rhouNext1, rhovNext1
  real(pr), dimension(size(rho,1),size(rho,2)) :: rhouNext2, rhovNext2
  zero = 0._pr
  un = 1._pr

  !/* DIVERGENCE NULL /*!
  call computeFunction2D(dt, rho, tn, rhoNext, scheme, BC, u, v)
  rho = rho + dt*rhoNext

end


!/*************************
! * RK4Iterator *
! *************************/

!! Main function that advances the scheme in time
subroutine AdvanceRK4(n, tn, scheme, BC, rho, u, v, rhou, rhov)
  integer, intent(in) :: n
  real(pr), intent(in) :: tn
  integer, intent(in) :: scheme, BC
  real(pr), dimension(:,:), intent(inout) :: rho, u, v, rhou, rhov
  integer :: i, j, k
  real(pr), dimension(size(rho,1),size(rho,2)) :: kRK, kRKu, kRKv, addedFunctionU,&
  addedFunctionV, addedFunction, rhoNext, rhouNext, rhovNext
  real(pr), dimension(size(rho,1),size(rho,2)) :: kRK1, kRK2, kRKu1, kRKu2, kRKv1,&
  kRKv2, rhoNext1, rhoNext2, rhouNext1, rhouNext2, rhovNext1, rhovNext2, un, zero
  
  un = 1._pr
  zero = 0._pr

  !/* NON NULL DIVERGENCE (WHATEVER) /*!
  rhoNext = rho
  rhouNext = rhou
  rhovNext = rhov
  do k=1,4
    addedFunction = rho + alp(k)*(kRK1 + kRK2)
    addedFunctionU = rhou + alp(k)*(kRKu1 + kRKu2)
    addedFunctionV = rhov + alp(k)*(kRKv1 + kRKv2)
    rhoNext = addedFunctionU
    call computeFunction2D(1._pr, rhoNext, tn, kRK1, scheme, BC, un, zero)
    rhoNext = addedFunctionV
    call computeFunction2D(1._pr, rhoNext, tn, kRK2, scheme, BC, zero, un)
    rhoNext = addedFunctionU*u
    call computeFunction2D(1._pr, rhoNext, tn, kRKu1, scheme, BC, un, zero)
    rhoNext = addedFunctionU*v
    call computeFunction2D(1._pr, rhoNext, tn, kRKu2, scheme, BC, zero, un)
    rhoNext = addedFunctionV*u
    call computeFunction2D(1._pr, rhoNext, tn, kRKv1, scheme, BC, un, zero)
    rhoNext = addedFunctionV*v
    call computeFunction2D(1._pr, rhoNext, tn, kRKv2, scheme, BC, zero, un)
    rhoNext = rhoNext + beta(k)*(kRK1 + kRK2)
    rhouNext = rhouNext + beta(k)*(kRKu1 + kRKu2)
    rhovNext = rhovNext + beta(k)*(kRKv1 + kRKv2)
    do i=1,nx
       do j=1,ny 
          if (rhoNext(i,j).gt.1e-6) then
             u(i,j) = rhouNext(i,j)/rhoNext(i,j) 
             v(i,j) = rhovNext(i,j)/rhoNext(i,j) 
          else
             rhoNext(i,j) = 0._pr
             rhouNext(i,j) = 0._pr
             rhovNext(i,j) = 0._pr
             u(i,j) = 0._pr
             v(i,j) = 0._pr
          endif
       enddo
    enddo
  enddo
  rho = rhoNext
  rhou = rhouNext
  rhov = rhovNext
end

!/*************************
! * RK3Iterator *
! *************************/

!! Main function that advances the scheme in time
subroutine AdvanceRK3(n, tn, scheme, BC, rho, u, v, rhou, rhov)
  integer, intent(in) :: n
  real(pr), intent(in) :: tn
  integer, intent(in) :: scheme, BC
  real(pr), dimension(:,:), intent(inout) :: rho, u, v, rhou, rhov
  integer :: i, j
  real(pr), dimension(size(rho,1),size(rho,2)) :: addedFunction, rhostar, rhostarstar, rhoNext
  real(pr), dimension(size(rho,1),size(rho,2)) :: addedFunctionU, rhoustar, rhoustarstar, rhouNext
  real(pr), dimension(size(rho,1),size(rho,2)) :: addedFunctionV, rhovstar, rhovstarstar, rhovNext
  real(pr), dimension(size(rho,1),size(rho,2)) :: ustar, vstar, ustarstar, vstarstar, un, zero
  real(pr), dimension(size(rho,1),size(rho,2)) :: rhoNext1, rhoNext2, rhouNext1, rhouNext2, rhovNext1, rhovNext2

  un = 1._pr
  zero = 0._pr

  call computeFunction2D(1._pr, rho, tn, addedFunction, scheme, BC, u, v)
  rhostar = rho + dt*addedFunction
  call computeFunction2D(1._pr, rhostar, tn+dt, addedFunction, scheme, BC, u, v) 
  rhostarstar = 0.75_pr*rho + 0.25_pr*rhostar + 0.25_pr*dt*addedFunction
  call computeFunction2D(1._pr, rhostarstar, tn+0.5_pr*dt, addedFunction, scheme, BC, u, v)
  rhoNext = rho/3._pr + 2._pr/3._pr*rhostarstar + 2._pr/3._pr*dt*addedFunction

  rho = rhoNext

end
end

