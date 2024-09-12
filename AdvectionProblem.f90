module AdvectionProblem
  use interpolation
  use variables
  use doubler

  implicit none

contains
real(pr) function interpLS(LS1,x1,y1,LS2,x2,y2,LS3,x3,y3,LS4,x4,y4,x,y) !bilinear interpolation of LS between X1,X2,X3,X4 (cartesian mesh)
    implicit none
    real(pr),intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4,LS1,LS2,lS3,LS4,x,y
    real(pr) :: deltaX, deltaY, mdX1, mdY1, mdX2, mdY2

    deltaX = x4-x1
    deltaY = y4-y1
    mdX1 = x-x1
    mdY1 = y-y1
    mdX2 = x4 - x
    mdY2 = y4 - y

    interpLS = mdX2*mdY2*LS1 + mdX1*mdY2*LS2 + mdX2*mdY1*LS3 + mdX1*mdY1*LS4
    interpLS = interpLS/(deltaX*deltaY)
end function
real(pr) function interpLS2D(LS1,x1,LS2,x2,x) !linear interpolation of LS between X1,X2 (cartesian mesh)
    implicit none
    real(pr),intent(in) :: x1,x2,LS1,LS2,x
    real(pr) :: deltaX, mdX1, mdX2

    deltaX = x2-x1
    mdX1 = x-x1
    mdX2 = x2 - x

    interpLS2D = mdX2*LS1 + mdX1*LS2
    interpLS2D = interpLS2D/deltaX
end function

real(pr) function dist(x1, y1, x2, y2)
    implicit none
    real(pr) :: x1,y1,x2,y2
    dist = sqrt((x2-x1)**2 + (y2-y1)**2)
end function dist

real(pr) function dist3D(x1, y1, z1, x2, y2, z2)
    implicit none
    real(pr) :: x1,y1,z1,x2,y2,z2
    dist3D = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
end function dist3D

real(pr) function det2(x1, y1, x2, y2)
    implicit none
    real(pr) :: x1,y1,x2,y2

    det2 = x1*y2 - y1*x2
end function det2

real(pr) function det(x1, y1, x2, y2, x3, y3) !det(A1A2,A1A3)
    implicit none
    real(pr) :: x1,y1,x2,y2,x3,y3

    det = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
end function det
real(pr) function dotProd(x1, y1, x2, y2, x3, y3) !A1A2.A1A3
    implicit none
    real(pr) :: x1,y1,x2,y2,x3,y3

    dotProd = (x2-x1)*(x3-x1) + (y2-y1)*(y3-y1)
end function dotProd
real(pr) function dotProd3D(x1, y1, z1, x2, y2, z2, x3, y3, z3) !A1A2.A1A3
    implicit none
    real(pr) :: x1,y1,z1,x2,y2,z2,x3,y3,z3

    dotProd3D = (x2-x1)*(x3-x1) + (y2-y1)*(y3-y1) + (z2-z1)*(z3-z1)
end function dotProd3D
real(pr) function norme(x1, y1, x2, y2) !||A1A2||
    implicit none
    real(pr) :: x1,y1,x2,y2

    norme = sqrt((x2-x1)**2 + (y2-y1)**2)
end function norme

character(len=3) function str(k)
    implicit none
    integer, intent(in) :: k
    write (str, '(I3.3)') k
    str = adjustl(str)
end function str
character(len=1) function stri(k)
    integer, intent(in) :: k
    write (stri, '(I1.1)') k
    stri = adjustl(stri)
end function stri


subroutine charge(ntot,i0,np,me,i1,iN)  !!np<=ntot
  implicit none
  integer :: ntot,i0,np,me,i1,iN,q,r

  r = mod(ntot,np)
  q = ntot/np

  if (me < r) then
    i1 = me*q + me + i0
    iN = i1 +q
  else
    i1 = me*q + r + i0
    iN = i1 + q -1
  endif
end subroutine charge

subroutine charge_step(ntot,i0,np,me,i1,iN,step)
  implicit none
  integer :: ntot,i0,np,me,i1,iN,q,r,step

  r = mod(ntot,np)
  q = ntot/np

  if (me < r) then
    i1 = me*(q+1)*step + i0
    iN = i1 + q*step
  else
    i1 = (me*q+r)*step + i0
    iN = i1 + (q-1)*step
  endif
end subroutine charge_step


subroutine gaussian(x, y, z, BC)
  real(pr), intent(in) :: x, y
  real(pr), intent(out) :: z
  integer, intent(in) :: BC
  real(pr) :: xt, yt

  xt = x
  yt = y

  if (BC == 1) then
          call TranslateX(x,xt)
          call TranslateX(y,yt)
  endif

  z = exp(-((xt+(xx(N)-xx(1))/4)*(xt+(xx(N)-xx(1))/4) + (yt+(xx(N)-xx(1))/4)*(yt+(xx(N)-xx(1))/4))/(xx(N)-xx(1)/2))

  if (BC == 2) then
          if ((xt<xx(1)).or.(xt>xx(N))) z = 0
          if ((yt<xx(1)).or.(yt>xx(N))) z = 0
  endif

end

subroutine creneau(x, y, z, BC)
  real(pr), intent(in) :: x, y
  real(pr), intent(out) :: z
  integer, intent(in) :: BC
  real(pr) :: xt, yt

  xt = x
  yt = y

  if (BC == 1) then
          call TranslateX(x,xt)
          call TranslateX(y,yt)
  endif

  if ((xt+(xx(N)-xx(1))/4)*(xt+(xx(N)-xx(1))/4)+(yt+(xx(N)-xx(1))/4)*(yt+(xx(N)-xx(1))/4)<=0.4*(xx(N)-xx(1))) then
    z = 1
  else
    z = 0
  endif

  if (BC == 2) then
          if ((xt<xx(1)).or.(xt>xx(N))) xt = 0
          if ((yt<xx(1)).or.(yt>xx(N))) yt = 0
  endif

end

subroutine initialFunction(x, y, z, scheme, BC)
  real(pr), intent(in) :: x, y
  real(pr), intent(out) :: z
  integer, intent(in) :: scheme, BC

  if (scheme == 1) then
          call creneau(x,y,z,BC)
  else
          call gaussian(x,y,z,BC)
  endif
end


!! Calculates the euclidean norm of a vector
subroutine Norm21D(vv,nn)
  real(pr), dimension(:), intent(in) :: vv
  real(pr), intent(out) :: nn
  integer :: i

  nn = 0
  do i=1,size(vv)
    nn = nn + vv(i)**2
  enddo
  nn = sqrt(nn)

end

!! Calculates the euclidean norm of a vector
subroutine Norm22D(vv,nn)
  real(pr), dimension(:,:), intent(in) :: vv
  real(pr), intent(out) :: nn
  integer :: i, j

  nn = 0
  do i=1,size(vv,1)
    do j=1,size(vv,2)
      nn = nn + vv(i,j)**2
    enddo
  enddo
  nn = sqrt(nn)

end

!! Calculates the euclidean norm of a vector
subroutine Norm23D(vv,nn)
  real(pr), dimension(:,:,:), intent(in) :: vv
  real(pr), intent(out) :: nn
  integer :: i, j, k

  nn = 0
  do i=1,size(vv,1)
    do j=1,size(vv,2)
      do k=1,size(vv,3)
        nn = nn + vv(i,j,k)**2
      enddo
    enddo
  enddo
  nn = sqrt(nn)

end

!! Returns a point in the interval [-L,L] from any point
subroutine TranslateX(x, y)
  real(pr), intent(in) :: x
  real(pr), intent(out) :: y
  real(pr) :: xMod, xMin, xMax
  integer :: Nmod
  xMin = xx(1)
  xMax = xx(N)
  xMod = x
  if (xMod < xMin) then
    Nmod = nint(-xMod/(xMax-xMin))
    xMod = xMod + Nmod*(xMax-xMin)
  endif
  if (xMod > xMax) then
    Nmod = nint(xMod/(xMax-Xmin))
    xMod = xMod - Nmod*(xMax-xMin)
  endif

  y = xMod

end

subroutine computeFunction(alpha, rho_, t, rhonext_, scheme)
  real(pr), intent(in) :: alpha, t
  integer, intent(in) :: scheme
  real(pr), dimension(:), intent(inout) :: rho_
  real(pr), dimension(:), intent(inout) :: rhonext_

  if (scheme == 1) then !Upwind
          call computeFunctionUpwind(alpha, rho_, t, rhonext_)
  else if (scheme == 5) then !WENO5
          call computeFunctionWENO5(alpha, rho_, t, rhonext_)
  endif
end

subroutine computeFunction2D(alpha, rho_, t, rhonext_, scheme, BC, u, v)
  real(pr), intent(in) :: alpha, t
  integer, intent(in) :: scheme, BC
  real(pr), dimension(:,:), intent(inout) :: rho_
  real(pr), dimension(:,:), intent(inout) :: u
  real(pr), dimension(:,:), intent(inout) :: v
  real(pr), dimension(:,:), intent(inout) :: rhonext_

  if ((scheme == 1).and.(BC == 1)) then !Upwind
          call computeFunctionUpwind2Dperiodic(alpha, rho_, t, rhonext_, u, v)
  else if ((scheme == 1).and.(BC == 2)) then !Upwind
          call computeFunctionUpwind2DdirichletH(alpha, rho_, t, rhonext_, u, v)
  else if ((scheme == 5).and.(BC == 1)) then !WENO5
          call computeFunctionWENO52Dperiodic(alpha, rho_, t, rhonext_, u, v)
  else if ((scheme == 5).and.(BC == 2)) then !WENO5
          call computeFunctionWENO52DdirichletH(alpha, rho_, t, rhonext_, u, v)
  endif
end

subroutine neumannBC2D(Vec, i, j, valueBC)
  real(pr), dimension(:,:), intent(in) :: Vec
  integer, intent(in) :: i, j
  integer :: tmpi, tmpj
  real(pr), intent(out) :: valueBC

  tmpi = i
  tmpj = j

  if (i < 1) tmpi = 1
  if (i > size(Vec,1)) tmpi = size(Vec,1)
  if (j < 1) tmpj = 1
  if (j > size(Vec,2)) tmpj = size(Vec,2)

  valueBC = Vec(tmpi, tmpj)
end
subroutine neumannBC(Vec, i, j, k, valueBC)
  real(pr), dimension(:,:,:), intent(in) :: Vec
  integer, intent(in) :: i, j, k
  integer :: tmpi, tmpj, tmpk
  real(pr), intent(out) :: valueBC

  tmpi = i
  tmpj = j
  tmpk = k

  if (i < 1) tmpi = 1
  !if (i > nx) tmpi = nx
  if (i > size(Vec,1)) tmpi = size(Vec,1)
  if (j < 1) tmpj = 1
  !if (j > ny) tmpj = ny
  if (j > size(Vec,2)) tmpj = size(Vec,2)
  if (k < 1) tmpk = 1
  !if (k > nz) tmpk = nz
  if (k > size(Vec,3)) tmpk = size(Vec,3)

  valueBC = Vec(tmpi, tmpj, tmpk)
end


!/*************************
! * Decentre Upwind Advection Scheme *
! *************************/

!! 1D
!! Implementation of the uncentered Upwind scheme
subroutine computeFunctionUpwind(alpha, rho_, t, u_)
  real(pr), intent(in) :: alpha, t
  real(pr), dimension(:), intent(inout) :: rho_
  real(pr), dimension(:), intent(inout) :: u_
  integer :: i

  if (a>0) then
      u_(1) = - a*(rho_(1)-rho_(N))/dx
      do i=2,N
        u_(i) = - a*(rho_(i)-rho_(i-1))/dx
      enddo
  else
      do i=1,N-1
        u_(i) = - a*(rho_(i+1)-rho_(i))/dx
      enddo
      u_(N) = - a*(rho_(1)-rho_(N))/dx
  endif
end

!! 2D
!! Implementation of the uncentered Upwind scheme
subroutine computeFunctionUpwind2Dperiodicold(alpha, rho_, t, u_, u, v)
  real(pr), intent(in) :: alpha, t
  real(pr), dimension(:,:), intent(inout) :: rho_
  real(pr), dimension(:,:), intent(inout) :: u_, u, v
  integer :: i, j
  if ((a>=0).and.(bbb>=0)) then
      u_(1,1) = - a*(rho_(1,1)-rho_(N,1))/dx - bbb*(rho_(1,1)-rho_(1,N))/dy
      do j=2,N
        u_(1,j) = - a*(rho_(1,j)-rho_(N,j))/dx - bbb*(rho_(1,j)-rho_(1,j-1))/dy
      enddo
      do i=2,N
        u_(i,1) = - a*(rho_(i,1)-rho_(i-1,1))/dx - bbb*(rho_(i,1)-rho_(i,N))/dy
        do j=2,N
          u_(i,j) = - a*(rho_(i,j)-rho_(i-1,j))/dx - bbb*(rho_(i,j)-rho_(i,j-1))/dy
        enddo
      enddo
  else if ((a<0).and.(bbb>=0)) then
      do i=1,N-1
        u_(i,1) = - a*(rho_(i+1,1)-rho_(i,1))/dx - bbb*(rho_(i,1)-rho_(i,N))/dy
        do j=2,N
          u_(i,j) = - a*(rho_(i+1,j)-rho_(i,j))/dx - bbb*(rho_(i,j)-rho_(i,j-1))/dy
        enddo
      enddo
      u_(N,1) = - a*(rho_(1,1)-rho_(N,1))/dx - bbb*(rho_(N,1)-rho_(N,N))/dy
      do j=2,N
        u_(N,j) = - a*(rho_(1,j)-rho_(N,j))/dx - bbb*(rho_(N,j)-rho_(N,j-1))/dy
      enddo
  else if ((a>=0).and.(bbb<0)) then
      do j=1,N-1
        u_(1,j) = - a*(rho_(1,j)-rho_(N,j))/dx - bbb*(rho_(1,j+1)-rho_(1,j))/dy
      enddo
      u_(1,N) = - a*(rho_(1,N)-rho_(N,N))/dx - bbb*(rho_(1,1)-rho_(1,N))/dy
      do i=2,N
        do j=1,N-1
          u_(i,j) = - a*(rho_(i,j)-rho_(i-1,j))/dx - bbb*(rho_(i,j+1)-rho_(i,j))/dy
        enddo
        u_(i,N) = - a*(rho_(i,N)-rho_(i-1,N))/dx - bbb*(rho_(i,1)-rho_(i,N))/dy
      enddo
  else
      do i=1,N-1
        do j=1,N-1
          u_(i,j) = - a*(rho_(i+1,j)-rho_(i,j))/dx - bbb*(rho_(i,j+1)-rho_(i,j))/dy
        enddo
        u_(i,N) = - a*(rho_(i+1,N)-rho_(i,N))/dx - bbb*(rho_(i,1)-rho_(i,N))/dy
      enddo
      do j=1,N-1
        u_(N,j) = - a*(rho_(1,j)-rho_(N,j))/dx - bbb*(rho_(N,j+1)-rho_(N,j))/dy
      enddo
      u_(N,N) = - a*(rho_(1,N)-rho_(N,N))/dx - bbb*(rho_(N,1)-rho_(N,N))/dy
  endif
end
subroutine computeFunctionUpwind2Dperiodic(alpha, rho_, t, u_, u, v)
  real(pr), intent(in) :: alpha, t
  real(pr), dimension(:,:), intent(inout) :: rho_
  real(pr), dimension(:,:), intent(inout) :: u_
  real(pr), dimension(:,:), intent(inout) :: u
  real(pr), dimension(:,:), intent(inout) :: v
  integer :: i, j, imod1m, imod1p, jmod1m, jmod1p

  do i=1,N
    do j=1,N
      imod1m = modulo(i-1-1,N)+1
      imod1p = modulo(i-1+1,N)+1
      jmod1m = modulo(j-1-1,N)+1
      jmod1p = modulo(j-1+1,N)+1

      if ((u(i,j)>=0).and.(v(i,j)>=0)) then
          u_(i,j) = - u(i,j)*(rho_(i,j)-rho_(imod1m,j))/dx - v(i,j)*(rho_(i,j)-rho_(i,jmod1m))/dy
      else if ((u(i,j)>=0).and.(v(i,j)<0)) then
          u_(i,j) = - u(i,j)*(rho_(i,j)-rho_(imod1m,j))/dx - v(i,j)*(rho_(i,jmod1p)-rho_(i,j))/dy
      else if ((u(i,j)<0).and.(v(i,j)>=0)) then
          u_(i,j) = - u(i,j)*(rho_(imod1p,j)-rho_(i,j))/dx - v(i,j)*(rho_(i,j)-rho_(i,jmod1m))/dy
      else
          u_(i,j) = - u(i,j)*(rho_(imod1p,j)-rho_(i,j))/dx - v(i,j)*(rho_(i,jmod1p)-rho_(i,j))/dy
      endif
    enddo
 enddo

end

subroutine computeFunctionUpwind2DdirichletH(alpha, rho_, t, u_, u, v)
  real(pr), intent(in) :: alpha, t
  real(pr), dimension(:,:), intent(inout) :: rho_
  real(pr), dimension(:,:), intent(inout) :: u_
  real(pr), dimension(:,:), intent(inout) :: u
  real(pr), dimension(:,:), intent(inout) :: v
  integer :: i, j, imod1m, imod1p, jmod1m, jmod1p

  do i=1,N
    do j=1,N
      imod1m = modulo(i-1-1,N)+1
      imod1p = modulo(i-1+1,N)+1
      jmod1m = modulo(j-1-1,N)+1
      jmod1p = modulo(j-1+1,N)+1

      if (i == 1) imod1m = 0._pr
      if (i == N) imod1p = 0._pr
      if (j == 1) jmod1m = 0._pr
      if (j == N) jmod1p = 0._pr

      if ((u(i,j)>0).and.(v(i,j)>0)) then
          u_(i,j) = - u(i,j)*(rho_(i,j)-rho_(imod1m,j))/dx - v(i,j)*(rho_(i,j)-rho_(i,jmod1m))/dy
      else if ((u(i,j)>0).and.(v(i,j)<0)) then
          u_(i,j) = - u(i,j)*(rho_(i,j)-rho_(imod1m,j))/dx - v(i,j)*(rho_(i,jmod1p)-rho_(i,j))/dy
      else if ((u(i,j)<0).and.(v(i,j)>0)) then
          u_(i,j) = - u(i,j)*(rho_(imod1p,j)-rho_(i,j))/dx - v(i,j)*(rho_(i,j)-rho_(i,jmod1m))/dy
      else
          u_(i,j) = - u(i,j)*(rho_(imod1p,j)-rho_(i,j))/dx - v(i,j)*(rho_(i,jmod1p)-rho_(i,j))/dy
      endif
    enddo
 enddo

end

!/*************************
! * Lax-Wendroff Advection Scheme *
! *************************/

!! Advances the scheme in time
subroutine AdvanceLWperiodic(nt, tn, rho, u, v)
  integer, intent(in) :: nt
  real(pr), intent(in) :: tn
  real(pr), dimension(:,:), intent(inout) :: rho, u, v
  real(pr), dimension(size(rho,1),size(rho,2)) :: rhoNext
  integer :: i, j


  !! 2D
  !! Lax-Wendroff : rho^n+1(i) = rho^n(i) - a*dt/(2*dx)*(rho^n(i+1) - rho^n(i-1)) + (a*dt)*(a*dt)/(2*dx*dx)*(rho^n(i+1) - 2*rho^n(i) + rho^n(i-1))
      rhoNext(1,1) = rho(1,1) - u(1,1)*dt/(2._pr*dx)*(rho(2,1) - rho(N,1)) + (u(1,1)*dt)*(u(1,1)*dt)/(2._pr*dx*dx)*(rho(2,1)&
 - 2._pr*rho(1,1) + rho(N,1)) - v(1,1)*dt/(2._pr*dy)*(rho(1,2) - rho(1,N)) + (v(1,1)*dt)*(v(1,1)*dt)/(2._pr*dy*dy)*(rho(1,2)&
 - 2._pr*rho(1,1) + rho(1,N)) + (u(1,1)*dt)*(v(1,1)*dt)*(rho(2,2)-rho(2,N)-rho(N,2)+rho(N,N))/(4._pr*dx*dy)
    do j=2,N-1
      rhoNext(1,j) = rho(1,j) - u(1,j)*dt/(2._pr*dx)*(rho(2,j) - rho(N,j)) + (u(1,j)*dt)*(u(1,j)*dt)/(2._pr*dx*dx)*(rho(2,j)&
 - 2._pr*rho(1,j) + rho(N,j)) - v(1,j)*dt/(2._pr*dy)*(rho(1,j+1) - rho(1,j-1)) + (v(1,j)*dt)*(v(1,j)*dt)/(2._pr*dy*dy)*(rho(1,j+1) &
- 2._pr*rho(1,j) + rho(1,j-1)) + (u(1,j)*dt)*(v(1,j)*dt)*(rho(2,j+1)-rho(2,j-1)-rho(1,j+1)+rho(1,j-1))/(4._pr*dx*dy)
    enddo
      rhoNext(1,N) = rho(1,N) - u(1,j)*dt/(2._pr*dx)*(rho(2,N) - rho(N,N)) + (u(1,j)*dt)*(u(1,j)*dt)/(2._pr*dx*dx)*(rho(2,N)&
 - 2._pr*rho(1,N) + rho(N,N)) - v(1,N)*dt/(2._pr*dy)*(rho(1,1) - rho(1,N-1)) + (v(1,N)*dt)*(v(1,N)*dt)/(2._pr*dy*dy)*(rho(1,1)&
 - 2._pr*rho(1,N) + rho(1,N-1)) + (u(1,j)*dt)*(v(1,N)*dt)*(rho(2,1)-rho(2,N-1)-rho(N,1)+rho(N,N-1))/(4._pr*dx*dy)


  do i=2,N-1
      rhoNext(i,1) = rho(i,1) - u(i,1)*dt/(2._pr*dx)*(rho(i+1,1) - rho(i-1,1)) + (u(i,1)*dt)*(u(i,1)*dt)/(2._pr*dx*dx)*(rho(i+1,1) &
- 2._pr*rho(i,1) + rho(i-1,1)) - v(i,1)*dt/(2._pr*dy)*(rho(i,2) - rho(i,N)) + (v(i,1)*dt)*(v(i,1)*dt)/(2._pr*dy*dy)*(rho(i,2)&
 - 2._pr*rho(i,1) + rho(i,N)) + (u(i,1)*dt)*(v(i,1)*dt)*(rho(i+1,2)-rho(i+1,N)-rho(i-1,2)+rho(i-1,N))/(4._pr*dx*dy)
    do j=2,N-1
      rhoNext(i,j) = rho(i,j) - u(i,j)*dt/(2._pr*dx)*(rho(i+1,j) - rho(i-1,j)) + (u(i,j)*dt)*(u(i,j)*dt)/(2._pr*dx*dx)*(rho(i+1,j) &
- 2._pr*rho(i,j) + rho(i-1,j)) - v(i,j)*dt/(2._pr*dy)*(rho(i,j+1) - rho(i,j-1)) + (v(i,j)*dt)*(v(i,j)*dt)/(2._pr*dy*dy)*(rho(i,j+1)&
 - 2._pr*rho(i,j) + rho(i,j-1)) + (u(i,j)*dt)*(v(i,j)*dt)*(rho(i+1,j+1)-rho(i+1,j-1)-rho(i-1,j+1)+rho(i-1,j-1))/(4._pr*dx*dy)
    enddo
      rhoNext(i,N) = rho(i,N) - u(i,N)*dt/(2._pr*dx)*(rho(i+1,N) - rho(i-1,N)) + (u(i,N)*dt)*(u(i,N)*dt)/(2._pr*dx*dx)*(rho(i+1,N) &
- 2._pr*rho(i,N) + rho(i-1,N)) - v(i,N)*dt/(2._pr*dy)*(rho(i,1) - rho(i,N-1)) + (v(i,N)*dt)*(v(i,N)*dt)/(2._pr*dy*dy)*(rho(i,1)&
 - 2._pr*rho(i,N) + rho(i,N-1)) + (u(i,N)*dt)*(v(i,N)*dt)*(rho(i+1,1)-rho(i+1,N-1)-rho(i-1,1)+rho(i-1,N-1))/(4._pr*dx*dy)
  enddo


      rhoNext(N,1) = rho(N,1) - u(N,1)*dt/(2._pr*dx)*(rho(1,1) - rho(N-1,1)) + (u(N,1)*dt)*(u(N,1)*dt)/(2._pr*dx*dx)*(rho(1,1)&
 - 2._pr*rho(N,1) + rho(N-1,1)) - v(N,1)*dt/(2._pr*dy)*(rho(N,2) - rho(N,N)) + (v(N,1)*dt)*(v(N,1)*dt)/(2._pr*dy*dy)*(rho(N,2)&
 - 2._pr*rho(N,1) + rho(N,N)) + (u(N,1)*dt)*(v(N,1)*dt)*(rho(1,2)-rho(1,N)-rho(N-1,2)+rho(N-1,N))/(4._pr*dx*dy)
    do j=2,N-1
      rhoNext(N,j) = rho(N,j) - u(N,j)*dt/(2._pr*dx)*(rho(1,j) - rho(N-1,j)) + (u(N,j)*dt)*(u(N,j)*dt)/(2._pr*dx*dx)*(rho(1,j)&
 - 2._pr*rho(N,j) + rho(N-1,j)) - v(N,j)*dt/(2._pr*dy)*(rho(N,j+1) - rho(N,j-1)) + (v(N,j)*dt)*(v(N,j)*dt)/(2._pr*dy*dy)*(rho(N,j+1&
) - 2._pr*rho(N,j) + rho(N,j-1)) + (u(N,j)*dt)*(v(N,j)*dt)*(rho(1,j+1)-rho(1,j-1)-rho(N-1,j+1)+rho(N-1,j-1))/(4._pr*dx*dy)
    enddo
      rhoNext(N,N) = rho(N,N) - u(N,N)*dt/(2._pr*dx)*(rho(1,N) - rho(N-1,N)) + (u(N,N)*dt)*(u(N,N)*dt)/(2._pr*dx*dx)*(rho(1,N)&
 - 2._pr*rho(N,N) + rho(N-1,N)) - v(N,N)*dt/(2._pr*dy)*(rho(N,1) - rho(N,N-1)) + (v(N,N)*dt)*(v(N,N)*dt)/(2._pr*dy*dy)*(rho(N,1)&
 - 2._pr*rho(N,N) + rho(N,N-1)) + (u(N,N)*dt)*(v(N,N)*dt)*(rho(1,1)-rho(1,N-1)-rho(N-1,1)+rho(N-1,N-1))/(4._pr*dx*dy)

  rho = rhoNext
end
subroutine AdvanceLWdirichletH(nt, tn, rho, u, v)
  integer, intent(in) :: nt
  real(pr), intent(in) :: tn
  real(pr), dimension(:,:), intent(inout) :: rho, u, v
  real(pr), dimension(size(rho,1),size(rho,2)) :: rhoNext
  integer :: i, j

  !! 2D
  !! Lax-Wendroff : rho^n+1(i) = rho^n(i) - a*dt/(2*dx)*(rho^n(i+1) - rho^n(i-1)) + (a*dt)*(a*dt)/(2*dx*dx)*(rho^n(i+1) - 2*rho^n(i) + rho^n(i-1))
      rhoNext(1,1) = rho(1,1) - u(1,1)*dt/(2._pr*dx)*rho(2,1) + (u(1,1)*dt)*(u(1,1)*dt)/(2._pr*dx*dx)*(rho(2,1) - 2._pr*rho(1,1)) -&
 v(1,1)*dt/(2._pr*dy)*(rho(1,2)) + (v(1,1)*dt)*(v(1,1)*dt)/(2._pr*dy*dy)*(rho(1,2) - 2._pr*rho(1,1)) + (u(1,1)*dt)*(v(1,1)*dt)*(rho&
(2,2))/(4._pr*dx*dy)
    do j=2,N-1
      rhoNext(1,j) = rho(1,j) - u(1,j)*dt/(2._pr*dx)*(rho(2,j)) + (u(1,j)*dt)*(u(1,j)*dt)/(2._pr*dx*dx)*(rho(2,j) - 2._pr*rho(1,j))&
 - v(1,j)*dt/(2._pr*dy)*(rho(1,j+1) - rho(1,j-1)) + (v(1,j)*dt)*(v(1,j)*dt)/(2._pr*dy*dy)*(rho(1,j+1) - 2._pr*rho(1,j) + rho(1,j-1)&
) + (u(1,j)*dt)*(v(1,j)*dt)*(rho(2,j+1)-rho(2,j-1)-rho(1,j+1)+rho(1,j-1))/(4._pr*dx*dy)
    enddo
      rhoNext(1,N) = rho(1,N) - u(1,j)*dt/(2._pr*dx)*(rho(2,N)) + (u(1,j)*dt)*(u(1,j)*dt)/(2._pr*dx*dx)*(rho(2,N) - 2._pr*rho(1,N))&
 - v(1,N)*dt/(2._pr*dy)*(- rho(1,N-1)) + (v(1,N)*dt)*(v(1,N)*dt)/(2._pr*dy*dy)*( - 2._pr*rho(1,N) + rho(1,N-1)) + (u(1,j)*dt)*(v(1,&
N)*dt)*(-rho(2,N-1))/(4._pr*dx*dy)


  do i=2,N-1
      rhoNext(i,1) = rho(i,1) - u(i,1)*dt/(2._pr*dx)*(rho(i+1,1) - rho(i-1,1)) + (u(i,1)*dt)*(u(i,1)*dt)/(2._pr*dx*dx)*(rho(i+1,1) &
- 2._pr*rho(i,1) + rho(i-1,1)) - v(i,1)*dt/(2._pr*dy)*rho(i,2) + (v(i,1)*dt)*(v(i,1)*dt)/(2._pr*dy*dy)*(rho(i,2) - 2._pr*rho(i,1)) &
+ (u(i,1)*dt)*(v(i,1)*dt)*(rho(i+1,2)-rho(i-1,2))/(4._pr*dx*dy)
    do j=2,N-1
      rhoNext(i,j) = rho(i,j) - u(i,j)*dt/(2._pr*dx)*(rho(i+1,j) - rho(i-1,j)) + (u(i,j)*dt)*(u(i,j)*dt)/(2._pr*dx*dx)*(rho(i+1,j) &
- 2._pr*rho(i,j) + rho(i-1,j)) - v(i,j)*dt/(2._pr*dy)*(rho(i,j+1) - rho(i,j-1)) + (v(i,j)*dt)*(v(i,j)*dt)/(2._pr*dy*dy)*(rho(i,j+1)&
 - 2._pr*rho(i,j) + rho(i,j-1)) + (u(i,j)*dt)*(v(i,j)*dt)*(rho(i+1,j+1)-rho(i+1,j-1)-rho(i-1,j+1)+rho(i-1,j-1))/(4._pr*dx*dy)
    enddo
      rhoNext(i,N) = rho(i,N) - u(i,N)*dt/(2._pr*dx)*(rho(i+1,N) - rho(i-1,N)) + (u(i,N)*dt)*(u(i,N)*dt)/(2._pr*dx*dx)*(rho(i+1,N) &
- 2._pr*rho(i,N) + rho(i-1,N)) - v(i,N)*dt/(2._pr*dy)*(- rho(i,N-1)) + (v(i,N)*dt)*(v(i,N)*dt)/(2._pr*dy*dy)*( - 2._pr*rho(i,N) + &
rho(i,N-1)) + (u(i,N)*dt)*(v(i,N)*dt)*(-rho(i+1,N-1)+rho(i-1,N-1))/(4._pr*dx*dy)
  enddo


      rhoNext(N,1) = rho(N,1) - u(N,1)*dt/(2._pr*dx)*( - rho(N-1,1)) + (u(N,1)*dt)*(u(N,1)*dt)/(2._pr*dx*dx)*( - 2._pr*rho(N,1) + &
rho(N-1,1)) - v(N,1)*dt/(2._pr*dy)*rho(N,2) + (v(N,1)*dt)*(v(N,1)*dt)/(2._pr*dy*dy)*(rho(N,2) - 2._pr*rho(N,1)) + (u(N,1)*dt)*(v(N,&
1)*dt)*(-rho(N-1,2))/(4._pr*dx*dy)
    do j=2,N-1
      rhoNext(N,j) = rho(N,j) - u(N,j)*dt/(2._pr*dx)*( - rho(N-1,j)) + (u(N,j)*dt)*(u(N,j)*dt)/(2._pr*dx*dx)*( - 2._pr*rho(N,j) + &
rho(N-1,j)) - v(N,j)*dt/(2._pr*dy)*(rho(N,j+1) - rho(N,j-1)) + (v(N,j)*dt)*(v(N,j)*dt)/(2._pr*dy*dy)*(rho(N,j+1) - 2._pr*rho(N,j) +&
 rho(N,j-1)) + (u(N,j)*dt)*(v(N,j)*dt)*(-rho(N-1,j+1)+rho(N-1,j-1))/(4._pr*dx*dy)
    enddo
      rhoNext(N,N) = rho(N,N) - u(N,N)*dt/(2._pr*dx)*( - rho(N-1,N)) + (u(N,N)*dt)*(u(N,N)*dt)/(2._pr*dx*dx)*( - 2._pr*rho(N,N) + &
rho(N-1,N)) - v(N,N)*dt/(2._pr*dy)*( - rho(N,N-1)) + (v(N,N)*dt)*(v(N,N)*dt)/(2._pr*dy*dy)*(- 2._pr*rho(N,N) + rho(N,N-1)) + (u(N,N&
)*dt)*(v(N,N)*dt)*(rho(N-1,N-1))/(4._pr*dx*dy)

  rho = rhoNext
end

!! 2D
!! Lax-Wendroff : rho^n+1(i) = rho^n(i) - a*dt/(2*dx)*(rho^n(i+1) - rho^n(i-1)) + (a*dt)*(a*dt)/(2*dx*dx)*(rho^n(i+1) - 2*rho^n(i) + rho^n(i-1))
subroutine addFunctionLW2Dold(alpha, rho_, t, u_)
  real(pr), intent(in) :: alpha, t
  real(pr), dimension(:,:), intent(inout) :: rho_
  real(pr), dimension(:,:), intent(inout) :: u_
  integer :: i, j

      u_(1,1) = u_(1,1) - a*alpha/(2._pr*dx)*(rho_(2,1) - rho_(N,1)) + (a*alpha)*(a*alpha)/(2._pr*dx*dx)*(rho_(2,1) - 2._pr*rho_(1,&
1) + rho_(N,1)) - bbb*alpha/(2._pr*dy)*(rho_(1,2) - rho_(1,N)) + (bbb*alpha)*(bbb*alpha)/(2._pr*dy*dy)*(rho_(1,2) - 2._pr*rho_(1,1)&
 + rho_(1,N))  + (a*alpha)*(bbb*alpha)*(rho_(2,2)-rho_(2,N)-rho_(N,2)+rho_(N,N))/(4._pr*dx*dy)
    do j=2,N-1
      u_(1,j) = u_(1,j) - a*alpha/(2._pr*dx)*(rho_(2,j) - rho_(N,j)) + (a*alpha)*(a*alpha)/(2._pr*dx*dx)*(rho_(2,j) - 2._pr*rho_(1,&
j) + rho_(N,j)) - bbb*alpha/(2._pr*dy)*(rho_(1,j+1) - rho_(1,j-1)) + (bbb*alpha)*(bbb*alpha)/(2._pr*dy*dy)*(rho_(1,j+1) - 2._pr*&
rho_(1,j) + rho_(1,j-1))  + (a*alpha)*(bbb*alpha)*(rho_(2,j+1)-rho_(2,j-1)-rho_(1,j+1)+rho_(1,j-1))/(4._pr*dx*dy)
    enddo
      u_(1,N) = u_(1,N) - a*alpha/(2._pr*dx)*(rho_(2,N) - rho_(N,N)) + (a*alpha)*(a*alpha)/(2._pr*dx*dx)*(rho_(2,N) - 2._pr*rho_(1,&
N) + rho_(N,N)) - bbb*alpha/(2._pr*dy)*(rho_(1,1) - rho_(1,N-1)) + (bbb*alpha)*(bbb*alpha)/(2._pr*dy*dy)*(rho_(1,1) - 2._pr*rho_(1,&
N) + rho_(1,N-1))  + (a*alpha)*(bbb*alpha)*(rho_(2,1)-rho_(2,N-1)-rho_(N,1)+rho_(N,N-1))/(4._pr*dx*dy)


  do i=2,N-1
      u_(i,1) = u_(i,1) - a*alpha/(2._pr*dx)*(rho_(i+1,1) - rho_(i-1,1)) + (a*alpha)*(a*alpha)/(2._pr*dx*dx)*(rho_(i+1,1) - 2._pr*&
rho_(i,1) + rho_(i-1,1)) - bbb*alpha/(2._pr*dy)*(rho_(i,2) - rho_(i,N)) + (bbb*alpha)*(bbb*alpha)/(2._pr*dy*dy)*(rho_(i,2) - 2._pr*&
rho_(i,1) + rho_(i,N))  + (a*alpha)*(bbb*alpha)*(rho_(i+1,2)-rho_(i+1,N)-rho_(i-1,2)+rho_(i-1,N))/(4._pr*dx*dy)
    do j=2,N-1
      u_(i,j) = u_(i,j) - a*alpha/(2._pr*dx)*(rho_(i+1,j) - rho_(i-1,j)) + (a*alpha)*(a*alpha)/(2._pr*dx*dx)*(rho_(i+1,j) - 2._pr*&
rho_(i,j) + rho_(i-1,j)) - bbb*alpha/(2._pr*dy)*(rho_(i,j+1) - rho_(i,j-1)) + (bbb*alpha)*(bbb*alpha)/(2._pr*dy*dy)*(rho_(i,j+1) - &
2._pr*rho_(i,j) + rho_(i,j-1))  + (a*alpha)*(bbb*alpha)*(rho_(i+1,j+1)-rho_(i+1,j-1)-rho_(i-1,j+1)+rho_(i-1,j-1))/(4._pr*dx*dy)
    enddo
      u_(i,N) = u_(i,N) - a*alpha/(2._pr*dx)*(rho_(i+1,N) - rho_(i-1,N)) + (a*alpha)*(a*alpha)/(2._pr*dx*dx)*(rho_(i+1,N) - 2._pr*&
rho_(i,N) + rho_(i-1,N)) - bbb*alpha/(2._pr*dy)*(rho_(i,1) - rho_(i,N-1)) + (bbb*alpha)*(bbb*alpha)/(2._pr*dy*dy)*(rho_(i,1)&
 - 2._pr*rho_(i,N) + rho_(i,N-1))  + (a*alpha)*(bbb*alpha)*(rho_(i+1,1)-rho_(i+1,N-1)-rho_(i-1,1)+rho_(i-1,N-1))/(4._pr*dx*dy)
  enddo


      u_(N,1) = u_(N,1) - a*alpha/(2._pr*dx)*(rho_(1,1) - rho_(N-1,1)) + (a*alpha)*(a*alpha)/(2._pr*dx*dx)*(rho_(1,1) - 2._pr*rho_(&
N,1) + rho_(N-1,1)) - bbb*alpha/(2._pr*dy)*(rho_(N,2) - rho_(N,N)) + (bbb*alpha)*(bbb*alpha)/(2._pr*dy*dy)*(rho_(N,2) - 2._pr*rho_(&
N,1) + rho_(N,N))  + (a*alpha)*(bbb*alpha)*(rho_(1,2)-rho_(1,N)-rho_(N-1,2)+rho_(N-1,N))/(4._pr*dx*dy)
    do j=2,N-1
      u_(N,j) = u_(N,j) - a*alpha/(2._pr*dx)*(rho_(1,j) - rho_(N-1,j)) + (a*alpha)*(a*alpha)/(2._pr*dx*dx)*(rho_(1,j) - 2._pr*rho_(&
N,j) + rho_(N-1,j)) - bbb*alpha/(2._pr*dy)*(rho_(N,j+1) - rho_(N,j-1)) + (bbb*alpha)*(bbb*alpha)/(2._pr*dy*dy)*(rho_(N,j+1) - 2._pr&
*rho_(N,j) + rho_(N,j-1))  + (a*alpha)*(bbb*alpha)*(rho_(1,j+1)-rho_(1,j-1)-rho_(N-1,j+1)+rho_(N-1,j-1))/(4._pr*dx*dy)
    enddo
      u_(N,N) = u_(N,N) - a*alpha/(2._pr*dx)*(rho_(1,N) - rho_(N-1,N)) + (a*alpha)*(a*alpha)/(2._pr*dx*dx)*(rho_(1,N) - 2._pr*rho_(&
N,N) + rho_(N-1,N)) - bbb*alpha/(2._pr*dy)*(rho_(N,1) - rho_(N,N-1)) + (bbb*alpha)*(bbb*alpha)/(2._pr*dy*dy)*(rho_(N,1) - 2._pr*&
rho_(N,N) + rho_(N,N-1))  + (a*alpha)*(bbb*alpha)*(rho_(1,1)-rho_(1,N-1)-rho_(N-1,1)+rho_(N-1,N-1))/(4._pr*dx*dy)
end


!/*************************
! * WENO5 Scheme *
! *************************/

!! 1D
subroutine computeFunctionWENO5(alpha, rho_, t, u_)
  real(pr), intent(in) :: alpha, t
  real(pr), dimension(:), intent(inout) :: rho_
  real(pr), dimension(:), intent(inout) :: u_
  real(pr) :: v1m, v2m, v3m, v4m, v5m
  real(pr) :: s1m, s2m, s3m
  real(pr) :: a1m, a2m, a3m
  real(pr) :: w1m, w2m, w3m
  real(pr) :: drhodxm
  real(pr) :: v1p, v2p, v3p, v4p, v5p
  real(pr) :: s1p, s2p, s3p
  real(pr) :: a1p, a2p, a3p
  real(pr) :: w1p, w2p, w3p
  real(pr) :: drhodxp
  integer :: i

  if (a>0) then
        v1m = (rho_(N-1) - rho_(N-2))/dx
        v2m = (rho_(N) - rho_(N-1))/dx
        v3m = (rho_(1) - rho_(N))/dx
        v4m = (rho_(1+1) - rho_(1))/dx
        v5m = (rho_(1+2) - rho_(1+1))/dx
        s1m = 13._pr/12._pr*(v1m - 2._pr*v2m + v3m)*(v1m - 2._pr*v2m + v3m) + 0.25_pr*(v1m - 4._pr*v2m + 3._pr*v3m)*(v1m - 4._pr*&
v2m + 3._pr*v3m)
        s2m = 13._pr/12._pr*(v2m - 2._pr*v3m + v4m)*(v2m - 2._pr*v3m + v4m) + 0.25_pr*(v2m - v4m)*(v2m - v4m)
        s3m = 13._pr/12._pr*(v3m - 2._pr*v4m + v5m)*(v3m - 2._pr*v4m + v5m) + 0.25_pr*(3._pr*v3m - 4._pr*v4m + v5m)*(3._pr*v3m&
 - 4._pr*v4m + v5m)
        a1m = 0.1_pr/((eepsilon + s1m)*(eepsilon + s1m))
        a2m = 0.6_pr/((eepsilon + s2m)*(eepsilon + s2m))
        a3m = 0.3_pr/((eepsilon + s3m)*(eepsilon + s3m))
        w1m = a1m/(a1m + a2m + a3m)
        w2m = a2m/(a1m + a2m + a3m)
        w3m = a3m/(a1m + a2m + a3m)
        drhodxm = w1m*(v1m/3._pr - 7._pr/6._pr*v2m + 11._pr/6._pr*v3m) + w2m*(-v2m/6._pr + 5._pr/6._pr*v3m + v4m/3._pr) + w3m*(v3m/&
3._pr + 5._pr/6._pr*v4m - v5m/6._pr)
        u_(1) = - a*drhodxm

        v1m = (rho_(N) - rho_(N-1))/dx
        v2m = (rho_(2-1) - rho_(N))/dx
        v3m = (rho_(2) - rho_(2-1))/dx
        v4m = (rho_(2+1) - rho_(2))/dx
        v5m = (rho_(2+2) - rho_(2+1))/dx
        s1m = 13._pr/12._pr*(v1m - 2._pr*v2m + v3m)*(v1m - 2._pr*v2m + v3m) + 0.25_pr*(v1m - 4._pr*v2m + 3._pr*v3m)*(v1m - 4._pr*&
v2m + 3._pr*v3m)
        s2m = 13._pr/12._pr*(v2m - 2._pr*v3m + v4m)*(v2m - 2._pr*v3m + v4m) + 0.25_pr*(v2m - v4m)*(v2m - v4m)
        s3m = 13._pr/12._pr*(v3m - 2._pr*v4m + v5m)*(v3m - 2._pr*v4m + v5m) + 0.25_pr*(3._pr*v3m - 4._pr*v4m + v5m)*(3._pr*v3m&
 - 4._pr*v4m + v5m)
        a1m = 0.1_pr/((eepsilon + s1m)*(eepsilon + s1m))
        a2m = 0.6_pr/((eepsilon + s2m)*(eepsilon + s2m))
        a3m = 0.3_pr/((eepsilon + s3m)*(eepsilon + s3m))
        w1m = a1m/(a1m + a2m + a3m)
        w2m = a2m/(a1m + a2m + a3m)
        w3m = a3m/(a1m + a2m + a3m)
        drhodxm = w1m*(v1m/3._pr - 7._pr/6._pr*v2m + 11._pr/6._pr*v3m) + w2m*(-v2m/6._pr + 5._pr/6._pr*v3m + v4m/3._pr) + w3m*(v3m/&
3._pr + 5._pr/6._pr*v4m - v5m/6._pr)
        u_(2) = - a*drhodxm

        v1m = (rho_(3-2) - rho_(N))/dx
        v2m = (rho_(3-1) - rho_(3-2))/dx
        v3m = (rho_(3) - rho_(3-1))/dx
        v4m = (rho_(3+1) - rho_(3))/dx
        v5m = (rho_(3+2) - rho_(3+1))/dx
        s1m = 13._pr/12._pr*(v1m - 2._pr*v2m + v3m)*(v1m - 2._pr*v2m + v3m) + 0.25_pr*(v1m - 4._pr*v2m + 3._pr*v3m)*(v1m - 4._pr*&
v2m + 3._pr*v3m)
        s2m = 13._pr/12._pr*(v2m - 2._pr*v3m + v4m)*(v2m - 2._pr*v3m + v4m) + 0.25_pr*(v2m - v4m)*(v2m - v4m)
        s3m = 13._pr/12._pr*(v3m - 2._pr*v4m + v5m)*(v3m - 2._pr*v4m + v5m) + 0.25_pr*(3._pr*v3m - 4._pr*v4m + v5m)*(3._pr*v3m&
 - 4._pr*v4m + v5m)
        a1m = 0.1_pr/((eepsilon + s1m)*(eepsilon + s1m))
        a2m = 0.6_pr/((eepsilon + s2m)*(eepsilon + s2m))
        a3m = 0.3_pr/((eepsilon + s3m)*(eepsilon + s3m))
        w1m = a1m/(a1m + a2m + a3m)
        w2m = a2m/(a1m + a2m + a3m)
        w3m = a3m/(a1m + a2m + a3m)
        drhodxm = w1m*(v1m/3._pr - 7._pr/6._pr*v2m + 11._pr/6._pr*v3m) + w2m*(-v2m/6._pr + 5._pr/6._pr*v3m + v4m/3._pr) + w3m*(v3m/&
3._pr + 5._pr/6._pr*v4m - v5m/6._pr)
        u_(3) = - a*drhodxm
      do i=4,N-2
        v1m = (rho_(i-2) - rho_(i-3))/dx
        v2m = (rho_(i-1) - rho_(i-2))/dx
        v3m = (rho_(i) - rho_(i-1))/dx
        v4m = (rho_(i+1) - rho_(i))/dx
        v5m = (rho_(i+2) - rho_(i+1))/dx

        s1m = 13._pr/12._pr*(v1m - 2._pr*v2m + v3m)*(v1m - 2._pr*v2m + v3m) + 0.25_pr*(v1m - 4._pr*v2m + 3._pr*v3m)*(v1m - 4._pr*&
v2m + 3._pr*v3m)
        s2m = 13._pr/12._pr*(v2m - 2._pr*v3m + v4m)*(v2m - 2._pr*v3m + v4m) + 0.25_pr*(v2m - v4m)*(v2m - v4m)
        s3m = 13._pr/12._pr*(v3m - 2._pr*v4m + v5m)*(v3m - 2._pr*v4m + v5m) + 0.25_pr*(3._pr*v3m - 4._pr*v4m + v5m)*(3._pr*v3m&
 - 4._pr*v4m + v5m)

        a1m = 0.1_pr/((eepsilon + s1m)*(eepsilon + s1m))
        a2m = 0.6_pr/((eepsilon + s2m)*(eepsilon + s2m))
        a3m = 0.3_pr/((eepsilon + s3m)*(eepsilon + s3m))

        w1m = a1m/(a1m + a2m + a3m)
        w2m = a2m/(a1m + a2m + a3m)
        w3m = a3m/(a1m + a2m + a3m)

        drhodxm = w1m*(v1m/3._pr - 7._pr/6._pr*v2m + 11._pr/6._pr*v3m) + w2m*(-v2m/6._pr + 5._pr/6._pr*v3m + v4m/3._pr) + w3m*(v3m/&
3._pr + 5._pr/6._pr*v4m - v5m/6._pr)

        u_(i) = - a*drhodxm
      enddo
        v1m = (rho_(N-1-2) - rho_(N-1-3))/dx
        v2m = (rho_(N-1-1) - rho_(N-1-2))/dx
        v3m = (rho_(N-1) - rho_(N-1-1))/dx
        v4m = (rho_(N-1+1) - rho_(N-1))/dx
        v5m = (rho_(1) - rho_(N-1+1))/dx
        s1m = 13._pr/12._pr*(v1m - 2._pr*v2m + v3m)*(v1m - 2._pr*v2m + v3m) + 0.25_pr*(v1m - 4._pr*v2m + 3._pr*v3m)*(v1m - 4._pr*&
v2m + 3._pr*v3m)
        s2m = 13._pr/12._pr*(v2m - 2._pr*v3m + v4m)*(v2m - 2._pr*v3m + v4m) + 0.25_pr*(v2m - v4m)*(v2m - v4m)
        s3m = 13._pr/12._pr*(v3m - 2._pr*v4m + v5m)*(v3m - 2._pr*v4m + v5m) + 0.25_pr*(3._pr*v3m - 4._pr*v4m + v5m)*(3._pr*v3m&
 - 4._pr*v4m + v5m)
        a1m = 0.1_pr/((eepsilon + s1m)*(eepsilon + s1m))
        a2m = 0.6_pr/((eepsilon + s2m)*(eepsilon + s2m))
        a3m = 0.3_pr/((eepsilon + s3m)*(eepsilon + s3m))
        w1m = a1m/(a1m + a2m + a3m)
        w2m = a2m/(a1m + a2m + a3m)
        w3m = a3m/(a1m + a2m + a3m)
        drhodxm = w1m*(v1m/3._pr - 7._pr/6._pr*v2m + 11._pr/6._pr*v3m) + w2m*(-v2m/6._pr + 5._pr/6._pr*v3m + v4m/3._pr) + w3m*(v3m/&
3._pr + 5._pr/6._pr*v4m - v5m/6._pr)
        u_(N-1) = - a*drhodxm

        v1m = (rho_(N-2) - rho_(N-3))/dx
        v2m = (rho_(N-1) - rho_(N-2))/dx
        v3m = (rho_(N) - rho_(N-1))/dx
        v4m = (rho_(1) - rho_(N))/dx
        v5m = (rho_(2) - rho_(1))/dx
        s1m = 13._pr/12._pr*(v1m - 2._pr*v2m + v3m)*(v1m - 2._pr*v2m + v3m) + 0.25_pr*(v1m - 4._pr*v2m + 3._pr*v3m)*(v1m - 4._pr*&
v2m + 3._pr*v3m)
        s2m = 13._pr/12._pr*(v2m - 2._pr*v3m + v4m)*(v2m - 2._pr*v3m + v4m) + 0.25_pr*(v2m - v4m)*(v2m - v4m)
        s3m = 13._pr/12._pr*(v3m - 2._pr*v4m + v5m)*(v3m - 2._pr*v4m + v5m) + 0.25_pr*(3._pr*v3m - 4._pr*v4m + v5m)*(3._pr*v3m&
 - 4._pr*v4m + v5m)
        a1m = 0.1_pr/((eepsilon + s1m)*(eepsilon + s1m))
        a2m = 0.6_pr/((eepsilon + s2m)*(eepsilon + s2m))
        a3m = 0.3_pr/((eepsilon + s3m)*(eepsilon + s3m))
        w1m = a1m/(a1m + a2m + a3m)
        w2m = a2m/(a1m + a2m + a3m)
        w3m = a3m/(a1m + a2m + a3m)
        drhodxm = w1m*(v1m/3._pr - 7._pr/6._pr*v2m + 11._pr/6._pr*v3m) + w2m*(-v2m/6._pr + 5._pr/6._pr*v3m + v4m/3._pr) + w3m*(v3m/&
3._pr + 5._pr/6._pr*v4m - v5m/6._pr)
        u_(N) = - a*drhodxm
  else
        v1p = (rho_(1+3) - rho_(1+2))/dx
        v2p = (rho_(1+2) - rho_(1+1))/dx
        v3p = (rho_(1+1) - rho_(1))/dx
        v4p = (rho_(1) - rho_(N))/dx
        v5p = (rho_(N) - rho_(N-1))/dx
        s1p = 13._pr/12._pr*(v1p - 2._pr*v2p + v3p)*(v1p - 2._pr*v2p + v3p) + 0.25_pr*(v1p - 4._pr*v2p + 3._pr*v3p)*(v1p - 4._pr*&
v2p + 3._pr*v3p)
        s2p = 13._pr/12._pr*(v2p - 2._pr*v3p + v4p)*(v2p - 2._pr*v3p + v4p) + 0.25_pr*(v2p - v4p)*(v2p - v4p)
        s3p = 13._pr/12._pr*(v3p - 2._pr*v4p + v5p)*(v3p - 2._pr*v4p + v5p) + 0.25_pr*(3._pr*v3p - 4._pr*v4p + v5p)*(3._pr*v3p&
 - 4._pr*v4p + v5p)
        a1p = 0.1_pr/((eepsilon + s1p)*(eepsilon + s1p))
        a2p = 0.6_pr/((eepsilon + s2p)*(eepsilon + s2p))
        a3p = 0.3_pr/((eepsilon + s3p)*(eepsilon + s3p))
        w1p = a1p/(a1p + a2p + a3p)
        w2p = a2p/(a1p + a2p + a3p)
        w3p = a3p/(a1p + a2p + a3p)
        drhodxp = w1p*(v1p/3._pr - 7._pr/6._pr*v2p + 11._pr/6._pr*v3p) + w2p*(-v2p/6._pr + 5._pr/6._pr*v3p + v4p/3._pr) + w3p*(v3p/&
3._pr + 5._pr/6._pr*v4p - v5p/6._pr)
        u_(1) = - a*drhodxp

        v1p = (rho_(2+3) - rho_(2+2))/dx
        v2p = (rho_(2+2) - rho_(2+1))/dx
        v3p = (rho_(2+1) - rho_(2))/dx
        v4p = (rho_(2) - rho_(2-1))/dx
        v5p = (rho_(2-1) - rho_(N))/dx
        s1p = 13._pr/12._pr*(v1p - 2._pr*v2p + v3p)*(v1p - 2._pr*v2p + v3p) + 0.25_pr*(v1p - 4._pr*v2p + 3._pr*v3p)*(v1p - 4._pr*&
v2p + 3._pr*v3p)
        s2p = 13._pr/12._pr*(v2p - 2._pr*v3p + v4p)*(v2p - 2._pr*v3p + v4p) + 0.25_pr*(v2p - v4p)*(v2p - v4p)
        s3p = 13._pr/12._pr*(v3p - 2._pr*v4p + v5p)*(v3p - 2._pr*v4p + v5p) + 0.25_pr*(3._pr*v3p - 4._pr*v4p + v5p)*(3._pr*v3p&
 - 4._pr*v4p + v5p)
        a1p = 0.1_pr/((eepsilon + s1p)*(eepsilon + s1p))
        a2p = 0.6_pr/((eepsilon + s2p)*(eepsilon + s2p))
        a3p = 0.3_pr/((eepsilon + s3p)*(eepsilon + s3p))
        w1p = a1p/(a1p + a2p + a3p)
        w2p = a2p/(a1p + a2p + a3p)
        w3p = a3p/(a1p + a2p + a3p)
        drhodxp = w1p*(v1p/3._pr - 7._pr/6._pr*v2p + 11._pr/6._pr*v3p) + w2p*(-v2p/6._pr + 5._pr/6._pr*v3p + v4p/3._pr) + w3p*(v3p/&
3._pr + 5._pr/6._pr*v4p - v5p/6._pr)
        u_(2) = - a*drhodxp
      do i=3,N-3
        v1p = (rho_(i+3) - rho_(i+2))/dx
        v2p = (rho_(i+2) - rho_(i+1))/dx
        v3p = (rho_(i+1) - rho_(i))/dx
        v4p = (rho_(i) - rho_(i-1))/dx
        v5p = (rho_(i-1) - rho_(i-2))/dx

        s1p = 13._pr/12._pr*(v1p - 2._pr*v2p + v3p)*(v1p - 2._pr*v2p + v3p) + 0.25_pr*(v1p - 4._pr*v2p + 3._pr*v3p)*(v1p - 4._pr*&
v2p + 3._pr*v3p)
        s2p = 13._pr/12._pr*(v2p - 2._pr*v3p + v4p)*(v2p - 2._pr*v3p + v4p) + 0.25_pr*(v2p - v4p)*(v2p - v4p)
        s3p = 13._pr/12._pr*(v3p - 2._pr*v4p + v5p)*(v3p - 2._pr*v4p + v5p) + 0.25_pr*(3._pr*v3p - 4._pr*v4p + v5p)*(3._pr*v3p&
 - 4._pr*v4p + v5p)

        a1p = 0.1_pr/((eepsilon + s1p)*(eepsilon + s1p))
        a2p = 0.6_pr/((eepsilon + s2p)*(eepsilon + s2p))
        a3p = 0.3_pr/((eepsilon + s3p)*(eepsilon + s3p))

        w1p = a1p/(a1p + a2p + a3p)
        w2p = a2p/(a1p + a2p + a3p)
        w3p = a3p/(a1p + a2p + a3p)

        drhodxp = w1p*(v1p/3._pr - 7._pr/6._pr*v2p + 11._pr/6._pr*v3p) + w2p*(-v2p/6._pr + 5._pr/6._pr*v3p + v4p/3._pr) + w3p*(v3p/&
3._pr + 5._pr/6._pr*v4p - v5p/6._pr)

        u_(i) = - a*drhodxp
      enddo
        v1p = (rho_(1) - rho_(N-2+2))/dx
        v2p = (rho_(N-2+2) - rho_(N-2+1))/dx
        v3p = (rho_(N-2+1) - rho_(N-2))/dx
        v4p = (rho_(N-2) - rho_(N-2-1))/dx
        v5p = (rho_(N-2-1) - rho_(N-2-2))/dx
        s1p = 13._pr/12._pr*(v1p - 2._pr*v2p + v3p)*(v1p - 2._pr*v2p + v3p) + 0.25_pr*(v1p - 4._pr*v2p + 3._pr*v3p)*(v1p - 4._pr*&
v2p + 3._pr*v3p)
        s2p = 13._pr/12._pr*(v2p - 2._pr*v3p + v4p)*(v2p - 2._pr*v3p + v4p) + 0.25_pr*(v2p - v4p)*(v2p - v4p)
        s3p = 13._pr/12._pr*(v3p - 2._pr*v4p + v5p)*(v3p - 2._pr*v4p + v5p) + 0.25_pr*(3._pr*v3p - 4._pr*v4p + v5p)*(3._pr*v3p&
 - 4._pr*v4p + v5p)
        a1p = 0.1_pr/((eepsilon + s1p)*(eepsilon + s1p))
        a2p = 0.6_pr/((eepsilon + s2p)*(eepsilon + s2p))
        a3p = 0.3_pr/((eepsilon + s3p)*(eepsilon + s3p))
        w1p = a1p/(a1p + a2p + a3p)
        w2p = a2p/(a1p + a2p + a3p)
        w3p = a3p/(a1p + a2p + a3p)
        drhodxp = w1p*(v1p/3._pr - 7._pr/6._pr*v2p + 11._pr/6._pr*v3p) + w2p*(-v2p/6._pr + 5._pr/6._pr*v3p + v4p/3._pr) + w3p*(v3p/&
3._pr + 5._pr/6._pr*v4p - v5p/6._pr)
        u_(N-2) = - a*drhodxp

        v1p = (rho_(2) - rho_(1))/dx
        v2p = (rho_(1) - rho_(N-1+1))/dx
        v3p = (rho_(N-1+1) - rho_(N-1))/dx
        v4p = (rho_(N-1) - rho_(N-1-1))/dx
        v5p = (rho_(N-1-1) - rho_(N-1-2))/dx
        s1p = 13._pr/12._pr*(v1p - 2._pr*v2p + v3p)*(v1p - 2._pr*v2p + v3p) + 0.25_pr*(v1p - 4._pr*v2p + 3._pr*v3p)*(v1p - 4._pr*&
v2p + 3._pr*v3p)
        s2p = 13._pr/12._pr*(v2p - 2._pr*v3p + v4p)*(v2p - 2._pr*v3p + v4p) + 0.25_pr*(v2p - v4p)*(v2p - v4p)
        s3p = 13._pr/12._pr*(v3p - 2._pr*v4p + v5p)*(v3p - 2._pr*v4p + v5p) + 0.25_pr*(3._pr*v3p - 4._pr*v4p + v5p)*(3._pr*v3p&
 - 4._pr*v4p + v5p)
        a1p = 0.1_pr/((eepsilon + s1p)*(eepsilon + s1p))
        a2p = 0.6_pr/((eepsilon + s2p)*(eepsilon + s2p))
        a3p = 0.3_pr/((eepsilon + s3p)*(eepsilon + s3p))
        w1p = a1p/(a1p + a2p + a3p)
        w2p = a2p/(a1p + a2p + a3p)
        w3p = a3p/(a1p + a2p + a3p)
        drhodxp = w1p*(v1p/3._pr - 7._pr/6._pr*v2p + 11._pr/6._pr*v3p) + w2p*(-v2p/6._pr + 5._pr/6._pr*v3p + v4p/3._pr) + w3p*(v3p/&
3._pr + 5._pr/6._pr*v4p - v5p/6._pr)
        u_(N-1) = - a*drhodxp

        v1p = (rho_(3) - rho_(2))/dx
        v2p = (rho_(2) - rho_(1))/dx
        v3p = (rho_(1) - rho_(N))/dx
        v4p = (rho_(N) - rho_(N-1))/dx
        v5p = (rho_(N-1) - rho_(N-2))/dx
        s1p = 13._pr/12._pr*(v1p - 2._pr*v2p + v3p)*(v1p - 2._pr*v2p + v3p) + 0.25_pr*(v1p - 4._pr*v2p + 3._pr*v3p)*(v1p - 4._pr*&
v2p + 3._pr*v3p)
        s2p = 13._pr/12._pr*(v2p - 2._pr*v3p + v4p)*(v2p - 2._pr*v3p + v4p) + 0.25_pr*(v2p - v4p)*(v2p - v4p)
        s3p = 13._pr/12._pr*(v3p - 2._pr*v4p + v5p)*(v3p - 2._pr*v4p + v5p) + 0.25_pr*(3._pr*v3p - 4._pr*v4p + v5p)*(3._pr*v3p&
 - 4._pr*v4p + v5p)
        a1p = 0.1_pr/((eepsilon + s1p)*(eepsilon + s1p))
        a2p = 0.6_pr/((eepsilon + s2p)*(eepsilon + s2p))
        a3p = 0.3_pr/((eepsilon + s3p)*(eepsilon + s3p))
        w1p = a1p/(a1p + a2p + a3p)
        w2p = a2p/(a1p + a2p + a3p)
        w3p = a3p/(a1p + a2p + a3p)
        drhodxp = w1p*(v1p/3._pr - 7._pr/6._pr*v2p + 11._pr/6._pr*v3p) + w2p*(-v2p/6._pr + 5._pr/6._pr*v3p + v4p/3._pr) + w3p*(v3p/&
3._pr + 5._pr/6._pr*v4p - v5p/6._pr)
        u_(N) = - a*drhodxp
  endif
end

!! 2D
subroutine computeFunctionWENO52Dperiodicoff(alpha, rho_, t, u_, u, v)
  real(pr), intent(in) :: alpha, t
  real(pr), dimension(:,:), intent(inout) :: rho_
  real(pr), dimension(:,:), intent(inout) :: u_, u, v
  real(pr) :: v1xm, v2xm, v3xm, v4xm, v5xm
  real(pr) :: s1xm, s2xm, s3xm
  real(pr) :: a1xm, a2xm, a3xm
  real(pr) :: w1xm, w2xm, w3xm
  real(pr) :: drhodxm
  real(pr) :: v1xp, v2xp, v3xp, v4xp, v5xp
  real(pr) :: s1xp, s2xp, s3xp
  real(pr) :: a1xp, a2xp, a3xp
  real(pr) :: w1xp, w2xp, w3xp
  real(pr) :: drhodxp
  real(pr) :: v1ym, v2ym, v3ym, v4ym, v5ym
  real(pr) :: s1ym, s2ym, s3ym
  real(pr) :: a1ym, a2ym, a3ym
  real(pr) :: w1ym, w2ym, w3ym
  real(pr) :: drhodym
  real(pr) :: v1yp, v2yp, v3yp, v4yp, v5yp
  real(pr) :: s1yp, s2yp, s3yp
  real(pr) :: a1yp, a2yp, a3yp
  real(pr) :: w1yp, w2yp, w3yp
  real(pr) :: drhodyp
  integer :: i, j, imod1m, imod1p, jmod1m, jmod1p, imod2m, imod2p, jmod2m, jmod2p, imod3m, imod3p, jmod3m, jmod3p

  if ((a>=0).and.(bbb>=0)) then
      do i=1,N
        do j=1,N
          imod1m = modulo(i-1-1,N)+1
          imod2m = modulo(i-1-2,N)+1
          imod3m = modulo(i-1-3,N)+1
          imod1p = modulo(i-1+1,N)+1
          imod2p = modulo(i-1+2,N)+1
          imod3p = modulo(i-1+3,N)+1
          v1xm = (rho_(imod2m,j) - rho_(imod3m,j))/dx
          v2xm = (rho_(imod1m,j) - rho_(imod2m,j))/dx
          v3xm = (rho_(i,j) - rho_(imod1m,j))/dx
          v4xm = (rho_(imod1p,j) - rho_(i,j))/dx
          v5xm = (rho_(imod2p,j) - rho_(imod1p,j))/dx
          s1xm = 13._pr/12._pr*(v1xm - 2._pr*v2xm + v3xm)*(v1xm - 2._pr*v2xm + v3xm) + 0.25_pr*(v1xm - 4._pr*v2xm + 3._pr*v3xm)*(&
v1xm - 4._pr*v2xm + 3._pr*v3xm)
          s2xm = 13._pr/12._pr*(v2xm - 2._pr*v3xm + v4xm)*(v2xm - 2._pr*v3xm + v4xm) + 0.25_pr*(v2xm - v4xm)*(v2xm - v4xm)
          s3xm = 13._pr/12._pr*(v3xm - 2._pr*v4xm + v5xm)*(v3xm - 2._pr*v4xm + v5xm) + 0.25_pr*(3._pr*v3xm - 4._pr*v4xm + v5xm)*&
(3._pr*v3xm - 4._pr*v4xm + v5xm)
          a1xm = 0.1_pr/((eepsilon + s1xm)*(eepsilon + s1xm))
          a2xm = 0.6_pr/((eepsilon + s2xm)*(eepsilon + s2xm))
          a3xm = 0.3_pr/((eepsilon + s3xm)*(eepsilon + s3xm))
          w1xm = a1xm/(a1xm + a2xm + a3xm)
          w2xm = a2xm/(a1xm + a2xm + a3xm)
          w3xm = a3xm/(a1xm + a2xm + a3xm)
          drhodxm = w1xm*(v1xm/3._pr - 7._pr/6._pr*v2xm + 11._pr/6._pr*v3xm) + w2xm*(-v2xm/6._pr + 5._pr/6._pr*v3xm + v4xm/3._pr) +&
 w3xm*(v3xm/3._pr + 5._pr/6._pr*v4xm - v5xm/6._pr)

          jmod1m = modulo(j-1-1,N)+1
          jmod2m = modulo(j-1-2,N)+1
          jmod3m = modulo(j-1-3,N)+1
          jmod1p = modulo(j-1+1,N)+1
          jmod2p = modulo(j-1+2,N)+1
          jmod3p = modulo(j-1+3,N)+1

          v1ym = (rho_(i,jmod2m) - rho_(i,jmod3m))/dy
          v2ym = (rho_(i,jmod1m) - rho_(i,jmod2m))/dy
          v3ym = (rho_(i,j) - rho_(i,jmod1m))/dy
          v4ym = (rho_(i,jmod1p) - rho_(i,j))/dy
          v5ym = (rho_(i,jmod2p) - rho_(i,jmod1p))/dy
          s1ym = 13._pr/12._pr*(v1ym - 2._pr*v2ym + v3ym)*(v1ym - 2._pr*v2ym + v3ym) + 0.25_pr*(v1ym - 4._pr*v2ym + 3._pr*v3ym)*(&
v1ym - 4._pr*v2ym + 3._pr*v3ym)
          s2ym = 13._pr/12._pr*(v2ym - 2._pr*v3ym + v4ym)*(v2ym - 2._pr*v3ym + v4ym) + 0.25_pr*(v2ym - v4ym)*(v2ym - v4ym)
          s3ym = 13._pr/12._pr*(v3ym - 2._pr*v4ym + v5ym)*(v3ym - 2._pr*v4ym + v5ym) + 0.25_pr*(3._pr*v3ym - 4._pr*v4ym + v5ym)*&
(3._pr*v3ym - 4._pr*v4ym + v5ym)
          a1ym = 0.1_pr/((eepsilon + s1ym)*(eepsilon + s1ym))
          a2ym = 0.6_pr/((eepsilon + s2ym)*(eepsilon + s2ym))
          a3ym = 0.3_pr/((eepsilon + s3ym)*(eepsilon + s3ym))
          w1ym = a1ym/(a1ym + a2ym + a3ym)
          w2ym = a2ym/(a1ym + a2ym + a3ym)
          w3ym = a3ym/(a1ym + a2ym + a3ym)
          drhodym = w1ym*(v1ym/3._pr - 7._pr/6._pr*v2ym + 11._pr/6._pr*v3ym) + w2ym*(-v2ym/6._pr + 5._pr/6._pr*v3ym + v4ym/3._pr) +&
 w3ym*(v3ym/3._pr + 5._pr/6._pr*v4ym - v5ym/6._pr)

          u_(i,j) = - a*drhodxm - bbb*drhodym
        enddo
      enddo
  else if ((a<0).and.(bbb>=0)) then
      do i=1,N
        do j=1,N
          imod1m = modulo(i-1-1,N)+1
          imod2m = modulo(i-1-2,N)+1
          imod3m = modulo(i-1-3,N)+1
          imod1p = modulo(i-1+1,N)+1
          imod2p = modulo(i-1+2,N)+1
          imod3p = modulo(i-1+3,N)+1

          v1xp = (rho_(imod3p,j) - rho_(imod2p,j))/dx
          v2xp = (rho_(imod2p,j) - rho_(imod1p,j))/dx
          v3xp = (rho_(imod1p,j) - rho_(i,j))/dx
          v4xp = (rho_(i,j) - rho_(imod1m,j))/dx
          v5xp = (rho_(imod1m,j) - rho_(imod2m,j))/dx
          s1xp = 13._pr/12._pr*(v1xp - 2._pr*v2xp + v3xp)*(v1xp - 2._pr*v2xp + v3xp) + 0.25_pr*(v1xp - 4._pr*v2xp + 3._pr*v3xp)*(&
v1xp - 4._pr*v2xp + 3._pr*v3xp)
          s2xp = 13._pr/12._pr*(v2xp - 2._pr*v3xp + v4xp)*(v2xp - 2._pr*v3xp + v4xp) + 0.25_pr*(v2xp - v4xp)*(v2xp - v4xp)
          s3xp = 13._pr/12._pr*(v3xp - 2._pr*v4xp + v5xp)*(v3xp - 2._pr*v4xp + v5xp) + 0.25_pr*(3._pr*v3xp - 4._pr*v4xp + v5xp)*&
(3._pr*v3xp - 4._pr*v4xp + v5xp)
          a1xp = 0.1_pr/((eepsilon + s1xp)*(eepsilon + s1xp))
          a2xp = 0.6_pr/((eepsilon + s2xp)*(eepsilon + s2xp))
          a3xp = 0.3_pr/((eepsilon + s3xp)*(eepsilon + s3xp))
          w1xp = a1xp/(a1xp + a2xp + a3xp)
          w2xp = a2xp/(a1xp + a2xp + a3xp)
          w3xp = a3xp/(a1xp + a2xp + a3xp)
          drhodxp = w1xp*(v1xp/3._pr - 7._pr/6._pr*v2xp + 11._pr/6._pr*v3xp) + w2xp*(-v2xp/6._pr + 5._pr/6._pr*v3xp + v4xp/3._pr) +&
 w3xp*(v3xp/3._pr + 5._pr/6._pr*v4xp - v5xp/6._pr)

          jmod1m = modulo(j-1-1,N)+1
          jmod2m = modulo(j-1-2,N)+1
          jmod3m = modulo(j-1-3,N)+1
          jmod1p = modulo(j-1+1,N)+1
          jmod2p = modulo(j-1+2,N)+1
          jmod3p = modulo(j-1+3,N)+1

          v1ym = (rho_(i,jmod2m) - rho_(i,jmod3m))/dy
          v2ym = (rho_(i,jmod1m) - rho_(i,jmod2m))/dy
          v3ym = (rho_(i,j) - rho_(i,jmod1m))/dy
          v4ym = (rho_(i,jmod1p) - rho_(i,j))/dy
          v5ym = (rho_(i,jmod2p) - rho_(i,jmod1p))/dy
          s1ym = 13._pr/12._pr*(v1ym - 2._pr*v2ym + v3ym)*(v1ym - 2._pr*v2ym + v3ym) + 0.25_pr*(v1ym - 4._pr*v2ym + 3._pr*v3ym)*(&
v1ym - 4._pr*v2ym + 3._pr*v3ym)
          s2ym = 13._pr/12._pr*(v2ym - 2._pr*v3ym + v4ym)*(v2ym - 2._pr*v3ym + v4ym) + 0.25_pr*(v2ym - v4ym)*(v2ym - v4ym)
          s3ym = 13._pr/12._pr*(v3ym - 2._pr*v4ym + v5ym)*(v3ym - 2._pr*v4ym + v5ym) + 0.25_pr*(3._pr*v3ym - 4._pr*v4ym + v5ym)*&
(3._pr*v3ym - 4._pr*v4ym + v5ym)
          a1ym = 0.1_pr/((eepsilon + s1ym)*(eepsilon + s1ym))
          a2ym = 0.6_pr/((eepsilon + s2ym)*(eepsilon + s2ym))
          a3ym = 0.3_pr/((eepsilon + s3ym)*(eepsilon + s3ym))
          w1ym = a1ym/(a1ym + a2ym + a3ym)
          w2ym = a2ym/(a1ym + a2ym + a3ym)
          w3ym = a3ym/(a1ym + a2ym + a3ym)
          drhodym = w1ym*(v1ym/3._pr - 7._pr/6._pr*v2ym + 11._pr/6._pr*v3ym) + w2ym*(-v2ym/6._pr + 5._pr/6._pr*v3ym + v4ym/3._pr) +&
 w3ym*(v3ym/3._pr + 5._pr/6._pr*v4ym - v5ym/6._pr)

          u_(i,j) = - a*drhodxp - bbb*drhodym
        enddo
      enddo
  else if ((a>=0).and.(bbb<0)) then
      do i=1,N
        do j=1,N
          imod1m = modulo(i-1-1,N)+1
          imod2m = modulo(i-1-2,N)+1
          imod3m = modulo(i-1-3,N)+1
          imod1p = modulo(i-1+1,N)+1
          imod2p = modulo(i-1+2,N)+1
          imod3p = modulo(i-1+3,N)+1

          v1xm = (rho_(imod2m,j) - rho_(imod3m,j))/dx
          v2xm = (rho_(imod1m,j) - rho_(imod2m,j))/dx
          v3xm = (rho_(i,j) - rho_(imod1m,j))/dx
          v4xm = (rho_(imod1p,j) - rho_(i,j))/dx
          v5xm = (rho_(imod2p,j) - rho_(imod1p,j))/dx
          s1xm = 13._pr/12._pr*(v1xm - 2._pr*v2xm + v3xm)*(v1xm - 2._pr*v2xm + v3xm) + 0.25_pr*(v1xm - 4._pr*v2xm + 3._pr*v3xm)*(&
v1xm - 4._pr*v2xm + 3._pr*v3xm)
          s2xm = 13._pr/12._pr*(v2xm - 2._pr*v3xm + v4xm)*(v2xm - 2._pr*v3xm + v4xm) + 0.25_pr*(v2xm - v4xm)*(v2xm - v4xm)
          s3xm = 13._pr/12._pr*(v3xm - 2._pr*v4xm + v5xm)*(v3xm - 2._pr*v4xm + v5xm) + 0.25_pr*(3._pr*v3xm - 4._pr*v4xm + v5xm)*&
(3._pr*v3xm - 4._pr*v4xm + v5xm)
          a1xm = 0.1_pr/((eepsilon + s1xm)*(eepsilon + s1xm))
          a2xm = 0.6_pr/((eepsilon + s2xm)*(eepsilon + s2xm))
          a3xm = 0.3_pr/((eepsilon + s3xm)*(eepsilon + s3xm))
          w1xm = a1xm/(a1xm + a2xm + a3xm)
          w2xm = a2xm/(a1xm + a2xm + a3xm)
          w3xm = a3xm/(a1xm + a2xm + a3xm)
          drhodxm = w1xm*(v1xm/3._pr - 7._pr/6._pr*v2xm + 11._pr/6._pr*v3xm) + w2xm*(-v2xm/6._pr + 5._pr/6._pr*v3xm + v4xm/3._pr) +&
 w3xm*(v3xm/3._pr + 5._pr/6._pr*v4xm - v5xm/6._pr)

          jmod1m = modulo(j-1-1,N)+1
          jmod2m = modulo(j-1-2,N)+1
          jmod3m = modulo(j-1-3,N)+1
          jmod1p = modulo(j-1+1,N)+1
          jmod2p = modulo(j-1+2,N)+1
          jmod3p = modulo(j-1+3,N)+1

          v1yp = (rho_(i,jmod3p) - rho_(i,jmod2p))/dy
          v2yp = (rho_(i,jmod2p) - rho_(i,jmod1p))/dy
          v3yp = (rho_(i,jmod1p) - rho_(i,j))/dy
          v4yp = (rho_(i,j) - rho_(i,jmod1m))/dy
          v5yp = (rho_(i,jmod1m) - rho_(i,jmod2m))/dy
          s1yp = 13._pr/12._pr*(v1yp - 2._pr*v2yp + v3yp)*(v1yp - 2._pr*v2yp + v3yp) + 0.25_pr*(v1yp - 4._pr*v2yp + 3._pr*v3yp)*(&
v1yp - 4._pr*v2yp + 3._pr*v3yp)
          s2yp = 13._pr/12._pr*(v2yp - 2._pr*v3yp + v4yp)*(v2yp - 2._pr*v3yp + v4yp) + 0.25_pr*(v2yp - v4yp)*(v2yp - v4yp)
          s3yp = 13._pr/12._pr*(v3yp - 2._pr*v4yp + v5yp)*(v3yp - 2._pr*v4yp + v5yp) + 0.25_pr*(3._pr*v3yp - 4._pr*v4yp + v5yp)*&
(3._pr*v3yp - 4._pr*v4yp + v5yp)
          a1yp = 0.1_pr/((eepsilon + s1yp)*(eepsilon + s1yp))
          a2yp = 0.6_pr/((eepsilon + s2yp)*(eepsilon + s2yp))
          a3yp = 0.3_pr/((eepsilon + s3yp)*(eepsilon + s3yp))
          w1yp = a1yp/(a1yp + a2yp + a3yp)
          w2yp = a2yp/(a1yp + a2yp + a3yp)
          w3yp = a3yp/(a1yp + a2yp + a3yp)
          drhodyp = w1yp*(v1yp/3._pr - 7._pr/6._pr*v2yp + 11._pr/6._pr*v3yp) + w2yp*(-v2yp/6._pr + 5._pr/6._pr*v3yp + v4yp/3._pr) +&
 w3yp*(v3yp/3._pr + 5._pr/6._pr*v4yp - v5yp/6._pr)

          u_(i,j) = - a*drhodxm - bbb*drhodyp
        enddo
      enddo
  else
      do i=1,N
        do j=1,N
          imod1m = modulo(i-1-1,N)+1
          imod2m = modulo(i-1-2,N)+1
          imod3m = modulo(i-1-3,N)+1
          imod1p = modulo(i-1+1,N)+1
          imod2p = modulo(i-1+2,N)+1
          imod3p = modulo(i-1+3,N)+1

          v1xp = (rho_(imod3p,j) - rho_(imod2p,j))/dx
          v2xp = (rho_(imod2p,j) - rho_(imod1p,j))/dx
          v3xp = (rho_(imod1p,j) - rho_(i,j))/dx
          v4xp = (rho_(i,j) - rho_(imod1m,j))/dx
          v5xp = (rho_(imod1m,j) - rho_(imod2m,j))/dx
          s1xp = 13._pr/12._pr*(v1xp - 2._pr*v2xp + v3xp)*(v1xp - 2._pr*v2xp + v3xp) + 0.25_pr*(v1xp - 4._pr*v2xp + 3._pr*v3xp)*(&
v1xp - 4._pr*v2xp + 3._pr*v3xp)
          s2xp = 13._pr/12._pr*(v2xp - 2._pr*v3xp + v4xp)*(v2xp - 2._pr*v3xp + v4xp) + 0.25_pr*(v2xp - v4xp)*(v2xp - v4xp)
          s3xp = 13._pr/12._pr*(v3xp - 2._pr*v4xp + v5xp)*(v3xp - 2._pr*v4xp + v5xp) + 0.25_pr*(3._pr*v3xp - 4._pr*v4xp + v5xp)*&
(3._pr*v3xp - 4._pr*v4xp + v5xp)
          a1xp = 0.1_pr/((eepsilon + s1xp)*(eepsilon + s1xp))
          a2xp = 0.6_pr/((eepsilon + s2xp)*(eepsilon + s2xp))
          a3xp = 0.3_pr/((eepsilon + s3xp)*(eepsilon + s3xp))
          w1xp = a1xp/(a1xp + a2xp + a3xp)
          w2xp = a2xp/(a1xp + a2xp + a3xp)
          w3xp = a3xp/(a1xp + a2xp + a3xp)
          drhodxp = w1xp*(v1xp/3._pr - 7._pr/6._pr*v2xp + 11._pr/6._pr*v3xp) + w2xp*(-v2xp/6._pr + 5._pr/6._pr*v3xp + v4xp/3._pr) +&
 w3xp*(v3xp/3._pr + 5._pr/6._pr*v4xp - v5xp/6._pr)

          jmod1m = modulo(j-1-1,N)+1
          jmod2m = modulo(j-1-2,N)+1
          jmod3m = modulo(j-1-3,N)+1
          jmod1p = modulo(j-1+1,N)+1
          jmod2p = modulo(j-1+2,N)+1
          jmod3p = modulo(j-1+3,N)+1

          v1yp = (rho_(i,jmod3p) - rho_(i,jmod2p))/dy
          v2yp = (rho_(i,jmod2p) - rho_(i,jmod1p))/dy
          v3yp = (rho_(i,jmod1p) - rho_(i,j))/dy
          v4yp = (rho_(i,j) - rho_(i,jmod1m))/dy
          v5yp = (rho_(i,jmod1m) - rho_(i,jmod2m))/dy
          s1yp = 13._pr/12._pr*(v1yp - 2._pr*v2yp + v3yp)*(v1yp - 2._pr*v2yp + v3yp) + 0.25_pr*(v1yp - 4._pr*v2yp + 3._pr*v3yp)*(&
v1yp - 4._pr*v2yp + 3._pr*v3yp)
          s2yp = 13._pr/12._pr*(v2yp - 2._pr*v3yp + v4yp)*(v2yp - 2._pr*v3yp + v4yp) + 0.25_pr*(v2yp - v4yp)*(v2yp - v4yp)
          s3yp = 13._pr/12._pr*(v3yp - 2._pr*v4yp + v5yp)*(v3yp - 2._pr*v4yp + v5yp) + 0.25_pr*(3._pr*v3yp - 4._pr*v4yp + v5yp)*&
(3._pr*v3yp - 4._pr*v4yp + v5yp)
          a1yp = 0.1_pr/((eepsilon + s1yp)*(eepsilon + s1yp))
          a2yp = 0.6_pr/((eepsilon + s2yp)*(eepsilon + s2yp))
          a3yp = 0.3_pr/((eepsilon + s3yp)*(eepsilon + s3yp))
          w1yp = a1yp/(a1yp + a2yp + a3yp)
          w2yp = a2yp/(a1yp + a2yp + a3yp)
          w3yp = a3yp/(a1yp + a2yp + a3yp)
          drhodyp = w1yp*(v1yp/3._pr - 7._pr/6._pr*v2yp + 11._pr/6._pr*v3yp) + w2yp*(-v2yp/6._pr + 5._pr/6._pr*v3yp + v4yp/3._pr) +&
 w3yp*(v3yp/3._pr + 5._pr/6._pr*v4yp - v5yp/6._pr)

          u_(i,j) = - a*drhodxp - bbb*drhodyp
        enddo
      enddo
  endif
end


!! 2D
subroutine computeFunctionWENO52Dperiodic(alpha, rho_, t, u_, u, v)
  real(pr), intent(in) :: alpha, t
  real(pr), dimension(:,:), intent(inout) :: rho_
  real(pr), dimension(:,:), intent(inout) :: u_
  real(pr), dimension(:,:), intent(inout) :: u
  real(pr), dimension(:,:), intent(inout) :: v
  real(pr) :: v1xm, v2xm, v3xm, v4xm, v5xm
  real(pr) :: s1xm, s2xm, s3xm
  real(pr) :: a1xm, a2xm, a3xm
  real(pr) :: w1xm, w2xm, w3xm
  real(pr) :: drhodxm
  real(pr) :: v1xp, v2xp, v3xp, v4xp, v5xp
  real(pr) :: s1xp, s2xp, s3xp
  real(pr) :: a1xp, a2xp, a3xp
  real(pr) :: w1xp, w2xp, w3xp
  real(pr) :: drhodxp
  real(pr) :: v1ym, v2ym, v3ym, v4ym, v5ym
  real(pr) :: s1ym, s2ym, s3ym
  real(pr) :: a1ym, a2ym, a3ym
  real(pr) :: w1ym, w2ym, w3ym
  real(pr) :: drhodym
  real(pr) :: v1yp, v2yp, v3yp, v4yp, v5yp
  real(pr) :: s1yp, s2yp, s3yp
  real(pr) :: a1yp, a2yp, a3yp
  real(pr) :: w1yp, w2yp, w3yp
  real(pr) :: drhodyp
  integer :: i, j, imod1m, imod1p, jmod1m, jmod1p, imod2m, imod2p, jmod2m, jmod2p, imod3m, imod3p, jmod3m, jmod3p

  do i=1,N
    do j=1,N
      imod1m = modulo(i-1-1,N)+1
      imod2m = modulo(i-1-2,N)+1
      imod3m = modulo(i-1-3,N)+1
      imod1p = modulo(i-1+1,N)+1
      imod2p = modulo(i-1+2,N)+1
      imod3p = modulo(i-1+3,N)+1

      jmod1m = modulo(j-1-1,N)+1
      jmod2m = modulo(j-1-2,N)+1
      jmod3m = modulo(j-1-3,N)+1
      jmod1p = modulo(j-1+1,N)+1
      jmod2p = modulo(j-1+2,N)+1
      jmod3p = modulo(j-1+3,N)+1

      v1xm = (rho_(imod2m,j) - rho_(imod3m,j))/dx
      v2xm = (rho_(imod1m,j) - rho_(imod2m,j))/dx
      v3xm = (rho_(i,j) - rho_(imod1m,j))/dx
      v4xm = (rho_(imod1p,j) - rho_(i,j))/dx
      v5xm = (rho_(imod2p,j) - rho_(imod1p,j))/dx
      s1xm = 13._pr/12._pr*(v1xm - 2._pr*v2xm + v3xm)*(v1xm - 2._pr*v2xm + v3xm) + 0.25_pr*(v1xm - 4._pr*v2xm + 3._pr*v3xm)*(v1xm -&
 4._pr*v2xm + 3._pr*v3xm)
      s2xm = 13._pr/12._pr*(v2xm - 2._pr*v3xm + v4xm)*(v2xm - 2._pr*v3xm + v4xm) + 0.25_pr*(v2xm - v4xm)*(v2xm - v4xm)
      s3xm = 13._pr/12._pr*(v3xm - 2._pr*v4xm + v5xm)*(v3xm - 2._pr*v4xm + v5xm) + 0.25_pr*(3._pr*v3xm - 4._pr*v4xm + v5xm)*(3._pr*&
v3xm - 4._pr*v4xm + v5xm)
      a1xm = 0.1_pr/((eepsilon + s1xm)*(eepsilon + s1xm))
      a2xm = 0.6_pr/((eepsilon + s2xm)*(eepsilon + s2xm))
      a3xm = 0.3_pr/((eepsilon + s3xm)*(eepsilon + s3xm))
      w1xm = a1xm/(a1xm + a2xm + a3xm)
      w2xm = a2xm/(a1xm + a2xm + a3xm)
      w3xm = a3xm/(a1xm + a2xm + a3xm)
      drhodxm = w1xm*(v1xm/3._pr - 7._pr/6._pr*v2xm + 11._pr/6._pr*v3xm) + w2xm*(-v2xm/6._pr + 5._pr/6._pr*v3xm + v4xm/3._pr) + &
w3xm*(v3xm/3._pr + 5._pr/6._pr*v4xm - v5xm/6._pr)

      v1ym = (rho_(i,jmod2m) - rho_(i,jmod3m))/dy
      v2ym = (rho_(i,jmod1m) - rho_(i,jmod2m))/dy
      v3ym = (rho_(i,j) - rho_(i,jmod1m))/dy
      v4ym = (rho_(i,jmod1p) - rho_(i,j))/dy
      v5ym = (rho_(i,jmod2p) - rho_(i,jmod1p))/dy
      s1ym = 13._pr/12._pr*(v1ym - 2._pr*v2ym + v3ym)*(v1ym - 2._pr*v2ym + v3ym) + 0.25_pr*(v1ym - 4._pr*v2ym + 3._pr*v3ym)*(v1ym -&
 4._pr*v2ym + 3._pr*v3ym)
      s2ym = 13._pr/12._pr*(v2ym - 2._pr*v3ym + v4ym)*(v2ym - 2._pr*v3ym + v4ym) + 0.25_pr*(v2ym - v4ym)*(v2ym - v4ym)
      s3ym = 13._pr/12._pr*(v3ym - 2._pr*v4ym + v5ym)*(v3ym - 2._pr*v4ym + v5ym) + 0.25_pr*(3._pr*v3ym - 4._pr*v4ym + v5ym)*(3._pr*&
v3ym - 4._pr*v4ym + v5ym)
      a1ym = 0.1_pr/((eepsilon + s1ym)*(eepsilon + s1ym))
      a2ym = 0.6_pr/((eepsilon + s2ym)*(eepsilon + s2ym))
      a3ym = 0.3_pr/((eepsilon + s3ym)*(eepsilon + s3ym))
      w1ym = a1ym/(a1ym + a2ym + a3ym)
      w2ym = a2ym/(a1ym + a2ym + a3ym)
      w3ym = a3ym/(a1ym + a2ym + a3ym)
      drhodym = w1ym*(v1ym/3._pr - 7._pr/6._pr*v2ym + 11._pr/6._pr*v3ym) + w2ym*(-v2ym/6._pr + 5._pr/6._pr*v3ym + v4ym/3._pr) + &
w3ym*(v3ym/3._pr + 5._pr/6._pr*v4ym - v5ym/6._pr)

      v1xp = (rho_(imod3p,j) - rho_(imod2p,j))/dx
      v2xp = (rho_(imod2p,j) - rho_(imod1p,j))/dx
      v3xp = (rho_(imod1p,j) - rho_(i,j))/dx
      v4xp = (rho_(i,j) - rho_(imod1m,j))/dx
      v5xp = (rho_(imod1m,j) - rho_(imod2m,j))/dx
      s1xp = 13._pr/12._pr*(v1xp - 2._pr*v2xp + v3xp)*(v1xp - 2._pr*v2xp + v3xp) + 0.25_pr*(v1xp - 4._pr*v2xp + 3._pr*v3xp)*(v1xp -&
 4._pr*v2xp + 3._pr*v3xp)
      s2xp = 13._pr/12._pr*(v2xp - 2._pr*v3xp + v4xp)*(v2xp - 2._pr*v3xp + v4xp) + 0.25_pr*(v2xp - v4xp)*(v2xp - v4xp)
      s3xp = 13._pr/12._pr*(v3xp - 2._pr*v4xp + v5xp)*(v3xp - 2._pr*v4xp + v5xp) + 0.25_pr*(3._pr*v3xp - 4._pr*v4xp + v5xp)*(3._pr*&
v3xp - 4._pr*v4xp + v5xp)
      a1xp = 0.1_pr/((eepsilon + s1xp)*(eepsilon + s1xp))
      a2xp = 0.6_pr/((eepsilon + s2xp)*(eepsilon + s2xp))
      a3xp = 0.3_pr/((eepsilon + s3xp)*(eepsilon + s3xp))
      w1xp = a1xp/(a1xp + a2xp + a3xp)
      w2xp = a2xp/(a1xp + a2xp + a3xp)
      w3xp = a3xp/(a1xp + a2xp + a3xp)
      drhodxp = w1xp*(v1xp/3._pr - 7._pr/6._pr*v2xp + 11._pr/6._pr*v3xp) + w2xp*(-v2xp/6._pr + 5._pr/6._pr*v3xp + v4xp/3._pr) + &
w3xp*(v3xp/3._pr + 5._pr/6._pr*v4xp - v5xp/6._pr)

      v1yp = (rho_(i,jmod3p) - rho_(i,jmod2p))/dy
      v2yp = (rho_(i,jmod2p) - rho_(i,jmod1p))/dy
      v3yp = (rho_(i,jmod1p) - rho_(i,j))/dy
      v4yp = (rho_(i,j) - rho_(i,jmod1m))/dy
      v5yp = (rho_(i,jmod1m) - rho_(i,jmod2m))/dy
      s1yp = 13._pr/12._pr*(v1yp - 2._pr*v2yp + v3yp)*(v1yp - 2._pr*v2yp + v3yp) + 0.25_pr*(v1yp - 4._pr*v2yp + 3._pr*v3yp)*(v1yp -&
 4._pr*v2yp + 3._pr*v3yp)
      s2yp = 13._pr/12._pr*(v2yp - 2._pr*v3yp + v4yp)*(v2yp - 2._pr*v3yp + v4yp) + 0.25_pr*(v2yp - v4yp)*(v2yp - v4yp)
      s3yp = 13._pr/12._pr*(v3yp - 2._pr*v4yp + v5yp)*(v3yp - 2._pr*v4yp + v5yp) + 0.25_pr*(3._pr*v3yp - 4._pr*v4yp + v5yp)*(3._pr*&
v3yp - 4._pr*v4yp + v5yp)
      a1yp = 0.1_pr/((eepsilon + s1yp)*(eepsilon + s1yp))
      a2yp = 0.6_pr/((eepsilon + s2yp)*(eepsilon + s2yp))
      a3yp = 0.3_pr/((eepsilon + s3yp)*(eepsilon + s3yp))
      w1yp = a1yp/(a1yp + a2yp + a3yp)
      w2yp = a2yp/(a1yp + a2yp + a3yp)
      w3yp = a3yp/(a1yp + a2yp + a3yp)
      drhodyp = w1yp*(v1yp/3._pr - 7._pr/6._pr*v2yp + 11._pr/6._pr*v3yp) + w2yp*(-v2yp/6._pr + 5._pr/6._pr*v3yp + v4yp/3._pr) + &
w3yp*(v3yp/3._pr + 5._pr/6._pr*v4yp - v5yp/6._pr)

      if ((u(i,j)>=0).and.(v(i,j)>=0)) then
          u_(i,j) = - u(i,j)*drhodxm - v(i,j)*drhodym
      else if ((u(i,j)<0).and.(v(i,j)>=0)) then
          u_(i,j) = - u(i,j)*drhodxp - v(i,j)*drhodym
      else if ((u(i,j)>=0).and.(v(i,j)<0)) then
          u_(i,j) = - u(i,j)*drhodxm - v(i,j)*drhodyp
      else
          u_(i,j) = - u(i,j)*drhodxp - v(i,j)*drhodyp
      endif

  enddo
enddo
end
!! 2D
subroutine computeFunctionWENO52DdirichletH(alpha, rho_, t, u_, u, v)
  real(pr), intent(in) :: alpha, t
  real(pr), dimension(:,:), intent(inout) :: rho_
  real(pr), dimension(:,:), intent(inout) :: u_
  real(pr), dimension(:,:), intent(inout) :: u
  real(pr), dimension(:,:), intent(inout) :: v
  real(pr) :: v1xm, v2xm, v3xm, v4xm, v5xm
  real(pr) :: s1xm, s2xm, s3xm
  real(pr) :: a1xm, a2xm, a3xm
  real(pr) :: w1xm, w2xm, w3xm
  real(pr) :: drhodxm
  real(pr) :: v1xp, v2xp, v3xp, v4xp, v5xp
  real(pr) :: s1xp, s2xp, s3xp
  real(pr) :: a1xp, a2xp, a3xp
  real(pr) :: w1xp, w2xp, w3xp
  real(pr) :: drhodxp
  real(pr) :: v1ym, v2ym, v3ym, v4ym, v5ym
  real(pr) :: s1ym, s2ym, s3ym
  real(pr) :: a1ym, a2ym, a3ym
  real(pr) :: w1ym, w2ym, w3ym
  real(pr) :: drhodym
  real(pr) :: v1yp, v2yp, v3yp, v4yp, v5yp
  real(pr) :: s1yp, s2yp, s3yp
  real(pr) :: a1yp, a2yp, a3yp
  real(pr) :: w1yp, w2yp, w3yp
  real(pr) :: drhodyp
  integer :: i, j, imod1m, imod1p, jmod1m, jmod1p, imod2m, imod2p, jmod2m, jmod2p, imod3m, imod3p, jmod3m, jmod3p

  do i=1,N
    do j=1,N
      imod1m = modulo(i-1-1,N)+1
      imod2m = modulo(i-1-2,N)+1
      imod3m = modulo(i-1-3,N)+1
      imod1p = modulo(i-1+1,N)+1
      imod2p = modulo(i-1+2,N)+1
      imod3p = modulo(i-1+3,N)+1

      jmod1m = modulo(j-1-1,N)+1
      jmod2m = modulo(j-1-2,N)+1
      jmod3m = modulo(j-1-3,N)+1
      jmod1p = modulo(j-1+1,N)+1
      jmod2p = modulo(j-1+2,N)+1
      jmod3p = modulo(j-1+3,N)+1

      if (i == 1) imod1m = 0._pr
      if (i == 2) then
              imod1m = 0._pr
              imod2m = 0._pr
      endif
      if (i == 3) then
              imod2m = 0._pr
              imod3m = 0._pr
      endif
      if (i == 4) imod3m = 0._pr
      if (i == N) imod1p = 0._pr
      if (i == N-1) then
              imod1p = 0._pr
              imod2p = 0._pr
      endif
      if (i == N-2) then
              imod2p = 0._pr
              imod3p = 0._pr
      endif
      if (i == N-3) imod3p = 0._pr
      if (j == 1) jmod1m = 0._pr
      if (j == 2) then
              jmod1m = 0._pr
              jmod2m = 0._pr
      endif
      if (j == 3) then
              jmod2m = 0._pr
              jmod3m = 0._pr
      endif
      if (j == 4) jmod3m = 0._pr
      if (j == N) jmod1p = 0._pr
      if (j == N-1) then
              jmod1p = 0._pr
              jmod2p = 0._pr
      endif
      if (j == N-2) then
              jmod2p = 0._pr
              jmod3p = 0._pr
      endif
      if (j == N-3) jmod3p = 0._pr

      v1xm = (rho_(imod2m,j) - rho_(imod3m,j))/dx
      v2xm = (rho_(imod1m,j) - rho_(imod2m,j))/dx
      v3xm = (rho_(i,j) - rho_(imod1m,j))/dx
      v4xm = (rho_(imod1p,j) - rho_(i,j))/dx
      v5xm = (rho_(imod2p,j) - rho_(imod1p,j))/dx
      s1xm = 13._pr/12._pr*(v1xm - 2._pr*v2xm + v3xm)*(v1xm - 2._pr*v2xm + v3xm) + 0.25_pr*(v1xm - 4._pr*v2xm + 3._pr*v3xm)*(v1xm -&
 4._pr*v2xm + 3._pr*v3xm)
      s2xm = 13._pr/12._pr*(v2xm - 2._pr*v3xm + v4xm)*(v2xm - 2._pr*v3xm + v4xm) + 0.25_pr*(v2xm - v4xm)*(v2xm - v4xm)
      s3xm = 13._pr/12._pr*(v3xm - 2._pr*v4xm + v5xm)*(v3xm - 2._pr*v4xm + v5xm) + 0.25_pr*(3._pr*v3xm - 4._pr*v4xm + v5xm)*(3._pr*&
v3xm - 4._pr*v4xm + v5xm)
      a1xm = 0.1_pr/((eepsilon + s1xm)*(eepsilon + s1xm))
      a2xm = 0.6_pr/((eepsilon + s2xm)*(eepsilon + s2xm))
      a3xm = 0.3_pr/((eepsilon + s3xm)*(eepsilon + s3xm))
      w1xm = a1xm/(a1xm + a2xm + a3xm)
      w2xm = a2xm/(a1xm + a2xm + a3xm)
      w3xm = a3xm/(a1xm + a2xm + a3xm)
      drhodxm = w1xm*(v1xm/3._pr - 7._pr/6._pr*v2xm + 11._pr/6._pr*v3xm) + w2xm*(-v2xm/6._pr + 5._pr/6._pr*v3xm + v4xm/3._pr) + &
w3xm*(v3xm/3._pr + 5._pr/6._pr*v4xm - v5xm/6._pr)

      v1ym = (rho_(i,jmod2m) - rho_(i,jmod3m))/dy
      v2ym = (rho_(i,jmod1m) - rho_(i,jmod2m))/dy
      v3ym = (rho_(i,j) - rho_(i,jmod1m))/dy
      v4ym = (rho_(i,jmod1p) - rho_(i,j))/dy
      v5ym = (rho_(i,jmod2p) - rho_(i,jmod1p))/dy
      s1ym = 13._pr/12._pr*(v1ym - 2._pr*v2ym + v3ym)*(v1ym - 2._pr*v2ym + v3ym) + 0.25_pr*(v1ym - 4._pr*v2ym + 3._pr*v3ym)*(v1ym -&
 4._pr*v2ym + 3._pr*v3ym)
      s2ym = 13._pr/12._pr*(v2ym - 2._pr*v3ym + v4ym)*(v2ym - 2._pr*v3ym + v4ym) + 0.25_pr*(v2ym - v4ym)*(v2ym - v4ym)
      s3ym = 13._pr/12._pr*(v3ym - 2._pr*v4ym + v5ym)*(v3ym - 2._pr*v4ym + v5ym) + 0.25_pr*(3._pr*v3ym - 4._pr*v4ym + v5ym)*(3._pr*&
v3ym - 4._pr*v4ym + v5ym)
      a1ym = 0.1_pr/((eepsilon + s1ym)*(eepsilon + s1ym))
      a2ym = 0.6_pr/((eepsilon + s2ym)*(eepsilon + s2ym))
      a3ym = 0.3_pr/((eepsilon + s3ym)*(eepsilon + s3ym))
      w1ym = a1ym/(a1ym + a2ym + a3ym)
      w2ym = a2ym/(a1ym + a2ym + a3ym)
      w3ym = a3ym/(a1ym + a2ym + a3ym)
      drhodym = w1ym*(v1ym/3._pr - 7._pr/6._pr*v2ym + 11._pr/6._pr*v3ym) + w2ym*(-v2ym/6._pr + 5._pr/6._pr*v3ym + v4ym/3._pr) + &
w3ym*(v3ym/3._pr + 5._pr/6._pr*v4ym - v5ym/6._pr)

      v1xp = (rho_(imod3p,j) - rho_(imod2p,j))/dx
      v2xp = (rho_(imod2p,j) - rho_(imod1p,j))/dx
      v3xp = (rho_(imod1p,j) - rho_(i,j))/dx
      v4xp = (rho_(i,j) - rho_(imod1m,j))/dx
      v5xp = (rho_(imod1m,j) - rho_(imod2m,j))/dx
      s1xp = 13._pr/12._pr*(v1xp - 2._pr*v2xp + v3xp)*(v1xp - 2._pr*v2xp + v3xp) + 0.25_pr*(v1xp - 4._pr*v2xp + 3._pr*v3xp)*(v1xp -&
 4._pr*v2xp + 3._pr*v3xp)
      s2xp = 13._pr/12._pr*(v2xp - 2._pr*v3xp + v4xp)*(v2xp - 2._pr*v3xp + v4xp) + 0.25_pr*(v2xp - v4xp)*(v2xp - v4xp)
      s3xp = 13._pr/12._pr*(v3xp - 2._pr*v4xp + v5xp)*(v3xp - 2._pr*v4xp + v5xp) + 0.25_pr*(3._pr*v3xp - 4._pr*v4xp + v5xp)*(3._pr*&
v3xp - 4._pr*v4xp + v5xp)
      a1xp = 0.1_pr/((eepsilon + s1xp)*(eepsilon + s1xp))
      a2xp = 0.6_pr/((eepsilon + s2xp)*(eepsilon + s2xp))
      a3xp = 0.3_pr/((eepsilon + s3xp)*(eepsilon + s3xp))
      w1xp = a1xp/(a1xp + a2xp + a3xp)
      w2xp = a2xp/(a1xp + a2xp + a3xp)
      w3xp = a3xp/(a1xp + a2xp + a3xp)
      drhodxp = w1xp*(v1xp/3._pr - 7._pr/6._pr*v2xp + 11._pr/6._pr*v3xp) + w2xp*(-v2xp/6._pr + 5._pr/6._pr*v3xp + v4xp/3._pr) + &
w3xp*(v3xp/3._pr + 5._pr/6._pr*v4xp - v5xp/6._pr)

      v1yp = (rho_(i,jmod3p) - rho_(i,jmod2p))/dy
      v2yp = (rho_(i,jmod2p) - rho_(i,jmod1p))/dy
      v3yp = (rho_(i,jmod1p) - rho_(i,j))/dy
      v4yp = (rho_(i,j) - rho_(i,jmod1m))/dy
      v5yp = (rho_(i,jmod1m) - rho_(i,jmod2m))/dy
      s1yp = 13._pr/12._pr*(v1yp - 2._pr*v2yp + v3yp)*(v1yp - 2._pr*v2yp + v3yp) + 0.25_pr*(v1yp - 4._pr*v2yp + 3._pr*v3yp)*(v1yp -&
 4._pr*v2yp + 3._pr*v3yp)
      s2yp = 13._pr/12._pr*(v2yp - 2._pr*v3yp + v4yp)*(v2yp - 2._pr*v3yp + v4yp) + 0.25_pr*(v2yp - v4yp)*(v2yp - v4yp)
      s3yp = 13._pr/12._pr*(v3yp - 2._pr*v4yp + v5yp)*(v3yp - 2._pr*v4yp + v5yp) + 0.25_pr*(3._pr*v3yp - 4._pr*v4yp + v5yp)*(3._pr*&
v3yp - 4._pr*v4yp + v5yp)
      a1yp = 0.1_pr/((eepsilon + s1yp)*(eepsilon + s1yp))
      a2yp = 0.6_pr/((eepsilon + s2yp)*(eepsilon + s2yp))
      a3yp = 0.3_pr/((eepsilon + s3yp)*(eepsilon + s3yp))
      w1yp = a1yp/(a1yp + a2yp + a3yp)
      w2yp = a2yp/(a1yp + a2yp + a3yp)
      w3yp = a3yp/(a1yp + a2yp + a3yp)
      drhodyp = w1yp*(v1yp/3._pr - 7._pr/6._pr*v2yp + 11._pr/6._pr*v3yp) + w2yp*(-v2yp/6._pr + 5._pr/6._pr*v3yp + v4yp/3._pr) + &
w3yp*(v3yp/3._pr + 5._pr/6._pr*v4yp - v5yp/6._pr)

      if ((u(i,j)>0).and.(v(i,j)>0)) then
          u_(i,j) = - u(i,j)*drhodxm - v(i,j)*drhodym
      else if ((u(i,j)<0).and.(v(i,j)>0)) then
          u_(i,j) = - u(i,j)*drhodxp - v(i,j)*drhodym
      else if ((u(i,j)>0).and.(v(i,j)<0)) then
          u_(i,j) = - u(i,j)*drhodxm - v(i,j)*drhodyp
      else
          u_(i,j) = - u(i,j)*drhodxp - v(i,j)*drhodyp
      endif

  enddo
enddo

end

subroutine updateDistance(rho,gradPhi)
  real(pr), dimension(:,:), intent(inout) :: rho
  real(pr),dimension(size(rho,1),size(rho,2)),intent(inout) :: gradPhi
  real(pr),dimension(size(gradPhi,1),size(gradPhi,2)) :: gradxPhi, gradyPhi
  real(pr) :: nr, nrtmp1, nrtmp2, nrprec, epsDist, dtDist, cflDist, h,a,b,c
  integer :: nIterMaxDist
  real(pr), dimension(size(rho,1),size(rho,2)) :: rhoDist, rho0Dist
  integer :: iter
  real(pr) :: xpp, xmp, ypp, ymp, xpm, xmm, ypm, ymm
  real(pr) :: v1xm, v2xm, v3xm, v4xm, v5xm
  real(pr) :: s1xm, s2xm, s3xm
  real(pr) :: a1xm, a2xm, a3xm
  real(pr) :: w1xm, w2xm, w3xm
  real(pr) :: drhodxm
  real(pr) :: v1xp, v2xp, v3xp, v4xp, v5xp
  real(pr) :: s1xp, s2xp, s3xp
  real(pr) :: a1xp, a2xp, a3xp
  real(pr) :: w1xp, w2xp, w3xp
  real(pr) :: drhodxp
  real(pr) :: drhodx
  real(pr) :: v1ym, v2ym, v3ym, v4ym, v5ym
  real(pr) :: s1ym, s2ym, s3ym
  real(pr) :: a1ym, a2ym, a3ym
  real(pr) :: w1ym, w2ym, w3ym
  real(pr) :: drhodym
  real(pr) :: v1yp, v2yp, v3yp, v4yp, v5yp
  real(pr) :: s1yp, s2yp, s3yp
  real(pr) :: a1yp, a2yp, a3yp
  real(pr) :: w1yp, w2yp, w3yp
  real(pr) :: drhodyp
  real(pr) :: drhody
  integer :: i, j
  real(pr) :: rhoi1mj, rhoi1pj, rhoij1m, rhoij1p, rhoi2mj, rhoi2pj, rhoij2m, rhoij2p, rhoi3mj, rhoi3pj,&
  rhoij3m, rhoij3p, rhoij

  cflDist = 0.5_pr
  epsDist = 0.01_pr
  nIterMaxDist = 4000
  dtDist = cflDist/(1._pr/dx + 1._pr/dy)
  rho0Dist = rho
  h = max(dx,dy)

  call Norm22D(rhoDist,nrtmp1)
  call Norm22D(rho,nrtmp2)

  call distanceIteration(rho0Dist, dtDist, rho, rhoDist)
  call Norm22D(rhoDist,nrtmp1)
  call Norm22D(rho,nrtmp2)

  rho = rho + rhoDist
  call Norm22D(rhoDist,nrtmp1)
  call Norm22D(rho,nrtmp2)
  nr = nrtmp1 / nrtmp2
  nrprec = 0._pr

  iter = 1

  do while ((nr > epsDist).and.(iter < nIterMaxDist))
    call distanceIteration(rho0Dist, dtDist, rho, rhoDist)
    nrprec = nr
    rho = rho + rhoDist
    call Norm22D(rhoDist,nrtmp1)
    call Norm22D(rho,nrtmp2)
    nr = nrtmp1 / nrtmp2
    iter = iter + 1
  enddo

  do j=1,size(rho,2)
    do i=1,size(rho,1)
              call neumannBC2D(rho,i-1,j,rhoi1mj)
              call neumannBC2D(rho,i+1,j,rhoi1pj)
              call neumannBC2D(rho,i,j-1,rhoij1m)
              call neumannBC2D(rho,i,j+1,rhoij1p)
              drhodx = (rhoi1pj-rhoi1mj)/(2*dx)
              drhody = (rhoij1p-rhoij1m)/(2*dy)
              a = sqrt(drhodx**2 + drhody**2)

              drhodx = (rhoi1pj-rho(i,j))/dx
              drhody = (rhoij1p-rho(i,j))/dy
              b = sqrt(drhodx**2 + drhody**2)

              drhodx = (rho(i,j)-rhoi1mj)/dx
              drhody = (rho(i,j)-rhoij1m)/dy
              c = sqrt(drhodx**2 + drhody**2)

              gradPhi(i,j) = min(min(a,b),c)
  enddo
  enddo

  if (iter >= nIterMaxDist) then
          write(*,*) "Distanciation did not converge, res norm = ", nr
  else
          write(*,*) "Distanciation nIter = ", iter, " nr = ", nr&
          ," max ",maxval(rho(:,:))," min ",minval(rho(:,:))&
          ," maxG minG ",maxval(gradPhi)," ",minval(gradPhi)
  endif
end
subroutine updateDistanceINI(rho,gradPhi)
  real(pr), dimension(:,:), intent(inout) :: rho
  real(pr),dimension(size(rho,1),size(rho,2)),intent(inout) :: gradPhi
  real(pr),dimension(size(gradPhi,1),size(gradPhi,2)) :: gradxPhi, gradyPhi
  real(pr) :: nr, nrtmp1, nrtmp2, nrprec, epsDist, dtDist, cflDist, h,a,b,c
  integer :: nIterMaxDist
  real(pr), dimension(size(rho,1),size(rho,2)) :: rhoDist, rho0Dist
  integer :: iter
  real(pr) :: xpp, xmp, ypp, ymp, xpm, xmm, ypm, ymm
  real(pr) :: v1xm, v2xm, v3xm, v4xm, v5xm
  real(pr) :: s1xm, s2xm, s3xm
  real(pr) :: a1xm, a2xm, a3xm
  real(pr) :: w1xm, w2xm, w3xm
  real(pr) :: drhodxm
  real(pr) :: v1xp, v2xp, v3xp, v4xp, v5xp
  real(pr) :: s1xp, s2xp, s3xp
  real(pr) :: a1xp, a2xp, a3xp
  real(pr) :: w1xp, w2xp, w3xp
  real(pr) :: drhodxp
  real(pr) :: drhodx
  real(pr) :: v1ym, v2ym, v3ym, v4ym, v5ym
  real(pr) :: s1ym, s2ym, s3ym
  real(pr) :: a1ym, a2ym, a3ym
  real(pr) :: w1ym, w2ym, w3ym
  real(pr) :: drhodym
  real(pr) :: v1yp, v2yp, v3yp, v4yp, v5yp
  real(pr) :: s1yp, s2yp, s3yp
  real(pr) :: a1yp, a2yp, a3yp
  real(pr) :: w1yp, w2yp, w3yp
  real(pr) :: drhodyp
  real(pr) :: drhody
  integer :: i, j
  real(pr) :: rhoi1mj, rhoi1pj, rhoij1m, rhoij1p, rhoi2mj, rhoi2pj, rhoij2m, rhoij2p, rhoi3mj, rhoi3pj,&
  rhoij3m, rhoij3p, rhoij

  cflDist = 0.5_pr
  epsDist = 0.001_pr
  nIterMaxDist = 4000
  dtDist = cflDist/(1._pr/dx + 1._pr/dy)
  rho0Dist = rho
  h = max(dx,dy)

  call Norm22D(rhoDist,nrtmp1)
  call Norm22D(rho,nrtmp2)

  call distanceIteration(rho0Dist, dtDist, rho, rhoDist)
  call Norm22D(rhoDist,nrtmp1)
  call Norm22D(rho,nrtmp2)

  rho = rho + rhoDist
  call Norm22D(rhoDist,nrtmp1)
  call Norm22D(rho,nrtmp2)
  nr = nrtmp1 / nrtmp2
  nrprec = 0._pr

  iter = 1

  do while ((nr > epsDist).and.(iter < nIterMaxDist))
    call distanceIteration(rho0Dist, dtDist, rho, rhoDist)
    nrprec = nr
    rho = rho + rhoDist
    call Norm22D(rhoDist,nrtmp1)
    call Norm22D(rho,nrtmp2)
    nr = nrtmp1 / nrtmp2
    iter = iter + 1
  enddo

  do j=1,size(rho,2)
    do i=1,size(rho,1)
              call neumannBC2D(rho,i-1,j,rhoi1mj)
              call neumannBC2D(rho,i+1,j,rhoi1pj)
              call neumannBC2D(rho,i,j-1,rhoij1m)
              call neumannBC2D(rho,i,j+1,rhoij1p)
              drhodx = (rhoi1pj-rhoi1mj)/(2*dx)
              drhody = (rhoij1p-rhoij1m)/(2*dy)
              a = sqrt(drhodx**2 + drhody**2)

              drhodx = (rhoi1pj-rho(i,j))/dx
              drhody = (rhoij1p-rho(i,j))/dy
              b = sqrt(drhodx**2 + drhody**2)

              drhodx = (rho(i,j)-rhoi1mj)/dx
              drhody = (rho(i,j)-rhoij1m)/dy
              c = sqrt(drhodx**2 + drhody**2)

              gradPhi(i,j) = a
  enddo
  enddo

  if (iter >= nIterMaxDist) then
          write(*,*) "Distanciation did not converge, res norm = ", nr
  else
          write(*,*) "Distanciation nIter = ", iter, " nr = ", nr&
          ," max ",maxval(rho(:,:))," min ",minval(rho(:,:))&
          ," maxG minG ",maxval(gradPhi)," ",minval(gradPhi)
  endif
end

subroutine distanceIteration(rho0Dist, dt_, rho_, u_)
  real(pr), intent(in) :: dt_
  real(pr), dimension(:,:), intent(in) :: rho0Dist
  real(pr), dimension(:,:), intent(inout) :: rho_
  real(pr), dimension(:,:), intent(inout) :: u_
  real(pr), dimension(size(u_,1),size(u_,2)) :: gradxPhi, gradyPhi
  real(pr) :: h, x, y, xpp, xmp, ypp, ymp, xpm, xmm, ypm, ymm, a, b, c
  real(pr) :: S, G, dUi0, U0_ij, Un_ij, U0_im1j, U0_ip1j, U0_ijm1, U0_ijp1
  real(pr) :: v1xm, v2xm, v3xm, v4xm, v5xm
  real(pr) :: s1xm, s2xm, s3xm
  real(pr) :: a1xm, a2xm, a3xm
  real(pr) :: w1xm, w2xm, w3xm
  real(pr) :: drhodxm
  real(pr) :: v1xp, v2xp, v3xp, v4xp, v5xp
  real(pr) :: s1xp, s2xp, s3xp
  real(pr) :: a1xp, a2xp, a3xp
  real(pr) :: w1xp, w2xp, w3xp
  real(pr) :: drhodxp
  real(pr) :: v1ym, v2ym, v3ym, v4ym, v5ym
  real(pr) :: s1ym, s2ym, s3ym
  real(pr) :: a1ym, a2ym, a3ym
  real(pr) :: w1ym, w2ym, w3ym
  real(pr) :: drhodym
  real(pr) :: v1yp, v2yp, v3yp, v4yp, v5yp
  real(pr) :: s1yp, s2yp, s3yp
  real(pr) :: a1yp, a2yp, a3yp
  real(pr) :: w1yp, w2yp, w3yp
  real(pr) :: drhodyp
  integer :: i, j
  real(pr) :: rho_i1mj, rho_i1pj, rho_ij1m, rho_ij1p, rho_i2mj, rho_i2pj, rho_ij2m, rho_ij2p, rho_i3mj, rho_i3pj,&
  rho_ij3m, rho_ij3p, rho_ij

  h = sqrt(2._pr)*dx*dy/(sqrt(dx*dx + dy*dy))

  do j=1,size(rho0Dist,2)
    do i=1,size(rho0Dist,1)
      U0_ij = rho0Dist(i,j)

      call neumannBC2D(rho0Dist,i-1,j,U0_im1j)
      call neumannBC2D(rho0Dist,i+1,j,U0_ip1j)
      call neumannBC2D(rho0Dist,i,j-1,U0_ijm1)
      call neumannBC2D(rho0Dist,i,j+1,U0_ijp1)

      Un_ij = rho_(i,j)
      S = sign(1._pr,U0_ij)

      if ((U0_ij*U0_im1j<0).or.(U0_ij*U0_ip1j<0).or.(U0_ij*U0_ijm1<0).or.(U0_ij*U0_ijp1<0)) then
              x = U0_ip1j - U0_im1j
              y = U0_ijp1 - U0_ijm1
              a = 0.5_pr*sqrt(x*x + y*y)

              x = U0_ip1j - U0_ij
              y = U0_ijp1 - U0_ij
              b = sqrt(x*x + y*y)

              x = U0_ij - U0_im1j
              y = U0_ij - U0_ijm1
              c = sqrt(x*x + y*y)

              dUi0 = max(max(a,b),max(c,eepsilon))
              G = abs(Un_ij)

              u_(i,j) = -dt_*(S*G/h - U0_ij/dUi0)
      else

                call neumannBC2D(rho_,i-3,j,rho_i3mj)
                call neumannBC2D(rho_,i-2,j,rho_i2mj)
                call neumannBC2D(rho_,i-1,j,rho_i1mj)
                call neumannBC2D(rho_,i+1,j,rho_i1pj)
                call neumannBC2D(rho_,i+2,j,rho_i2pj)
                call neumannBC2D(rho_,i+3,j,rho_i3pj)
                call neumannBC2D(rho_,i,j-3,rho_ij3m)
                call neumannBC2D(rho_,i,j-2,rho_ij2m)
                call neumannBC2D(rho_,i,j-1,rho_ij1m)
                call neumannBC2D(rho_,i,j+1,rho_ij1p)
                call neumannBC2D(rho_,i,j+2,rho_ij2p)
                call neumannBC2D(rho_,i,j+3,rho_ij3p)
                rho_ij = rho_(i,j)

              v1xm = (rho_i2mj - rho_i3mj)/dx
              v2xm = (rho_i1mj - rho_i2mj)/dx
              v3xm = (rho_ij - rho_i1mj)/dx
              v4xm = (rho_i1pj - rho_ij)/dx
              v5xm = (rho_i2pj - rho_i1pj)/dx

              s1xm = 13._pr/12._pr*(v1xm - 2._pr*v2xm + v3xm)*(v1xm - 2._pr*v2xm + v3xm) + 0.25_pr*(v1xm - 4._pr*v2xm + 3._pr*v3xm)&
*(v1xm - 4._pr*v2xm + 3._pr*v3xm)
              s2xm = 13._pr/12._pr*(v2xm - 2._pr*v3xm + v4xm)*(v2xm - 2._pr*v3xm + v4xm) + 0.25_pr*(v2xm - v4xm)*(v2xm - v4xm)
              s3xm = 13._pr/12._pr*(v3xm - 2._pr*v4xm + v5xm)*(v3xm - 2._pr*v4xm + v5xm) + 0.25_pr*(3._pr*v3xm - 4._pr*v4xm + v5xm)&
*(3._pr*v3xm - 4._pr*v4xm + v5xm)
              a1xm = 0.1_pr/((eepsilon + s1xm)*(eepsilon + s1xm))
              a2xm = 0.6_pr/((eepsilon + s2xm)*(eepsilon + s2xm))
              a3xm = 0.3_pr/((eepsilon + s3xm)*(eepsilon + s3xm))
              w1xm = a1xm/(a1xm + a2xm + a3xm)
              w2xm = a2xm/(a1xm + a2xm + a3xm)
              w3xm = a3xm/(a1xm + a2xm + a3xm)
              drhodxm = w1xm*(v1xm/3._pr - 7._pr/6._pr*v2xm + 11._pr/6._pr*v3xm) + w2xm*(-v2xm/6._pr + 5._pr/6._pr*v3xm + v4xm/&
3._pr) + w3xm*(v3xm/3._pr + 5._pr/6._pr*v4xm - v5xm/6._pr)

              v1ym = (rho_ij2m - rho_ij3m)/dy
              v2ym = (rho_ij1m - rho_ij2m)/dy
              v3ym = (rho_ij - rho_ij1m)/dy
              v4ym = (rho_ij1p - rho_ij)/dy
              v5ym = (rho_ij2p - rho_ij1p)/dy

              s1ym = 13._pr/12._pr*(v1ym - 2._pr*v2ym + v3ym)*(v1ym - 2._pr*v2ym + v3ym) + 0.25_pr*(v1ym - 4._pr*v2ym + 3._pr*v3ym)&
*(v1ym - 4._pr*v2ym + 3._pr*v3ym)
              s2ym = 13._pr/12._pr*(v2ym - 2._pr*v3ym + v4ym)*(v2ym - 2._pr*v3ym + v4ym) + 0.25_pr*(v2ym - v4ym)*(v2ym - v4ym)
              s3ym = 13._pr/12._pr*(v3ym - 2._pr*v4ym + v5ym)*(v3ym - 2._pr*v4ym + v5ym) + 0.25_pr*(3._pr*v3ym - 4._pr*v4ym + v5ym)&
*(3._pr*v3ym - 4._pr*v4ym + v5ym)
              a1ym = 0.1_pr/((eepsilon + s1ym)*(eepsilon + s1ym))
              a2ym = 0.6_pr/((eepsilon + s2ym)*(eepsilon + s2ym))
              a3ym = 0.3_pr/((eepsilon + s3ym)*(eepsilon + s3ym))
              w1ym = a1ym/(a1ym + a2ym + a3ym)
              w2ym = a2ym/(a1ym + a2ym + a3ym)
              w3ym = a3ym/(a1ym + a2ym + a3ym)
              drhodym = w1ym*(v1ym/3._pr - 7._pr/6._pr*v2ym + 11._pr/6._pr*v3ym) + w2ym*(-v2ym/6._pr + 5._pr/6._pr*v3ym + v4ym/&
3._pr) + w3ym*(v3ym/3._pr + 5._pr/6._pr*v4ym - v5ym/6._pr)

              v1xp = (rho_i3pj - rho_i2pj)/dx
              v2xp = (rho_i2pj - rho_i1pj)/dx
              v3xp = (rho_i1pj - rho_ij)/dx
              v4xp = (rho_ij - rho_i1mj)/dx
              v5xp = (rho_i1mj - rho_i2mj)/dx

              s1xp = 13._pr/12._pr*(v1xp - 2._pr*v2xp + v3xp)*(v1xp - 2._pr*v2xp + v3xp) + 0.25_pr*(v1xp - 4._pr*v2xp + 3._pr*v3xp)&
*(v1xp - 4._pr*v2xp + 3._pr*v3xp)
              s2xp = 13._pr/12._pr*(v2xp - 2._pr*v3xp + v4xp)*(v2xp - 2._pr*v3xp + v4xp) + 0.25_pr*(v2xp - v4xp)*(v2xp - v4xp)
              s3xp = 13._pr/12._pr*(v3xp - 2._pr*v4xp + v5xp)*(v3xp - 2._pr*v4xp + v5xp) + 0.25_pr*(3._pr*v3xp - 4._pr*v4xp + v5xp)&
*(3._pr*v3xp - 4._pr*v4xp + v5xp)
              a1xp = 0.1_pr/((eepsilon + s1xp)*(eepsilon + s1xp))
              a2xp = 0.6_pr/((eepsilon + s2xp)*(eepsilon + s2xp))
              a3xp = 0.3_pr/((eepsilon + s3xp)*(eepsilon + s3xp))
              w1xp = a1xp/(a1xp + a2xp + a3xp)
              w2xp = a2xp/(a1xp + a2xp + a3xp)
              w3xp = a3xp/(a1xp + a2xp + a3xp)
              drhodxp = w1xp*(v1xp/3._pr - 7._pr/6._pr*v2xp + 11._pr/6._pr*v3xp) + w2xp*(-v2xp/6._pr + 5._pr/6._pr*v3xp + v4xp/&
3._pr) + w3xp*(v3xp/3._pr + 5._pr/6._pr*v4xp - v5xp/6._pr)

              v1yp = (rho_ij3p - rho_ij2p)/dy
              v2yp = (rho_ij2p - rho_ij1p)/dy
              v3yp = (rho_ij1p - rho_ij)/dy
              v4yp = (rho_ij - rho_ij1m)/dy
              v5yp = (rho_ij1m - rho_ij2m)/dy

              s1yp = 13._pr/12._pr*(v1yp - 2._pr*v2yp + v3yp)*(v1yp - 2._pr*v2yp + v3yp) + 0.25_pr*(v1yp - 4._pr*v2yp + 3._pr*v3yp)&
*(v1yp - 4._pr*v2yp + 3._pr*v3yp)
              s2yp = 13._pr/12._pr*(v2yp - 2._pr*v3yp + v4yp)*(v2yp - 2._pr*v3yp + v4yp) + 0.25_pr*(v2yp - v4yp)*(v2yp - v4yp)
              s3yp = 13._pr/12._pr*(v3yp - 2._pr*v4yp + v5yp)*(v3yp - 2._pr*v4yp + v5yp) + 0.25_pr*(3._pr*v3yp - 4._pr*v4yp + v5yp)&
*(3._pr*v3yp - 4._pr*v4yp + v5yp)
              a1yp = 0.1_pr/((eepsilon + s1yp)*(eepsilon + s1yp))
              a2yp = 0.6_pr/((eepsilon + s2yp)*(eepsilon + s2yp))
              a3yp = 0.3_pr/((eepsilon + s3yp)*(eepsilon + s3yp))
              w1yp = a1yp/(a1yp + a2yp + a3yp)
              w2yp = a2yp/(a1yp + a2yp + a3yp)
              w3yp = a3yp/(a1yp + a2yp + a3yp)
              drhodyp = w1yp*(v1yp/3._pr - 7._pr/6._pr*v2yp + 11._pr/6._pr*v3yp) + w2yp*(-v2yp/6._pr + 5._pr/6._pr*v3yp + v4yp/&
3._pr) + w3yp*(v3yp/3._pr + 5._pr/6._pr*v4yp - v5yp/6._pr)

              xpp = max(drhodxp,0._pr)
              xpm = -min(drhodxp,0._pr)
              xmp = max(drhodxm,0._pr)
              xmm = -min(drhodxm,0._pr)

              ypp = max(drhodyp,0._pr)
              ypm = -min(drhodyp,0._pr)
              ymp = max(drhodym,0._pr)
              ymm = -min(drhodym,0._pr)


              if (U0_ij > 0._pr) then
                      G = sqrt(max(xmp*xmp,xpm*xpm) + max(ymp*ymp,ypm*ypm)) - 1._pr
              else
                      G = sqrt(max(xmm*xmm,xpp*xpp) + max(ymm*ymm,ypp*ypp)) - 1._pr
              endif

              u_(i,j) = -dt_*S*G
      endif
    enddo
  enddo
end


subroutine updateDistance3D(rhoSlices,gradPhi,zslice)
  real(pr),dimension(:,:),intent(out) :: gradPhi
  real(pr), dimension(:,:,:), intent(inout) :: rhoSlices
  real(pr), intent(in) :: zslice
  real(pr) :: nr, nrtmp1, nrtmp2, nrprec, epsDist, dtDist, cflDist, h
  integer :: nIterMaxDist
  real(pr), dimension(size(rhoSlices,1),size(rhoSlices,2),size(rhoSlices,3)) :: rhoDist, rho0Dist
  real(pr) :: drhodx,drhody,drhodz
  integer :: iter,i,j
  real(pr) :: rhoki1mj, rhoki1pj, rhokij1m, rhokij1p, rhok1pij, rhok1mij

  cflDist = 0.5_pr
  epsDist = 0.01_pr
  nIterMaxDist = 4000
  dtDist = cflDist/(1._pr/dx + 1._pr/dy + 1._pr/dz)
  rho0Dist = rhoSlices
  h = max(dx,max(dy,dz))

  call distanceIteration3D(rho0Dist, dtDist, rhoSlices, rhoDist)
  call Norm23D(rhoDist,nrtmp1)
  call Norm23D(rhoSlices,nrtmp2)
  rhoSlices = rhoSlices + rhoDist
  nr = nrtmp1 / nrtmp2
  nrprec = 0._pr
  iter = 1

  do while ((nr > epsDist).and.(iter < nIterMaxDist))
    call distanceIteration3D(rho0Dist, dtDist, rhoSlices, rhoDist)
    nrprec = nr
    rhoSlices = rhoSlices + rhoDist
    call Norm23D(rhoDist,nrtmp1)
    call Norm23D(rhoSlices,nrtmp2)
    nr = nrtmp1 / nrtmp2
    iter = iter + 1
  enddo
  do j=1,N
    do i=1,N
                call neumannBC(rhoSlices,nint(zslice),i-1,j,rhoki1mj)
                call neumannBC(rhoSlices,nint(zslice),i+1,j,rhoki1pj)
                call neumannBC(rhoSlices,nint(zslice),i,j-1,rhokij1m)
                call neumannBC(rhoSlices,nint(zslice),i,j+1,rhokij1p)
                call neumannBC(rhoSlices,nint(zslice)+1,i,j,rhok1pij)
                call neumannBC(rhoSlices,nint(zslice)-1,i,j,rhok1mij)
                drhodx = (rhoki1pj-rhoki1mj)/(2*dx)
                drhody = (rhokij1p-rhokij1m)/(2*dy)
                drhodz = (rhok1pij-rhok1mij)/(2*dz)
              gradPhi(i,j) = sqrt(drhodx**2 + drhody**2 + drhodz**2)
  enddo
  enddo

  if (iter >= nIterMaxDist) then
          write(*,*) "Distanciation did not converge, res norm = ", nr
  else
          write(*,*) "Distanciation nIter = ", iter, " nr = ", nr
  endif
end


subroutine distanceIteration3D(rho0Dist, dt_, rhoSlices_, u_)
  real(pr), intent(in) :: dt_
  real(pr), dimension(:,:,:), intent(in) :: rho0Dist
  real(pr), dimension(:,:,:), intent(inout) :: rhoSlices_
  real(pr), dimension(:,:,:), intent(inout) :: u_
  real(pr) :: h, x, y, z, xpp, xmp, ypp, ymp, zpp, zmp, xpm, xmm, ypm, ymm, zpm, zmm, a, b, c
  real(pr) :: S, G, dUi0, U0_ijk, Un_ijk, U0_im1jk, U0_ip1jk, U0_ijm1k, U0_ijp1k, U0_ijkm1, U0_ijkp1
  real(pr) :: v1xm, v2xm, v3xm, v4xm, v5xm
  real(pr) :: s1xm, s2xm, s3xm
  real(pr) :: a1xm, a2xm, a3xm
  real(pr) :: w1xm, w2xm, w3xm
  real(pr) :: drhodxm
  real(pr) :: v1xp, v2xp, v3xp, v4xp, v5xp
  real(pr) :: s1xp, s2xp, s3xp
  real(pr) :: a1xp, a2xp, a3xp
  real(pr) :: w1xp, w2xp, w3xp
  real(pr) :: drhodxp
  real(pr) :: v1ym, v2ym, v3ym, v4ym, v5ym
  real(pr) :: s1ym, s2ym, s3ym
  real(pr) :: a1ym, a2ym, a3ym
  real(pr) :: w1ym, w2ym, w3ym
  real(pr) :: drhodym
  real(pr) :: v1yp, v2yp, v3yp, v4yp, v5yp
  real(pr) :: s1yp, s2yp, s3yp
  real(pr) :: a1yp, a2yp, a3yp
  real(pr) :: w1yp, w2yp, w3yp
  real(pr) :: drhodyp
  real(pr) :: v1zm, v2zm, v3zm, v4zm, v5zm
  real(pr) :: s1zm, s2zm, s3zm
  real(pr) :: a1zm, a2zm, a3zm
  real(pr) :: w1zm, w2zm, w3zm
  real(pr) :: drhodzm
  real(pr) :: v1zp, v2zp, v3zp, v4zp, v5zp
  real(pr) :: s1zp, s2zp, s3zp
  real(pr) :: a1zp, a2zp, a3zp
  real(pr) :: w1zp, w2zp, w3zp
  real(pr) :: drhodzp

  integer :: i, j, k
  real(pr) :: rho_i1mjk, rho_i1pjk, rho_ij1mk, rho_ij1pk, rho_ijk1m, rho_ijk1p, rho_i2mjk, rho_i2pjk, rho_ij2mk, rho_ij2pk,&
 rho_ijk2m, rho_ijk2p, rho_i3mjk, rho_i3pjk, rho_ij3mk, rho_ij3pk, rho_ijk3m, rho_ijk3p, rho_ijk

  h = sqrt(3._pr)*dx*dy*dz/(sqrt(dx*dx*dz*dz + dy*dy*dz*dz + dx*dx*dy*dy))

  do k=1,N
    do j=1,N
      do i=1,N
        U0_ijk = rho0Dist(i,j,k)

        call neumannBC(rho0Dist,i-1,j,k,U0_im1jk)
        call neumannBC(rho0Dist,i+1,j,k,U0_ip1jk)
        call neumannBC(rho0Dist,i,j-1,k,U0_ijm1k)
        call neumannBC(rho0Dist,i,j+1,k,U0_ijp1k)
        call neumannBC(rho0Dist,i,j,k-1,U0_ijkm1)
        call neumannBC(rho0Dist,i,j,k+1,U0_ijkp1)


        Un_ijk = rhoSlices_(i,j,k)
        S = sign(1._pr,U0_ijk)

        if ((U0_ijk*U0_im1jk<0).or.(U0_ijk*U0_ip1jk<0).or.(U0_ijk*U0_ijm1k<0).or.(U0_ijk*U0_ijp1k<0).or.(U0_ijk*U0_ijkm1<0).or.&
(U0_ijk*U0_ijkp1<0)) then
                x = U0_ip1jk - U0_im1jk
                y = U0_ijp1k - U0_ijm1k
                z = U0_ijkp1 - U0_ijkm1
                a = 0.5_pr*sqrt(x*x + y*y + z*z)

                x = U0_ip1jk - U0_ijk
                y = U0_ijp1k - U0_ijk
                z = U0_ijkp1 - U0_ijk
                b = sqrt(x*x + y*y + z*z)

                x = U0_ijk - U0_im1jk
                y = U0_ijk - U0_ijm1k
                z = U0_ijk - U0_ijkm1
                c = sqrt(x*x + y*y + z*z)

                dUi0 = max(max(a,b),max(c,eepsilon))
                G = abs(Un_ijk)

                u_(i,j,k) = -dt_*(S*G/h - U0_ijk/dUi0)
        else
                call neumannBC(rhoSlices_,i-3,j,k,rho_i3mjk)
                call neumannBC(rhoSlices_,i-2,j,k,rho_i2mjk)
                call neumannBC(rhoSlices_,i-1,j,k,rho_i1mjk)
                call neumannBC(rhoSlices_,i+1,j,k,rho_i1pjk)
                call neumannBC(rhoSlices_,i+2,j,k,rho_i2pjk)
                call neumannBC(rhoSlices_,i+3,j,k,rho_i3pjk)
                call neumannBC(rhoSlices_,i,j-3,k,rho_ij3mk)
                call neumannBC(rhoSlices_,i,j-2,k,rho_ij2mk)
                call neumannBC(rhoSlices_,i,j-1,k,rho_ij1mk)
                call neumannBC(rhoSlices_,i,j+1,k,rho_ij1pk)
                call neumannBC(rhoSlices_,i,j+2,k,rho_ij2pk)
                call neumannBC(rhoSlices_,i,j+3,k,rho_ij3pk)
                call neumannBC(rhoSlices_,i,j,k-3,rho_ijk3m)
                call neumannBC(rhoSlices_,i,j,k-2,rho_ijk2m)
                call neumannBC(rhoSlices_,i,j,k-1,rho_ijk1m)
                call neumannBC(rhoSlices_,i,j,k+1,rho_ijk1p)
                call neumannBC(rhoSlices_,i,j,k+2,rho_ijk2p)
                call neumannBC(rhoSlices_,i,j,k+3,rho_ijk3p)
                rho_ijk = rhoSlices_(i,j,k)

                v1xm = (rho_i2mjk - rho_i3mjk)/dx
                v2xm = (rho_i1mjk - rho_i2mjk)/dx
                v3xm = (rho_ijk - rho_i1mjk)/dx
                v4xm = (rho_i1pjk - rho_ijk)/dx
                v5xm = (rho_i2pjk - rho_i1pjk)/dx

                s1xm = 13._pr/12._pr*(v1xm - 2._pr*v2xm + v3xm)*(v1xm - 2._pr*v2xm + v3xm) + 0.25_pr*(v1xm - 4._pr*v2xm + 3._pr*&
v3xm)&
*(  v1xm - 4._pr*v2xm + 3._pr*v3xm)
                s2xm = 13._pr/12._pr*(v2xm - 2._pr*v3xm + v4xm)*(v2xm - 2._pr*v3xm + v4xm) + 0.25_pr*(v2xm - v4xm)*(v2xm - v4xm)
                s3xm = 13._pr/12._pr*(v3xm - 2._pr*v4xm + v5xm)*(v3xm - 2._pr*v4xm + v5xm) + 0.25_pr*(3._pr*v3xm - 4._pr*v4xm +&
 v5xm)*(  3._pr*v3xm - 4._pr*v4xm + v5xm)
                a1xm = 0.1_pr/((eepsilon + s1xm)*(eepsilon + s1xm))
                a2xm = 0.6_pr/((eepsilon + s2xm)*(eepsilon + s2xm))
                a3xm = 0.3_pr/((eepsilon + s3xm)*(eepsilon + s3xm))
                w1xm = a1xm/(a1xm + a2xm + a3xm)
                w2xm = a2xm/(a1xm + a2xm + a3xm)
                w3xm = a3xm/(a1xm + a2xm + a3xm)
                drhodxm = w1xm*(v1xm/3._pr - 7._pr/6._pr*v2xm + 11._pr/6._pr*v3xm) + w2xm*(-v2xm/6._pr + 5._pr/6._pr*v3xm + v4xm/&
3._pr) + w3xm*(v3xm/3._pr + 5._pr/6._pr*v4xm - v5xm/6._pr)

                v1zm = (rho_ij2mk - rho_ij3mk)/dz
                v2zm = (rho_ij1mk - rho_ij2mk)/dz
                v3zm = (rho_ijk - rho_ij1mk)/dz
                v4zm = (rho_ij1pk - rho_ijk)/dz
                v5zm = (rho_ij2pk - rho_ij1pk)/dz

                s1zm = 13._pr/12._pr*(v1zm - 2._pr*v2zm + v3zm)*(v1zm - 2._pr*v2zm + v3zm) + 0.25_pr*(v1zm - 4._pr*v2zm + 3._pr*&
v3zm)*(  v1zm - 4._pr*v2zm + 3._pr*v3zm)
                s2zm = 13._pr/12._pr*(v2zm - 2._pr*v3zm + v4zm)*(v2zm - 2._pr*v3zm + v4zm) + 0.25_pr*(v2zm - v4zm)*(v2zm - v4zm)
                s3zm = 13._pr/12._pr*(v3zm - 2._pr*v4zm + v5zm)*(v3zm - 2._pr*v4zm + v5zm) + 0.25_pr*(3._pr*v3zm - 4._pr*v4zm +&
 v5zm)*(  3._pr*v3zm - 4._pr*v4zm + v5zm)
                a1zm = 0.1_pr/((eepsilon + s1zm)*(eepsilon + s1zm))
                a2zm = 0.6_pr/((eepsilon + s2zm)*(eepsilon + s2zm))
                a3zm = 0.3_pr/((eepsilon + s3zm)*(eepsilon + s3zm))
                w1zm = a1zm/(a1zm + a2zm + a3zm)
                w2zm = a2zm/(a1zm + a2zm + a3zm)
                w3zm = a3zm/(a1zm + a2zm + a3zm)
                drhodzm = w1zm*(v1zm/3._pr - 7._pr/6._pr*v2zm + 11._pr/6._pr*v3zm) + w2zm*(-v2zm/6._pr + 5._pr/6._pr*v3zm + v4zm/&
3._pr) + w3zm*(v3zm/3._pr + 5._pr/6._pr*v4zm - v5zm/6._pr)

                v1ym = (rho_ijk2m - rho_ijk3m)/dy
                v2ym = (rho_ijk1m - rho_ijk2m)/dy
                v3ym = (rho_ijk - rho_ijk1m)/dy
                v4ym = (rho_ijk1p - rho_ijk)/dy
                v5ym = (rho_ijk2p - rho_ijk1p)/dy

                s1ym = 13._pr/12._pr*(v1ym - 2._pr*v2ym + v3ym)*(v1ym - 2._pr*v2ym + v3ym) + 0.25_pr*(v1ym - 4._pr*v2ym + 3._pr*&
v3ym)*(  v1ym - 4._pr*v2ym + 3._pr*v3ym)
                s2ym = 13._pr/12._pr*(v2ym - 2._pr*v3ym + v4ym)*(v2ym - 2._pr*v3ym + v4ym) + 0.25_pr*(v2ym - v4ym)*(v2ym - v4ym)
                s3ym = 13._pr/12._pr*(v3ym - 2._pr*v4ym + v5ym)*(v3ym - 2._pr*v4ym + v5ym) + 0.25_pr*(3._pr*v3ym - 4._pr*v4ym +&
 v5ym)*(  3._pr*v3ym - 4._pr*v4ym + v5ym)
                a1ym = 0.1_pr/((eepsilon + s1ym)*(eepsilon + s1ym))
                a2ym = 0.6_pr/((eepsilon + s2ym)*(eepsilon + s2ym))
                a3ym = 0.3_pr/((eepsilon + s3ym)*(eepsilon + s3ym))
                w1ym = a1ym/(a1ym + a2ym + a3ym)
                w2ym = a2ym/(a1ym + a2ym + a3ym)
                w3ym = a3ym/(a1ym + a2ym + a3ym)
                drhodym = w1ym*(v1ym/3._pr - 7._pr/6._pr*v2ym + 11._pr/6._pr*v3ym) + w2ym*(-v2ym/6._pr + 5._pr/6._pr*v3ym + v4ym/&
3._pr) + w3ym*(v3ym/3._pr + 5._pr/6._pr*v4ym - v5ym/6._pr)

                v1xp = (rho_i3pjk - rho_i2pjk)/dx
                v2xp = (rho_i2pjk - rho_i1pjk)/dx
                v3xp = (rho_i1pjk - rho_ijk)/dx
                v4xp = (rho_ijk - rho_i1mjk)/dx
                v5xp = (rho_i1mjk - rho_i2mjk)/dx

                s1xp = 13._pr/12._pr*(v1xp - 2._pr*v2xp + v3xp)*(v1xp - 2._pr*v2xp + v3xp) + 0.25_pr*(v1xp - 4._pr*v2xp + 3._pr*&
v3xp)*(  v1xp - 4._pr*v2xp + 3._pr*v3xp)
                s2xp = 13._pr/12._pr*(v2xp - 2._pr*v3xp + v4xp)*(v2xp - 2._pr*v3xp + v4xp) + 0.25_pr*(v2xp - v4xp)*(v2xp - v4xp)
                s3xp = 13._pr/12._pr*(v3xp - 2._pr*v4xp + v5xp)*(v3xp - 2._pr*v4xp + v5xp) + 0.25_pr*(3._pr*v3xp - 4._pr*v4xp +&
 v5xp)*(  3._pr*v3xp - 4._pr*v4xp + v5xp)
                a1xp = 0.1_pr/((eepsilon + s1xp)*(eepsilon + s1xp))
                a2xp = 0.6_pr/((eepsilon + s2xp)*(eepsilon + s2xp))
                a3xp = 0.3_pr/((eepsilon + s3xp)*(eepsilon + s3xp))
                w1xp = a1xp/(a1xp + a2xp + a3xp)
                w2xp = a2xp/(a1xp + a2xp + a3xp)
                w3xp = a3xp/(a1xp + a2xp + a3xp)
                drhodxp = w1xp*(v1xp/3._pr - 7._pr/6._pr*v2xp + 11._pr/6._pr*v3xp) + w2xp*(-v2xp/6._pr + 5._pr/6._pr*v3xp + v4xp/&
3._pr) + w3xp*(v3xp/3._pr + 5._pr/6._pr*v4xp - v5xp/6._pr)

                v1zp = (rho_ij3pk - rho_ij2pk)/dz
                v2zp = (rho_ij2pk - rho_ij1pk)/dz
                v3zp = (rho_ij1pk - rho_ijk)/dz
                v4zp = (rho_ijk - rho_ij1mk)/dz
                v5zp = (rho_ij1mk - rho_ij2mk)/dz

                s1zp = 13._pr/12._pr*(v1zp - 2._pr*v2zp + v3zp)*(v1zp - 2._pr*v2zp + v3zp) + 0.25_pr*(v1zp - 4._pr*v2zp + 3._pr*&
v3zp)*(  v1zp - 4._pr*v2zp + 3._pr*v3zp)
                s2zp = 13._pr/12._pr*(v2zp - 2._pr*v3zp + v4zp)*(v2zp - 2._pr*v3zp + v4zp) + 0.25_pr*(v2zp - v4zp)*(v2zp - v4zp)
                s3zp = 13._pr/12._pr*(v3zp - 2._pr*v4zp + v5zp)*(v3zp - 2._pr*v4zp + v5zp) + 0.25_pr*(3._pr*v3zp - 4._pr*v4zp +&
 v5zp)*(  3._pr*v3zp - 4._pr*v4zp + v5zp)
                a1zp = 0.1_pr/((eepsilon + s1zp)*(eepsilon + s1zp))
                a2zp = 0.6_pr/((eepsilon + s2zp)*(eepsilon + s2zp))
                a3zp = 0.3_pr/((eepsilon + s3zp)*(eepsilon + s3zp))
                w1zp = a1zp/(a1zp + a2zp + a3zp)
                w2zp = a2zp/(a1zp + a2zp + a3zp)
                w3zp = a3zp/(a1zp + a2zp + a3zp)
                drhodzp = w1zp*(v1zp/3._pr - 7._pr/6._pr*v2zp + 11._pr/6._pr*v3zp) + w2zp*(-v2zp/6._pr + 5._pr/6._pr*v3zp + v4zp/&
3._pr) + w3zp*(v3zp/3._pr + 5._pr/6._pr*v4zp - v5zp/6._pr)

                v1yp = (rho_ijk3p - rho_ijk2p)/dy
                v2yp = (rho_ijk2p - rho_ijk1p)/dy
                v3yp = (rho_ijk1p - rho_ijk)/dy
                v4yp = (rho_ijk - rho_ijk1m)/dy
                v5yp = (rho_ijk1m - rho_ijk2m)/dy

                s1yp = 13._pr/12._pr*(v1yp - 2._pr*v2yp + v3yp)*(v1yp - 2._pr*v2yp + v3yp) + 0.25_pr*(v1yp - 4._pr*v2yp + 3._pr*&
v3yp)*(  v1yp - 4._pr*v2yp + 3._pr*v3yp)
                s2yp = 13._pr/12._pr*(v2yp - 2._pr*v3yp + v4yp)*(v2yp - 2._pr*v3yp + v4yp) + 0.25_pr*(v2yp - v4yp)*(v2yp - v4yp)
                s3yp = 13._pr/12._pr*(v3yp - 2._pr*v4yp + v5yp)*(v3yp - 2._pr*v4yp + v5yp) + 0.25_pr*(3._pr*v3yp - 4._pr*v4yp +&
 v5yp)*(  3._pr*v3yp - 4._pr*v4yp + v5yp)
                a1yp = 0.1_pr/((eepsilon + s1yp)*(eepsilon + s1yp))
                a2yp = 0.6_pr/((eepsilon + s2yp)*(eepsilon + s2yp))
                a3yp = 0.3_pr/((eepsilon + s3yp)*(eepsilon + s3yp))
                w1yp = a1yp/(a1yp + a2yp + a3yp)
                w2yp = a2yp/(a1yp + a2yp + a3yp)
                w3yp = a3yp/(a1yp + a2yp + a3yp)
                drhodyp = w1yp*(v1yp/3._pr - 7._pr/6._pr*v2yp + 11._pr/6._pr*v3yp) + w2yp*(-v2yp/6._pr + 5._pr/6._pr*v3yp + v4yp/&
3._pr) + w3yp*(v3yp/3._pr + 5._pr/6._pr*v4yp - v5yp/6._pr)

                xpp = max(drhodxp,0._pr)
                xpm = -min(drhodxp,0._pr)
                xmp = max(drhodxm,0._pr)
                xmm = -min(drhodxm,0._pr)

                ypp = max(drhodyp,0._pr)
                ypm = -min(drhodyp,0._pr)
                ymp = max(drhodym,0._pr)
                ymm = -min(drhodym,0._pr)

                zpp = max(drhodzp,0._pr)
                zpm = -min(drhodzp,0._pr)
                zmp = max(drhodzm,0._pr)
                zmm = -min(drhodzm,0._pr)

                if (U0_ijk > 0._pr) then
                        G = sqrt(max(xmp*xmp,xpm*xpm) + max(ymp*ymp,ypm*ypm) + max(zmp*zmp,zpm*zpm)) - 1._pr
                else
                        G = sqrt(max(xmm*xmm,xpp*xpp) + max(ymm*ymm,ypp*ypp) + max(zmm*zmm,zpp*zpp)) - 1._pr
                endif

                u_(i,j,k) = -dt_*S*G
        endif
      enddo
    enddo
  enddo
end

end
