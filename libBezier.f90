module libBezier
  use doubler
  use AdvectionProblem

  implicit none

contains  

  subroutine intersection(x1,x2,y1,y2,xc,yc,d1l,d2l,d1r,d2r,xl,yl,xr,yr,errorl,errorr,deltal,deltar)
    implicit none
    real(pr),intent(in) :: x1,x2,y1,y2,d1l,d2l,d1r,d2r
    real(pr),intent(in) :: xc,yc
    real(pr),intent(out) :: xl,yl,xr,yr,deltal,deltar
    integer,intent(out) :: errorl,errorr
    real(pr) :: A,Bl,Br,Cl,Cr,xal,xar,xbl,xbr,yal,yar,ybl,ybr,cstl,cstr

    errorr = 0

    if (abs(y2-y1)/dx>10*eepsilon) then
       cstl = ((d1l*d1l - x1*x1 - y1*y1) - (d2l*d2l - x2*x2 - y2*y2))/(2*(y2-y1))
       A = 1._pr + (x2-x1)*(x2-x1)/((y2-y1)*(y2-y1))
       Bl = -2*(cstl-y1)*(x2-x1)/(y2-y1) - 2*x1
       Cl = x1*x1 - d1l*d1l + (cstl-y1)*(cstl-y1)

       call racine(A,Bl,Cl,xal,xbl,errorl,deltal)

       yal = cstl - xal*(x2-x1)/(y2-y1)
       yar = cstr - xar*(x2-x1)/(y2-y1)
       ybl = cstl - xbl*(x2-x1)/(y2-y1)
       ybr = cstr - xbr*(x2-x1)/(y2-y1)
       if (det(xal,yal,x2,y2,x1,y1) <= 0) then
          xl = xal
          yl = yal
       else
          xl = xbl
          yl = ybl
       endif
       if (det(xar,yar,x2,y2,x1,y1) >= 0) then
          xr = xar
          yr = yar
       else
          xr = xbr
          yr = ybr
       endif
    else
       xl = ((d1l*d1l - x1*x1) - (d2l*d2l - x2*x2))/(2*(x2-x1))
       xr = ((d1r*d1r - x1*x1) - (d2r*d2r - x2*x2))/(2*(x2-x1))

       A = 1._pr
       Bl = -2*y1
       Br = -2*y1
       Cl = (xl-x1)*(xl-x1) + y1*y1 - d1l*d1l 
       Cr = (xr-x1)*(xr-x1) + y1*y1 - d1r*d1r 
       call racine(A,Bl,Cl,yal,ybl,errorl,deltal)

       if (det(xl,yal,x2,y2,x1,y1) <= 0) then
          yl = yal
       else
          yl = ybl
       endif
       if (det(xr,yar,x2,y2,x1,y1) >= 0) then
          yr = yar
       else
          yr = ybr
       endif
    endif
  end subroutine intersection


  subroutine intersectionR(x1,x2,y1,y2,d1l,d2l,d1r,d2r,xl,yl,xr,yr,errorl,errorr,deltal,deltar)
    implicit none
    real(pr),intent(in) :: x1,x2,y1,y2,d1l,d2l,d1r,d2r
    real(pr),intent(out) :: xl,yl,xr,yr,deltal,deltar
    integer,intent(out) :: errorl,errorr
    real(pr) :: A,B,Cl,Cr,xal,xar,xbl,xbr,yal,yar,ybl,ybr,cst

    cst = (x2*x1 + y2*y1 - x2*x2 - x1*x1 - y2*y2 - y1*y1)/(y2-y1)

    A = 1._pr + (x2-x1)*(x2-x1)/((y2-y1)*(y2-y1))
    B = -2*(cst-y1)*(x2-x1)/(y2-y1) - 2*x1
    Cl = x1*x1 - d1l*d1l + (cst-y1)*(cst-y1)
    Cr = x1*x1 - d1r*d1r + (cst-y1)*(cst-y1)
    call racine(A,B,Cl,xal,xbl,errorl,deltal)
    call racine(A,B,Cr,xar,xbr,errorr,deltar)
    if (abs(y2-y1)/dx>eepsilon) then
       yal = ((x2*x1 + y2*y1 - x2*x2 - x1*x1 - y2*y2 - y1*y1) - xal*(x2-x1))/(y2-y1)
       ybl = ((x2*x1 + y2*y1 - x2*x2 - x1*x1 - y2*y2 - y1*y1) - xbl*(x2-x1))/(y2-y1)
       yar = ((x2*x1 + y2*y1 - x2*x2 - x1*x1 - y2*y2 - y1*y1) - xar*(x2-x1))/(y2-y1)
       ybr = ((x2*x1 + y2*y1 - x2*x2 - x1*x1 - y2*y2 - y1*y1) - xbr*(x2-x1))/(y2-y1)
       if (det(xal,yal,x2,y2,x1,y1) <= 0) then
          xl = xal
          yl = yal
       else
          xl = xbl
          yl = ybl
       endif
       if (det(xar,yar,x2,y2,x1,y1) >= 0) then
          xr = xar
          yr = yar
       else
          xr = xbr
          yr = ybr
       endif
    else
       xl = x1
       yl = y1 + d1l
       xr = x1
       yr = y1 - d1r
    endif
  end subroutine intersectionR


  subroutine racine(a,b,c,x1,x2,error,delta)
    implicit none
    real(pr),intent(in) :: a,b,c
    real(pr),intent(out) :: x1,x2
    integer,intent(out) :: error
    real(pr),intent(out) :: delta

    delta = b*b - 4*a*c
    if (delta>0._pr) then
       x1 = (-b-sqrt(abs(delta)))/(2*a)
       x2 = (-b+sqrt(abs(delta)))/(2*a)
       error = 0

    else if (abs(delta)/(dx*dx)<10*eepsilon) then
       x1 = -b/(2*a)
       x2 = x1
       error = 0
    else 
       error = 1
       write(*,*) "error ROOT FINDING : delta<0    ",a," ",b," ",c," ",delta

    endif
  end subroutine racine


  logical function appartient(tab,i,j,lmax)
    implicit none
    integer,dimension(:,:),intent(in) :: tab
    integer,intent(in) :: i,j,lmax
    integer :: l
    logical :: bool

    l = 1
    bool = .false.
    do while ((bool.eqv..false.).and.(l<lmax).and.(l<size(tab,1))) 
       if ((tab(l,1)==i).and.(tab(l,2)==j)) bool = .true.
       l = l+1
    enddo

    appartient = bool
  end function appartient

  logical function appartientElmtTh(tab,theta,l)
    implicit none
    integer,dimension(:,:),intent(in) :: tab
    integer,intent(in) :: theta,l
    integer :: k
    logical :: bool

    k = 1
    bool = .false.
    do while ((bool.eqv..false.).and.(k<size(tab,2)+1)) 
       if (tab(theta,k)==l) bool = .true.
       k = k+1
    enddo

    appartientElmtTh = bool
  end function appartientElmtTh

  logical function appartientElmtL(tab,theta,l)
    implicit none
    integer,dimension(:,:),intent(in) :: tab
    integer,intent(in) :: theta,l
    integer :: k
    logical :: bool

    k = 1
    bool = .false.
    do while ((bool.eqv..false.).and.(k<size(tab,1)+1)) 
       if (tab(k,l)==theta) bool = .true.
       k = k+1
    enddo

    appartientElmtL = bool
  end function appartientElmtL



  subroutine initialize4(points_control)
    implicit none
    real(pr),dimension(:,:),intent(inout) :: points_control
    integer :: l

    points_control(1,1) = 0._pr
    points_control(1,2) = 0._pr

    points_control(2,1) = 1._pr
    points_control(2,2) = 1._pr

    points_control(3,1) = 2._pr
    points_control(3,2) = 2._pr

    points_control(4,1) = 4._pr
    points_control(4,2) = 4._pr
  end subroutine initialize4

  subroutine initialize7(points_control)
    implicit none
    real(pr),dimension(:,:),intent(inout) :: points_control
    integer :: l

    points_control(1,1) = 0._pr
    points_control(1,2) = 0._pr

    points_control(2,1) = 0.3_pr
    points_control(2,2) = 1._pr

    points_control(3,1) = 1._pr
    points_control(3,2) = 1._pr

    points_control(4,1) = 1._pr
    points_control(4,2) = 0._pr

    points_control(5,1) = 1._pr
    points_control(5,2) = -1._pr

    points_control(6,1) = 1.5_pr
    points_control(6,2) = -0.5_pr

    points_control(7,1) = 2._pr
    points_control(7,2) = 0._pr

    points_control(:,1) = points_control(:,1) + 1
    points_control(:,2) = points_control(:,2) + 1
  end subroutine initialize7


  subroutine pointBezier3(points_control,t,px,py)
    implicit none
    real(pr),dimension(:,:),intent(in) :: points_control
    real(pr),intent(in) :: t
    real(pr),intent(out) :: px, py
    real(pr) :: x,y
    real(pr) :: Ax,Ay,Bx,By

    x = (1-t)*(1-t)
    y = t*t

    Ax = (1-t)*x*points_control(1,1) + 3*t*x*points_control(2,1)
    Ay = (1-t)*x*points_control(1,2) + 3*t*x*points_control(2,2)

    Bx = 3*y*(1-t)*points_control(3,1) + y*t*points_control(4,1)
    By = 3*y*(1-t)*points_control(3,2) + y*t*points_control(4,2)

    px = Ax + Bx
    py = Ay + By
  end subroutine pointBezier3


  subroutine reduction1D(points_control,t,N,points_sortie)
    implicit none
    real(pr),intent(in) :: t
    integer,intent(in) :: N
    real(pr),dimension(:),intent(in) :: points_control
    real(pr),dimension(:),intent(out) :: points_sortie
    integer :: l

    do l=1,N-1
       points_sortie(l) = (1-t)*points_control(l) + t*points_control(l+1)
    enddo
  end subroutine reduction1D

  subroutine reduction(points_control,t,N,points_sortie)
    implicit none
    real(pr),intent(in) :: t
    integer,intent(in) :: N
    real(pr),dimension(:,:),intent(in) :: points_control
    real(pr),dimension(:,:),intent(out) :: points_sortie
    integer :: l

    do l=1,N-1
       points_sortie(l,1) = (1-t)*points_control(l,1) + t*points_control(l+1,1)
       points_sortie(l,2) = (1-t)*points_control(l,2) + t*points_control(l+1,2)
    enddo
  end subroutine reduction

  subroutine reduction3D(points_control,t,N,points_sortie)
    implicit none
    real(pr),intent(in) :: t
    integer,intent(in) :: N
    real(pr),dimension(:,:),intent(in) :: points_control
    real(pr),dimension(:,:),intent(out) :: points_sortie
    integer :: l

    do l=1,N-1
       points_sortie(l,1) = (1-t)*points_control(l,1) + t*points_control(l+1,1)
       points_sortie(l,2) = (1-t)*points_control(l,2) + t*points_control(l+1,2)
       points_sortie(l,3) = (1-t)*points_control(l,3) + t*points_control(l+1,3)
    enddo
  end subroutine reduction3D


  subroutine pointsBezierN1D(points_control,t,P)
    implicit none
    real(pr),intent(in) :: t
    real(pr),dimension(:),intent(in) :: points_control
    real(pr),intent(out) :: P
    real(pr),dimension(size(points_control)) :: tab
    integer :: N,l

    N = size(points_control)
    do l=1,N
       tab(l) = points_control(l)
    enddo
    do while (N>1)
       call reduction1D(tab,t,N,tab)
       N = N-1
    enddo

    P = tab(1)
  end subroutine pointsBezierN1D

  subroutine pointsBezierN(points_control,t,Px,Py)
    implicit none
    real(pr),intent(in) :: t
    real(pr),dimension(:,:),intent(in) :: points_control
    real(pr),intent(out) :: Px,Py
    real(pr),dimension(size(points_control,1),2) :: tab
    integer :: N,l

    N = size(points_control,1)
    do l=1,N
       tab(l,1) = points_control(l,1)
       tab(l,2) = points_control(l,2)
    enddo
    do while (N>1)
       call reduction(tab,t,N,tab)
       N = N-1
    enddo

    Px = tab(1,1)
    Py = tab(1,2)
  end subroutine pointsBezierN

  subroutine pointsBezierN3D(points_control,t,Px,Py,Pz)
    implicit none
    real(pr),intent(in) :: t
    real(pr),dimension(:,:),intent(in) :: points_control
    real(pr),intent(out) :: Px,Py,Pz
    real(pr),dimension(size(points_control,1),3) :: tab
    integer :: N,l

    N = size(points_control,1)
    do l=1,N
       tab(l,1) = points_control(l,1)
       tab(l,2) = points_control(l,2)
       tab(l,3) = points_control(l,3)
    enddo
    do while (N>1)
       call reduction3D(tab,t,N,tab)
       N = N-1
    enddo

    Px = tab(1,1)
    Py = tab(1,2)
    Pz = tab(1,3)
  end subroutine pointsBezierN3D


  subroutine courbeBezierN(points_courbe,t,Px,Py)
    implicit none
    real(pr),intent(in) :: t
    real(pr),dimension(:,:),intent(in) :: points_courbe
    real(pr),intent(out) :: Px,Py

    if (t<0._pr) then
       Px = (1-t)*points_courbe(1,1) + t*points_courbe(2,1)
       Py = (1-t)*points_courbe(1,2) + t*points_courbe(2,2)
    else if (t>1._pr) then

       Px = (1-t)*points_courbe(size(points_courbe,1)-1,1) + t*points_courbe(size(points_courbe,1),1)
       Py = (1-t)*points_courbe(size(points_courbe,1)-1,2) + t*points_courbe(size(points_courbe,1),2)
    else
       write(*,*) "error problem with t-value ",t
    endif
  end subroutine courbeBezierN

  subroutine courbeBezierN3D(points_courbe,t,Px,Py,Pz)
    implicit none
    real(pr),intent(in) :: t
    real(pr),dimension(:,:),intent(in) :: points_courbe
    real(pr),intent(out) :: Px,Py,Pz

    if (t<0._pr) then
       Px = (1-t)*points_courbe(1,1) + t*points_courbe(2,1)
       Py = (1-t)*points_courbe(1,2) + t*points_courbe(2,2)
       Pz = (1-t)*points_courbe(1,3) + t*points_courbe(2,3)
    else if (t>1._pr) then
       Px = (1-t)*points_courbe(size(points_courbe,1)-1,1) + t*points_courbe(size(points_courbe,1),1)
       Py = (1-t)*points_courbe(size(points_courbe,1)-1,2) + t*points_courbe(size(points_courbe,1),2)
       Pz = (1-t)*points_courbe(size(points_courbe,1)-1,3) + t*points_courbe(size(points_courbe,1),3)

    else
       write(*,*) "error problem with t-value ",t
    endif
  end subroutine courbeBezierN3D


  subroutine poly2(pointA,pointB,pointC,t,tB,outputP)
    implicit none
    real(pr),intent(in) :: pointA,pointB,pointC,t,tB
    real(pr),intent(out) :: outputP
    real(pr) :: a

    a = (pointC-pointA)/(1._pr-tB) - (pointB-pointA)/((1._pr-tB)*tB)

    outputP = (1._pr-t)*pointA + t*pointC - t*(1._pr-t)*a
  end subroutine poly2

  subroutine splineN(points_control,t,Px,Py)
    implicit none
    real(pr),intent(in) :: t
    real(pr),dimension(:,:),intent(in) :: points_control
    real(pr),intent(out) :: Px,Py
    integer :: N,i,j,k
    real(pr),dimension(size(points_control,1)) :: u
    real(pr),dimension(size(points_control,1)) :: a,c,l,z
    real(pr),dimension(size(points_control,1)-1) :: b,d,mu,h
    real(pr),dimension(size(points_control,1)-1-1) :: alpha

    N = size(points_control,1)
    do i=1,N
       u(i) = (i-1)*1._pr/(1._pr*(N-1))
    enddo

    do i=1,N
       a(i) = points_control(i,2)
    enddo
    do i=1,N-1
       h(i) = u(i+1)-u(i)
    enddo
    do i=2,N-1
       alpha(i) = 3._pr*(a(i+1)-a(i))/h(i) - 3._pr*(a(i)-a(i-1))/(h(i-1))
    enddo
    l(1) = 1._pr
    mu(1) = 0._pr
    z(1) = 0._pr
    do i=2,N-1
       l(i) = 2._pr*(u(i+1)-u(i-1)) - h(i-1)*mu(i-1)
       mu(i) = h(i)/l(i)
       z(i) = (alpha(i) - h(i-1)*z(i-1))/(l(i))
    enddo
    l(N) = 1._pr
    z(N) = 0._pr
    c(N) = 0._pr
    do i=N-1,1,-1
       c(i) = z(i) - mu(i)*c(i+1)
       b(i) = (a(i+1)-a(i))/h(i) - (h(i)*(c(i+1)+2._pr*c(i)))/3._pr
       d(i) = (c(i+1)-c(i))/(3._pr*h(i))
    enddo
    do j=1,N
       if ((t>=u(j)).and.(t<u(j+1))) k=j
    enddo
    Py = a(k) + b(k)*(t-u(k)) + c(k)*(t-u(k))**2 + d(k)*(t-u(k))**3

    do i=1,N
       a(i) = points_control(i,1)
    enddo
    do i=1,N-1
       h(i) = u(i+1)-u(i)
    enddo
    do i=2,N-1
       alpha(i) = 3._pr*(a(i+1)-a(i))/h(i) - 3._pr*(a(i)-a(i-1))/(h(i-1))
    enddo
    l(1) = 1._pr
    mu(1) = 0._pr
    z(1) = 0._pr
    do i=2,N-1
       l(i) = 2._pr*(u(i+1)-u(i-1)) - h(i-1)*mu(i-1)
       mu(i) = h(i)/l(i)
       z(i) = (alpha(i) - h(i-1)*z(i-1))/(l(i))
    enddo
    l(N) = 1._pr
    z(N) = 0._pr
    c(N) = 0._pr
    do i=N-1,1,-1
       c(i) = z(i) - mu(i)*c(i+1)
       b(i) = (a(i+1)-a(i))/h(i) - (h(i)*(c(i+1)+2._pr*c(i)))/3._pr
       d(i) = (c(i+1)-c(i))/(3._pr*h(i))
    enddo
    do j=1,N
       if ((t>=u(j)).and.(t<u(j+1))) k=j
    enddo
    Px = a(k) + b(k)*(t-u(k)) + c(k)*(t-u(k))**2 + d(k)*(t-u(k))**3
  end subroutine splineN


  subroutine splineN3D(points_control,t,Px,Py,Pz)
    implicit none
    real(pr),intent(in) :: t
    real(pr),dimension(:,:),intent(in) :: points_control
    real(pr),intent(out) :: Px,Py,Pz
    integer :: N,i,j,k
    real(pr),dimension(size(points_control,1)) :: u
    real(pr),dimension(size(points_control,1)) :: a,c,l,z
    real(pr),dimension(size(points_control,1)-1) :: b,d,mu,h
    real(pr),dimension(size(points_control,1)-1-1) :: alpha

    N = size(points_control,1)
    do i=1,N
       u(i) = (i-1)*1._pr/(1._pr*(N-1))
    enddo

    do i=1,N
       a(i) = points_control(i,2)
    enddo
    do i=1,N-1
       h(i) = u(i+1)-u(i)
    enddo
    do i=2,N-1
       alpha(i) = 3._pr*(a(i+1)-a(i))/h(i) - 3._pr*(a(i)-a(i-1))/(h(i-1))
    enddo
    l(1) = 1._pr
    mu(1) = 0._pr
    z(1) = 0._pr
    do i=2,N-1
       l(i) = 2._pr*(u(i+1)-u(i-1)) - h(i-1)*mu(i-1)
       mu(i) = h(i)/l(i)
       z(i) = (alpha(i) - h(i-1)*z(i-1))/(l(i))
    enddo
    l(N) = 1._pr
    z(N) = 0._pr
    c(N) = 0._pr
    do i=N-1,1,-1
       c(i) = z(i) - mu(i)*c(i+1)
       b(i) = (a(i+1)-a(i))/h(i) - (h(i)*(c(i+1)+2._pr*c(i)))/3._pr
       d(i) = (c(i+1)-c(i))/(3._pr*h(i))
    enddo
    do j=1,N
       if ((t>=u(j)).and.(t<u(j+1))) k=j
    enddo
    Py = a(k) + b(k)*(t-u(k)) + c(k)*(t-u(k))**2 + d(k)*(t-u(k))**3

    do i=1,N
       a(i) = points_control(i,1)
    enddo
    do i=1,N-1
       h(i) = u(i+1)-u(i)
    enddo
    do i=2,N-1
       alpha(i) = 3._pr*(a(i+1)-a(i))/h(i) - 3._pr*(a(i)-a(i-1))/(h(i-1))
    enddo
    l(1) = 1._pr
    mu(1) = 0._pr
    z(1) = 0._pr
    do i=2,N-1
       l(i) = 2._pr*(u(i+1)-u(i-1)) - h(i-1)*mu(i-1)
       mu(i) = h(i)/l(i)
       z(i) = (alpha(i) - h(i-1)*z(i-1))/(l(i))
    enddo
    l(N) = 1._pr
    z(N) = 0._pr
    c(N) = 0._pr
    do i=N-1,1,-1
       c(i) = z(i) - mu(i)*c(i+1)
       b(i) = (a(i+1)-a(i))/h(i) - (h(i)*(c(i+1)+2._pr*c(i)))/3._pr
       d(i) = (c(i+1)-c(i))/(3._pr*h(i))
    enddo
    do j=1,N
       if ((t>=u(j)).and.(t<u(j+1))) k=j
    enddo
    Px = a(k) + b(k)*(t-u(k)) + c(k)*(t-u(k))**2 + d(k)*(t-u(k))**3

    do i=1,N
       a(i) = points_control(i,3)
    enddo
    do i=1,N-1
       h(i) = u(i+1)-u(i)
    enddo
    do i=2,N-1
       alpha(i) = 3._pr*(a(i+1)-a(i))/h(i) - 3._pr*(a(i)-a(i-1))/(h(i-1))
    enddo
    l(1) = 1._pr
    mu(1) = 0._pr
    z(1) = 0._pr
    do i=2,N-1
       l(i) = 2._pr*(u(i+1)-u(i-1)) - h(i-1)*mu(i-1)
       mu(i) = h(i)/l(i)
       z(i) = (alpha(i) - h(i-1)*z(i-1))/(l(i))
    enddo
    l(N) = 1._pr
    z(N) = 0._pr
    c(N) = 0._pr
    do i=N-1,1,-1
       c(i) = z(i) - mu(i)*c(i+1)
       b(i) = (a(i+1)-a(i))/h(i) - (h(i)*(c(i+1)+2._pr*c(i)))/3._pr
       d(i) = (c(i+1)-c(i))/(3._pr*h(i))
    enddo
    do j=1,N
       if ((t>=u(j)).and.(t<u(j+1))) k=j
    enddo
    Pz = a(k) + b(k)*(t-u(k)) + c(k)*(t-u(k))**2 + d(k)*(t-u(k))**3
  end subroutine splineN3D

  subroutine filterskel(skel,rhoSlices)
    implicit none
    real(pr),dimension(:,:),intent(in) :: rhoSlices
    integer,dimension(:,:),intent(inout) :: skel
    integer,dimension(size(skel,1),size(skel,2)) :: skeltmp,bool
    integer,dimension(3,3) :: filt,filt2,filt3,filt4,filt5,filt6,filt7,filt8,filt9,filt1,filt10,filt14
    integer :: i,j,ll

    bool = 0

    filt = -1
    filt(2,2) = 1
    filt(2,3) = 1
    filt(3,1) = 1
    filt(3,2) = 1
    filt(1,2) = 0
    filt(1,3) = 0
    filt(3,3) = 0

    filt2 = -1
    filt2(2,2) = 1
    filt2(2,3) = 1
    filt2(3,1) = 1
    filt2(3,2) = 1
    filt2(1,2) = 0
    filt2(1,3) = 0
    filt2(3,3) = 1

    filt3 = -1
    filt3(2,2) = 1
    filt3(2,3) = 1
    filt3(2,1) = 1
    filt3(3,2) = 1
    filt3(1,2) = 0
    filt3(1,3) = 0
    filt3(3,3) = 1

    filt4(1,2) = 0
    filt4(1,1) = 0
    filt4(1,3) = 0
    filt4(2,2) = 1
    filt4(2,3) = 0
    filt4(2,1) = 0
    filt4(3,2) = 1
    filt4(3,1) = 1
    filt4(3,3) = 1
    filt14(1,2) = 0
    filt14(1,1) = 0
    filt14(1,3) = 1
    filt14(2,2) = 1
    filt14(2,3) = 1
    filt14(2,1) = 1
    filt14(3,2) = 0
    filt14(3,1) = 0
    filt14(3,3) = 1

    filt5 = -1
    filt5(2,2) = 1
    filt5(2,3) = 1
    filt5(3,1) = 1
    filt5(3,2) = 1
    filt5(1,3) = 1
    filt5(1,2) = 0
    filt5(3,3) = 0

    filt6 = -1
    filt6(2,2) = 1
    filt6(2,3) = 1
    filt6(1,3) = 1
    filt6(1,2) = 0
    filt6(3,2) = 0
    filt6(3,3) = 0

    filt7 = -1
    filt7(2,2) = 1
    filt7(2,3) = 1
    filt7(1,3) = 1
    filt7(1,2) = 1

    filt8 = -1
    filt8(1,1) = 1
    filt8(2,2) = 1
    filt8(2,3) = 1
    filt8(3,2) = 1
    filt8(3,3) = 1
    filt8(1,2) = 0
    filt8(1,3) = 0

    filt9 = -1
    filt9(1,1) = 1
    filt9(2,2) = 1
    filt9(3,3) = 1
    filt9(1,3) = 1
    filt9(1,2) = 0
    filt9(2,3) = 0
    filt9(2,1) = 0
    filt9(3,1) = 0
    filt9(3,2) = 0

    filt1 = -1
    filt1(1,3) = 1
    filt1(2,3) = 1
    filt1(3,1) = 1
    filt1(3,2) = 1
    filt1(2,2) = 0
    filt1(1,2) = 0
    filt1(3,3) = 0

    filt10 = -1
    filt10(2,2) = 1
    filt10(2,3) = 1
    filt10(3,2) = 1
    filt10(1,2) = 0
    filt10(1,3) = 0
    filt10(3,3) = 0
    filt10(1,1) = 1
    filt10(2,1) = 0
    filt10(3,1) = 0


    skeltmp = skel
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt2(2,2)).and.(skel(i,j+1)==filt2(2,3)).and.(skel(i-1,j-1)==filt2(3,1)).and.(skel(i-1,j)==filt2(3,2))&
        .and.(skel(i-1,j+1)==filt2(3,3)).and.(skel(i+1,j)==filt2(1,2)).and.(skel(i+1,j+1)==filt2(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i-1,j)) then
            skel(i,j) = 0
            bool(i-1,j) = 1
          else
            skel(i-1,j) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i,j+1)<rhoSlices(i-1,j+1)) then
              skel(i,j+1) = 0
              bool(i-1,j+1) = 1
            else
              skel(i-1,j+1) = 0
              bool(i,j+1) = 1
            endif
          if (skel(i+1,j+2)==1) skel(i,j+1) = 1
          if (skel(i-2,j+2)==1) skel(i-1,j+1) = 1
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt2(2,2)).and.(skel(i,j-1)==filt2(2,3)).and.(skel(i-1,j+1)==filt2(3,1)).and.(skel(i-1,j)==filt2(3,2))&
        .and.(skel(i-1,j-1)==filt2(3,3)).and.(skel(i+1,j)==filt2(1,2)).and.(skel(i+1,j-1)==filt2(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i-1,j)) then
            skel(i,j) = 0
            bool(i-1,j) = 1
          else
            skel(i-1,j) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i,j-1)<rhoSlices(i-1,j-1)) then
              skel(i,j-1) = 0
              bool(i-1,j-1) = 1
            else
              skel(i-1,j-1) = 0
              bool(i,j-1) = 1
            endif
          if (skel(i+1,j-2)==1) skel(i,j-1) = 1
          if (skel(i-2,j-2)==1) skel(i-1,j-1) = 1
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt2(2,2)).and.(skel(i,j+1)==filt2(2,3)).and.(skel(i+1,j-1)==filt2(3,1)).and.(skel(i+1,j)==filt2(3,2))&
        .and.(skel(i+1,j+1)==filt2(3,3)).and.(skel(i-1,j)==filt2(1,2)).and.(skel(i-1,j+1)==filt2(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i+1,j)) then
            skel(i,j) = 0
            bool(i+1,j) = 1
          else
            skel(i+1,j) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i,j+1)<rhoSlices(i+1,j+1)) then
              skel(i,j+1) = 0
              bool(i+1,j+1) = 1
            else
              skel(i+1,j+1) = 0
              bool(i,j+1) = 1
            endif
          if (skel(i-1,j+2)==1) skel(i,j+1) = 1
          if (skel(i+2,j+2)==1) skel(i+1,j+1) = 1
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt2(2,2)).and.(skel(i,j-1)==filt2(2,3)).and.(skel(i+1,j+1)==filt2(3,1)).and.(skel(i+1,j)==filt2(3,2))&
        .and.(skel(i+1,j-1)==filt2(3,3)).and.(skel(i-1,j)==filt2(1,2)).and.(skel(i-1,j-1)==filt2(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i+1,j)) then
            skel(i,j) = 0
            bool(i+1,j) = 1
          else
            skel(i+1,j) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i,j-1)<rhoSlices(i+1,j-1)) then
              skel(i,j-1) = 0
              bool(i+1,j-1) = 1
            else
              skel(i+1,j-1) = 0
              bool(i,j-1) = 1
            endif
          if (skel(i-1,j-2)==1) skel(i,j-1) = 1
          if (skel(i+2,j-2)==1) skel(i+1,j-1) = 1
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2

        if ((skel(i,j)==filt2(2,2)).and.(skel(i+1,j)==filt2(2,3)).and.(skel(i-1,j-1)==filt2(3,1)).and.(skel(i,j-1)==filt2(3,2))&
        .and.(skel(i+1,j-1)==filt2(3,3)).and.(skel(i,j+1)==filt2(1,2)).and.(skel(i+1,j+1)==filt2(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j-1)) then
            skel(i,j) = 0
            bool(i,j-1) = 1
          else
            skel(i,j-1) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i+1,j)<rhoSlices(i+1,j-1)) then
              skel(i+1,j) = 0
              bool(i+1,j-1) = 1
            else
              skel(i+1,j-1) = 0
              bool(i+1,j) = 1
            endif
          if (skel(i+2,j+1)==1) skel(i+1,j) = 1
          if (skel(i+2,j-2)==1) skel(i+1,j-1) = 1
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt2(2,2)).and.(skel(i-1,j)==filt2(2,3)).and.(skel(i+1,j-1)==filt2(3,1)).and.(skel(i,j-1)==filt2(3,2))&
        .and.(skel(i-1,j-1)==filt2(3,3)).and.(skel(i,j+1)==filt2(1,2)).and.(skel(i-1,j+1)==filt2(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j-1)) then
            skel(i,j) = 0
            bool(i,j-1) = 1
          else
            skel(i,j-1) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i-1,j)<rhoSlices(i-1,j-1)) then
              skel(i-1,j) = 0
              bool(i-1,j-1) = 1
            else
              skel(i-1,j-1) = 0
              bool(i-1,j) = 1
            endif
          if (skel(i-2,j+1)==1) skel(i-1,j) = 1
          if (skel(i-2,j-2)==1) skel(i-1,j-1) = 1
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt2(2,2)).and.(skel(i+1,j)==filt2(2,3)).and.(skel(i-1,j+1)==filt2(3,1)).and.(skel(i,j-1)==filt2(3,2))&
        .and.(skel(i+1,j+1)==filt2(3,3)).and.(skel(i,j+1)==filt2(1,2)).and.(skel(i+1,j-1)==filt2(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j+1)) then
            skel(i,j) = 0
            bool(i,j+1) = 1
          else
            skel(i,j+1) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i+1,j)<rhoSlices(i+1,j+1)) then
              skel(i+1,j) = 0
              bool(i+1,j+1) = 1
            else
              skel(i+1,j+1) = 0
              bool(i+1,j) = 1
            endif
          if (skel(i+2,j-1)==1) skel(i+1,j) = 1
          if (skel(i+2,j+2)==1) skel(i+1,j+1) = 1
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt2(2,2)).and.(skel(i-1,j)==filt2(2,3)).and.(skel(i+1,j+1)==filt2(3,1)).and.(skel(i,j-1)==filt2(3,2))&
        .and.(skel(i-1,j+1)==filt2(3,3)).and.(skel(i,j+1)==filt2(1,2)).and.(skel(i-1,j-1)==filt2(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j+1)) then
            skel(i,j) = 0
            bool(i,j+1) = 1
          else
            skel(i,j+1) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i-1,j)<rhoSlices(i-1,j+1)) then
              skel(i-1,j) = 0
              bool(i-1,j+1) = 1
            else
              skel(i-1,j+1) = 0
              bool(i-1,j) = 1
            endif
          if (skel(i-2,j-1)==1) skel(i-1,j) = 1
          if (skel(i-2,j+2)==1) skel(i-1,j+1) = 1
        endif
      enddo
    enddo



    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt3(2,2)).and.(skel(i,j+1)==filt3(2,3)).and.(skel(i,j-1)==filt3(2,1)).and.(skel(i-1,j)==filt3(3,2))&
        .and.(skel(i-1,j+1)==filt3(3,3)).and.(skel(i+1,j)==filt3(1,2)).and.(skel(i+1,j+1)==filt3(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i-1,j)) then
            skel(i,j) = 0
            bool(i-1,j) = 1
          else
            skel(i-1,j) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i,j+1)<rhoSlices(i-1,j+1)) then
              skel(i,j+1) = 0
              bool(i-1,j+1) = 1
            else
              skel(i-1,j+1) = 0
              bool(i,j+1) = 1
            endif
          if (skel(i+1,j+2)==1) skel(i,j+1) = 1
          if (skel(i-2,j+2)==1) skel(i-1,j+1) = 1
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt3(2,2)).and.(skel(i,j-1)==filt3(2,3)).and.(skel(i,j+1)==filt3(2,1)).and.(skel(i-1,j)==filt3(3,2))&
        .and.(skel(i-1,j-1)==filt3(3,3)).and.(skel(i+1,j)==filt3(1,2)).and.(skel(i+1,j-1)==filt3(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i-1,j)) then
            skel(i,j) = 0
            bool(i-1,j) = 1
          else
            skel(i-1,j) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i,j-1)<rhoSlices(i-1,j-1)) then
              skel(i,j-1) = 0
              bool(i-1,j-1) = 1
            else
              skel(i-1,j-1) = 0
              bool(i,j-1) = 1
            endif
          if (skel(i+1,j-2)==1) skel(i,j-1) = 1
          if (skel(i-2,j-2)==1) skel(i-1,j-1) = 1
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt3(2,2)).and.(skel(i,j+1)==filt3(2,3)).and.(skel(i,j-1)==filt3(2,1)).and.(skel(i+1,j)==filt3(3,2))&
        .and.(skel(i+1,j+1)==filt3(3,3)).and.(skel(i-1,j)==filt3(1,2)).and.(skel(i-1,j+1)==filt3(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i+1,j)) then
            skel(i,j) = 0
            bool(i+1,j) = 1
          else
            skel(i+1,j) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i,j+1)<rhoSlices(i+1,j+1)) then
              skel(i,j+1) = 0
              bool(i+1,j+1) = 1
            else
              skel(i+1,j+1) = 0
              bool(i,j+1) = 1
            endif
          if (skel(i-1,j+2)==1) skel(i,j+1) = 1
          if (skel(i+2,j+2)==1) skel(i+1,j+1) = 1
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt3(2,2)).and.(skel(i,j-1)==filt3(2,3)).and.(skel(i,j+1)==filt3(2,1)).and.(skel(i+1,j)==filt3(3,2))&
        .and.(skel(i+1,j-1)==filt3(3,3)).and.(skel(i-1,j)==filt3(1,2)).and.(skel(i-1,j-1)==filt3(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i+1,j)) then
            skel(i,j) = 0
            bool(i+1,j) = 1
          else
            skel(i+1,j) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i,j-1)<rhoSlices(i+1,j-1)) then
              skel(i,j-1) = 0
              bool(i+1,j-1) = 1
            else
              skel(i+1,j-1) = 0
              bool(i,j-1) = 1
            endif
          if (skel(i-1,j-2)==1) skel(i,j-1) = 1
          if (skel(i+2,j-2)==1) skel(i+1,j-1) = 1
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2

        if ((skel(i,j)==filt3(2,2)).and.(skel(i+1,j)==filt3(2,3)).and.(skel(i-1,j)==filt3(2,1)).and.(skel(i,j-1)==filt3(3,2))&
        .and.(skel(i+1,j-1)==filt3(3,3)).and.(skel(i,j+1)==filt3(1,2)).and.(skel(i+1,j+1)==filt3(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j-1)) then
            skel(i,j) = 0
            bool(i,j-1) = 1
          else
            skel(i,j-1) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i+1,j)<rhoSlices(i+1,j-1)) then
              skel(i+1,j) = 0
              bool(i+1,j-1) = 1
            else
              skel(i+1,j-1) = 0
              bool(i+1,j) = 1
            endif
          if (skel(i+2,j+1)==1) skel(i+1,j) = 1
          if (skel(i+2,j-2)==1) skel(i+1,j-1) = 1
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt3(2,2)).and.(skel(i-1,j)==filt3(2,3)).and.(skel(i+1,j)==filt3(2,1)).and.(skel(i,j-1)==filt3(3,2))&
        .and.(skel(i-1,j-1)==filt3(3,3)).and.(skel(i,j+1)==filt3(1,2)).and.(skel(i-1,j+1)==filt3(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j-1)) then
            skel(i,j) = 0
            bool(i,j-1) = 1
          else
            skel(i,j-1) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i-1,j)<rhoSlices(i-1,j-1)) then
              skel(i-1,j) = 0
              bool(i-1,j-1) = 1
            else
              skel(i-1,j-1) = 0
              bool(i-1,j) = 1
            endif
          if (skel(i-2,j+1)==1) skel(i-1,j) = 1
          if (skel(i-2,j-2)==1) skel(i-1,j-1) = 1
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt3(2,2)).and.(skel(i+1,j)==filt3(2,3)).and.(skel(i-1,j)==filt3(2,1)).and.(skel(i,j-1)==filt3(3,2))&
        .and.(skel(i+1,j+1)==filt3(3,3)).and.(skel(i,j+1)==filt3(1,2)).and.(skel(i+1,j-1)==filt3(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j+1)) then
            skel(i,j) = 0
            bool(i,j+1) = 1
          else
            skel(i,j+1) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i+1,j)<rhoSlices(i+1,j+1)) then
              skel(i+1,j) = 0
              bool(i+1,j+1) = 1
            else
              skel(i+1,j+1) = 0
              bool(i+1,j) = 1
            endif
          if (skel(i+2,j-1)==1) skel(i+1,j) = 1
          if (skel(i+2,j+2)==1) skel(i+1,j+1) = 1
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt3(2,2)).and.(skel(i-1,j)==filt3(2,3)).and.(skel(i+1,j)==filt3(2,1)).and.(skel(i,j-1)==filt3(3,2))&
        .and.(skel(i-1,j+1)==filt3(3,3)).and.(skel(i,j+1)==filt3(1,2)).and.(skel(i-1,j-1)==filt3(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j+1)) then
            skel(i,j) = 0
            bool(i,j+1) = 1
          else
            skel(i,j+1) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i-1,j)<rhoSlices(i-1,j+1)) then
              skel(i-1,j) = 0
              bool(i-1,j+1) = 1
            else
              skel(i-1,j+1) = 0
              bool(i-1,j) = 1
            endif
          if (skel(i-2,j-1)==1) skel(i-1,j) = 1
          if (skel(i-2,j+2)==1) skel(i-1,j+1) = 1
        endif
      enddo
    enddo


    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt8(2,2)).and.(skel(i,j+1)==filt8(2,3)).and.(skel(i+1,j-1)==filt8(1,1)).and.(skel(i-1,j)==filt8(3,2))&
        .and.(skel(i-1,j+1)==filt8(3,3)).and.(skel(i+1,j)==filt8(1,2)).and.(skel(i+1,j+1)==filt8(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i-1,j)) then
            skel(i,j) = 0
            bool(i-1,j) = 1
          else
            skel(i-1,j) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i,j+1)<rhoSlices(i-1,j+1)) then
              skel(i,j+1) = 0
              bool(i-1,j+1) = 1
            else
              skel(i-1,j+1) = 0
              bool(i,j+1) = 1
            endif
          if (skel(i+1,j+2)==1) skel(i,j+1) = 1
          if (skel(i-2,j+2)==1) skel(i-1,j+1) = 1

        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt8(2,2)).and.(skel(i,j-1)==filt8(2,3)).and.(skel(i+1,j+1)==filt8(1,1)).and.(skel(i-1,j)==filt8(3,2))&
        .and.(skel(i-1,j-1)==filt8(3,3)).and.(skel(i+1,j)==filt8(1,2)).and.(skel(i+1,j-1)==filt8(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i-1,j)) then
            skel(i,j) = 0
            bool(i-1,j) = 1
          else
            skel(i-1,j) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i,j-1)<rhoSlices(i-1,j-1)) then
              skel(i,j-1) = 0
              bool(i-1,j-1) = 1
            else
              skel(i-1,j-1) = 0
              bool(i,j-1) = 1
            endif
          if (skel(i+1,j-2)==1) skel(i,j-1) = 1
          if (skel(i-2,j-2)==1) skel(i-1,j-1) = 1

        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt8(2,2)).and.(skel(i,j+1)==filt8(2,3)).and.(skel(i-1,j-1)==filt8(1,1)).and.(skel(i+1,j)==filt8(3,2))&
        .and.(skel(i+1,j+1)==filt8(3,3)).and.(skel(i-1,j)==filt8(1,2)).and.(skel(i-1,j+1)==filt8(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i+1,j)) then
            skel(i,j) = 0
            bool(i+1,j) = 1
          else
            skel(i+1,j) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i,j+1)<rhoSlices(i+1,j+1)) then
              skel(i,j+1) = 0
              bool(i+1,j+1) = 1
            else
              skel(i+1,j+1) = 0
              bool(i,j+1) = 1
            endif
          if (skel(i-1,j+2)==1) skel(i,j+1) = 1
          if (skel(i+2,j+2)==1) skel(i+1,j+1) = 1

        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt8(2,2)).and.(skel(i,j-1)==filt8(2,3)).and.(skel(i-1,j+1)==filt8(1,1)).and.(skel(i+1,j)==filt8(3,2))&
        .and.(skel(i+1,j-1)==filt8(3,3)).and.(skel(i-1,j)==filt8(1,2)).and.(skel(i-1,j-1)==filt8(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i+1,j)) then
            skel(i,j) = 0
            bool(i+1,j) = 1
          else
            skel(i+1,j) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i,j-1)<rhoSlices(i+1,j-1)) then
              skel(i,j-1) = 0
              bool(i+1,j-1) = 1
            else
              skel(i+1,j-1) = 0
              bool(i,j-1) = 1
            endif
          if (skel(i-1,j-2)==1) skel(i,j-1) = 1
          if (skel(i+2,j-2)==1) skel(i+1,j-1) = 1

        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2

        if ((skel(i,j)==filt8(2,2)).and.(skel(i+1,j)==filt8(2,3)).and.(skel(i-1,j+1)==filt8(1,1)).and.(skel(i,j-1)==filt8(3,2))&
        .and.(skel(i+1,j-1)==filt8(3,3)).and.(skel(i,j+1)==filt8(1,2)).and.(skel(i+1,j+1)==filt8(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j-1)) then
            skel(i,j) = 0
            bool(i,j-1) = 1
          else
            skel(i,j-1) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i+1,j)<rhoSlices(i+1,j-1)) then
              skel(i+1,j) = 0
              bool(i+1,j-1) = 1
            else
              skel(i+1,j-1) = 0
              bool(i+1,j) = 1
            endif
          if (skel(i+2,j+1)==1) skel(i+1,j) = 1
          if (skel(i+2,j-2)==1) skel(i+1,j-1) = 1

        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt8(2,2)).and.(skel(i-1,j)==filt8(2,3)).and.(skel(i+1,j+1)==filt8(1,1)).and.(skel(i,j-1)==filt8(3,2))&
        .and.(skel(i-1,j-1)==filt8(3,3)).and.(skel(i,j+1)==filt8(1,2)).and.(skel(i-1,j+1)==filt8(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j-1)) then
            skel(i,j) = 0
            bool(i,j-1) = 1
          else
            skel(i,j-1) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i-1,j)<rhoSlices(i-1,j-1)) then
              skel(i-1,j) = 0
              bool(i-1,j-1) = 1
            else
              skel(i-1,j-1) = 0
              bool(i-1,j) = 1
            endif
          if (skel(i-2,j+1)==1) skel(i-1,j) = 1
          if (skel(i-2,j-2)==1) skel(i-1,j-1) = 1

        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt8(2,2)).and.(skel(i+1,j)==filt8(2,3)).and.(skel(i-1,j-1)==filt8(1,1)).and.(skel(i,j+1)==filt8(3,2))&
        .and.(skel(i+1,j+1)==filt8(3,3)).and.(skel(i,j-1)==filt8(1,2)).and.(skel(i+1,j-1)==filt8(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j+1)) then
            skel(i,j) = 0
            bool(i,j+1) = 1
          else
            skel(i,j+1) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i+1,j)<rhoSlices(i+1,j+1)) then
              skel(i+1,j) = 0
              bool(i+1,j+1) = 1
            else
              skel(i+1,j+1) = 0
              bool(i+1,j) = 1
            endif
          if (skel(i+2,j-1)==1) skel(i+1,j) = 1
          if (skel(i+2,j+2)==1) skel(i+1,j+1) = 1

        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt8(2,2)).and.(skel(i-1,j)==filt8(2,3)).and.(skel(i+1,j-1)==filt8(1,1)).and.(skel(i,j+1)==filt8(3,2))&
        .and.(skel(i-1,j+1)==filt8(3,3)).and.(skel(i,j-1)==filt8(1,2)).and.(skel(i-1,j-1)==filt8(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j+1)) then
            skel(i,j) = 0
            bool(i,j+1) = 1
          else
            skel(i,j+1) = 0
            bool(i,j) = 1
          endif
            if (rhoSlices(i-1,j)<rhoSlices(i-1,j+1)) then
              skel(i-1,j) = 0
              bool(i-1,j+1) = 1
            else
              skel(i-1,j+1) = 0
              bool(i-1,j) = 1
            endif
          if (skel(i-2,j-1)==1) skel(i-1,j) = 1
          if (skel(i-2,j+2)==1) skel(i-1,j+1) = 1

        endif
      enddo
    enddo


    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt14(2,2)).and.(skel(i,j+1)==filt14(2,3)).and.(skel(i,j-1)==filt14(2,1)).and.(skel(i+1,j)==filt14(1,2))&
        .and.(skel(i-1,j+1)==filt14(3,3)).and.(skel(i-1,j-1)==filt14(3,1)).and.(skel(i-1,j)==filt14(3,2))&
        .and.(skel(i+1,j-1)==filt14(1,1)).and.(skel(i+1,j+1)==filt14(1,3))) then
            skel(i,j+1) = 0
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt14(2,2)).and.(skel(i,j-1)==filt14(2,3)).and.(skel(i,j+1)==filt14(2,1)).and.(skel(i+1,j)==filt14(1,2))&
        .and.(skel(i-1,j-1)==filt14(3,3)).and.(skel(i-1,j+1)==filt14(3,1)).and.(skel(i-1,j)==filt14(3,2))&
        .and.(skel(i+1,j+1)==filt14(1,1)).and.(skel(i+1,j-1)==filt14(1,3))) then

            skel(i,j-1) = 0
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt14(2,2)).and.(skel(i,j+1)==filt14(2,3)).and.(skel(i,j-1)==filt14(2,1)).and.(skel(i-1,j)==filt14(1,2))&
        .and.(skel(i+1,j+1)==filt14(3,3)).and.(skel(i+1,j-1)==filt14(3,1)).and.(skel(i+1,j)==filt14(3,2))&
        .and.(skel(i-1,j-1)==filt14(1,1)).and.(skel(i-1,j+1)==filt14(1,3))) then

            skel(i,j+1) = 0
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt14(2,2)).and.(skel(i,j-1)==filt14(2,3)).and.(skel(i,j+1)==filt14(2,1)).and.(skel(i-1,j)==filt14(1,2))&
        .and.(skel(i+1,j-1)==filt14(3,3)).and.(skel(i+1,j+1)==filt14(3,1)).and.(skel(i+1,j)==filt14(3,2))&
        .and.(skel(i-1,j+1)==filt14(1,1)).and.(skel(i-1,j-1)==filt14(1,3))) then
            skel(i,j-1) = 0
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1

        if ((skel(i,j)==filt14(2,2)).and.(skel(i+1,j)==filt14(2,3)).and.(skel(i-1,j)==filt14(2,1)).and.(skel(i,j+1)==filt14(1,2))&
        .and.(skel(i+1,j-1)==filt14(3,3)).and.(skel(i-1,j-1)==filt14(3,1)).and.(skel(i,j-1)==filt14(3,2))&
        .and.(skel(i-1,j+1)==filt14(1,1)).and.(skel(i+1,j+1)==filt14(1,3))) then
            skel(i+1,j) = 0
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt14(2,2)).and.(skel(i+1,j)==filt14(2,3)).and.(skel(i-1,j)==filt14(2,1)).and.(skel(i,j+1)==filt14(1,2))&
        .and.(skel(i-1,j-1)==filt14(3,3)).and.(skel(i+1,j-1)==filt14(3,1)).and.(skel(i,j-1)==filt14(3,2))&
        .and.(skel(i+1,j+1)==filt14(1,1)).and.(skel(i-1,j+1)==filt14(1,3))) then
            skel(i-1,j) = 0
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt14(2,2)).and.(skel(i+1,j)==filt14(2,3)).and.(skel(i-1,j)==filt14(2,1)).and.(skel(i,j-1)==filt14(1,2))&
        .and.(skel(i+1,j+1)==filt14(3,3)).and.(skel(i-1,j+1)==filt14(3,1)).and.(skel(i,j+1)==filt14(3,2))&
        .and.(skel(i-1,j-1)==filt14(1,1)).and.(skel(i+1,j-1)==filt14(1,3))) then
            skel(i+1,j) = 0
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt14(2,2)).and.(skel(i+1,j)==filt14(2,3)).and.(skel(i-1,j)==filt14(2,1)).and.(skel(i,j-1)==filt14(1,2))&
        .and.(skel(i-1,j+1)==filt14(3,3)).and.(skel(i+1,j+1)==filt14(3,1)).and.(skel(i,j+1)==filt14(3,2))&
        .and.(skel(i+1,j-1)==filt14(1,1)).and.(skel(i-1,j-1)==filt14(1,3))) then
            skel(i-1,j) = 0
        endif
      enddo
    enddo


    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt4(2,2)).and.(skel(i,j+1)==filt4(2,3)).and.(skel(i,j-1)==filt4(2,1)).and.(skel(i+1,j)==filt4(1,2))&
        .and.(skel(i-1,j+1)==filt4(3,3)).and.(skel(i-1,j-1)==filt4(3,1)).and.(skel(i-1,j)==filt4(3,2))&
        .and.(skel(i+1,j-1)==filt4(1,1)).and.(skel(i+1,j+1)==filt4(1,3))) then
            skel(i,j) = 0
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt4(2,2)).and.(skel(i,j-1)==filt4(2,3)).and.(skel(i,j+1)==filt4(2,1)).and.(skel(i+1,j)==filt4(1,2))&
        .and.(skel(i-1,j-1)==filt4(3,3)).and.(skel(i-1,j+1)==filt4(3,1)).and.(skel(i-1,j)==filt4(3,2))&
        .and.(skel(i+1,j+1)==filt4(1,1)).and.(skel(i+1,j-1)==filt4(1,3))) then
            skel(i,j) = 0
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt4(2,2)).and.(skel(i,j+1)==filt4(2,3)).and.(skel(i,j-1)==filt4(2,1)).and.(skel(i-1,j)==filt4(1,2))&
        .and.(skel(i+1,j+1)==filt4(3,3)).and.(skel(i+1,j-1)==filt4(3,1)).and.(skel(i+1,j)==filt4(3,2))&
        .and.(skel(i-1,j-1)==filt4(1,1)).and.(skel(i-1,j+1)==filt4(1,3))) then
            skel(i,j) = 0
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt4(2,2)).and.(skel(i,j-1)==filt4(2,3)).and.(skel(i,j+1)==filt4(2,1)).and.(skel(i-1,j)==filt4(1,2))&
        .and.(skel(i+1,j-1)==filt4(3,3)).and.(skel(i+1,j+1)==filt4(3,1)).and.(skel(i+1,j)==filt4(3,2))&
        .and.(skel(i-1,j+1)==filt4(1,1)).and.(skel(i-1,j-1)==filt4(1,3))) then
            skel(i,j) = 0
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1

        if ((skel(i,j)==filt4(2,2)).and.(skel(i+1,j)==filt4(2,3)).and.(skel(i-1,j)==filt4(2,1)).and.(skel(i,j+1)==filt4(1,2))&
        .and.(skel(i+1,j-1)==filt4(3,3)).and.(skel(i-1,j-1)==filt4(3,1)).and.(skel(i,j-1)==filt4(3,2))&
        .and.(skel(i-1,j+1)==filt4(1,1)).and.(skel(i+1,j+1)==filt4(1,3))) then
            skel(i,j) = 0
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt4(2,2)).and.(skel(i+1,j)==filt4(2,3)).and.(skel(i-1,j)==filt4(2,1)).and.(skel(i,j+1)==filt4(1,2))&
        .and.(skel(i-1,j-1)==filt4(3,3)).and.(skel(i+1,j-1)==filt4(3,1)).and.(skel(i,j-1)==filt4(3,2))&
        .and.(skel(i+1,j+1)==filt4(1,1)).and.(skel(i-1,j+1)==filt4(1,3))) then
            skel(i,j) = 0
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt4(2,2)).and.(skel(i+1,j)==filt4(2,3)).and.(skel(i-1,j)==filt4(2,1)).and.(skel(i,j-1)==filt4(1,2))&
        .and.(skel(i+1,j+1)==filt4(3,3)).and.(skel(i-1,j+1)==filt4(3,1)).and.(skel(i,j+1)==filt4(3,2))&
        .and.(skel(i-1,j-1)==filt4(1,1)).and.(skel(i+1,j-1)==filt4(1,3))) then
            skel(i,j) = 0
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt4(2,2)).and.(skel(i+1,j)==filt4(2,3)).and.(skel(i-1,j)==filt4(2,1)).and.(skel(i,j-1)==filt4(1,2))&
        .and.(skel(i-1,j+1)==filt4(3,3)).and.(skel(i+1,j+1)==filt4(3,1)).and.(skel(i,j+1)==filt4(3,2))&
        .and.(skel(i+1,j-1)==filt4(1,1)).and.(skel(i-1,j-1)==filt4(1,3))) then
            skel(i,j) = 0
        endif
      enddo
    enddo


    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt5(2,2)).and.(skel(i,j+1)==filt5(2,3)).and.(skel(i-1,j-1)==filt5(3,1)).and.(skel(i-1,j)==filt5(3,2))&
        .and.(skel(i-1,j+1)==filt5(3,3)).and.(skel(i+1,j)==filt5(1,2)).and.(skel(i+1,j+1)==filt5(1,3))) then
          if (rhoSlices(i-1,j)<rhoSlices(i,j+1)) then
             skel(i-1,j) = 0
             bool(i,j+1) = 1
          else
             skel(i,j+1) = 0
             bool(i-1,j) = 1
          endif
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt5(2,2)).and.(skel(i,j-1)==filt5(2,3)).and.(skel(i-1,j+1)==filt5(3,1)).and.(skel(i-1,j)==filt5(3,2))&
        .and.(skel(i-1,j-1)==filt5(3,3)).and.(skel(i+1,j)==filt5(1,2)).and.(skel(i+1,j-1)==filt5(1,3))) then

          if (rhoSlices(i-1,j)<rhoSlices(i,j-1)) then
             skel(i-1,j) = 0
             bool(i,j-1) = 1
          else
             skel(i,j-1) = 0
             bool(i-1,j) = 1
          endif
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt5(2,2)).and.(skel(i,j+1)==filt5(2,3)).and.(skel(i+1,j-1)==filt5(3,1)).and.(skel(i+1,j)==filt5(3,2))&
        .and.(skel(i+1,j+1)==filt5(3,3)).and.(skel(i-1,j)==filt5(1,2)).and.(skel(i-1,j+1)==filt5(1,3))) then

          if (rhoSlices(i+1,j)<rhoSlices(i,j+1)) then
             skel(i+1,j) = 0
             bool(i,j+1) = 1
          else
             skel(i,j+1) = 0
             bool(i+1,j) = 1
          endif
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt5(2,2)).and.(skel(i,j-1)==filt5(2,3)).and.(skel(i+1,j+1)==filt5(3,1)).and.(skel(i+1,j)==filt5(3,2))&
        .and.(skel(i+1,j-1)==filt5(3,3)).and.(skel(i-1,j)==filt5(1,2)).and.(skel(i-1,j-1)==filt5(1,3))) then

          if (rhoSlices(i+1,j)<rhoSlices(i,j-1)) then
             skel(i+1,j) = 0
             bool(i,j-1) = 1
          else
             skel(i,j-1) = 0
             bool(i+1,j) = 1
          endif
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1

        if ((skel(i,j)==filt5(2,2)).and.(skel(i+1,j)==filt5(2,3)).and.(skel(i-1,j-1)==filt5(3,1)).and.(skel(i,j-1)==filt5(3,2))&
        .and.(skel(i+1,j-1)==filt5(3,3)).and.(skel(i,j+1)==filt5(1,2)).and.(skel(i+1,j+1)==filt5(1,3))) then

          if (rhoSlices(i,j-1)<rhoSlices(i+1,j)) then
             skel(i,j-1) = 0
             bool(i+1,j) = 1
          else
             skel(i+1,j) = 0
             bool(i,j-1) = 1
          endif
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt5(2,2)).and.(skel(i-1,j)==filt5(2,3)).and.(skel(i+1,j-1)==filt5(3,1)).and.(skel(i,j-1)==filt5(3,2))&
        .and.(skel(i-1,j-1)==filt5(3,3)).and.(skel(i,j+1)==filt5(1,2)).and.(skel(i-1,j+1)==filt5(1,3))) then

          if (rhoSlices(i,j-1)<rhoSlices(i-1,j)) then
             skel(i,j-1) = 0
             bool(i-1,j) = 1
          else
             skel(i-1,j) = 0
             bool(i,j-1) = 1
          endif
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt5(2,2)).and.(skel(i+1,j)==filt5(2,3)).and.(skel(i-1,j+1)==filt5(3,1)).and.(skel(i,j+1)==filt5(3,2))&
        .and.(skel(i+1,j+1)==filt5(3,3)).and.(skel(i,j-1)==filt5(1,2)).and.(skel(i+1,j-1)==filt5(1,3))) then

          if (rhoSlices(i,j+1)<rhoSlices(i+1,j)) then
             skel(i,j+1) = 0
             bool(i+1,j) = 1
          else
             skel(i+1,j) = 0
             bool(i,j+1) = 1
          endif
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt5(2,2)).and.(skel(i-1,j)==filt5(2,3)).and.(skel(i+1,j+1)==filt5(3,1)).and.(skel(i,j+1)==filt5(3,2))&
        .and.(skel(i-1,j+1)==filt5(3,3)).and.(skel(i,j-1)==filt5(1,2)).and.(skel(i-1,j-1)==filt5(1,3))) then

          if (rhoSlices(i,j+1)<rhoSlices(i-1,j)) then
             skel(i,j+1) = 0
             bool(i-1,j) = 1
          else
             skel(i-1,j) = 0
             bool(i,j+1) = 1
          endif
        endif
      enddo
    enddo


    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt(2,2)).and.(skel(i,j+1)==filt(2,3)).and.(skel(i-1,j-1)==filt(3,1)).and.(skel(i-1,j)==filt(3,2))&
        .and.(skel(i-1,j+1)==filt(3,3)).and.(skel(i+1,j)==filt(1,2)).and.(skel(i+1,j+1)==filt(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i-1,j)) then
            skel(i,j) = 0
            bool(i-1,j) = 1
          else
            skel(i-1,j) = 0
            bool(i,j) = 1
          endif
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt(2,2)).and.(skel(i,j-1)==filt(2,3)).and.(skel(i-1,j+1)==filt(3,1)).and.(skel(i-1,j)==filt(3,2))&
        .and.(skel(i-1,j-1)==filt(3,3)).and.(skel(i+1,j)==filt(1,2)).and.(skel(i+1,j-1)==filt(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i-1,j)) then
            skel(i,j) = 0
            bool(i-1,j) = 1
          else
            skel(i-1,j) = 0
            bool(i,j) = 1
          endif
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt(2,2)).and.(skel(i,j+1)==filt(2,3)).and.(skel(i+1,j-1)==filt(3,1)).and.(skel(i+1,j)==filt(3,2))&
        .and.(skel(i+1,j+1)==filt(3,3)).and.(skel(i-1,j)==filt(1,2)).and.(skel(i-1,j+1)==filt(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i+1,j)) then
            skel(i,j) = 0
            bool(i+1,j) = 1
          else
            skel(i+1,j) = 0
            bool(i,j) = 1
          endif
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt(2,2)).and.(skel(i,j-1)==filt(2,3)).and.(skel(i+1,j+1)==filt(3,1)).and.(skel(i+1,j)==filt(3,2))&
        .and.(skel(i+1,j-1)==filt(3,3)).and.(skel(i-1,j)==filt(1,2)).and.(skel(i-1,j-1)==filt(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i+1,j)) then
            skel(i,j) = 0
            bool(i+1,j) = 1
          else
            skel(i+1,j) = 0
            bool(i,j) = 1
          endif
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1

        if ((skel(i,j)==filt(2,2)).and.(skel(i+1,j)==filt(2,3)).and.(skel(i-1,j-1)==filt(3,1)).and.(skel(i,j-1)==filt(3,2))&
        .and.(skel(i+1,j-1)==filt(3,3)).and.(skel(i,j+1)==filt(1,2)).and.(skel(i+1,j+1)==filt(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j-1)) then
            skel(i,j) = 0
            bool(i,j-1) = 1
          else
            skel(i,j-1) = 0
            bool(i,j) = 1
          endif
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt(2,2)).and.(skel(i-1,j)==filt(2,3)).and.(skel(i+1,j-1)==filt(3,1)).and.(skel(i,j-1)==filt(3,2))&
        .and.(skel(i-1,j-1)==filt(3,3)).and.(skel(i,j+1)==filt(1,2)).and.(skel(i-1,j+1)==filt(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j-1)) then
            skel(i,j) = 0
            bool(i,j-1) = 1
          else
            skel(i,j-1) = 0
            bool(i,j) = 1
          endif
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt(2,2)).and.(skel(i+1,j)==filt(2,3)).and.(skel(i-1,j+1)==filt(3,1)).and.(skel(i,j-1)==filt(3,2))&
        .and.(skel(i+1,j+1)==filt(3,3)).and.(skel(i,j+1)==filt(1,2)).and.(skel(i+1,j-1)==filt(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j+1)) then
            skel(i,j) = 0
            bool(i,j+1) = 1
          else
            skel(i,j+1) = 0
            bool(i,j) = 1
          endif
        endif
      enddo
    enddo
    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1
        if ((skel(i,j)==filt(2,2)).and.(skel(i-1,j)==filt(2,3)).and.(skel(i+1,j+1)==filt(3,1)).and.(skel(i,j-1)==filt(3,2))&
        .and.(skel(i-1,j+1)==filt(3,3)).and.(skel(i,j+1)==filt(1,2)).and.(skel(i-1,j-1)==filt(1,3))) then
          if (rhoSlices(i,j)<rhoSlices(i,j+1)) then
            skel(i,j) = 0
            bool(i,j+1) = 1
          else
            skel(i,j+1) = 0
            bool(i,j) = 1
          endif
        endif
      enddo
    enddo


    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt1(2,2)).and.(skel(i,j+1)==filt1(2,3)).and.(skel(i-1,j-1)==filt1(3,1)).and.(skel(i-1,j)==filt1(3,2))&
        .and.((skel(i+1,j+2)+skel(i+2,j+2)+skel(i+2,j+1)+skel(i,j+1)+skel(i+1,j)+skel(i,j)+skel(i,j+2)+skel(i+2,j)==1).or.&
        (skel(i-1,j)+skel(i,j)+skel(i,j-1)+skel(i-2,j-1)+skel(i-1,j-2)+skel(i-2,j-2)+skel(i-2,j)+skel(i,j-2)==1))&
        .and.(skel(i-1,j+1)==filt1(3,3)).and.(skel(i+1,j)==filt1(1,2)).and.(skel(i+1,j+1)==filt1(1,3))) then
          if (skel(i+1,j+2)+skel(i+2,j+2)+skel(i+2,j+1)+skel(i,j+1)+skel(i+1,j)+skel(i,j)+skel(i,j+2)+skel(i+2,j)==1)&
          skel(i+1,j+1) = 0
          if (skel(i-1,j)+skel(i,j)+skel(i,j-1)+skel(i-2,j-1)+skel(i-1,j-2)+skel(i-2,j-2)+skel(i-2,j)+skel(i,j-2)==1)&
          skel(i-1,j-1) = 0

        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt1(2,2)).and.(skel(i,j-1)==filt1(2,3)).and.(skel(i-1,j+1)==filt1(3,1)).and.(skel(i-1,j)==filt1(3,2))&
        .and.((skel(i+1,j-2)+skel(i+2,j-2)+skel(i+2,j-1)+skel(i,j-1)+skel(i+1,j)+skel(i,j)+skel(i,j-2)+skel(i+2,j)==1).or.&
        (skel(i-1,j)+skel(i,j)+skel(i,j+1)+skel(i-2,j+1)+skel(i-1,j+2)+skel(i-2,j+2)+skel(i-2,j)+skel(i,j+2)==1))&
        .and.(skel(i-1,j-1)==filt1(3,3)).and.(skel(i+1,j)==filt1(1,2)).and.(skel(i+1,j-1)==filt1(1,3))) then
          if (skel(i+1,j-2)+skel(i+2,j-2)+skel(i+2,j-1)+skel(i,j-1)+skel(i+1,j)+skel(i,j)+skel(i,j-2)+skel(i+2,j)==1)&
          skel (i+1,j-1) = 0
          if (skel(i-1,j)+skel(i,j)+skel(i,j+1)+skel(i-2,j+1)+skel(i-1,j+2)+skel(i-2,j+2)+skel(i-2,j)+skel(i,j+2)==1)&
          skel(i-1,j+1) = 0

        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt1(2,2)).and.(skel(i,j+1)==filt1(2,3)).and.(skel(i+1,j-1)==filt1(3,1)).and.(skel(i+1,j)==filt1(3,2))&
        .and.((skel(i-1,j+2)+skel(i-2,j+2)+skel(i-2,j+1)+skel(i,j+1)+skel(i-1,j)+skel(i,j)+skel(i,j+2)+skel(i-2,j)==1).or.&
        (skel(i+1,j)+skel(i,j)+skel(i,j-1)+skel(i+2,j-1)+skel(i+1,j-2)+skel(i+2,j-2)+skel(i+2,j)+skel(i,j-2)==1))&
        .and.(skel(i+1,j+1)==filt1(3,3)).and.(skel(i-1,j)==filt1(1,2)).and.(skel(i-1,j+1)==filt1(1,3))) then
          if (skel(i-1,j+2)+skel(i-2,j+2)+skel(i-2,j+1)+skel(i,j+1)+skel(i-1,j)+skel(i,j)+skel(i,j+2)+skel(i-2,j)==1)&
          skel(i-1,j+1) = 0
          if (skel(i+1,j)+skel(i,j)+skel(i,j-1)+skel(i+2,j-1)+skel(i+1,j-2)+skel(i+2,j-2)+skel(i+2,j)+skel(i,j-2)==1)&
          skel(i+1,j-1) = 0

        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt1(2,2)).and.(skel(i,j-1)==filt1(2,3)).and.(skel(i+1,j+1)==filt1(3,1)).and.(skel(i+1,j)==filt1(3,2))&
        .and.((skel(i-1,j-2)+skel(i-2,j-2)+skel(i-2,j-1)+skel(i,j-1)+skel(i-1,j)+skel(i,j)+skel(i,j-2)+skel(i-2,j)==1).or.&
        (skel(i+1,j)+skel(i,j)+skel(i,j+1)+skel(i+2,j+1)+skel(i+1,j+2)+skel(i+2,j+2)+skel(i+2,j)+skel(i,j+2)==1))&
        .and.(skel(i+1,j-1)==filt1(3,3)).and.(skel(i-1,j)==filt1(1,2)).and.(skel(i-1,j-1)==filt1(1,3))) then
          if (skel(i-1,j-2)+skel(i-2,j-2)+skel(i-2,j-1)+skel(i,j-1)+skel(i-1,j)+skel(i,j)+skel(i,j-2)+skel(i-2,j)==1)&
          skel(i-1,j-1) = 0
          if (skel(i+1,j)+skel(i,j)+skel(i,j+1)+skel(i+2,j+1)+skel(i+1,j+2)+skel(i+2,j+2)+skel(i+2,j)+skel(i,j+2)==1)&
          skel(i+1,j+1) = 0

        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2

        if ((skel(i,j)==filt1(2,2)).and.(skel(i+1,j)==filt1(2,3)).and.(skel(i-1,j-1)==filt1(3,1)).and.(skel(i,j-1)==filt1(3,2))&
        .and.((skel(i+2,j+1)+skel(i+2,j+2)+skel(i+1,j+2)+skel(i+1,j)+skel(i,j+1)+skel(i,j)+skel(i+2,j)+skel(i,j+2)==1).or.&
        (skel(i,j-1)+skel(i,j)+skel(i-1,j)+skel(i-1,j-2)+skel(i-2,j-1)+skel(i-2,j-2)+skel(i,j-2)+skel(i-2,j)==1))&
        .and.(skel(i+1,j-1)==filt1(3,3)).and.(skel(i,j+1)==filt1(1,2)).and.(skel(i+1,j+1)==filt1(1,3))) then
          if (skel(i+2,j+1)+skel(i+2,j+2)+skel(i+1,j+2)+skel(i+1,j)+skel(i,j+1)+skel(i,j)+skel(i+2,j)+skel(i,j+2)==1)&
          skel(i+1,j+1) = 0
          if (skel(i,j-1)+skel(i,j)+skel(i-1,j)+skel(i-1,j-2)+skel(i-2,j-1)+skel(i-2,j-2)+skel(i,j-2)+skel(i-2,j)==1)&
          skel(i-1,j-1) = 0

        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt1(2,2)).and.(skel(i-1,j)==filt1(2,3)).and.(skel(i+1,j-1)==filt1(3,1)).and.(skel(i,j-1)==filt1(3,2))&
        .and.((skel(i-2,j+1)+skel(i-2,j+2)+skel(i-1,j+2)+skel(i-1,j)+skel(i,j+1)+skel(i,j)+skel(i-2,j)+skel(i,j+2)==1).or.&
        (skel(i,j-1)+skel(i,j)+skel(i+1,j)+skel(i+1,j-2)+skel(i+2,j-1)+skel(i+2,j-2)+skel(i,j-2)+skel(i+2,j)==1))&
        .and.(skel(i-1,j-1)==filt1(3,3)).and.(skel(i,j+1)==filt1(1,2)).and.(skel(i-1,j+1)==filt1(1,3))) then
          if (skel(i-2,j+1)+skel(i-2,j+2)+skel(i-1,j+2)+skel(i-1,j)+skel(i,j+1)+skel(i,j)+skel(i-2,j)+skel(i,j+2)==1)&
          skel(i-1,j+1) = 0
          if (skel(i,j-1)+skel(i,j)+skel(i+1,j)+skel(i+1,j-2)+skel(i+2,j-1)+skel(i+2,j-2)+skel(i,j-2)+skel(i+2,j)==1)&
          skel(i+1,j-1) = 0

        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt1(2,2)).and.(skel(i+1,j)==filt1(2,3)).and.(skel(i-1,j+1)==filt1(3,1)).and.(skel(i,j+1)==filt1(3,2))&
        .and.((skel(i+2,j-1)+skel(i+2,j-2)+skel(i+1,j-2)+skel(i+1,j)+skel(i,j-1)+skel(i,j)+skel(i+2,j)+skel(i,j-2)==1).or.&
        (skel(i,j+1)+skel(i,j)+skel(i-1,j)+skel(i-1,j+2)+skel(i-2,j+1)+skel(i-2,j+2)+skel(i,j+2)+skel(i-2,j)==1))&
        .and.(skel(i+1,j+1)==filt1(3,3)).and.(skel(i,j-1)==filt1(1,2)).and.(skel(i+1,j-1)==filt1(1,3))) then
          if (skel(i+2,j-1)+skel(i+2,j-2)+skel(i+1,j-2)+skel(i+1,j)+skel(i,j-1)+skel(i,j)+skel(i+2,j)+skel(i,j-2)==1)&
          skel(i+1,j-1) = 0
          if (skel(i,j+1)+skel(i,j)+skel(i-1,j)+skel(i-1,j+2)+skel(i-2,j+1)+skel(i-2,j+2)+skel(i,j+2)+skel(i-2,j)==1)&
          skel(i-1,j+1) = 0

        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt1(2,2)).and.(skel(i-1,j)==filt1(2,3)).and.(skel(i+1,j+1)==filt1(3,1)).and.(skel(i,j+1)==filt1(3,2))&
        .and.((skel(i-2,j-1)+skel(i-2,j-2)+skel(i-1,j-2)+skel(i-1,j)+skel(i,j-1)+skel(i,j)+skel(i-2,j)+skel(i,j-2)==1).or.&
        (skel(i,j+1)+skel(i,j)+skel(i+1,j)+skel(i+1,j+2)+skel(i+2,j+1)+skel(i+2,j+2)+skel(i,j+2)+skel(i+2,j)==1))&
        .and.(skel(i-1,j+1)==filt1(3,3)).and.(skel(i,j-1)==filt1(1,2)).and.(skel(i-1,j-1)==filt1(1,3))) then
          if (skel(i-2,j-1)+skel(i-2,j-2)+skel(i-1,j-2)+skel(i-1,j)+skel(i,j-1)+skel(i,j)+skel(i-2,j)+skel(i,j-2)==1)&
          skel(i-1,j-1) = 0
          if (skel(i,j+1)+skel(i,j)+skel(i+1,j)+skel(i+1,j+2)+skel(i+2,j+1)+skel(i+2,j+2)+skel(i,j+2)+skel(i+2,j)==1)&
          skel(i+1,j+1) = 0

        endif
      enddo
    enddo

    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt10(2,2)).and.(skel(i,j+1)==filt10(2,3)).and.(skel(i+1,j-1)==filt10(1,1)).and.(skel(i-1,j)==filt10(3,2))&
        .and.(skel(i-1,j-1)==filt10(3,1)).and.(skel(i,j-1)==filt10(2,1))&
        .and.(skel(i-1,j+1)==filt10(3,3)).and.(skel(i+1,j)==filt10(1,2)).and.(skel(i+1,j+1)==filt10(1,3))) then
          skel(i-1,j+1) = 1  
          skel(i,j+1) = 0
          skel(i-1,j) = 0
          if (skel(i+1,j+2)+skel(i,j+2)+skel(i-1,j+2)>0) then
                  skel(i-1,j+1) = 0
                  skel(i,j+1) = 1
          endif
          if (skel(i-2,j-1)+skel(i-2,j)+skel(i-2,j+1)>0) then
                  skel(i-1,j+1) = 0
                  skel(i-1,j) = 1
          endif
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt10(2,2)).and.(skel(i,j-1)==filt10(2,3)).and.(skel(i+1,j+1)==filt10(1,1)).and.(skel(i-1,j)==filt10(3,2))&
        .and.(skel(i-1,j+1)==filt10(3,1)).and.(skel(i,j+1)==filt10(2,1))&
        .and.(skel(i-1,j-1)==filt10(3,3)).and.(skel(i+1,j)==filt10(1,2)).and.(skel(i+1,j-1)==filt10(1,3))) then
          skel(i-1,j-1) = 1  
          skel(i,j-1) = 0
          skel(i-1,j) = 0
          if (skel(i+1,j-2)+skel(i,j-2)+skel(i-1,j-2)>0) then
                  skel(i-1,j-1) = 0
                  skel(i,j-1) = 1
          endif
          if (skel(i-2,j+1)+skel(i-2,j)+skel(i-2,j-1)>0) then
                  skel(i-1,j-1) = 0
                  skel(i-1,j) = 1
          endif
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt10(2,2)).and.(skel(i,j+1)==filt10(2,3)).and.(skel(i-1,j-1)==filt10(1,1)).and.(skel(i+1,j)==filt10(3,2))&
        .and.(skel(i+1,j-1)==filt10(3,1)).and.(skel(i,j-1)==filt10(2,1))&
        .and.(skel(i+1,j+1)==filt10(3,3)).and.(skel(i-1,j)==filt10(1,2)).and.(skel(i-1,j+1)==filt10(1,3))) then
          skel(i+1,j+1) = 1  
          skel(i,j+1) = 0
          skel(i+1,j) = 0
          if (skel(i-1,j+2)+skel(i,j+2)+skel(i+1,j+2)>0) then
                  skel(i+1,j+1) = 0
                  skel(i,j+1) = 1
          endif
          if (skel(i+2,j-1)+skel(i+2,j)+skel(i+2,j+1)>0) then
                  skel(i+1,j+1) = 0
                  skel(i+1,j) = 1
          endif
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt10(2,2)).and.(skel(i,j-1)==filt10(2,3)).and.(skel(i-1,j+1)==filt10(1,1)).and.(skel(i+1,j)==filt10(3,2))&
        .and.(skel(i+1,j+1)==filt10(3,1)).and.(skel(i,j+1)==filt10(2,1))&
        .and.(skel(i+1,j-1)==filt10(3,3)).and.(skel(i-1,j)==filt10(1,2)).and.(skel(i-1,j-1)==filt10(1,3))) then
          skel(i+1,j-1) = 1  
          skel(i,j-1) = 0
          skel(i+1,j) = 0
          if (skel(i-1,j-2)+skel(i,j-2)+skel(i+1,j-2)>0) then
                  skel(i+1,j-1) = 0
                  skel(i,j-1) = 1
          endif
          if (skel(i+2,j+1)+skel(i+2,j)+skel(i+2,j-1)>0) then
                  skel(i+1,j-1) = 0
                  skel(i+1,j) = 1
          endif
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2

        if ((skel(i,j)==filt10(2,2)).and.(skel(i+1,j)==filt10(2,3)).and.(skel(i-1,j+1)==filt10(1,1)).and.(skel(i,j-1)==filt10(3,2))&
        .and.(skel(i-1,j-1)==filt10(3,1)).and.(skel(i-1,j)==filt10(2,1))&
        .and.(skel(i+1,j-1)==filt10(3,3)).and.(skel(i,j+1)==filt10(1,2)).and.(skel(i+1,j+1)==filt10(1,3))) then
          skel(i+1,j-1) = 1  
          skel(i+1,j) = 0
          skel(i,j-1) = 0
          if (skel(i+2,j+1)+skel(i+2,j)+skel(i+2,j-1)>0) then
                  skel(i+1,j-1) = 0
                  skel(i+1,j) = 1
          endif
          if (skel(i-1,j-2)+skel(i,j-2)+skel(i+1,j-2)>0) then
                  skel(i+1,j-1) = 0
                  skel(i,j-1) = 1
          endif
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt10(2,2)).and.(skel(i-1,j)==filt10(2,3)).and.(skel(i+1,j+1)==filt10(1,1)).and.(skel(i,j-1)==filt10(3,2))&
        .and.(skel(i+1,j-1)==filt10(3,1)).and.(skel(i+1,j)==filt10(2,1))&
        .and.(skel(i-1,j-1)==filt10(3,3)).and.(skel(i,j+1)==filt10(1,2)).and.(skel(i-1,j+1)==filt10(1,3))) then
          skel(i-1,j-1) = 1  
          skel(i-1,j) = 0
          skel(i,j-1) = 0
          if (skel(i-2,j+1)+skel(i-2,j)+skel(i-2,j-1)>0) then
                  skel(i-1,j-1) = 0
                  skel(i-1,j) = 1
          endif
          if (skel(i+1,j-2)+skel(i,j-2)+skel(i-1,j-2)>0) then
                  skel(i-1,j-1) = 0
                  skel(i,j-1) = 1
          endif
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt10(2,2)).and.(skel(i+1,j)==filt10(2,3)).and.(skel(i-1,j-1)==filt10(1,1)).and.(skel(i,j-1)==filt10(3,2))&
        .and.(skel(i-1,j+1)==filt10(3,1)).and.(skel(i-1,j)==filt10(2,1))&
        .and.(skel(i+1,j+1)==filt10(3,3)).and.(skel(i,j+1)==filt10(1,2)).and.(skel(i+1,j-1)==filt10(1,3))) then
          skel(i+1,j+1) = 1  
          skel(i+1,j) = 0
          skel(i,j+1) = 0
          if (skel(i+2,j-1)+skel(i+2,j)+skel(i+2,j+1)>0) then
                  skel(i+1,j+1) = 0
                  skel(i+1,j) = 1
          endif
          if (skel(i-1,j+2)+skel(i,j+2)+skel(i+1,j+2)>0) then
                  skel(i+1,j+1) = 0
                  skel(i,j+1) = 1
          endif
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt10(2,2)).and.(skel(i-1,j)==filt10(2,3)).and.(skel(i+1,j-1)==filt10(1,1)).and.(skel(i,j-1)==filt10(3,2))&
        .and.(skel(i+1,j+1)==filt10(3,1)).and.(skel(i+1,j)==filt10(2,1))&
        .and.(skel(i-1,j+1)==filt10(3,3)).and.(skel(i,j+1)==filt10(1,2)).and.(skel(i-1,j-1)==filt10(1,3))) then
          skel(i-1,j+1) = 1  
          skel(i-1,j) = 0
          skel(i,j+1) = 0
          if (skel(i-2,j-1)+skel(i-2,j)+skel(i-2,j+1)>0) then
                  skel(i-1,j+1) = 0
                  skel(i-1,j) = 1
          endif
          if (skel(i+1,j+2)+skel(i,j+2)+skel(i-1,j+2)>0) then
                  skel(i-1,j+1) = 0
                  skel(i,j+1) = 1
          endif
        endif
      enddo
    enddo


    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt6(2,2)).and.(skel(i,j+1)==filt6(2,3)).and.(skel(i+1,j+1)==filt6(1,3))&
        .and.(skel(i+1,j)==filt6(1,2)).and.(skel(i-1,j)==filt6(3,2)).and.(skel(i-1,j+1)==filt6(3,3))) then
          if ((skel(i+1,j-1)==1).or.(skel(i,j-1)==1)) then
            if (rhoSlices(i,j+1)<rhoSlices(i+1,j+1)) then
              skel(i,j+1) = 0
              bool(i+1,j+1) = 1
            else
              skel(i+1,j+1) = 0
              bool(i,j+1) = 1
            endif
          else
          skel(i,j+1)=0
          endif

          if (skel(i-1,j+2)==1) then
            skel(i,j+1) = 1
            skel(i+1,j+1) = 0
          else if (skel(i+2,j)+skel(i+2,j+1)+skel(i+2,j+2)>0) then
            skel(i+1,j+1) = 1
            skel(i,j+1) = 0
          endif

          if (skel(i+1,j+2)+skel(i+2,j+2)+skel(i+2,j+1)+skel(i,j+1)+skel(i+1,j)+skel(i,j)+skel(i,j+2)+skel(i+2,j)==2) then
            skel(i+1,j+1) = 1
            skel(i,j+1) = 0
          endif
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt6(2,2)).and.(skel(i,j-1)==filt6(2,3)).and.(skel(i+1,j-1)==filt6(1,3))&
        .and.(skel(i+1,j)==filt6(1,2)).and.(skel(i-1,j)==filt6(3,2)).and.(skel(i-1,j-1)==filt6(3,3))) then
          if ((skel(i+1,j+1)==1).or.(skel(i,j+1)==1)) then
            if (rhoSlices(i,j-1)<rhoSlices(i+1,j-1)) then
              skel(i,j-1) = 0
              bool(i+1,j-1) = 1
            else
              skel(i+1,j-1) = 0
              bool(i,j-1) = 1
            endif
            else
            skel(i,j-1)=0
            endif

          if (skel(i-1,j-2)==1) then
            skel(i,j-1) = 1
            skel(i+1,j-1) = 0
          else if (skel(i+2,j)+skel(i+2,j-1)+skel(i+2,j-2)>0) then
            skel(i+1,j-1) = 1
            skel(i,j-1) = 0
          endif

          if (skel(i+1,j-2)+skel(i+2,j-2)+skel(i+2,j-1)+skel(i,j-1)+skel(i+1,j)+skel(i,j)+skel(i,j-2)+skel(i+2,j)==2) then
            skel(i+1,j-1) = 1
            skel(i,j-1) = 0
          endif
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt6(2,2)).and.(skel(i,j+1)==filt6(2,3)).and.(skel(i-1,j+1)==filt6(1,3))&
        .and.(skel(i-1,j)==filt6(1,2)).and.(skel(i+1,j)==filt6(3,2)).and.(skel(i+1,j+1)==filt6(3,3))) then
          if ((skel(i-1,j-1)==1).or.(skel(i,j-1)==1)) then
            if (rhoSlices(i,j+1)<rhoSlices(i-1,j+1)) then
              skel(i,j+1) = 0
              bool(i-1,j+1) = 1
            else
              skel(i-1,j+1) = 0
              bool(i,j+1) = 1
            endif
            else
          skel(i,j+1)=0
          endif

          if (skel(i+1,j+2)==1) then
            skel(i,j+1) = 1
            skel(i-1,j+1) = 0
          else if (skel(i-2,j)+skel(i-2,j+1)+skel(i-2,j+2)>0) then
            skel(i-1,j+1) = 1
            skel(i,j+1) = 0
          endif

          if (skel(i-1,j+2)+skel(i-2,j+2)+skel(i-2,j+1)+skel(i,j+1)+skel(i-1,j)+skel(i,j)+skel(i,j+2)+skel(i-2,j)==2) then
            skel(i-1,j+1) = 1
            skel(i,j+1) = 0
          endif
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt6(2,2)).and.(skel(i,j-1)==filt6(2,3)).and.(skel(i-1,j-1)==filt6(1,3))&
        .and.(skel(i-1,j)==filt6(1,2)).and.(skel(i+1,j)==filt6(3,2)).and.(skel(i+1,j-1)==filt6(3,3))) then
          if ((skel(i-1,j+1)==1).or.(skel(i,j+1)==1)) then
            if (rhoSlices(i,j-1)<rhoSlices(i-1,j-1)) then
              skel(i,j-1) = 0
              bool(i-1,j-1) = 1
            else
              skel(i-1,j-1) = 0
              bool(i,j-1) = 1
            endif
            else
          skel(i,j-1)=0
          endif

          if (skel(i+1,j-2)==1) then
            skel(i,j-1) = 1
            skel(i-1,j-1) = 0
          else if (skel(i-2,j)+skel(i-2,j-1)+skel(i-2,j-2)>0) then
            skel(i-1,j-1) = 1
            skel(i,j-1) = 0
          endif

          if (skel(i-1,j-2)+skel(i-2,j-2)+skel(i-2,j-1)+skel(i,j-1)+skel(i-1,j)+skel(i,j)+skel(i,j-2)+skel(i-2,j)==2) then
            skel(i-1,j-1) = 1
            skel(i,j-1) = 0
          endif
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2

        if ((skel(i,j)==filt6(2,2)).and.(skel(i+1,j)==filt6(2,3)).and.(skel(i+1,j+1)==filt6(1,3))&
        .and.(skel(i,j+1)==filt6(1,2)).and.(skel(i,j-1)==filt6(3,2)).and.(skel(i+1,j-1)==filt6(3,3))) then
          if ((skel(i-1,j+1)==1).or.(skel(i-1,j)==1)) then
            if (rhoSlices(i+1,j)<rhoSlices(i+1,j+1)) then
              skel(i+1,j) = 0
              bool(i+1,j+1) = 1
            else
              skel(i+1,j+1) = 0
              bool(i+1,j) = 1
            endif
            else
          skel(i+1,j)=0
          endif

          if (skel(i+2,j-1)==1) then
            skel(i+1,j) = 1
            skel(i+1,j+1) = 0
          else if (skel(i,j+2)+skel(i+1,j+2)+skel(i+2,j+2)>0) then
            skel(i+1,j+1) = 1
            skel(i+1,j) = 0
          endif

          if (skel(i+2,j+1)+skel(i+2,j+2)+skel(i+1,j+2)+skel(i+1,j)+skel(i,j+1)+skel(i,j)+skel(i+2,j)+skel(i,j+2)==2) then
            skel(i+1,j+1) = 1
            skel(i+1,j) = 0
          endif
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt6(2,2)).and.(skel(i-1,j)==filt6(2,3)).and.(skel(i-1,j+1)==filt6(1,3))&
        .and.(skel(i,j+1)==filt6(1,2)).and.(skel(i,j-1)==filt6(3,2)).and.(skel(i-1,j-1)==filt6(3,3))) then
          if ((skel(i+1,j+1)==1).or.(skel(i+1,j)==1)) then
            if (rhoSlices(i-1,j)<rhoSlices(i-1,j+1)) then
              skel(i-1,j) = 0
              bool(i-1,j+1) = 1
            else
              skel(i-1,j+1) = 0
              bool(i-1,j) = 1
            endif
           else
          skel(i-1,j)=0
          endif

          if (skel(i-2,j-1)==1) then
            skel(i-1,j) = 1
            skel(i-1,j+1) = 0
          else if (skel(i,j+2)+skel(i-1,j+2)+skel(i-2,j+2)>0) then
            skel(i-1,j+1) = 1
            skel(i-1,j) = 0
          endif

          if (skel(i-2,j+1)+skel(i-2,j+2)+skel(i-1,j+2)+skel(i-1,j)+skel(i,j+1)+skel(i,j)+skel(i-2,j)+skel(i,j+2)==2) then
            skel(i-1,j+1) = 1
            skel(i-1,j) = 0
          endif
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt6(2,2)).and.(skel(i+1,j)==filt6(2,3)).and.(skel(i+1,j-1)==filt6(1,3))&
        .and.(skel(i,j+1)==filt6(1,2)).and.(skel(i,j-1)==filt6(3,2)).and.(skel(i+1,j+1)==filt6(3,3))) then
          if ((skel(i-1,j-1)==1).or.(skel(i-1,j)==1)) then
            if (rhoSlices(i+1,j)<rhoSlices(i+1,j-1)) then
              skel(i+1,j) = 0
              bool(i+1,j-1) = 1
            else
              skel(i+1,j-1) = 0
              bool(i+1,j) = 1
            endif
            else
          skel(i+1,j)=0
          endif

          if (skel(i+2,j+1)==1) then
            skel(i+1,j) = 1
            skel(i+1,j-1) = 0
          else if (skel(i,j-2)+skel(i+1,j-2)+skel(i+2,j-2)>0) then
            skel(i+1,j-1) = 1
            skel(i+1,j) = 0
          endif

          if (skel(i+2,j-1)+skel(i+2,j-2)+skel(i+1,j-2)+skel(i+1,j)+skel(i,j-1)+skel(i,j)+skel(i+2,j)+skel(i,j-2)==2) then
            skel(i+1,j-1) = 1
            skel(i+1,j) = 0
          endif
        endif
      enddo
    enddo
    do i=3,size(skel,1)-2
      do j=3,size(skel,2)-2
        if ((skel(i,j)==filt6(2,2)).and.(skel(i-1,j)==filt6(2,3)).and.(skel(i-1,j-1)==filt6(1,3))&
        .and.(skel(i,j+1)==filt6(1,2)).and.(skel(i,j-1)==filt6(3,2)).and.(skel(i-1,j+1)==filt6(3,3))) then
          if ((skel(i+1,j-1)==1).or.(skel(i+1,j)==1)) then
            if (rhoSlices(i-1,j)<rhoSlices(i-1,j-1)) then
              skel(i-1,j) = 0
              bool(i-1,j-1) = 1
            else
              skel(i-1,j-1) = 0
              bool(i-1,j) = 1
            endif
            else
          skel(i-1,j)=0
          endif

          if (skel(i-2,j+1)==1) then
            skel(i-1,j) = 1
            skel(i-1,j-1) = 0
          else if (skel(i,j-2)+skel(i-1,j-2)+skel(i-2,j-2)>0) then
            skel(i-1,j-1) = 1
            skel(i-1,j) = 0
          endif

          if (skel(i-2,j-1)+skel(i-2,j-2)+skel(i-1,j-2)+skel(i-1,j)+skel(i,j-1)+skel(i,j)+skel(i-2,j)+skel(i,j-2)==2) then
            skel(i-1,j-1) = 1
            skel(i-1,j) = 0
          endif
        endif



      enddo
    enddo

    do i=2,size(skel,1)-1
      do j=2,size(skel,2)-1

        if ((skel(i,j)==1).and.(skel(i,j+1)==1).and.(skel(i+1,j+1)==1)) write(*,*) "error NAN1"
        if ((skel(i,j)==1).and.(skel(i,j+1)==1).and.(skel(i-1,j+1)==1)) write(*,*) "error NAN2"
        if ((skel(i,j)==1).and.(skel(i,j-1)==1).and.(skel(i+1,j-1)==1)) write(*,*) "error NAN3"
        if ((skel(i,j)==1).and.(skel(i,j-1)==1).and.(skel(i-1,j-1)==1)) write(*,*) "error NAN4"
        if ((skel(i,j)==1).and.(skel(i+1,j)==1).and.(skel(i+1,j+1)==1)) write(*,*) "error NAN5"
        if ((skel(i,j)==1).and.(skel(i+1,j)==1).and.(skel(i+1,j-1)==1)) write(*,*) "error NAN6"
        if ((skel(i,j)==1).and.(skel(i-1,j)==1).and.(skel(i-1,j+1)==1)) write(*,*) "error NAN7"
        if ((skel(i,j)==1).and.(skel(i-1,j)==1).and.(skel(i-1,j-1)==1)) write(*,*) "error NAN8"
      enddo
    enddo

  end subroutine filterskel

  subroutine filtermidline(midline,rhoSlices)
    implicit none
    real(pr),dimension(:,:),intent(in) :: rhoSlices
    integer,dimension(:,:),intent(inout) :: midline
    integer,dimension(size(midline,1),size(midline,2)) :: midlinetmp
    integer :: l

    do l=1,size(midline,1)
       midlinetmp(l,1) = midline(size(midline,1)-l+1,1)
       midlinetmp(l,2) = midline(size(midline,1)-l+1,2)
    enddo
    midline = midlinetmp

    midlinetmp = midline
    do l=2,size(midline,1)-1
       if ((midline(l+1,1)-midline(l,1)==+2)) then
          midlinetmp(l+1,1) = midline(l,1)+1
       else if ((midline(l+1,1)-midline(l,1)==-2)) then
          midlinetmp(l+1,1) = midline(l,1)-1
       else if ((midline(l+1,2)-midline(l,2)==+2)) then
          midlinetmp(l+1,2) = midline(l,2)+1
       else if ((midline(l+1,2)-midline(l,2)==-2)) then
          midlinetmp(l+1,2) = midline(l,2)-1
       else if (midline(l-1,1)==midline(l,1)) then !2pts aligned
          
          if ((midline(l-1,2)+midline(l+1,2)==2*midline(l,2)).and.(abs(midline(l+1,1)-midline(l,1))==1)) then
             midlinetmp(l+1,1) = midline(l,1)
        endif

       else if (midline(l-1,2)==midline(l,2)) then !2pts aligned
          if ((midline(l-1,1)+midline(l+1,1)==2*midline(l,1)).and.(abs(midline(l+1,2)-midline(l,2))==1)) then
             midlinetmp(l+1,2) = midline(l,2)
          endif

       else if ((abs(midline(l-1,1)-midline(l,1))==1).and.(abs(midline(l-1,2)-midline(l,2))==1)) then !2 pts diagonal
          if (midline(l+1,1)==midline(l-1,1)) midlinetmp(l,1) = midline(l-1,1)
          if (midline(l+1,2)==midline(l-1,2)) midlinetmp(l,2) = midline(l-1,2)
          if (midline(l+1,1)==midline(l,1)) then !2pts aligned encore
            midlinetmp(l,1) = midline(l-1,1)
             midlinetmp(l+1,1) = midline(l-1,1)
          endif
          if (midline(l+1,2)==midline(l,2)) then !2pts aligned encore
             midlinetmp(l,2) = midline(l-1,2)
             midlinetmp(l+1,2) = midline(l-1,2)
          endif
       endif
       midline = midlinetmp
    enddo

    do l=1,size(midline,1)
       midlinetmp(l,1) = midline(size(midline,1)-l+1,1)
       midlinetmp(l,2) = midline(size(midline,1)-l+1,2)
    enddo
    midline = midlinetmp

  end subroutine filtermidline

  subroutine findL(midline,ii,jj,ll)
    implicit none
    integer,intent(in) :: ii,jj
    integer,intent(inout) :: ll
    integer,dimension(:,:),intent(in) :: midline

    ll=1
    do while (((.not.(ii==midline(ll,1))).or.(.not.(jj==midline(ll,2)))).and.(ll<size(midline,1)))
      ll = ll+1
    enddo
  end subroutine


  subroutine cutheadskel(midline,skel,skel2,nl,ii,jj,kt)
    implicit none
    integer,intent(in) :: nl,kt
    integer,intent(inout) :: ii,jj
    integer,dimension(:,:),intent(inout) :: skel2
    integer,dimension(:,:),intent(inout) :: skel,midline
    integer,dimension(size(skel2,1),size(skel2,2)) :: skeltmp
    integer,dimension(size(midline,1),size(midline,2)) :: midlinetmp
    integer :: l,bool,ll,lm,lp,bool22,iinew,jjnew
    integer :: ll1,ll2,ll3,ll4,ll5,ll6,ll7,ll8,ll9,i,j,nnz
    integer,dimension(:),allocatable :: tab
    integer,dimension(9) :: tabll
    real(pr) :: mintab
    ll1 = 0
    ll2 = 0
    ll3 = 0
    ll4 = 0
    ll5 = 0
    ll6 = 0
    ll7 = 0
    ll8 = 0
    ll9 = 0
    ll = 0

    skeltmp = 0
    do i=1,size(skel2,1)
      do j=1,size(skel2,2)
        if (skel2(i,j)==1) then
          skeltmp(i,j) = 1
          skeltmp(i,j+1) = 1
          skeltmp(i,j-1) = 1
          skeltmp(i-1,j) = 1
          skeltmp(i+1,j) = 1
          skeltmp(i+1,j+1) = 1
          skeltmp(i-1,j-1) = 1
          skeltmp(i+1,j-1) = 1
          skeltmp(i-1,j+1) = 1
        endif
      enddo
    enddo
    skel2 = skeltmp

    do l=1,size(midline,1)
       midlinetmp(l,1) = midline(size(midline,1)-l+1,1)
       midlinetmp(l,2) = midline(size(midline,1)-l+1,2)
    enddo
    midline = midlinetmp

    skeltmp = 0
    if ((ii==-1).and.(jj==-1)) then
      skeltmp = skel2

    else
    if (skel(ii-1,jj-1)==1) &
      call findL(midline,ii-1,jj-1,ll1)
    if ((skel(ii-1,jj)==1).or.(skel(ii-2,jj)==1)) then
            if (skel(ii-1,jj)==1) then
                    call findL(midline,ii-1,jj,ll2)
            else
                    call findL(midline,ii-2,jj,ll2)
            endif
    endif
    if (skel(ii-1,jj+1)==1) &
      call findL(midline,ii-1,jj+1,ll3)
    if ((skel(ii,jj-1)==1).or.(skel(ii,jj-2)==1)) then
            if (skel(ii,jj-1)==1) then
                    call findL(midline,ii,jj-1,ll4)
            else
                    call findL(midline,ii,jj-2,ll4)
            endif
    endif
    if (skel(ii,jj)==1) &
      call findL(midline,ii,jj,ll5)
    if ((skel(ii,jj+1)==1).or.(skel(ii,jj+2)==1)) then
            if (skel(ii,jj+1)==1) then
                    call findL(midline,ii,jj+1,ll6)
            else
                    call findL(midline,ii,jj+2,ll6)
            endif
    endif
    if (skel(ii+1,jj-1)==1) &
      call findL(midline,ii+1,jj-1,ll7)
    if ((skel(ii+1,jj)==1).or.(skel(ii+2,jj)==1)) then
            if (skel(ii+1,jj)==1) then
                    call findL(midline,ii+1,jj,ll8)
            else
                    call findL(midline,ii+2,jj,ll8)
            endif
     endif
    if (skel(ii+1,jj+1)==1) &
      call findL(midline,ii+1,jj+1,ll9)
      tabll(1) = ll1
      tabll(2) = ll2
      tabll(3) = ll3
      tabll(4) = ll4
      tabll(5) = ll5
      tabll(6) = ll6
      tabll(7) = ll7
      tabll(8) = ll8
      tabll(9) = ll9
      nnz = 9
      if (ll1==0) nnz = nnz-1
      if (ll2==0) nnz = nnz-1
      if (ll3==0) nnz = nnz-1
      if (ll4==0) nnz = nnz-1
      if (ll5==0) nnz = nnz-1
      if (ll6==0) nnz = nnz-1
      if (ll7==0) nnz = nnz-1
      if (ll8==0) nnz = nnz-1
      if (ll9==0) nnz = nnz-1
      if (nnz>0) then
              allocate(tab(nnz))
              mintab = 1E6
              j = 0
              do i=1,9
                if (tabll(i)>0) then
                        j=j+1
                        tab(j) = tabll(i)
                        if (sqrt((ii*1._pr-midline(tabll(i),1))**2+(jj*1._pr-midline(tabll(i),2))**2)<mintab) then
                                mintab = sqrt((ii*1._pr-midline(tabll(i),1))**2+(jj*1._pr-midline(tabll(i),2))**2)
                                ll = tabll(i)
                        endif
                endif
              enddo
              if ((ii==-1).and.(jj==-1)) ll = minval(tab)
              ll = minval(tab)
              if (ll==ll1) write(*,*) "CASE 1 ",mintab
              if (ll==ll2) write(*,*) "CASE 2 ",mintab
              if (ll==ll3) write(*,*) "CASE 3 ",mintab
              if (ll==ll4) write(*,*) "CASE 4 ",mintab
              if (ll==ll5) write(*,*) "CASE 5 ",mintab
              if (ll==ll6) write(*,*) "CASE 6 ",mintab
              if (ll==ll7) write(*,*) "CASE 7 ",mintab
              if (ll==ll8) write(*,*) "CASE 8 ",mintab
              if (ll==ll9) write(*,*) "CASE 9 ",mintab

              bool = 0
              lp = ll
              do while ((bool==0).and.(lp<size(midline,1)))
               
                if (skel2(midline(lp,1),midline(lp,2))==1) skeltmp(midline(lp,1),midline(lp,2)) = 1
                if (skel2(midline(lp,1),midline(lp,2))==1) write(*,*) "ppp1  ",midline(lp,1)," ",midline(lp,2)
                if ((skeltmp(midline(lp,1),midline(lp,2))==1).and.(skel2(midline(lp+1,1),midline(lp+1,2))==0)) bool = 1   
                lp=lp+1
              enddo
              if ((lp==size(midline,1)).and.(skel2(midline(size(midline,1),1),midline(size(midline,1),2))==1))&
                skeltmp(midline(size(midline,1),1),midline(size(midline,1),2)) = 1
              if ((lp==size(midline,1)).and.(skel2(midline(size(midline,1),1),midline(size(midline,1),2))==1))&
                write(*,*) "P1 ",midline(lp,1)," ",midline(lp,2)
              bool = 0
              lm = ll
              do while ((bool==0).and.(lm>size(midline,1)-nl+1))
                write(*,*) "mm0  ",midline(lm,1)," ",midline(lm,2)
               
                if (skel2(midline(lm,1),midline(lm,2))==1) skeltmp(midline(lm,1),midline(lm,2)) = 1
                if (skel2(midline(lm,1),midline(lm,2))==1) write(*,*) "mm1  ",midline(lm,1)," ",midline(lm,2)
                if ((skeltmp(midline(lm,1),midline(lm,2))==1).and.(skel2(midline(lm-1,1),midline(lm-1,2))==0)) bool = 1   
                lm=lm-1
              enddo
              if ((lm==size(midline,1)-nl+1).and.(skel2(midline(lm,1),midline(lm,2))==1)) skeltmp(midline(lm,1),midline(lm,2)) = 1
              if ((lm==size(midline,1)-nl+1).and.(skel2(midline(lm,1),midline(lm,2))==1)) &
              write(*,*) "M1 ",midline(lm,1)," ",midline(lm,2)
              deallocate(tab)
      else
              write(*,*) "error STOP"
      endif



    endif

    skel2 = skeltmp


    bool = 0
    bool22 = 0
    skeltmp = skel
    l=size(midline,1)-nl+1
    do while ((bool==0).and.(l<size(midline,1)).and.(nl>1))
       write(*,*) "writeL (i,j) ",midline(l,1)," ",midline(l,2)," ",l
      if (((skel(midline(l,1),midline(l,2))==1).and.(skel2(midline(l+1,1),midline(l+1,2))==1)).and.&
      (sqrt((ii*1._pr-midline(l+1,1))**2+(jj*1._pr-midline(l+1,2))**2)<&
      sqrt((ii*1._pr-midline(l,1))**2+(jj*1._pr-midline(l,2))**2))) then
              skeltmp(midline(l,1),midline(l,2)) = 0   
      write(*,*) "AIE 11"
      elseif (((skel(midline(l,1),midline(l,2))==1).and.(skel2(midline(l+1,1),midline(l+1,2))==1)).and.&
      (sqrt((ii*1._pr-midline(l+1,1))**2+(jj*1._pr-midline(l+1,2))**2)>=&
      sqrt((ii*1._pr-midline(l,1))**2+(jj*1._pr-midline(l,2))**2))) then
      if (skel2(midline(l,1),midline(l,2))==0) then
              skeltmp(midline(l,1),midline(l,2)) = 0
      write(*,*) "AIE 22"
      else
              iinew = midline(l,1)
              jjnew = midline(l,2)
              bool22 = 1
              bool = 1
      endif
      endif
      if ((skel2(midline(l,1),midline(l,2))==1).and.(skel2(midline(l+1,1),midline(l+1,2))==0)) bool = 1   
      l=l+1
    enddo
    if (bool==0) write(*,*) "error ouch"

    if (l>size(midline,1)-nl+1) then
      if ((sqrt((ii*1._pr-midline(l-1,1))**2+(jj*1._pr-midline(l-1,2))**2)<2.001).or.((ii==-1).and.(jj==-1))) then
      ii = midline(l,1)
      jj = midline(l,2)
      if (bool22==1) then
              ii = iinew
              jj = jjnew
      endif
      skel = skeltmp
      else
      skel = skeltmp
      write(*,*) "OUCH  ",sqrt((ii*1._pr-midline(l-1,1))**2+(jj*1._pr-midline(l-1,2))**2)
      ii = midline(l,1)
      jj = midline(l,2)
      endif
    else
      if ((sqrt((ii*1._pr-midline(l,1))**2+(jj*1._pr-midline(l,2))**2)<2.001).or.((ii==-1).and.(jj==-1))) then
      ii = midline(l,1)
      jj = midline(l,2)
      skel = skeltmp
      else
      write(*,*) "OUCH22  ",sqrt((ii*1._pr-midline(l,1))**2+(jj*1._pr-midline(l,2))**2)
      endif
    endif
      write(*,*) "Last tracked midline point ",ii," ",jj


    do l=1,size(midline,1)
       midlinetmp(l,1) = midline(size(midline,1)-l+1,1)
       midlinetmp(l,2) = midline(size(midline,1)-l+1,2)
    end do
    midline = midlinetmp

  end subroutine cutheadskel




subroutine body_rotationdef(courbe,xg,yg,alpha)
  implicit none
  real(pr),dimension(:,:),intent(inout) :: courbe
  real(pr),intent(inout) :: xg,yg
  real(pr),intent(out) :: alpha
  integer :: l,nb
  real(pr) :: dt,eps,dalpha,d,sinm,cosm,PI

  PI=acos(-1._pr)
  alpha=0._pr
  eps=1e-12
  dt=0._pr
  nb=0._pr
  do l=1,size(courbe,1)
    dalpha=atan((courbe(l,2)-yg)/(courbe(l,1)-xg+eps))
    d=sqrt((courbe(l,2)-yg)**2+(courbe(l,1)-xg+eps)**2)
    if ((abs(courbe(l,2)-yg)<eps).and.(abs(courbe(l,1)-xg)<eps)) then
       dalpha=0
    elseif (abs(abs(dalpha)-pi/2)<0.1*PI/2) then
        dalpha=atan((courbe(l,1)-xg)/(courbe(l,2)-yg+eps))
        d=sqrt((courbe(l,2)-yg+eps)**2+(courbe(l,1)-xg)**2)
        
        if ((courbe(l,1)-xg)*(courbe(l,2)-yg+eps)>0) then
            dalpha = PI/2 - dalpha
        else
            dalpha = -PI/2 - dalpha
        endif
        nb=nb+1
    endif
    dt=dt+d
    alpha=alpha+dalpha*d
  enddo
  alpha=alpha/dt


end subroutine body_rotationdef

subroutine body_rotationdefPS(courberef,courbe,xg,yg,alphadef,dx)
  implicit none
  real(pr),dimension(:,:),intent(inout) :: courbe,courberef
  real(pr),intent(inout) :: xg,yg,dx
  real(pr),intent(out) :: alphadef
  integer :: l,nb,iter,maxiter
  real(pr) :: alpha,dt,eps,dalpha,d,sinm,cosm,PI
  real(pr) :: PS,PSprev,sgn

  PI=acos(-1._pr)
  alpha=0._pr
  eps=1e-12
  dt=0._pr
  nb=0._pr
  do l=1,size(courbe,1)
    dalpha=atan((courbe(l,2)/dx-yg/dx)/(courbe(l,1)/dx-xg/dx+eps))
    d=sqrt((courbe(l,2)/dx-yg/dx)**2+(courbe(l,1)/dx-xg/dx+eps)**2)
    if ((abs(courbe(l,2)/dx-yg/dx)<eps).and.(abs(courbe(l,1)/dx-xg/dx)<eps)) then
       dalpha=0
    elseif (abs(abs(dalpha)-PI/2)<0.1*PI/2) then
        dalpha=atan((courbe(l,1)/dx-xg/dx)/(courbe(l,2)/dx-yg/dx+eps))
        d=sqrt((courbe(l,2)/dx-yg/dx+eps)**2+(courbe(l,1)/dx-xg/dx)**2)
        
        if ((courbe(l,1)/dx-xg/dx)*(courbe(l,2)/dx-yg/dx+eps)>0) then
            dalpha = PI/2 - dalpha
        else
            dalpha = -PI/2 - dalpha
        endif
        nb=nb+1
    endif
    dt=dt+d
    alpha=alpha+dalpha*d
  enddo
  alpha=alpha/dt
  alphadef=alpha

  PS=0._pr
  do l=1,size(courbe,1)-1
    PS = PS + (courbe(l+1,1)-courbe(l,1))/dx*(courberef(l+1,1)-courberef(l,1))/dx&
    + (courbe(l+1,2)-courbe(l,2))/dx*(courberef(l+1,2)-courberef(l,2))/dx
  enddo
  PS = PS/(size(courbe,1)-1)
  write(*,*) "PS11  0 ",PS
  if (PS<0._pr) then
          PS = abs(PS)
          alpha = alpha+PI
  endif
  alpha=10._pr*PI/180
  alphadef=alpha

  iter = 0
  PSprev = 0._pr
  alphadef=0._pr
  maxiter = 100
  do while(((abs(PSprev-PS)>eps*PSprev).or.(iter==0)).and.(iter<=maxiter))

    PSprev = PS
    call body_rotating_theta(courbe,xg,yg,-alpha)
    PS=0._pr
    do l=1,size(courbe,1)-1
      PS = PS + (courbe(l+1,1)-courbe(l,1))/dx*(courberef(l+1,1)-courberef(l,1))/dx&
      + (courbe(l+1,2)-courbe(l,2))/dx*(courberef(l+1,2)-courberef(l,2))/dx
    enddo
    PS = PS/(size(courbe,1)-1)
    alphadef=alphadef+alpha
    if (abs(PSprev-PS)<eps*PSprev) alphadef=alphadef-alpha*0.5_pr
    write(*,*) "ITER  ",iter," ",alpha*180/PI," ",PS," ",PSprev," ",abs(PSprev-PS)," ",alphadef*180/PI
    if (PS > PSprev) then
      alpha = alpha
    else
      alpha = -alpha*0.5_pr
    endif

    iter = iter+1
  enddo
    

  write(*,*) "ERREUR  ",modulo(alphadef,2*PI)-PI," ",PS," ",PSprev," ",abs(PSprev-PS)

end subroutine body_rotationdefps

subroutine body_rotating_theta(courbe,xg,yg,ang)
  implicit none
  real(pr),dimension(:,:),intent(inout) :: courbe
  real(pr),intent(inout) :: xg,yg
  real(pr),intent(in) :: ang
  integer :: l,method
  real(pr) :: xy,yy,x0,y0
 
  method=2
 
  do l=1,size(courbe,1)
       x0=courbe(l,1)
       y0=courbe(l,2)
       xy=(x0-xg)*cos(ang)+(y0-yg)*sin(ang)+xg
       yy=-(x0-xg)*sin(ang)+(y0-yg)*cos(ang)+yg
       courbe(l,1)=xy
       courbe(l,2)=yy
  enddo
end subroutine body_rotating_theta

subroutine body_masscenter(courbe,xg,yg)
  implicit none
  real(pr),dimension(:,:),intent(inout) :: courbe
  real(pr),intent(inout) :: xg,yg
  integer :: l

  xg=0._pr
  yg=0._pr
  do l=1,size(courbe,1)
    xg=xg+courbe(l,1)
    yg=yg+courbe(l,2)
  enddo
  xg=xg/size(courbe,1)
  yg=yg/size(courbe,1)

end subroutine body_masscenter

subroutine body_masscentering(courbe,xg,yg,xgref,ygref)
  implicit none
  real(pr),dimension(:,:),intent(inout) :: courbe
  real(pr),intent(inout) :: xg,yg,xgref,ygref
  integer :: l

  do l=1,size(courbe,1)
    courbe(l,1)=courbe(l,1)+xgref-xg
    courbe(l,2)=courbe(l,2)+ygref-yg
  enddo

end subroutine body_masscentering

subroutine compute_thetadef(courbe,alpha,calpha,salpha,kt)
  implicit none
  real(pr),dimension(:,:),intent(inout) :: courbe
  real(pr),intent(inout) :: alpha
  real(pr),intent(inout) :: calpha,salpha
  integer,intent(in) :: kt
  integer :: l,nb
  real(pr) :: PI,PS1,PS2,PV1,PV2,beta1,beta2,eps,norml,normlm,normlp,dx

  PI=acos(-1._pr)
  nb=0
  eps=1e-12
  dx=0.0256*1e-3

  do l=2,size(courbe,1)-1
    PS1 = (courbe(l+1,1)-courbe(l,1))*(courbe(l-1,1)-courbe(l,1)) + (courbe(l+1,2)-courbe(l,2))*(courbe(l-1,2)-courbe(l,2))
    PV1 = (courbe(l+1,1)-courbe(l,1))*(courbe(l-1,2)-courbe(l,2)) - (courbe(l+1,2)-courbe(l,2))*(courbe(l-1,1)-courbe(l,1))
    norml = sqrt((courbe(l+1,1)-courbe(l,1))*(courbe(l+1,1)-courbe(l,1))&
    + (courbe(l+1,2)-courbe(l,2))*(courbe(l+1,2)-courbe(l,2)))
    normlm = sqrt((courbe(l-1,1)-courbe(l,1))*(courbe(l-1,1)-courbe(l,1))&
    + (courbe(l-1,2)-courbe(l,2))*(courbe(l-1,2)-courbe(l,2)))
    beta1 = atan(PV1/(PS1+eps))


    if (abs(abs(beta1)-PI/2)<0.1*PI/2) then
        beta1=0
        if (PS1*PV1>0) then
        beta1=PI/2 - atan(PS1/(PV1+eps))
        else
        beta1=-PI/2 - atan(PS1/(PV1+eps))
        endif
        write(*,*) "TANGENTEpb1 ",abs(abs(beta1)-PI/2)
    endif





    alpha=alpha - beta1
    nb=nb+1
  enddo
  nb=nb+1


  calpha = cos(alpha)
  salpha = sin(alpha)



end subroutine compute_thetadef

end module libBezier




