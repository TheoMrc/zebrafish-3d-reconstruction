program def3D
   use variables
   use interpolation
   use AdvectionProblem
   use libBezier
   use midline2D_modules_theo

   implicit none

   !**
   !! Declaration of all variables (to sort)
   !**
   integer :: i, j, k, kt, idisplay, iUpdateDist, iter, picNum, l, ll, pix1, pix2, l0, pos, sizeSkel, &
         nl, lp, bool, booltmp, boolPhi, loopbool, xskelL, yskelL, xskelR, yskelR, ii, lph, lpt, lpf, &
         hh, ic, jc, bool1, bool2, n_dim, file_id
   real(pr) :: Tps, threshold, surf, maxv, minimL, minimL_m, minimL_p, minimR, long, tb, dtb, long2, long3, &
         dti, dtf, dsi, dsf, xr, yr, xl, yl, long00, tp, longexp, longratio, tbb
   real(pr), dimension(:), allocatable :: longslice, longslicei, longslicef, longTh, longTheta
   real(pr) :: t, tPic, xi, yj, distW, distWb, xLeft, yLeft, xRight, yRight
   real(pr), dimension(:, :), allocatable :: rhoSlices2, gradPhi
   integer, dimension(:, :), allocatable :: midline, midlinebis
   real(pr), dimension(:, :), allocatable :: tmp, tmp1, tmp2, un, zero, distslice, voisin
   integer, dimension(:, :), allocatable :: dir1, dir2, dir3, dir4, Nseed, skel, skel2, tmpbool, skel3
   real(pr) :: px, py, pz, s0, sinit, tbm, tbp, s, ds, rr, pxx, pyy, pzz
   real(pr), dimension(:), allocatable :: sslice, dsslice, stheta, dstheta, sth, dsth
   real(pr), dimension(:, :), allocatable :: points_control, points_courbe, points_courbe_equal, tail_courbe, head_courbe, &
         points_courbe_equal_ref
   real(pr), dimension(:, :, :), allocatable :: slicecontrolm, slicecontrolf, slicecontrol, slicemid, slicecontroli, &
         slicecontroltmpf, slicecontroltmpi, slicecontroltmp
   real(pr), dimension(:, :, :), allocatable :: slice_courbe, slice_courbe_equal
   integer :: Ns, Ni, Nf, nt, errorl, errorr, itail, ihead
   real(pr) :: deltal, deltar, disttail, disthead, rhead, rtail, disttaillY, disttailrY, distheadlY, distheadrY
   real(pr), dimension(:, :, :), allocatable :: slice, slice2, vect, slicetmp
   integer, dimension(:), allocatable :: thetatab, indextheta
   integer, dimension(:, :), allocatable :: indextab
   real(pr), dimension(:), allocatable :: xTheta, yTheta
   real(pr), dimension(:, :), allocatable :: valDist, valTheta, valTh, valThtmp
   integer :: th, theta, nbtheta, boolskel
   real(pr) :: PI, sigma, sigma2 
   real(pr) :: cosPhi, sinPhi, cosTheta, sinTheta, oldS, oldC, cosAlpha, sinAlpha, cosPhil, cosPhir, sinPhil, sinPhir
   real(pr), dimension(:, :), allocatable :: cosTheta_tab, sinTheta_tab
   real(pr), dimension(:, :, :), allocatable :: rhoSlices
   real(pr) :: zslice, area, meshRatio
   real(pr) :: x1, x2, y1, y2, z1, z2, xc, yc, zc, delta, rt, xt, yt, xp, yp, x0, xg, yg, xgref, ygref
   real(pr) :: LS1, LS2, LS3, LS4, LSp, LSrr, rp, alpha, alphadef
   integer :: ip1, ip2, jp1, jp2
   real(pr), dimension(:, :, :), allocatable :: slice_courbemidLarge, slicecontrolLarge
   real(pr), dimension(:, :), allocatable :: valThLarge, valThetaLarge
   integer :: nbthetaLarge, NsLarge, n_frames
   real :: start_time, end_time
   integer, dimension(3) :: file_ids
   character(len=256) :: target_folder, subfolder, results_folder

   !**
   !! Main Initialisation
   !**

   call cpu_time(start_time)

   target_folder = "D:/data_modelisation/results3D"
   subfolder = "second_test"      ! CHANGE DESTINATION SUBFOLDER HERE

   target_folder = trim(trim(adjustl(target_folder))//"/"//trim(adjustl(subfolder)))
   results_folder = trim(trim(target_folder)//'/surf')

   call create_directory(result_folder)

   n_frames = 580
   n_dim = 3
   meshRatio = 0.24_pr
   meshRatio = 1._pr 
   area = 1._pr 
   PI = acos(-1.0_pr)
   sigma = 1._pr
   sigma2 = 1._pr
   !N =  200 
   !nx = 200 
   !ny = 200 
   !nz = 200
   nx = 1602 
   ny = 300 
   nz = 300
   zslice = 150._pr
   x0 = 1._pr
   eepsilon = 1.e-6_pr
   dt = 1._pr
   threshold = 0.001
   !dx = 1._pr 
   !dy = 1._pr 
   !dz = 1._pr
   dx = 2.4 !24 !2.4!*0.000001_pr 
   dy = 2.25171 !24.7687 !2.25171!*0.000001_pr 
   dz = 2.25171 !24.7687 !2.25171!*0.000001_pr
   Tps = 1._pr
   !  Ns = 290 
   !  Ni = 22  
   !  Nf = 30  
   !Ns = 180 
   Ns = 263 
   Ni = 4  
   !  Nf = 54  
   !Nf = 28  
   Nf = 35  
   dtb = 1._pr/(Ns-1)
   nbtheta = 180 
   !nbtheta = 270 
   NsLarge = 900
   nbthetaLarge = 180!200 !3600

   !**
   !! Main Allocation of arrays
   !**
   allocate(cosTheta_tab(Ns+Ni+Nf-2, nbtheta), sinTheta_tab(Ns+Ni+Nf-2, nbtheta))
   allocate(slicemid(nbtheta, Ns, 3))
   allocate(longslice(nbtheta), longslicei(nbtheta), longslicef(nbtheta), sslice(nbtheta), dsslice(nbtheta), stheta(Ns+Ni+Nf-2)&
         , dstheta(Ns+Ni+Nf-2), sth(nx), dsth(nx))
   allocate(thetatab(nbtheta), indextheta(Ns+Ni+Nf-2), indextab(nbtheta, Ns+Ni+Nf-2))
   allocate(xTheta(nbtheta))
   allocate(yTheta(nbtheta))
   allocate(valTheta(Ns+Ni+Nf-2, nbtheta), longTheta(Ns+Ni+Nf-2), longTh(nx))
   allocate(valDist(Ns+Ni+Nf-2, nbTheta))
   allocate(head_courbe(Nf, 3), tail_courbe(Ni, 3))
   allocate(points_courbe(Ns, 3))
   allocate(points_courbe_equal(Ns, 3))
   allocate(points_courbe_equal_ref(Ns, 3))
   allocate(rhoSlices2(nx, ny), gradPhi(nx, ny))
   allocate(tmp(nx, ny), tmp1(nx, ny), tmp2(nx, ny), tmpbool(nx, ny))
   allocate(rhoSlices(nz, nx, ny))
   allocate(xx(nx))
   allocate(yy(ny)) 
   allocate(zz(nz)) 
   allocate(dir1(nx, ny), dir2(nx, ny), dir3(nx, ny), dir4(nx, ny), Nseed(nx, ny), skel(nx, ny), skel2(nx, ny), skel3(nx, ny))
   allocate(distslice(2*(Ns+Ni+Nf), 3))
   allocate(slice_courbe(nbtheta, Ns+Ni+Nf-2, 3), slice_courbemidLarge(nbtheta, NsLarge+Ni+Nf-2, 3))
   allocate(slice_courbe_equal(nbtheta, Ns+Ni+Nf-2, 3))
   allocate(slice(nbtheta, Ns+Ni+Nf-2, 3), slicetmp(nbtheta, Ns+Ni+Nf-2, 3))
   allocate(slice2(nbtheta, Ns+Ni+Nf-2, 3), vect(nbtheta, Ns+Ni+Nf-2, 3))
   allocate(valThetaLarge(Ns+Ni+Nf-2, nbthetaLarge))

   x0 = 0._pr
   do i=1, nx
      xx(i) = x0 + (float(i)-1)*dx
   enddo
   do j=1, ny
      yy(j) = x0 + (float(j)-1)*dy
   enddo
   do k=1, nz
      zz(k) = x0 + (float(k)-1)*dz
   enddo
   do l=1, size(valTheta, 1)
      do theta=1, nbtheta/4+1
         t = (theta-1)/(nbtheta*0.25_pr)
         valTheta(l, theta) = (nbtheta/4+1-1)*2*PI/nbtheta*t
      enddo
      do theta=nbtheta/4+1, nbtheta/2+1
         t = (theta-(nbtheta/4+1))/(nbtheta*0.5_pr-nbtheta*0.25_pr)
         valTheta(l, theta) = (nbtheta/4+1-1)*2*PI/nbtheta*(1-t) + (nbtheta/2+1-1)*2*PI/nbtheta*t
      enddo
      do theta=nbtheta/2+1, nbtheta*3/4+1
         t = (theta-(nbtheta/2+1))/(nbtheta*0.75_pr-nbtheta*0.5_pr)
         valTheta(l, theta) = (nbtheta/2+1-1)*2*PI/nbtheta*(1-t) + (nbtheta*3/4+1-1)*2*PI/nbtheta*t
      enddo
      do theta=nbtheta*3/4+1, nbtheta
         t = (theta-(nbtheta*3/4+1))/(nbtheta-(nbtheta*0.75_pr+1._pr))
         valTheta(l, theta) = (nbtheta*3/4+1-1)*2*PI/nbtheta*(1-t) + (nbtheta-1)*2*PI/nbtheta*t
      enddo
   enddo
   open(unit=78, file='D:\Users\Théo\GitHub\Repositories\deformation3d\fish_shapes\3Dshape.dat', status='unknown')
   do k=1, nx
      do j=1, ny
         do i=1, nz
            read(78, *) rhoSlices(i, k, j)
         enddo
      enddo
   enddo
   close(78)
   do k=1, nx
      rhoSlices(:, k, :) = rhoSlices(:, k, :)/maxval(rhoSlices(:, k, :))
   enddo

   rhoSlices = 10*rhoSlices
   do i=1, nz
      do k=1, nx
         do j=1, ny
            if (rhoSlices(i, k, j)>1._pr) then
               rhoSlices(i, k, j) = 1._pr
            endif
         enddo
      enddo
   enddo
   rhoSlices = rhoSlices - 0.5_pr*(maxval(rhoSlices)+minval(rhoSlices))

   !! FILM
   picNum = 1
   rhoSlices2(:, :) = rhoSlices(nint(zslice), :, :)
   rhoSlices2 = rhoSlices2 - 0.5_pr*(maxval(rhoSlices2) + minval(rhoSlices2))
 
   dy = 2.4 
   dx = 2.25171 
   dz = 2.25171 

   call cpu_time(end_time)
   write(*, *) "******Initialisation and shape loading time = ", end_time-start_time, " sec******"  ! Around unknown sec of exec

   !**
   !! Now, the level-set 3D is computed in rhoSlices
   !**
   call cpu_time(start_time)

   call updateDistance3D(rhoSlices, gradPhi, zslice)
   !**
   !! Now, the level-set 2D is computed in rhoSlices2
   !**
   call updateDistanceINI(rhoSlices2, gradPhi)
   dx = 2.4*0.000001_pr 
   dy = 2.25171*0.000001_pr 
   dz = 2.25171*0.000001_pr
   do i=1, nx
      xx(i) = x0 + (float(i)-1)*dx
   enddo
   do j=1, ny
      yy(j) = x0 + (float(j)-1)*dy
   enddo
   do k=1, nz
      zz(k) = x0 + (float(k)-1)*dz
   enddo

   open(unit=79, file=trim(target_folder)//'/sol00.vtk', status='unknown')
   file_id = 79
   n_dim = 2
   call write_vtk_header(file_id, nx, ny, nz, dx, dy, dz, n_dim)


   !  do k=1, nx
      do j=1, ny
         do i=1, nz
            !write(79, *) rhoSlices(i, 1596, j)
            write(79, *) rhoSlices(i, 6, j)
         enddo
      enddo
   !  enddo
   close(79)   
   
   call cpu_time(end_time)
   write(*, *) "******Shape processing through updateDistance3D and updateDistanceINI and sol00 saving time= ",&
         end_time-start_time, " sec******" ! Around unknown sec of exec

   !**
   !! Definition of the zebrafish length
   !**

   call cpu_time(start_time)

   long00 = 3.8647606990128360E-003 + 1.3658546530224907E-005
   longexp = 164*0.001_pr*0.0256
   longratio = longexp/long00
   write(*, *) "longueurs initiales ", long00, " ", longexp, " ", longratio

   !**
   !! Construction of the initial midline
   !**

   dir1 = 0
   dir2 = 0
   dir3 = 0
   dir4 = 0
   Nseed = 0
   skel = 0
   do i=2, ny-1
      do j=2, nx-1
         if ((gradPhi(j, i)<0.74).and.(rhoSlices2(j, i)>0._pr)) skel(j, i) = 1
      enddo
   enddo

   tmpbool = 0
   nl = sum(skel)
   sizeSkel = sum(skel)
   allocate(midlinebis(sizeSkel, 3))
   l=1
   midlinebis = 0
   do j=100, nx
      do k=1, ny
         if ((skel(j, k)==1).and.(l==1)) then
            midlinebis(l, 2) = k
            if (j==nx/2) l0 = l
            l = l+1
         endif
      enddo
   enddo
   l=1
   do j=1, nx
      do k=1, ny
         if ((skel(j, k)==1).and.(l==1)) then
            midlinebis(l, 1) = j
            if (j==nx/2) l0 = l
            l = l+1
         endif
      enddo
   enddo
   midlinebis(l, 3) = nint(zslice)
   !**
   !! first point of the midline (tail)
   !**
   xskelL = midlinebis(1, 1)
   yskelL = midlinebis(1, 2)
   lp = 1
   boolskel=0

   !do i=2, nx-1
   do i=midlinebis(1, 1)+1, nx-1
      !**
      !! horizontal midline
      !**
      if (lp<=1462-midlinebis(1, 1)+1) then
         if (lp<20) write(*, *) "MIDLINEBIS : ", midlinebis(lp, 1), " ", midlinebis(lp, 2)
         lp=lp+1
         midlinebis(lp, 1) = i
         midlinebis(lp, 2) = midlinebis(1, 2)
      endif
   enddo
   midlinebis(lp, 3) = nint(zslice)
   nl = lp
   !**
   !! Now we have the midline coordinates
   !**
   allocate(midline(nl, 3))
   midline(1:nl, 1) = midlinebis(1:nl, 1)
   midline(1:nl, 2) = midlinebis(1:nl, 2)
   midline(1:nl, 3) = midlinebis(1:nl, 3)

   !**
   !! Definition of the control points of the midline
   !**
   allocate(points_control(nl, 3))
   do l=1, nl
      points_control(l, 1) = xx(midline(l, 1))
      points_control(l, 2) = yy(midline(1, 2))
      points_control(l, 3) = zz(nint(zslice))
   enddo


   !**
   !! Spline approximation of the midline
   !**
   points_courbe(1, 1) = points_control(1, 1)
   points_courbe(1, 2) = points_control(1, 2)
   points_courbe(1, 3) = points_control(1, 3)

   tb = dtb
   l = 1
   long = 0._pr
   do while ((tb<1._pr).and.(l+1<Ns+1))
      l = l+1
      call pointsBezierN3D(points_control, tb, px, py, pz)
      points_courbe(l, 1) = px
      points_courbe(l, 2) = py
      points_courbe(l, 3) = pz
      long = long + sqrt((points_courbe(l, 1)-points_courbe(l-1, 1))**2 + (points_courbe(l, 2)-points_courbe(l-1, 2))**2 +&
            (points_courbe(l, 3)-points_courbe(l-1, 3))**2)
      tb = tb+dtb
   enddo
   if (l==Ns-1) then
      l = l+1
      write(*, *) "ET VOILA"
      points_courbe(Ns, 1) = points_control(size(points_control, 1), 1)
      points_courbe(Ns, 2) = points_control(size(points_control, 1), 2)
      points_courbe(Ns, 3) = points_control(size(points_control, 1), 3)
      long = long + sqrt((points_courbe(l, 1)-points_courbe(l-1, 1))**2 + (points_courbe(l, 2)-points_courbe(l-1, 2))**2 +&
            (points_courbe(l, 3)-points_courbe(l-1, 3))**2)
   endif
   if (.not.(l==Ns)) write(*, *) "WARNIINNGG ", l
   write(*, *) "longueur ==  ", (points_control(size(points_control, 1), 1)-points_control(1, 1)), " ", long, " ", &
         (points_courbe(Ns, 1) - points_courbe(1, 1)), " ", points_courbe(1, 1), " ", points_courbe(Ns, 1)
   long2 = long
   ds = long/(Ns-1)
   !ds = 0.1*long/(0.2*Ns-1)
   !ds = area*long/(meshRatio*Ns-1)

   !**
   !! Uniform spline approximation of the midline
   !**
   points_courbe_equal(1, 1) = points_control(1, 1)
   points_courbe_equal(1, 2) = points_control(1, 2)
   points_courbe_equal(1, 3) = points_control(1, 3)
   l = 1
   long = 0._pr
   dtb = 1._pr/(Ns-1)
   tb = 0._pr
   !do while ((tb<1._pr).and.(l+1<meshRatio*Ns+1))
   !do while ((tb<1._pr).and.(l+1<Ns+1))
   do while ((tb<1._pr).and.(l+1<Ns))
      l = l+1
      nt = 1
      s = 0._pr
      do while ((l-1)*ds-s>0._pr) 
         nt = nt+1
         s = s + sqrt((points_courbe(nt, 1)-points_courbe(nt-1, 1))**2 + (points_courbe(nt, 2)-points_courbe(nt-1, 2))**2 +&
               (points_courbe(nt, 3)-points_courbe(nt-1, 3))**2)
      enddo
      tbm = (nt-2)*dtb
      tbp = (nt-1)*dtb
      tb = tbm !(tbm + tbp)*0.5_pr
      s0 = s
      sinit = s - sqrt((points_courbe(nt, 1)-points_courbe(nt-1, 1))**2 + (points_courbe(nt, 2)-points_courbe(nt-1, 2))**2 +&
            (points_courbe(nt, 3)-points_courbe(nt-1, 3))**2)
      s = sinit
      bool = 0
      !do while ((abs((l-1)*ds-s)>eepsilon*dx).and.(bool==0)) !.and.((tbm+tbp)*0.5_pr<1._pr))
      do while ((abs((l-1)*ds-s)/dx>eepsilon).and.(bool==0)) !.and.((tbm+tbp)*0.5_pr<1._pr))
         tb = (tbm + tbp)*0.5_pr
         !        call pointsBezierN(points_control, tb, px, py)
         !call pointsBezierN3D(points_control(1:nint(0.1*size(points_control, 1)), :), 0.1*tb, px, py, pz)
         call pointsBezierN3D(points_control, tb, px, py, pz)
         s = sinit + sqrt((px-points_courbe(nt-1, 1))**2 + (py-points_courbe(nt-1, 2))**2 + (pz-points_courbe(nt-1, 3))**2)
         if ((l-1)*ds-s>0._pr) then
            tbm = tb
         else
            tbp = tb
            if (tbp>1._pr) tbp = 1._pr
         endif
         if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then
            bool = 1
            !           write(*, *) "aie aie aie"
         endif
      enddo
      !     call pointsBezierN(points_control, tb, px, py)
      write(*, *) "TBpoints  ", l, " ", tb, " ", ds, " ", tbm, " ", tbp
      !call pointsBezierN3D(points_control(1:nint(0.1*size(points_control, 1)), :), 0.1*tb, px, py, pz)
      call pointsBezierN3D(points_control, tb, px, py, pz)
      points_courbe_equal(l, 1) = px
      points_courbe_equal(l, 2) = py
      points_courbe_equal(l, 3) = pz
      long = long +&
            sqrt((points_courbe_equal(l, 1)-points_courbe_equal(l-1, 1))**2 + &
            (points_courbe_equal(l, 2)-points_courbe_equal(l-1, 2))**2 &
            + (points_courbe_equal(l, 3)-points_courbe_equal(l-1, 3))**2)
      write(*, *) "LLfirst  ", l, "       ", px, " ", py, " ", pz, " ", ds
   enddo
   lp = l

      if (l==Ns-1) then
         l = l+1
         write(*, *) "ET VOILA1"
   points_courbe_equal(Ns, 1) = points_control(size(points_control, 1), 1)
   points_courbe_equal(Ns, 2) = points_control(size(points_control, 1), 2)
   points_courbe_equal(Ns, 3) = points_control(size(points_control, 1), 3)
         long = long + sqrt((points_courbe(l, 1)-points_courbe(l-1, 1))**2 + (points_courbe(l, 2)-points_courbe(l-1, 2))**2 +&
            (points_courbe(l, 3)-points_courbe(l-1, 3))**2)
      write(*, *) "LLfirst  ", l, " ", points_courbe_equal(Ns, 1), " ", points_courbe_equal(Ns, 2), " "&
      , points_courbe_equal(Ns, 3), " "&
      , sqrt((points_courbe(l, 1)-points_courbe(l-1, 1))**2 + (points_courbe(l, 2)-points_courbe(l-1, 2))**2 +&
                  (points_courbe(l, 3)-points_courbe(l-1, 3))**2), " ", ds
      endif
   if (.not.(l==Ns)) write(*, *) "WARNIINNGG2 ", l
   !if (.not.(l+lp-1==Ns)) write(*, *) "WARNIINNGG2 ", l
   write(*, *) "longueur =  ", (points_control(size(points_control, 1), 1)-points_control(1, 1)), " ", long, " ", &
         (points_courbe_equal(Ns, 1)-points_courbe_equal(1, 1)), "  LL  ", l
   long3 = 0._pr
   do l=1, nl-1
      long3 = long3 +sqrt((xx(midline(l+1, 1))-xx(midline(l, 1)))**2 + (yy(midline(l+1, 2))-yy(midline(l, 2)))**2)
   enddo
   write(*, *) "longueur =  ", long3

   !**
   !! Interpolation of endpoints of the real midline 
   !**
   rtail = xx(1)
   itail = floor((rtail-x0)/dx+1)
   rhead = xx(nx)
   ihead = ceiling((rhead-x0)/dx+1)
   bool1 = 0
   bool2 = 0
   if (bool1==0) write(*, *) rhoSlices2(1, midline(1, 2))
   do i=2, nx-1
      if (bool1==0) write(*, *) rhoSlices2(i, midline(1, 2))
      if ((rhoSlices2(i, midline(1, 2))>0._pr).and.(rhoSlices2(i-1, midline(1, 2))<0._pr).and.(bool1==0)) then 
         itail = i
         rtail = xx(i-1) - rhoSlices2(i, midline(1, 2))*dx/(rhoSlices2(i-1, midline(1, 2)) - rhoSlices2(i, midline(1, 2)))
         itail = i
         bool1 = 1
         write(*, *) "benvoila  ", rtail, " ", itail, " ", i
      endif
      if ((rhoSlices2(i, midline(1, 2))>0._pr).and.(rhoSlices2(i+1, midline(1, 2))<0._pr)) then
         ihead = i
         rhead = xx(i) - rhoSlices2(i, midline(1, 2))*dx/(rhoSlices2(i+1, midline(1, 2)) - rhoSlices2(i, midline(1, 2)))
         ihead = ceiling((rhead-x0)/dx+1)
         bool2 = 1
      endif
   enddo
   rhead = rhead + 0.0000425/longratio
   ihead = ihead-1
   rtail = rtail - 0.000025/longratio
   disttail = abs(rtail - points_courbe_equal(1, 1))

   disthead = abs(rhead-points_courbe_equal(size(points_courbe_equal, 1), 1))
   write(*, *) "longueur =  ", rhead-rtail, " ", rhead-rtail-disttail, " ", rhead-rtail-disthead, " ", rhead-rtail-disttail-disthead
   long00 = rhead-rtail-disttail-disthead
   write(*, *) "longueur ", long3, " ", long2, " ", long, " ", long00
   write(*, *) "DIST : TAIL ", disttail, " HEAD ", disthead, " NINT ", nint(disttail/dx), " ", nint(disthead/dx), &
          " ", ceiling(disthead/dx), " ", floor(disthead/dx), " ", rtail, " ", rhead&
         , " ", itail, " ", ihead, " ", ceiling(points_courbe_equal(1, 1)/dx), " ", floor(points_courbe_equal(Ns, 1)/dx)&
         , " ", points_courbe_equal(1, 1), " ", points_courbe_equal(1, 1)/dx!, " nbtheta ", nbtheta

   !**
   !! Allocation of control arrays 
   !**

   allocate(slicecontrolm(nbtheta, nint(points_courbe_equal(Ns, 1)/dx)-nint(points_courbe_equal(1, 1)/dx)+1, 3))
   allocate(slicecontrolf(nbtheta, nint(disthead/dx)+1, 3))
   allocate(slicecontroltmpf(nbtheta, nint(&
   abs(rhead-0.000045/longratio-points_courbe_equal(size(points_courbe_equal, 1), 1))/dx)+1, 3))
   allocate(slicecontroli(nbtheta, nint(disttail/dx)+2, 3))
   allocate(slicecontroltmpi(nbtheta, nint(abs(rtail + 0.000025/longratio - points_courbe_equal(1, 1))/dx)+2, 3))
   allocate(slicecontrol(nbtheta, ihead-itail+1, 3), slicecontroltmp(nbtheta, ihead-itail+1, 3)) 
   allocate(slicecontrolLarge(nbthetaLarge, ihead-itail+1, 3)) !nbthetaLarge, nint(points_courbe_equal(Ns, 1)+disthead-points_courbe_equal(1, 1)-disttail), 3))
   allocate(valTh(ihead-itail+1, nbtheta), valThLarge(ihead-itail+1, nbthetaLarge))
   allocate(valThtmp(nbtheta, ihead-itail+1))
   write(*, *) "SIZE CONTROLLARGE  ", size(slicecontrolLarge, 2), " ", points_courbe_equal(Ns, 1)+&
         disthead-points_courbe_equal(1, 1) - disttail, " ", rhead, " ", rtail

   do l=1, size(valTh, 1)
      do theta=1, nbtheta
         valTh(l, theta) = (theta-1)*2*PI/nbtheta
         valThtmp(theta, l) = (theta-1)*2*PI/nbtheta
      enddo
   enddo

   !**
   !! Filling control arrays (based on level-set zero searching - dichotomy) 
   !**
   lp = 0
   lph = 0
   lpf = 0
   lpt = 0
   lpt = size(slicecontroli, 2)-size(slicecontroltmpi, 2)
   write(*, *) "SIZEC  ", size(slicecontroltmpi, 2), " ", size(slicecontroli, 2), " ", itail
   longTh = 0._pr
   do l=itail, ihead
      do theta=1, nbtheta
         rr = 0._pr
         bool = 0
         do while (((LSrr<0._pr).and.(rr<yy(ny))).or.(bool==0))
            rp = rr
            if (bool==1) then
               LSp = LSrr
            else
               if (cos(valTh(l-itail+1, theta))<0._pr) then
                  xp = (-yy(1)-rp+points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
               else
                  xp = (yy(ny)-rp-points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
               endif
               if (sin(valTh(l-itail+1, theta))<0._pr) then
                  yp = (-zz(1)-rp+zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
               else
                  yp = (zz(nz)-rp-zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
               endif

               ip1 = int((xp-x0)/dy)+1
               jp1 = int((yp-x0)/dz)+1              
               if (cos(valTh(l-itail+1, theta))<0._pr) then
                  ip2 = ip1+1
               else
                  ip2 = ip1-1
               endif
               if (sin(valTh(l-itail+1, theta))<0._pr) then
                  jp2 = jp1+1
               else
                  jp2 = jp1-1
               endif

               if (jp2==0) write(*, *) "ATTENTIONSTOPZERO  ", l, " ", theta, " ", jp1, " "&
               , cos(valTh(l-itail+1, theta)), " ", sin(valTh(l-itail+1, theta)), " ", (xp-x0)/dy, " ", (yp-x0)/dz, " ", rp, " ", rr
               if (ip2==301) write(*, *) "ATTENTIONSTOPMAX  ", l, " ", theta, " ", jp1, " "&
               , cos(valTh(l-itail+1, theta)), " ", sin(valTh(l-itail+1, theta)), " ", (xp-x0)/dy, " ", (yp-x0)/dz, " ", rp, " ", rr

               LS1 = rhoSlices(jp1, l, ip1)
               LS2 = rhoSlices(jp1, l, ip2)
               LS3 = rhoSlices(jp2, l, ip1)
               LS4 = rhoSlices(jp2, l, ip2)

               x1 = yy(ip1)
               y1 = zz(jp1)

               if (cos(valTh(l-itail+1, theta))<0._pr) then
                  x2 = x1 + dy
               else
                  x2 = x1 - dy
               endif
               if (sin(valTh(l-itail+1, theta))<0._pr) then
                  y2 = y1 + dz
               else
                  y2 = y1 - dz
               endif

               LSp = interpLS(LS1, x1, y1, LS2, x2, y1, LS3, x1, y2, LS4, x2, y2, xp, yp) 
               if (l==itail+1) write(*, *) "TEST checkAVANT  ::    ", l, " ", theta, " ", rp, " ", LSp, " ", LSrr, "  ", &
                  xp, " ", yp, " ", (xp-x0)/dy, " ", (yp-x0)/dz, " ", ip1, " ", jp1&
                  , "  ", cos(valTh(l-itail+1, theta)), " ", sin(valTh(l-itail+1, theta)), "   ", yy(ny-2)-yy(ny/2)
            endif

            bool = 1
            rr=rr+dy

            if (cos(valTh(l-itail+1, theta))<0._pr) then
               xp = (-yy(1)-rr+points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
            else
               xp = (yy(ny)-rr-points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
            endif
            if (sin(valTh(l-itail+1, theta))<0._pr) then
               yp = (-zz(1)-rr+zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
            else
               yp = (zz(nz)-rr-zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
            endif

            ip1 = int((xp-x0)/dy)+1
            jp1 = int((yp-x0)/dz)+1
               if (l==itail+1) write(*, *) "TEST checkAVANT  ::    ", l, " ", theta, " ", rr, " ", LSp, " ", LSrr, "  ", & 
               xp, " ", yp, " ", (xp-x0)/dy, " ", (yp-x0)/dz, " ", ip1, " ", jp1&
               , "  ", cos(valTh(l-itail+1, theta)), " ", sin(valTh(l-itail+1, theta)), "   ", yy(ny-2)-yy(ny/2)&
               , " ", int((points_courbe_equal(1, 2)-x0)/dy+1), " ", int((zz(nint(zslice))-x0)/dz+1)

            if (cos(valTh(l-itail+1, theta))<0._pr) then
               ip2 = ip1+1
            else
               ip2 = ip1-1
            endif
            if (sin(valTh(l-itail+1, theta))<0._pr) then
               jp2 = jp1+1
            else
               jp2 = jp1-1
            endif
               if (jp2==0) write(*, *) "ATTENTIONSTOPZERO  ", l, " ", theta, " ", jp1, " "&
               , cos(valTh(l-itail+1, theta)), " ", sin(valTh(l-itail+1, theta)), " ", (xp-x0)/dy, " ", (yp-x0)/dz, " ", rp, " ", rr

            LS1 = rhoSlices(jp1, l, ip1)
            LS2 = rhoSlices(jp1, l, ip2)
            LS3 = rhoSlices(jp2, l, ip1)
            LS4 = rhoSlices(jp2, l, ip2)

            x1 = yy(ip1)
            y1 = zz(jp1)

            if (cos(valTh(l-itail+1, theta))<0._pr) then
               x2 = x1 + dy
            else
               x2 = x1 - dy
            endif
            if (sin(valTh(l-itail+1, theta))<0._pr) then
               y2 = y1 + dz
            else
               y2 = y1 - dz
            endif

            LSrr = interpLS(LS1, x1, y1, LS2, x2, y1, LS3, x1, y2, LS4, x2, y2, xp, yp) 
         enddo
         rr = rp
         LSrr = LSp
         do while ((LSrr<0._pr).and.(rr<yy(ny)))
            rp = rr
            LSp = LSrr

            rr=rr+0.1*dy

            if (cos(valTh(l-itail+1, theta))<0._pr) then
               xp = (-yy(1)-rr+points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
            else
               xp = (yy(ny)-rr-points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
            endif
            if (sin(valTh(l-itail+1, theta))<0._pr) then
               yp = (-zz(1)-rr+zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
            else
               yp = (zz(nz)-rr-zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
            endif
            ip1 = int((xp-x0)/dy)+1
            jp1 = int((yp-x0)/dz)+1

            if (cos(valTh(l-itail+1, theta))<0._pr) then
               ip2 = ip1+1
            else
               ip2 = ip1-1
            endif
            if (sin(valTh(l-itail+1, theta))<0._pr) then
               jp2 = jp1+1
            else
               jp2 = jp1-1
            endif

            LS1 = rhoSlices(jp1, l, ip1)
            LS2 = rhoSlices(jp1, l, ip2)
            LS3 = rhoSlices(jp2, l, ip1)
            LS4 = rhoSlices(jp2, l, ip2)
            x1 = yy(ip1)
            y1 = zz(jp1)

            if (cos(valTh(l-itail+1, theta))<0._pr) then
               x2 = x1 + dy
            else
               x2 = x1 - dy
            endif
            if (sin(valTh(l-itail+1, theta))<0._pr) then
               y2 = y1 + dz
            else
               y2 = y1 - dz
            endif

            LSrr = interpLS(LS1, x1, y1, LS2, x2, y1, LS3, x1, y2, LS4, x2, y2, xp, yp) 
         enddo
         px = xx(l)

         if (cos(valTh(l-itail+1, theta))<0._pr) then
            py = (-yy(1)-rp+points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2) +&
                  -cos(valTh(l-itail+1, theta))*LSP*(rr-rp)/abs(LSrr - LSp)  !rr+1
         else
            py = (yy(ny)-rp-points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2) +&
                  -cos(valTh(l-itail+1, theta))*LSP*(rr-rp)/abs(LSrr - LSp)  !rr+1
         endif
         if (sin(valTh(l-itail+1, theta))<0._pr) then
            pz = (-zz(1)-rp+zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice)) +&
                  -sin(valTh(l-itail+1, theta))*LSP*(rr-rp)/abs(LSrr - LSp)  !rr+1
         else
            pz = (zz(nz)-rp-zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))  +&
                  -sin(valTh(l-itail+1, theta))*LSP*(rr-rp)/abs(LSrr - LSp)  !rr+1
         endif

         if (l<=nint((points_courbe_equal(1, 1)-x0)/dx+1)) then
            if (lpt==size(slicecontroli, 2)-size(slicecontroltmpi, 2)) lpt = l-lpt

            slicecontroli(theta, l-lpt+2, 1) = px
            slicecontroli(theta, l-lpt+2, 2) = py
            slicecontroli(theta, l-lpt+2, 3) = pz
            if (theta==1+nbtheta/2) write(*, *) "totototo00  ", l, " ", px, " ", py, " ", pz&
            , " ", ceiling((points_courbe_equal(1, 1)-x0)/dx+1), " ", lpt
         endif
         if ((l>=nint((points_courbe_equal(1, 1)-x0)/dx+1)).and.(l<=nint((points_courbe_equal(Ns, 1)-x0)/dx+1))) then
            if (lp==0) lpf = l
            if (lp==0) lp = l
            slicecontrolm(theta, l-lp+1, 1) = px
            slicecontrolm(theta, l-lp+1, 2) = py                
            slicecontrolm(theta, l-lp+1, 3) = pz
         endif
         if (l>=nint((points_courbe_equal(Ns, 1)-x0)/dx+1)) then
            if (lph==0) lph = l
            slicecontrolf(theta, l-lph+1, 1) = px
            slicecontrolf(theta, l-lph+1, 2) = py
            slicecontrolf(theta, l-lph+1, 3) = pz
         endif
         slicecontrol(theta, l-itail+1, 1) = px
         slicecontrol(theta, l-itail+1, 2) = py
         slicecontrol(theta, l-itail+1, 3) = pz
         if (theta>1) longTh(l) = longTh(l) + dist(slicecontrol(theta, l-itail+1, 2), slicecontrol(theta, l-itail+1, 3)&
               , slicecontrol(theta-1, l-itail+1, 2), slicecontrol(theta-1, l-itail+1, 3))
         if (theta==nbtheta) longTh(l) = longTh(l) + dist(slicecontrol(theta, l-itail+1, 2)&
               , slicecontrol(theta, l-itail+1, 3), slicecontrol(1, l-itail+1, 2), slicecontrol(1, l-itail+1, 3))
                  if (theta==1) write(*, *) "LongTH  ", l, " ", longTh(l), " ", LSrr, " ", LSp
      enddo
      if (lp==0) then
      else
         lpf = lpf+1
      endif

   enddo
   slicecontroltmp = slicecontrol
   do l=itail, ihead
      dsth(l) = longTh(l)/nbtheta
   enddo
   dtb = 2*PI/nbtheta

   do l=itail, ihead
      ll = 1
      tb = 0._pr
      do while ((tb<2*PI).and.(ll+1<nbtheta+1))
         ll = ll+1
         nt = 1
         sth(l) = 0._pr
         do while (((ll-1)*dsth(l)-sth(l)>0._pr).and.(nt<size(slicecontroltmp, 1)))
            nt = nt+1
            sth(l) = sth(l) +&
                  sqrt((slicecontroltmp(nt, l-itail+1, 1)-slicecontroltmp(nt-1, l-itail+1, 1))**2 + &
                  (slicecontroltmp(nt, l-itail+1, 2)-slicecontroltmp(nt-1, l-itail+1, 2))**2 + &
                  (slicecontroltmp(nt, l-itail+1, 3)-slicecontroltmp(nt-1, l-itail+1, 3))**2)
         enddo

         tbm = (nt-2)*dtb
         tbp = (nt-1)*dtb

         tb = tbm 

         s0 = sth(l)
         if (nt>1) then
         sinit = sth(l) -&
               sqrt((slicecontroltmp(nt, l-itail+1, 1)-slicecontroltmp(nt-1, l-itail+1, 1))**2 + &
               (slicecontroltmp(nt, l-itail+1, 2)-slicecontroltmp(nt-1, l-itail+1, 2))**2 + &
               (slicecontroltmp(nt, l-itail+1, 3)-slicecontroltmp(nt-1, l-itail+1, 3))**2)
         else
         sinit = sth(l)
         endif
         sth(l) = sinit
         bool = 0
         do while ((abs((ll-1)*dsth(l)-sth(l))/dx>eepsilon).and.(bool==0))
            tb = (tbm + tbp)*0.5_pr
            rr = 0._pr
            booltmp = 0
            do while (((LSrr<0._pr).and.(rr<yy(ny))).or.(booltmp==0))
               rp = rr
               if (booltmp==1) then
                  LSp = LSrr
               else
                  if (cos(tb)<0._pr) then
                     xp = (-yy(1)-rp+points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2)
                  else
                     xp = (yy(ny)-rp-points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2)
                  endif
                  if (sin(tb)<0._pr) then
                     yp = (-zz(1)-rp+zz(nint(zslice)))*sin(tb)+zz(nint(zslice))
                  else
                     yp = (zz(nz)-rp-zz(nint(zslice)))*sin(tb)+zz(nint(zslice))
                  endif
                  ip1 = int((xp-x0)/dy)+1
                  jp1 = int((yp-x0)/dz)+1
                  if (cos(tb)<0._pr) then
                     ip2 = ip1+1
                  else
                     ip2 = ip1-1
                  endif
                  if (sin(tb)<0._pr) then
                     jp2 = jp1+1
                  else
                     jp2 = jp1-1
                  endif
                  LS1 = rhoSlices(jp1, l, ip1)
                  LS2 = rhoSlices(jp1, l, ip2)
                  LS3 = rhoSlices(jp2, l, ip1)
                  LS4 = rhoSlices(jp2, l, ip2)
                  x1 = yy(ip1)
                  y1 = zz(jp1)

                  if (cos(tb)<0._pr) then
                     x2 = x1 + dy
                  else
                     x2 = x1 - dy
                  endif
                  if (sin(tb)<0._pr) then
                     y2 = y1 + dz
                  else
                     y2 = y1 - dz
                  endif
                  LSp = interpLS(LS1, x1, y1, LS2, x2, y1, LS3, x1, y2, LS4, x2, y2, xp, yp) 

               endif
               booltmp = 1
               rr=rr+dy

               if (cos(tb)<0._pr) then
                  xp = (-yy(1)-rr+points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2)
               else
                  xp = (yy(ny)-rr-points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2)
               endif
               if (sin(tb)<0._pr) then
                  yp = (-zz(1)-rr+zz(nint(zslice)))*sin(tb)+zz(nint(zslice))
               else
                  yp = (zz(nz)-rr-zz(nint(zslice)))*sin(tb)+zz(nint(zslice))
               endif
               ip1 = int((xp-x0)/dy)+1
               jp1 = int((yp-x0)/dz)+1

               if (cos(tb)<0._pr) then
                  ip2 = ip1+1
               else
                  ip2 = ip1-1
               endif
               if (sin(tb)<0._pr) then
                  jp2 = jp1+1
               else
                  jp2 = jp1-1
               endif
               LS1 = rhoSlices(jp1, l, ip1)
               LS2 = rhoSlices(jp1, l, ip2)
               LS3 = rhoSlices(jp2, l, ip1)
               LS4 = rhoSlices(jp2, l, ip2)
               x1 = yy(ip1)
               y1 = zz(jp1)

               if (cos(tb)<0._pr) then
                  x2 = x1 + dy
               else
                  x2 = x1 - dy
               endif
               if (sin(tb)<0._pr) then
                  y2 = y1 + dz
               else
                  y2 = y1 - dz
               endif
               LSrr = interpLS(LS1, x1, y1, LS2, x2, y1, LS3, x1, y2, LS4, x2, y2, xp, yp) 
            enddo
            rr = rp
            LSrr = LSp
            do while ((LSrr<0._pr).and.(rr<yy(ny)))
               rp = rr
               LSp = LSrr
               rr=rr+0.1*dy

               if (cos(tb)<0._pr) then
                  xp = (-yy(1)-rr+points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2)
               else
                  xp = (yy(ny)-rr-points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2)
               endif
               if (sin(tb)<0._pr) then
                  yp = (-zz(1)-rr+zz(nint(zslice)))*sin(tb)+zz(nint(zslice))
               else
                  yp = (zz(nz)-rr-zz(nint(zslice)))*sin(tb)+zz(nint(zslice))
               endif
               ip1 = int((xp-x0)/dy)+1
               jp1 = int((yp-x0)/dz)+1

               if (cos(tb)<0._pr) then
                  ip2 = ip1+1
               else
                  ip2 = ip1-1
               endif
               if (sin(tb)<0._pr) then
                  jp2 = jp1+1
               else
                  jp2 = jp1-1
               endif
               LS1 = rhoSlices(jp1, l, ip1)
               LS2 = rhoSlices(jp1, l, ip2)
               LS3 = rhoSlices(jp2, l, ip1)
               LS4 = rhoSlices(jp2, l, ip2)
               x1 = yy(ip1)
               y1 = zz(jp1)
               if (cos(tb)<0._pr) then
                  x2 = x1 + dy
               else
                  x2 = x1 - dy
               endif
               if (sin(tb)<0._pr) then
                  y2 = y1 + dz
               else
                  y2 = y1 - dz
               endif
               LSrr = interpLS(LS1, x1, y1, LS2, x2, y1, LS3, x1, y2, LS4, x2, y2, xp, yp) 
            enddo
            px = xx(l)

            if (cos(tb)<0._pr) then
               py = (-yy(1)-rp+points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2) +&
                     -cos(tb)*LSP*(rr-rp)/abs(LSrr - LSp)
            else
               py = (yy(ny)-rp-points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2) +&
                     -cos(tb)*LSP*(rr-rp)/abs(LSrr - LSp)
            endif
            if (sin(tb)<0._pr) then
               pz = (-zz(1)-rp+zz(nint(zslice)))*sin(tb)+zz(nint(zslice)) +&
                     -sin(tb)*LSP*(rr-rp)/abs(LSrr - LSp)
            else
               pz = (zz(nz)-rp-zz(nint(zslice)))*sin(tb)+zz(nint(zslice))  +&
                     -sin(tb)*LSP*(rr-rp)/abs(LSrr - LSp)
            endif

            sth(l) = sinit + sqrt((px-slicecontroltmp(nt-1, l-itail+1, 1))**2 + (py-slicecontroltmp(nt-1, l-itail+1, 2))**2 +&
                  (pz-slicecontroltmp(nt-1, l-itail+1, 3))**2)
            if ((ll-1)*dsth(l)-sth(l)>0._pr) then
               tbm = tb
            else
               tbp = tb
            endif
            if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then
               bool = 1
            endif
         enddo
         rr = 0._pr
         booltmp = 0
         do while (((LSrr<0._pr).and.(rr<yy(ny))).or.(booltmp==0))
            rp = rr
            if (booltmp==1) then
               LSp = LSrr
            else
               if (cos(tb)<0._pr) then
                  xp = (-yy(1)-rp+points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2)
               else
                  xp = (yy(ny)-rp-points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2)
               endif
               if (sin(tb)<0._pr) then
                  yp = (-zz(1)-rp+zz(nint(zslice)))*sin(tb)+zz(nint(zslice))
               else
                  yp = (zz(nz)-rp-zz(nint(zslice)))*sin(tb)+zz(nint(zslice))
               endif
               ip1 = int((xp-x0)/dy)+1
               jp1 = int((yp-x0)/dz)+1

               if (cos(tb)<0._pr) then
                  ip2 = ip1+1
               else
                  ip2 = ip1-1
               endif
               if (sin(tb)<0._pr) then
                  jp2 = jp1+1
               else
                  jp2 = jp1-1
               endif
               LS1 = rhoSlices(jp1, l, ip1)
               LS2 = rhoSlices(jp1, l, ip2)
               LS3 = rhoSlices(jp2, l, ip1)
               LS4 = rhoSlices(jp2, l, ip2)
               x1 = yy(ip1)
               y1 = zz(jp1)

               if (cos(tb)<0._pr) then
                  x2 = x1 + dy
               else
                  x2 = x1 - dy
               endif
               if (sin(tb)<0._pr) then
                  y2 = y1 + dz
               else
                  y2 = y1 - dz
               endif
               LSp = interpLS(LS1, x1, y1, LS2, x2, y1, LS3, x1, y2, LS4, x2, y2, xp, yp) 
            endif
            booltmp = 1
            rr=rr+dy

            if (cos(tb)<0._pr) then
               xp = (-yy(1)-rr+points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2)
            else
               xp = (yy(ny)-rr-points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2)
            endif
            if (sin(tb)<0._pr) then
               yp = (-zz(1)-rr+zz(nint(zslice)))*sin(tb)+zz(nint(zslice))
            else
               yp = (zz(nz)-rr-zz(nint(zslice)))*sin(tb)+zz(nint(zslice))
            endif
            ip1 = int((xp-x0)/dy)+1
            jp1 = int((yp-x0)/dz)+1

            if (cos(tb)<0._pr) then
               ip2 = ip1+1
            else
               ip2 = ip1-1
            endif
            if (sin(tb)<0._pr) then
               jp2 = jp1+1
            else
               jp2 = jp1-1
            endif
            LS1 = rhoSlices(jp1, l, ip1)
            LS2 = rhoSlices(jp1, l, ip2)
            LS3 = rhoSlices(jp2, l, ip1)
            LS4 = rhoSlices(jp2, l, ip2)
            x1 = yy(ip1)
            y1 = zz(jp1)

            if (cos(tb)<0._pr) then
               x2 = x1 + dy
            else
               x2 = x1 - dy
            endif
            if (sin(tb)<0._pr) then
               y2 = y1 + dz
            else
               y2 = y1 - dz
            endif
            LSrr = interpLS(LS1, x1, y1, LS2, x2, y1, LS3, x1, y2, LS4, x2, y2, xp, yp) 
         enddo
         rr = rp
         LSrr = LSp
         do while ((LSrr<0._pr).and.(rr<yy(ny)))
            rp = rr
            LSp = LSrr
            rr=rr+0.1*dy

            if (cos(tb)<0._pr) then
               xp = (-yy(1)-rr+points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2)
            else
               xp = (yy(ny)-rr-points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2)
            endif
            if (sin(tb)<0._pr) then
               yp = (-zz(1)-rr+zz(nint(zslice)))*sin(tb)+zz(nint(zslice))
            else
               yp = (zz(nz)-rr-zz(nint(zslice)))*sin(tb)+zz(nint(zslice))
            endif
            ip1 = int((xp-x0)/dy)+1
            jp1 = int((yp-x0)/dz)+1

            if (cos(tb)<0._pr) then
               ip2 = ip1+1
            else
               ip2 = ip1-1
            endif
            if (sin(tb)<0._pr) then
               jp2 = jp1+1
            else
               jp2 = jp1-1
            endif
            LS1 = rhoSlices(jp1, l, ip1)
            LS2 = rhoSlices(jp1, l, ip2)
            LS3 = rhoSlices(jp2, l, ip1)
            LS4 = rhoSlices(jp2, l, ip2)
            x1 = yy(ip1)
            y1 = zz(jp1)

            if (cos(tb)<0._pr) then
               x2 = x1 + dy
            else
               x2 = x1 - dy
            endif
            if (sin(tb)<0._pr) then
               y2 = y1 + dz
            else
               y2 = y1 - dz
            endif
            LSrr = interpLS(LS1, x1, y1, LS2, x2, y1, LS3, x1, y2, LS4, x2, y2, xp, yp) 
            if (ll==1+nbtheta/2)   write(*, *) "totototo00  ", l, " ", px, " ", py, " ", pz, " ", LSrr, &
            " ", LSp, " ", rr, " ", yy(ny)
         enddo
         px = xx(l)

         if (cos(tb)<0._pr) then
            py = (-yy(1)-rp+points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2) +&
                  -cos(tb)*LSP*(rr-rp)/abs(LSrr - LSp) 
         else
            py = (yy(ny)-rp-points_courbe_equal(1, 2))*cos(tb)+points_courbe_equal(1, 2) +&
                  -cos(tb)*LSP*(rr-rp)/abs(LSrr - LSp)
         endif
         if (sin(tb)<0._pr) then
            pz = (-zz(1)-rp+zz(nint(zslice)))*sin(tb)+zz(nint(zslice)) +&
                  -sin(tb)*LSP*(rr-rp)/abs(LSrr - LSp)
         else
            pz = (zz(nz)-rp-zz(nint(zslice)))*sin(tb)+zz(nint(zslice))  +&
                  -sin(tb)*LSP*(rr-rp)/abs(LSrr - LSp)
         endif

         if (l<=nint((points_courbe_equal(1, 1)-x0)/dx+1)) then
            slicecontroli(ll, l-lpt+2, 1) = px 
            slicecontroli(ll, l-lpt+2, 2) = py 
            slicecontroli(ll, l-lpt+2, 3) = pz 
            if (ll==1+nbtheta/2) write(*, *) "totototo1  ", l, " ", ll, " ", px, " ", py, " ", pz&
            , " ", ceiling((points_courbe_equal(1, 1)-x0)/dx+1), " ", l-lpt+1

         endif
         if ((l>=nint((points_courbe_equal(1, 1)-x0)/dx+1)).and.(l<=nint((points_courbe_equal(Ns, 1)-x0)/dx+1))) then

            slicecontrolm(ll, l-lp+1, 1) = px 
            slicecontrolm(ll, l-lp+1, 2) = py 
            slicecontrolm(ll, l-lp+1, 3) = pz 
            if (ll==1+nbtheta/2) write(*, *) "totototo2  ", l, " ", ll, " ", px, " ", py, " ", pz&
            , " ", ceiling((points_courbe_equal(1, 1)-x0)/dx+1), " ", l-lp+1
         endif
         if (l>=nint((points_courbe_equal(Ns, 1)-x0)/dx+1)) then

            slicecontrolf(ll, l-lph+1, 1) = px 
            slicecontrolf(ll, l-lph+1, 2) = py 
            slicecontrolf(ll, l-lph+1, 3) = pz 
         endif

         slicecontrol(ll, l-itail+1, 1) = px 
         slicecontrol(ll, l-itail+1, 2) = py 
         slicecontrol(ll, l-itail+1, 3) = pz 
         valTh(l-itail+1, ll) = tb
         valThtmp(ll, l-itail+1) = tb

      enddo

      if (.not.(ll==nbtheta)) write(*, *) "WARNIINNGG c est moiTHETA  ", ll, " ", l, " ", tb
   enddo
   do theta=1, nbtheta
      l=itail+1
      rr = 0._pr
      bool = 0
      do while (((LSrr<0._pr).and.(rr<float(ny))).or.(bool==0))
         rp = rr
         if (bool==1) then
            LSp = LSrr
         else

            if (cos(valTh(l-itail+1, theta))<0._pr) then
               xp = (-yy(1)-rp+points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
            else
               xp = (yy(ny)-rp-points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
            endif
            if (sin(valTh(l-itail+1, theta))<0._pr) then
               yp = (-zz(1)-rp+zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
            else
               yp = (zz(nz)-rp-zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
            endif
            ip1 = int((xp-x0)/dy)+1
            jp1 = int((yp-x0)/dz)+1

            if (cos(valTh(l-itail+1, theta))<0._pr) then
               ip2 = ip1+1
            else
               ip2 = ip1-1
            endif
            if (sin(valTh(l-itail+1, theta))<0._pr) then
               jp2 = jp1+1
            else
               jp2 = jp1-1
            endif
            LS1 = rhoSlices(jp1, l, ip1)
            LS2 = rhoSlices(jp1, l, ip2)
            LS3 = rhoSlices(jp2, l, ip1)
            LS4 = rhoSlices(jp2, l, ip2)
            x1 = yy(ip1)
            y1 = zz(jp1)

            if (cos(valTh(l-itail+1, theta))<0._pr) then
               x2 = x1 + dy
            else
               x2 = x1 - dy
            endif
            if (sin(valTh(l-itail+1, theta))<0._pr) then
               y2 = y1 + dz
            else
               y2 = y1 - dz
            endif
            LSp = interpLS(LS1, x1, y1, LS2, x2, y1, LS3, x1, y2, LS4, x2, y2, xp, yp) 
            if (theta==1) write(*, *) "checknantail  ", x1, " ", y1, " ", LSrr, " ", LSp
         endif

         bool = 1
         rr=rr+dy

         if (cos(valTh(l-itail+1, theta))<0._pr) then
            xp = (-yy(1)-rr+points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
         else
            xp = (yy(ny)-rr-points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
         endif
         if (sin(valTh(l-itail+1, theta))<0._pr) then
            yp = (-zz(1)-rr+zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
         else
            yp = (zz(nz)-rr-zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
         endif
         ip1 = int((xp-x0)/dy)+1
         jp1 = int((yp-x0)/dz)+1

         if (cos(valTh(l-itail+1, theta))<0._pr) then
            ip2 = ip1+1
         else
            ip2 = ip1-1
         endif
         if (sin(valTh(l-itail+1, theta))<0._pr) then
            jp2 = jp1+1
         else
            jp2 = jp1-1
         endif
         LS1 = rhoSlices(jp1, l, ip1)
         LS2 = rhoSlices(jp1, l, ip2)
         LS3 = rhoSlices(jp2, l, ip1)
         LS4 = rhoSlices(jp2, l, ip2)
         x1 = yy(ip1)
         y1 = zz(jp1)

         if (cos(valTh(l-itail+1, theta))<0._pr) then
            x2 = x1 + dy
         else
            x2 = x1 - dy
         endif
         if (sin(valTh(l-itail+1, theta))<0._pr) then
            y2 = y1 + dz
         else
            y2 = y1 - dz
         endif
         LSrr = interpLS(LS1, x1, y1, LS2, x2, y1, LS3, x1, y2, LS4, x2, y2, xp, yp) 
      enddo
      rr = rp
      LSrr = LSp
      do while ((LSrr<0._pr).and.(rr<yy(ny)))
         rp = rr
         LSp = LSrr

         rr=rr+0.1*dy

         if (cos(valTh(l-itail+1, theta))<0._pr) then
            xp = (-yy(1)-rr+points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
         else
            xp = (yy(ny)-rr-points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
         endif
         if (sin(valTh(l-itail+1, theta))<0._pr) then
            yp = (-zz(1)-rr+zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
         else
            yp = (zz(nz)-rr-zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
         endif
         ip1 = int((xp-x0)/dy)+1
         jp1 = int((yp-x0)/dz)+1

         if (cos(valTh(l-itail+1, theta))<0._pr) then
            ip2 = ip1+1
         else
            ip2 = ip1-1
         endif
         if (sin(valTh(l-itail+1, theta))<0._pr) then
            jp2 = jp1+1
         else
            jp2 = jp1-1
         endif
         LS1 = rhoSlices(jp1, l, ip1)
         LS2 = rhoSlices(jp1, l, ip2)
         LS3 = rhoSlices(jp2, l, ip1)
         LS4 = rhoSlices(jp2, l, ip2)
         x1 = yy(ip1)
         y1 = zz(jp1)

         if (cos(valTh(l-itail+1, theta))<0._pr) then
            x2 = x1 + dy
         else
            x2 = x1 - dy
         endif
         if (sin(valTh(l-itail+1, theta))<0._pr) then
            y2 = y1 + dz
         else
            y2 = y1 - dz
         endif
         LSrr = interpLS(LS1, x1, y1, LS2, x2, y1, LS3, x1, y2, LS4, x2, y2, xp, yp) 
      enddo

      px = rtail

      py = 0.01_pr*rp*cos(valTh(l-itail+1, theta))+(slicecontroltmp(1, 2, 2) + slicecontroltmp(nint(1+0.5_pr*nbtheta), &
            2, 2))*0.5_pr + 0.01_pr*cos(valTh(l-itail+1, theta))*LSP*(rr-rp)/abs(LSrr - LSp) !rr+1

      pz = 0.01_pr*rp*sin(valTh(l-itail+1, theta)) + (slicecontroltmp(nint(1+0.25_pr*nbtheta), 2, 3)+&
            slicecontroltmp(nint(1+0.75_pr*nbtheta), 2, 3))*0.5_pr&
            + 0.01_pr*sin(valTh(l-itail+1, theta))*LSp*(rr-rp)/abs(LSrr - LSp) !rr+1

      slicecontroli(theta, 1, 1) = px 
      slicecontroli(theta, 1, 2) = py 
      slicecontroli(theta, 1, 3) = pz 
      slicecontrol(theta, 1, 1) = px 
      slicecontrol(theta, 1, 2) = py 
      slicecontrol(theta, 1, 3) = pz 
      if (theta==1) write(*, *) "checknantail  ", slicecontrol(theta, 1, 2), " ", slicecontrol(theta, 1, 3), " ", LSrr, " ", LSp
      l=ihead-1
      rr = 0._pr
      bool = 0
      do while (((LSrr<0._pr).and.(rr<yy(ny))).or.(bool==0))
         rp = rr
         if (bool==1) then
            LSp = LSrr
         else

            if (cos(valTh(l-itail+1, theta))<0._pr) then
               xp = (-yy(1)-rp+points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
            else
               xp = (yy(ny)-rp-points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
            endif
            if (sin(valTh(l-itail+1, theta))<0._pr) then
               yp = (-zz(1)-rp+zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
            else
               yp = (zz(nz)-rp-zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
            endif
            ip1 = int((xp-x0)/dy)+1
            jp1 = int((yp-x0)/dz)+1

            if (cos(valTh(l-itail+1, theta))<0._pr) then
               ip2 = ip1+1
            else
               ip2 = ip1-1
            endif
            if (sin(valTh(l-itail+1, theta))<0._pr) then
               jp2 = jp1+1
            else
               jp2 = jp1-1
            endif
            LS1 = rhoSlices(jp1, l, ip1)
            LS2 = rhoSlices(jp1, l, ip2)
            LS3 = rhoSlices(jp2, l, ip1)
            LS4 = rhoSlices(jp2, l, ip2)
            x1 = yy(ip1)
            y1 = zz(jp1)

            if (cos(valTh(l-itail+1, theta))<0._pr) then
               x2 = x1 + dy
            else
               x2 = x1 - dy
            endif
            if (sin(valTh(l-itail+1, theta))<0._pr) then
               y2 = y1 + dz
            else
               y2 = y1 - dz
            endif
            LSp = interpLS(LS1, x1, y1, LS2, x2, y1, LS3, x1, y2, LS4, x2, y2, xp, yp) 
         endif

         bool = 1
         rr=rr+dy

         if (cos(valTh(l-itail+1, theta))<0._pr) then
            xp = (-yy(1)-rr+points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
         else
            xp = (yy(ny)-rr-points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
         endif
         if (sin(valTh(l-itail+1, theta))<0._pr) then
            yp = (-zz(1)-rr+zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
         else
            yp = (zz(nz)-rr-zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
         endif
         ip1 = int((xp-x0)/dy)+1
         jp1 = int((yp-x0)/dz)+1

         if (cos(valTh(l-itail+1, theta))<0._pr) then
            ip2 = ip1+1
         else
            ip2 = ip1-1
         endif
         if (sin(valTh(l-itail+1, theta))<0._pr) then
            jp2 = jp1+1
         else
            jp2 = jp1-1
         endif
         LS1 = rhoSlices(jp1, l, ip1)
         LS2 = rhoSlices(jp1, l, ip2)
         LS3 = rhoSlices(jp2, l, ip1)
         LS4 = rhoSlices(jp2, l, ip2)
         x1 = yy(ip1)
         y1 = zz(jp1)

         if (cos(valTh(l-itail+1, theta))<0._pr) then
            x2 = x1 + dy
         else
            x2 = x1 - dy
         endif
         if (sin(valTh(l-itail+1, theta))<0._pr) then
            y2 = y1 + dz
         else
            y2 = y1 - dz
         endif
         LSrr = interpLS(LS1, x1, y1, LS2, x2, y1, LS3, x1, y2, LS4, x2, y2, xp, yp) 
      enddo
      rr = rp
      LSrr = LSp
      do while ((LSrr<0._pr).and.(rr<yy(ny)))
         rp = rr
         LSp = LSrr

         rr=rr+0.1*dy

         if (cos(valTh(l-itail+1, theta))<0._pr) then
            xp = (-yy(1)-rr+points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
         else
            xp = (yy(ny)-rr-points_courbe_equal(1, 2))*cos(valTh(l-itail+1, theta))+points_courbe_equal(1, 2)
         endif
         if (sin(valTh(l-itail+1, theta))<0._pr) then
            yp = (-zz(1)-rr+zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
         else
            yp = (zz(nz)-rr-zz(nint(zslice)))*sin(valTh(l-itail+1, theta))+zz(nint(zslice))
         endif
         ip1 = int((xp-x0)/dy)+1
         jp1 = int((yp-x0)/dz)+1

         if (cos(valTh(l-itail+1, theta))<0._pr) then
            ip2 = ip1+1
         else
            ip2 = ip1-1
         endif
         if (sin(valTh(l-itail+1, theta))<0._pr) then
            jp2 = jp1+1
         else
            jp2 = jp1-1
         endif
         LS1 = rhoSlices(jp1, l, ip1)
         LS2 = rhoSlices(jp1, l, ip2)
         LS3 = rhoSlices(jp2, l, ip1)
         LS4 = rhoSlices(jp2, l, ip2)
         x1 = yy(ip1)
         y1 = zz(jp1)

         if (cos(valTh(l-itail+1, theta))<0._pr) then
            x2 = x1 + dy
         else
            x2 = x1 - dy
         endif
         if (sin(valTh(l-itail+1, theta))<0._pr) then
            y2 = y1 + dz
         else
            y2 = y1 - dz
         endif
         LSrr = interpLS(LS1, x1, y1, LS2, x2, y1, LS3, x1, y2, LS4, x2, y2, xp, yp) 
      enddo

      px = rhead

      py = 0.01_pr*rp*cos(valTh(l-itail+1, theta)) + (slicecontroltmp(1, size(slicecontroltmp, 2)-1, 2)+&
            slicecontroltmp(nint(1+0.5_pr*nbtheta), size(slicecontroltmp, 2)-1, 2))*0.5_pr&
            + 0.01_pr*cos(valTh(l-itail+1, theta))*LSP*(rr-rp)/abs(LSrr - LSp) 
      pz = 0.01_pr*rp*sin(valTh(l-itail+1, theta)) + (slicecontroltmp(nint(1+0.25_pr*nbtheta), size(slicecontroltmp, 2)-1, 3)+&
            slicecontroltmp(nint(1+0.75_pr*nbtheta), size(slicecontroltmp, 2)-1, 3))*0.5_pr&
            + 0.01_pr*sin(valTh(l-itail+1, theta))*LSp*(rr-rp)/abs(LSrr - LSp) 
      slicecontrolf(theta, size(slicecontrolf, 2), 1) = px 
      slicecontrolf(theta, size(slicecontrolf, 2), 2) = py 
      slicecontrolf(theta, size(slicecontrolf, 2), 3) = pz 
      slicecontrol(theta, size(slicecontrol, 2), 1) = px 
      slicecontrol(theta, size(slicecontrol, 2), 2) = py 
      slicecontrol(theta, size(slicecontrol, 2), 3) = pz 
   enddo

   !**
   !! Spline approximation of contour/surface control points
   !**
   dtb = 1._pr/(Ni-1)
   tb = dtb
   do theta=1, nbtheta

      slice_courbe(theta, 1, 1) = slicecontroli(theta, 1, 1)
      slice_courbe(theta, 1, 2) = slicecontroli(theta, 1, 2)
      slice_courbe(theta, 1, 3) = slicecontroli(theta, 1, 3)
      valTheta(1, theta) = valTh(1, theta)
   enddo

   dsi = disttail/(Ni-1)
   dsf = disthead/(Nf-1)
   do theta=1, nbtheta
      tb = 0._pr
      if (theta==1) write(*, *) "XX  ", nint((points_courbe_equal(1, 1)-x0)/dx+1)&
      , " ", floor((points_courbe_equal(1, 1)-x0)/dx+1)&
      , " ", ceiling((points_courbe_equal(1, 1)-x0)/dx+1)
      if (theta==1) write(*, *) "XX  ", size(slicecontroltmpi, 2)&
      , " ", size(slicecontroli, 2)&
      , " ", itail-lpt+2&
      , " ", nint((points_courbe_equal(1, 1)-x0)/dx+1)-lpt+2
      slicecontroltmpi(theta, 2:size(slicecontroltmpi, 2), :) = slicecontroli(theta, itail-lpt+2:&
      nint((points_courbe_equal(1, 1)-x0)/dx+1)-lpt+2, :)
      slicecontroltmpi(theta, 1, :) = slicecontroli(theta, 1, :)
      if (theta==1) write(*, *) "XX  ", slicecontroli(theta, :, 1)
      if (theta==1) write(*, *) "XXTMP  ", slicecontroltmpi(theta, :, 1)
      if (theta==1) write(*, *) "YY  ", slicecontroli(theta, :, 2)
      if (theta==1) write(*, *) "YYTMP  ", slicecontroltmpi(theta, :, 2)
      do l=Ni, 1, -1
         nt = 1
         tbm = 0._pr
         tbp = 1._pr
         tb = tbm 
         bool = 0
         do while (((abs(sslice(theta))/(dx*dx)>eepsilon).and.(bool==0)).or.(nt==1))
            tb = (tbm + tbp)*0.5_pr
            nt = nt+1
            call pointsBezierN3D(slicecontroltmpi(theta, :, :), tb, px, py, pz)
            sslice(theta) = dotProd3D(points_courbe_equal(1, 1)-(Ni-l)*dsi, points_courbe_equal(1, 2), points_courbe_equal(1, 3)&
                  , points_courbe_equal(1, 1)-(Ni-l-1)*dsi, points_courbe_equal(1, 2), points_courbe_equal(1, 3), px, py, pz)

            if (sslice(theta)<0._pr) then
               tbm = tb
            else
               tbp = tb
            endif
            if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then
               bool = 1
            endif
         enddo
         call pointsBezierN3D(slicecontroltmpi(theta, :, :), tb, px, py, pz)
         slice_courbe(theta, l, 1) = px
         slice_courbe(theta, l, 2) = py
         slice_courbe(theta, l, 3) = pz

         if (l<=Ni) then
            call pointsBezierN3D(slicecontrolm(theta, :, :), 5.0258636474609375E-003_pr, pxx, pyy, pzz)
            tbb = (px-rtail)*1._pr/(pxx-rtail)
            tp = (slice_courbe(theta, Ni-1, 1)-slicecontroli(theta, 1, 1))/(pxx-slicecontroli(theta, 1, 1))
            call poly2(slicecontroli(theta, 1, 1), slice_courbe(theta, Ni-1, 1), pxx, tbb, tp, px) 
            call poly2(slicecontroli(theta, 1, 2), slice_courbe(theta, Ni-1, 2), pyy, tbb, tp, py) 
            call poly2(slicecontroli(theta, 1, 3), slice_courbe(theta, Ni-1, 3), pzz, tbb, tp, pz) 
            slice_courbe(theta, l, 1) = px
            slice_courbe(theta, l, 2) = py
            slice_courbe(theta, l, 3) = pz
         endif
      enddo
      nt = 1
      tbm = 0._pr
      tbp = 1._pr
      tb = tbm 
      bool = 0
      do while (((abs(sslice(theta))/(dx*dx)>eepsilon).and.(bool==0)).or.(nt==1))
         tb = (tbm + tbp)*0.5_pr
         nt = nt+1
         call pointsBezierN3D(slicecontroltmpi(theta, :, :), tb, px, py, pz)
         sslice(theta) = dotProd3D(points_courbe_equal(1, 1), points_courbe_equal(1, 2), points_courbe_equal(1, 3)&
               , points_courbe_equal(1, 1)+dsi, points_courbe_equal(1, 2), points_courbe_equal(1, 3), px, py, pz)
         if (theta==1) write(*, *) "SSLICE  ", sslice(theta), "   ", slicecontroltmpi(1, :, 1)

         if (sslice(theta)<0._pr) then
            tbm = tb
         else
            tbp = tb
         endif
         if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then
            bool = 1
         endif
      enddo
      call pointsBezierN3D(slicecontroltmpi(theta, :, :), tb, px, py, pz)
      if (theta==1) write(*, *) "TBsliceI  ", Ni, " ", tb, " ", px, " ", py, " ", pz
            call pointsBezierN3D(slicecontrolm(theta, :, :), 5.0258636474609375E-003_pr, pxx, pyy, pzz)
            tbb = (px-rtail)*1._pr/(pxx-rtail)
            tp = (slice_courbe(theta, Ni-1, 1)-slicecontroli(theta, 1, 1))/(pxx-slicecontroli(theta, 1, 1))
            call poly2(slicecontroli(theta, 1, 1), slice_courbe(theta, Ni-1, 1), pxx, tbb, tp, px) 
            call poly2(slicecontroli(theta, 1, 2), slice_courbe(theta, Ni-1, 2), pyy, tbb, tp, py) 
            call poly2(slicecontroli(theta, 1, 3), slice_courbe(theta, Ni-1, 3), pzz, tbb, tp, pz) 
      slice_courbe(theta, Ni, 1) = px
      slice_courbe(theta, Ni, 2) = py
      slice_courbe(theta, Ni, 3) = pz
   enddo


   dtb = 1._pr/(Ns-1)
   tb = dtb
   longslice = 0._pr

   do theta=1, nbtheta
      tb = 0._pr

      do l=2, Ns-1
         nt = 1
         tbm = 0._pr
         tbp = 1._pr
         tb = tbm 
         bool = 0
         do while (((abs(sslice(theta))/(dx*dx)>eepsilon).and.(bool==0)).or.(nt==1))
            tb = (tbm + tbp)*0.5_pr
            nt = nt+1
            call pointsBezierN3D(slicecontrolm(theta, :, :), tb, px, py, pz)

            if (theta==1) write(*, *) "TBslicem  ", l, " ", tb, " ", px, " ", py, " ", pz
            sslice(theta) = dotProd3D(points_courbe_equal(l, 1), points_courbe_equal(l, 2), points_courbe_equal(l, 3)&
                  , points_courbe_equal(l+1, 1), points_courbe_equal(l+1, 2), points_courbe_equal(l+1, 3), px, py, pz)

            if (sslice(theta)<0._pr) then
               tbm = tb
            else
               tbp = tb
            endif
            if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then
               bool = 1
            endif
         enddo
         call pointsBezierN3D(slicecontrolm(theta, :, :), tb, px, py, pz)
         if (theta==1) write(*, *) "TBslicem  ", l, " ", tb, " ", px, " ", py, " ", pz
         slice_courbe(theta, l+Ni-1, 1) = px
         slice_courbe(theta, l+Ni-1, 2) = py
         slice_courbe(theta, l+Ni-1, 3) = pz
      enddo
      nt = 1
      tbm = 0._pr
      tbp = 1._pr
      tb = tbm 
      bool = 0
      do while (((abs(sslice(theta))/(dx*dx)>eepsilon).and.(bool==0)).or.(nt==1))
         tb = (tbm + tbp)*0.5_pr
         nt = nt+1
         call pointsBezierN3D(slicecontrolm(theta, :, :), tb, px, py, pz)
         sslice(theta) = dotProd3D(points_courbe_equal(Ns, 1), points_courbe_equal(Ns, 2), points_courbe_equal(Ns, 3)&
               , points_courbe_equal(Ns, 1)+dsf, points_courbe_equal(1, 2), points_courbe_equal(1, 3), px, py, pz)

         if (sslice(theta)<0._pr) then
            tbm = tb
         else
            tbp = tb
         endif
         if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then
            bool = 1
         endif
      enddo
      call pointsBezierN3D(slicecontrolm(theta, :, :), tb, px, py, pz)
      if (theta==1) write(*, *) "TBslice  ", Ns, " ", tb, " ", px, " ", py, " ", pz
      slice_courbe(theta, Ns+Ni-1, 1) = px
      slice_courbe(theta, Ns+Ni-1, 2) = py
      slice_courbe(theta, Ns+Ni-1, 3) = pz
   enddo

   dtb = 1._pr/(Nf-1)
   tb = dtb
   longslicef = 0._pr
   do theta=1, nbtheta
      slice_courbe(theta, size(slice_courbe, 2), 1) = slicecontrolf(theta, size(slicecontrolf, 2), 1)
      slice_courbe(theta, size(slice_courbe, 2), 2) = slicecontrolf(theta, size(slicecontrolf, 2), 2)
      slice_courbe(theta, size(slice_courbe, 2), 3) = slicecontrolf(theta, size(slicecontrolf, 2), 3)
   enddo

   do theta=1, nbtheta
      tb = 0._pr
      slicecontroltmpf(theta, 1:ihead-lph, :) = slicecontrolf(theta, 1:ihead-lph, :)
      slicecontroltmpf(theta, size(slicecontroltmpf, 2), :) = slicecontrolf(theta, size(slicecontrolf, 2), :)
      do l=1, Nf-1
         nt = 1
         tbm = 0._pr
         tbp = 1._pr
         tb = tbm 
         bool = 0
         do while (((abs(sslice(theta))/(dx*dx)>eepsilon).and.(bool==0)).or.(nt==1))
            tb = (tbm + tbp)*0.5_pr
            nt = nt+1
            call pointsBezierN3D(slicecontroltmpf(theta, :, :), tb, px, py, pz)
            if (theta==1) write(*, *) "TBsliceF  ", l, " ", tb, " ", px, " ", py, " ", pz
            sslice(theta) = dotProd3D(points_courbe_equal(Ns, 1)+(l-1)*dsf, points_courbe_equal(1, 2), points_courbe_equal(1, 3)&
                  , points_courbe_equal(Ns, 1)+l*dsf, points_courbe_equal(1, 2), points_courbe_equal(1, 3), px, py, pz)

            if (sslice(theta)<0._pr) then
               tbm = tb
            else
               tbp = tb
            endif
            if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then
               bool = 1
            endif
         enddo
         if (theta==1) write(*, *) "TBsliceF  ", l, " ", tb
         call pointsBezierN3D(slicecontroltmpf(theta, :, :), tb, px, py, pz)
         slice_courbe(theta, l+Ni-1+Ns-1, 1) = px
         slice_courbe(theta, l+Ni-1+Ns-1, 2) = py
         slice_courbe(theta, l+Ni-1+Ns-1, 3) = pz
         if (l>Nf-5) then
            pxx = slice_courbe(theta, Nf-7-1+Ni-1+Ns-1, 1)
            pyy = slice_courbe(theta, Nf-7-1+Ni-1+Ns-1, 2)
            pzz = slice_courbe(theta, Nf-7-1+Ni-1+Ns-1, 3)
            tbb = (px-pxx)*1._pr/(rhead-pxx)
            tp = (slice_courbe(theta, Nf-4+Ni-1+Ns-1, 1)-pxx)/(slicecontrolf(theta, size(slicecontrolf, 2), 1)-pxx)
            call poly2(pxx, slice_courbe(theta, Nf-4+Ni-1+Ns-1, 1), slicecontrolf(theta, size(slicecontrolf, 2), 1), tbb, tp, px) 
            call poly2(pyy, slice_courbe(theta, Nf-4+Ni-1+Ns-1, 2), slicecontrolf(theta, size(slicecontrolf, 2), 2), tbb, tp, py) 
            call poly2(pzz, slice_courbe(theta, Nf-4+Ni-1+Ns-1, 3), slicecontrolf(theta, size(slicecontrolf, 2), 3), tbb, tp, pz) 
            slice_courbe(theta, l+Ni-1+Ns-1, 1) = px 
            slice_courbe(theta, l+Ni-1+Ns-1, 2) = py 
            slice_courbe(theta, l+Ni-1+Ns-1, 3) = pz 
         if (theta==1) write(*, *) "AH      ", floor((points_courbe_equal(Ns, 1)-x0)/dx+1), " "&
         , ceiling((rhead-x0)/dx+1), " ", floor((points_courbe_equal(Ns, 1)-x0)/dx+1)
         if (theta==1) write(*, *) "AH      ", nint((points_courbe_equal(Ns, 1)-x0)/dx+1), " "&
         , nint((rhead-x0)/dx+1), " ", nint((points_courbe_equal(Ns, 1)-x0)/dx+1)
         if (theta==1) write(*, *) "TBsliceFIN  ", l, " ", tbb, " ", tbb*tp, " ", tp, " ", px, " ", py, " ", pz, " ", pxx&
         , " ", slice_courbe(theta, Nf-10+Ni-1+Ns-1, 1), " ", slicecontrolf(theta, size(slicecontrolf, 2), 1)
         if (theta==1) write(*, *) "TBsliceFIN  ", l, " ", tbb, " ", tbb*tp, " ", tp, " ", px, " ", py, " ", pz, " ", pyy&
         , " ", slice_courbe(theta, Nf-10+Ni-1+Ns-1, 2), " ", slicecontrolf(theta, size(slicecontrolf, 2), 2)
         if (theta==1) write(*, *) "TBsliceFIN  ", l, " ", tbb, " ", tbb*tp, " ", tp, " ", px, " ", py, " ", pz, " ", pzz&
         , " ", slice_courbe(theta, Nf-10+Ni-1+Ns-1, 3), " ", slicecontrolf(theta, size(slicecontrolf, 2), 3)
         endif
      enddo
   enddo
   l = Nf-1+Ns-1+Ni-1
   if (l==Ni-1+Ns-1+Nf-1) then
      l = l+1
      do theta=1, nbtheta
         slice_courbe(theta, l, 1) = slicecontrolf(theta, size(slicecontrolf, 2), 1)
         slice_courbe(theta, l, 2) = slicecontrolf(theta, size(slicecontrolf, 2), 2)
         slice_courbe(theta, l, 3) = slicecontrolf(theta, size(slicecontrolf, 2), 3)
         call pointsBezierN1D(valThtmp(theta, lph-itail+1:ihead-itail+1)&
               , tb, alpha)
         valTheta(l, theta) = alpha
         longslicef(theta) = longslicef(theta) +&
               sqrt((slice_courbe(theta, l, 1)-slice_courbe(theta, l-1, 1))**2 + (slice_courbe(theta, l, 2) - &
               slice_courbe(theta, l-1, 2))**2 + (slice_courbe(theta, l, 3)-slice_courbe(theta, l-1, 3))**2)
      enddo
      write(*, *) "ET VOILA F "
   endif
   if (.not.(l==Ni-1+Ns-1+Nf)) write(*, *) "WARNIINNGG  F ", l

   call cpu_time(end_time)
   write(*, *) "******Frame 1 processing part 1 time = ", end_time-start_time, " sec******" ! Around 25 sec
 
   !**
   !! Scaling the contour/surface Lagrangian markers
   !**

   call cpu_time(start_time)

   do l=1, Ni-1+Ns-1+Nf
      do theta=1, nbtheta
         slice(theta, l, 1) = (slice_courbe(theta, l, 1)-x0)*longratio+x0
         slice(theta, l, 2) = (slice_courbe(theta, l, 2)-x0)*longratio+x0
         slice(theta, l, 3) = (slice_courbe(theta, l, 3)-x0)*longratio+x0
      enddo
   enddo
   do l=1, size(points_courbe_equal, 1)
      points_courbe_equal(l, 1) = (points_courbe_equal(l, 1)-x0)*longratio+x0
      points_courbe_equal(l, 2) = (points_courbe_equal(l, 2)-x0)*longratio+x0
      points_courbe_equal(l, 3) = (points_courbe_equal(l, 3)-x0)*longratio+x0
   enddo
   do l=1, size(slicecontrol, 2)
      do theta=1, nbtheta
         slicecontrol(theta, l, 1) = (slicecontrol(theta, l, 1)-x0)*longratio+x0
         slicecontrol(theta, l, 2) = (slicecontrol(theta, l, 2)-x0)*longratio+x0
         slicecontrol(theta, l, 3) = (slicecontrol(theta, l, 3)-x0)*longratio+x0
      enddo
   enddo
   long00 = long00*longratio
   disthead = disthead*longratio
   disttail = disttail*longratio
   dsi = disttail/(Ni-1)
   dsf = disthead/(Nf-1)
   dtb = 1._pr/(Ns-1)

   do l =1, Ni+Ns+Nf-2
      do theta=1, nbtheta
         if (theta==nbtheta/2) write(*, *) "FINAL THETA  ", l, " ", theta, " ", valTheta(l, theta)
      enddo
   enddo
   slicetmp = slice

   !**
   !! Computing the distance and angles between each Lagrangian markers and associated midline points
   !**
   distslice(1+Ni, 1) = dist3D(slice(1, 1+Ni-1, 1), slice(1, 1+Ni-1, 2), slice(1, 1+Ni-1, 3)&
         , points_courbe_equal(1, 1), points_courbe_equal(1, 2), points_courbe_equal(1, 3))
   distslice(1+Ni, 2) = dist3D(slice(1, 1+Ni-1, 1), slice(1, 1+Ni-1, 2), slice(1, 1+Ni-1, 3)&
         , points_courbe_equal(2, 1), points_courbe_equal(2, 2), points_courbe_equal(1, 3))
   do theta=1, nbtheta
      valDist(1+Ni-1, theta) = dist3D(slice(theta, 1+Ni-1, 1), slice(theta, 1+Ni-1, 2), slice(theta, 1+Ni-1, 3)&
            , points_courbe_equal(1, 1), points_courbe_equal(1, 2), points_courbe_equal(1, 3))
      sinTheta_tab(1+Ni-1, theta) = (slice(theta, 1+Ni-1, 3)-points_courbe_equal(1, 3))/valDist(1+Ni-1, theta) !(sqrt((slice(theta, 1+Ni-1, 2)-points_courbe_equal(1, 2))**2+(slice(theta, 1+Ni-1, 3)-points_courbe_equal(1, 3))**2))
      cosTheta_tab(1+Ni-1, theta) = (slice(theta, 1+Ni-1, 2)-points_courbe_equal(1, 2))/valDist(1+Ni-1, theta) !(sqrt((slice(theta, 1+Ni-1, 2)-points_courbe_equal(1, 2))**2+(slice(theta, 1+Ni-1, 3)-points_courbe_equal(1, 3))**2)) 
   enddo
   do l=2, Ns-1
      distslice(l+Ni, 1) = dist3D(slice(1, l+Ni-1, 1), slice(1, l+Ni-1, 2), slice(1, l+Ni-1, 3)&
            , points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2), points_courbe_equal(l-1, 3))
      distslice(l+Ni, 2) = dist3D(slice(1, l+Ni-1, 1), slice(1, l+Ni-1, 2), slice(1, l+Ni-1, 3)&
            , points_courbe_equal(l+1, 1), points_courbe_equal(l+1, 2), points_courbe_equal(l+1, 3))
      write(*, *) "DOTPROD  ", l, " ", dotProd(points_courbe_equal(l, 1), points_courbe_equal(l, 2), slice(1, l+Ni-1, 1)&
            , slice(1, l+Ni-1, 2), points_courbe_equal(l+1, 1), points_courbe_equal(l+1, 2))
      do theta=1, nbtheta
         valDist(l+Ni-1, theta) = dist3D(slice(theta, l+Ni-1, 1), slice(theta, l+Ni-1, 2), slice(theta, l+Ni-1, 3)&
               , points_courbe_equal(l, 1), points_courbe_equal(l, 2), points_courbe_equal(l, 3))
         sinTheta_tab(l+Ni-1, theta) = (slice(theta, l+Ni-1, 3)-points_courbe_equal(l, 3))/valDist(l+Ni-1, theta) !(sqrt((slice(theta, l+Ni-1, 2)-points_courbe_equal(l, 2))**2+(slice(theta, l+Ni-1, 3)-points_courbe_equal(l, 3))**2))
         cosTheta_tab(l+Ni-1, theta) = (slice(theta, l+Ni-1, 2)-points_courbe_equal(l, 2))/valDist(l+Ni-1, theta) !(sqrt((slice(theta, l+Ni-1, 2)-points_courbe_equal(l, 2))**2+(slice(theta, l+Ni-1, 3)-points_courbe_equal(l, 3))**2)) 
      enddo
   enddo
   distslice(Ns+Ni, 1) = dist3D(slice(1, Ns+Ni-1, 1), slice(1, Ns+Ni-1, 2), slice(1, Ns+Ni-1, 3)&
         , points_courbe_equal(Ns-1, 1), points_courbe_equal(Ns-1, 2), points_courbe_equal(Ns-1, 3))
   distslice(Ns+Ni, 2) = dist3D(slice(1, Ns+Ni-1, 1), slice(1, Ns+Ni-1, 2), slice(1, Ns+Ni-1, 3)&
         , points_courbe_equal(Ns, 1), points_courbe_equal(Ns, 2), points_courbe_equal(Ns, 3))
   do theta=1, nbtheta
      valDist(Ns+Ni-1, theta) =&
            dist3D(slice(theta, Ns+Ni-1, 1), slice(theta, Ns+Ni-1, 2), slice(theta, Ns+Ni-1, 3)&
            , points_courbe_equal(Ns, 1), points_courbe_equal(Ns, 2), points_courbe_equal(Ns, 3))
      sinTheta_tab(Ns+Ni-1, theta) = (slice(theta, Ns+Ni-1, 3)-points_courbe_equal(Ns, 3))/valDist(Ns+Ni-1, theta) !(sqrt((slice(theta, Ns+Ni-1, 2)-points_courbe_equal(Ns, 2))**2+(slice(theta, Ns+Ni-1, 3)-points_courbe_equal(Ns, 3))**2))
      cosTheta_tab(Ns+Ni-1, theta) = (slice(theta, Ns+Ni-1, 2)-points_courbe_equal(Ns, 2))/valDist(Ns+Ni-1, theta) !(sqrt((slice(theta, Ns+Ni-1, 2)-points_courbe_equal(Ns, 2))**2+(slice(theta, Ns+Ni-1, 3)-points_courbe_equal(Ns, 3))**2)) 
   enddo

   dsi = disttail/(Ni-1)
   distslice(1, 1) = dist3D(slice(1, 1, 1), slice(1, 1, 2), slice(1, 1, 3)&
         , points_courbe_equal(1, 1)-(Ni-1)*dsi, points_courbe_equal(1, 2), points_courbe_equal(1, 3))
   distslice(1, 2) = dist3D(slice(1, 1, 1), slice(1, 1, 2), slice(1, 1, 3)&
         , points_courbe_equal(1, 1)-(Ni-2)*dsi, points_courbe_equal(1, 2), points_courbe_equal(1, 3))
   do theta=1, nbtheta 
      valDist(1, theta) = dist3D(slice(theta, 1, 1), slice(theta, 1, 2), slice(theta, 1, 3)&
            , points_courbe_equal(1, 1)-(Ni-1)*dsi, points_courbe_equal(1, 2), points_courbe_equal(1, 3))
      sinTheta_tab(1, theta) = (slice(theta, 1, 3)-points_courbe_equal(1, 3))/valDist(1, theta) !(sqrt((slice(theta, 1, 2)-points_courbe_equal(1, 2))**2+(slice(theta, 1, 3)-points_courbe_equal(1, 3))**2))
      cosTheta_tab(1, theta) = (slice(theta, 1, 2)-points_courbe_equal(1, 2))/valDist(1, theta) !(sqrt((slice(theta, l+Ni-1, 2)-points_courbe_equal(l, 2))**2+(slice(theta, l+Ni-1, 3)-points_courbe_equal(l, 3))**2)) 
   enddo
   do l=2, Ni-1
      distslice(l, 1) = dist3D(slice(1, l, 1), slice(1, l, 2), slice(1, l, 3)&
            , points_courbe_equal(1, 1)-(Ni-l+1)*dsi, points_courbe_equal(1, 2), points_courbe_equal(1, 3))
      distslice(l, 2) = dist3D(slice(1, l, 1), slice(1, l, 2), slice(1, l, 3)&
            , points_courbe_equal(1, 1)-(Ni-l-1)*dsi, points_courbe_equal(1, 2), points_courbe_equal(1, 3))
      write(*, *) "DOTPRODi  ", l, " ", dotProd(points_courbe_equal(1, 1)-(Ni-l)*dsi, points_courbe_equal(1, 2), slice(1, l, 1)&
            , slice(1, l, 2), points_courbe_equal(1, 1)-(Ni-l-1)*dsi, points_courbe_equal(1, 2))
      do theta=1, nbtheta
         valDist(l, theta) = dist3D(slice(theta, l, 1), slice(theta, l, 2), slice(theta, l, 3)&
               , points_courbe_equal(1, 1)-(Ni-l)*dsi, points_courbe_equal(1, 2), points_courbe_equal(1, 3))
         sinTheta_tab(l, theta) = (slice(theta, l, 3)-points_courbe_equal(1, 3))/valDist(l, theta) 
         cosTheta_tab(l, theta) = (slice(theta, l, 2)-points_courbe_equal(1, 2))/valDist(l, theta)  
      enddo
   enddo
   distslice(Ni, 1) = dist3D(slice(1, Ni, 1), slice(1, Ni, 2), slice(1, Ni, 3)&
         , points_courbe_equal(1, 1)-dsi, points_courbe_equal(1, 2), points_courbe_equal(1, 3))
   distslice(Ni, 2) = dist3D(slice(1, Ni, 1), slice(1, Ni, 2), slice(1, Ni, 3)&
         , points_courbe_equal(1, 1), points_courbe_equal(1, 2), points_courbe_equal(1, 3))
   do theta=1, nbtheta
      valDist(Ni, theta) = dist3D(slice(theta, Ni, 1), slice(theta, Ni, 2), slice(theta, Ni, 3)&
            , points_courbe_equal(1, 1), points_courbe_equal(1, 2), points_courbe_equal(1, 3))
      sinTheta_tab(Ni, theta) = (slice(theta, Ni, 3)-points_courbe_equal(1, 3))/valDist(Ni, theta) 
      cosTheta_tab(Ni, theta) = (slice(theta, Ni, 2)-points_courbe_equal(1, 2))/valDist(Ni, theta)  
   enddo

   dsf = disthead/(Nf-1)
   distslice(1+Ni+Ns, 1) = dist3D(slice(1, 1+Ni-1+Ns-1, 1), slice(1, 1+Ni-1+Ns-1, 2), slice(1, 1+Ni-1+Ns-1, 3)&
         , points_courbe_equal(Ns, 1), points_courbe_equal(1, 2), points_courbe_equal(1, 3))
   distslice(1+Ni+Ns, 2) = dist3D(slice(1, 1+Ni-1+Ns-1, 1), slice(1, 1+Ni-1+Ns-1, 2), slice(1, 1+Ni-1+Ns-1, 3)&
         , points_courbe_equal(Ns, 1)+dsf, points_courbe_equal(1, 2), points_courbe_equal(1, 3))
   do theta=1, nbtheta
      valDist(1+Ni-1+Ns-1, theta) =&
            dist3D(slice(theta, 1+Ni-1+Ns-1, 1), slice(theta, 1+Ni-1+Ns-1, 2), slice(theta, 1+Ni-1+Ns-1, 3)&
            , points_courbe_equal(Ns, 1), points_courbe_equal(1, 2), points_courbe_equal(1, 3))
      sinTheta_tab(1+Ni-1+Ns-1, theta) = (slice(theta, 1+Ni-1+Ns-1, 3)-points_courbe_equal(1, 3))/valDist(1+Ni-1+Ns-1, theta) 
      cosTheta_tab(1+Ni-1+Ns-1, theta) = (slice(theta, 1+Ni-1+Ns-1, 2)-points_courbe_equal(1, 2))/valDist(1+Ni-1+Ns-1, theta)  
   enddo
   do l=2, Nf-1
      distslice(l+Ni+Ns, 1) =&
            dist3D(slice(1, l+Ni-1+Ns-1, 1), slice(1, l+Ni-1+Ns-1, 2), slice(1, l+Ni-1+Ns-1, 3)&
            , points_courbe_equal(Ns, 1)+(l-2)*dsf, points_courbe_equal(1, 2), points_courbe_equal(1, 3))
      distslice(l+Ni+Ns, 2) =&
            dist3D(slice(1, l+Ni-1+Ns-1, 1), slice(1, l+Ni-1+Ns-1, 2), slice(1, l+Ni-1+Ns-1, 3)&
            , points_courbe_equal(Ns, 1)+l*dsf, points_courbe_equal(1, 2), points_courbe_equal(1, 3))
      write(*, *) "DOTPRODf  ", l, " ", dotProd(points_courbe_equal(Ns, 1)+(l-1)*dsf, points_courbe_equal(1, 2), &
            slice(1, l+Ni-1+Ns-1, 1), slice(1, l+Ni-1+Ns-1, 2), points_courbe_equal(Ns, 1) + &
            l*dsf, points_courbe_equal(1, 2))
      do theta=1, nbtheta
         valDist(l+Ni-1+Ns-1, theta) =&
               dist3D(slice(theta, l+Ni-1+Ns-1, 1), slice(theta, l+Ni-1+Ns-1, 2), slice(theta, l+Ni-1+Ns-1, 3)&
               , points_courbe_equal(Ns, 1)+(l-1)*dsf, points_courbe_equal(1, 2), points_courbe_equal(1, 3))
         sinTheta_tab(l+Ni-1+Ns-1, theta) = (slice(theta, l+Ni-1+Ns-1, 3)-points_courbe_equal(1, 3))/valDist(l+Ni-1+Ns-1, theta) 
         cosTheta_tab(l+Ni-1+Ns-1, theta) = (slice(theta, l+Ni-1+Ns-1, 2)-points_courbe_equal(1, 2))/valDist(l+Ni-1+Ns-1, theta)  
      enddo
   enddo
   distslice(Nf+Ni+Ns, 1) =&
         dist3D(slice(1, Nf+Ni-1+Ns-1, 1), slice(1, Nf+Ni-1+Ns-1, 2), slice(1, Nf+NI-1+Ns-1, 3)&
         , points_courbe_equal(Ns, 1)+(Nf-2)*dsf, points_courbe_equal(1, 2), points_courbe_equal(1, 3))
   distslice(Nf+Ni+Ns, 2) =&
         dist3D(slice(1, Nf+Ni-1+Ns-1, 1), slice(1, Nf+Ni-1+Ns-1, 2), slice(1, Nf+Ni-1+Ns-1, 3)&
         , points_courbe_equal(Ns, 1)+(Nf-1)*dsf, points_courbe_equal(1, 2), points_courbe_equal(1, 3))
   do theta=1, nbtheta
      valDist(Nf+Ni-1+Ns-1, theta) =&
            dist3D(slice(theta, Nf+Ni-1+Ns-1, 1), slice(theta, Nf+Ni-1+Ns-1, 2), slice(theta, Nf+Ni-1+Ns-1, 3)&
            , points_courbe_equal(Ns, 1)+(Nf-1)*dsf, points_courbe_equal(1, 2), points_courbe_equal(1, 3))
      sinTheta_tab(Nf+Ni-1+Ns-1, theta) = (slice(theta, Nf+Ni-1+Ns-1, 3)-points_courbe_equal(1, 3))/valDist(Nf+Ni-1+Ns-1, theta) 
      cosTheta_tab(Nf+Ni-1+Ns-1, theta) = (slice(theta, Nf+Ni-1+Ns-1, 2)-points_courbe_equal(1, 2))/valDist(Nf+Ni-1+Ns-1, theta)  
   enddo

   call cpu_time(end_time)
   write(*, *) "Frame 1 Langrangian markers processing time = ", end_time-start_time, " sec" ! unknown

   !**
   !! Now, we can write in output the initial shape of the zebrafish (Lagrangian markers / midline)
   !**


   open(unit=91, file=trim(target_folder)//'/theta00.txt', status='unknown')
   open(unit=90, file=trim(target_folder)//'/controltheta00.txt', status='unknown')
   open(unit=88, file=trim(target_folder)//'/coupeLarge00.txt', status='unknown')
   open(unit=89, file=trim(target_folder)//'/coupetmp00.txt', status='unknown')
   open(unit=87, file=trim(target_folder)//'/coupe00.txt', status='unknown')
   open(unit=84, file=trim(target_folder)//'/right00.txt', status='unknown')
   open(unit=83, file=trim(target_folder)//'/left00.txt', status='unknown')
   open(unit=86, file=trim(target_folder)//'/controlright00.txt', status='unknown')
   open(unit=85, file=trim(target_folder)//'/controlleft00.txt', status='unknown')
   open(unit=79, file=trim(target_folder)//'/skelet00.vtk', status='unknown')
   open(unit=81, file=trim(target_folder)//'/skelett00.txt', status='unknown')
   open(unit=82, file=trim(target_folder)//'/skeletteq00.txt', status='unknown')
   file_id = 79
   n_dim = 2
   call write_vtk_header(file_id, nx, ny, nz, dx, dy, dz, n_dim)

   open(unit=78, file=trim(target_folder)//'/surf00.vts', status='unknown')
   open(unit=92, file=trim(target_folder)//'/surf00.dat', status='unknown')
   write(78, '(a)') "<?xml version=""1.0""?>"
   write(78, '(a)') "<VTKFile type=""StructuredGrid"" version=""0.1"" byte_order&
   &     =""LittleEndian"" compressor=""vtkZLibDataCompressor"">"
   write(78, '(a, I3, a, I3, a)') "<StructuredGrid WholeExtent=""0 ", nbtheta, " 0 ", Ni-1+Ns+Nf-1-1, " 0 0"">"
   write(78, '(a, I3, a, I3, a)') "<Piece Extent=""0 ", nbtheta, " 0 ", Ni-1+Ns+Nf-1-1, " 0 0"">"
   write(78, '(a)') "<PointData >"
   write(78, '(a)') "</PointData>"
   write(78, '(a)') "<CellData>"
   write(78, '(a)') "</CellData>"
   write(78, '(a)') "<Points>"
   write(78, '(a)') "<DataArray NumberOfComponents=""3"" type=""Float64"" format=""ascii"" >"  
   do l=1, Ni+Ns+Nf-2
      do theta=1, nbtheta 
         write(78, *) slice(theta, l, 1), " ", slice(theta, l, 2), " ", slice(theta, l, 3)
         write(92, *) slice(theta, l, 1), " ", slice(theta, l, 2), " ", slice(theta, l, 3)
         if (theta==88) write(*, *) "cmptheta75  ", l, " ", dist3D(slice(theta, l, 1), slice(theta, l, 2), slice(theta, l, 3)&
               , slicetmp(theta, l, 1), slicetmp(theta, l, 2), slicetmp(theta, l, 3))
         if (l==Ni+1+10) write(*, *) "cmpslice10  ", theta, " ", dist3D(slice(theta, l, 1), slice(theta, l, 2), slice(theta, l, 3)&
               , slicetmp(theta, l, 1), slicetmp(theta, l, 2), slicetmp(theta, l, 3))
         if (l==Ni+1+5) write(*, *) "cmpslice5  ", theta, " ", dist3D(slice(theta, l, 1), slice(theta, l, 2), slice(theta, l, 3)&
               , slicetmp(theta, l, 1), slicetmp(theta, l, 2), slicetmp(theta, l, 3))
         if (l==Ni+1+2) write(*, *) "cmpslice2  ", theta, " ", dist3D(slice(theta, l, 1), slice(theta, l, 2), slice(theta, l, 3)&
               , slicetmp(theta, l, 1), slicetmp(theta, l, 2), slicetmp(theta, l, 3))
         if (l==Ni+1) write(*, *) "cmpsliceNI  ", theta, " ", dist3D(slice(theta, l, 1), slice(theta, l, 2), slice(theta, l, 3)&
               , slicetmp(theta, l, 1), slicetmp(theta, l, 2), slicetmp(theta, l, 3))
      enddo
         write(78, *) slice(1, l, 1), " ", slice(1, l, 2), " ", slice(1, l, 3)
         write(92, *) slice(1, l, 1), " ", slice(1, l, 2), " ", slice(1, l, 3)
   enddo
   write(78, '(a)') "</DataArray>"
   write(78, '(a)') "</Points>"
   write(78, '(a)') "</Piece>"
   write(78, '(a)') "</StructuredGrid>"
   write(78, '(a)') "</VTKFile>"
   close(78)    
   close(92)    

   do k=1, ny
      do j=1, nx
         write(79, *) rhoSlices2(j, k)
      enddo
   enddo

   do l=1, Ns
      write(81, *) points_courbe(l, 1), " ", points_courbe(l, 2), " ", points_courbe(l, 3)
      write(82, *) points_courbe_equal(l, 1), " ", points_courbe_equal(l, 2), " ", points_courbe_equal(l, 3)
   enddo

   do l=1, size(slice, 2)

      theta = 1
      write(83, *) slice(theta, l, 1), " ", slice(theta, l, 2), " ", slice(theta, l, 3)
      theta = 1+nbtheta/2
      write(84, *) slice(theta, l, 1), " ", slice(theta, l, 2), " ", slice(theta, l, 3)
   enddo

   do l=1, size(slicecontrol, 2) 
      theta = 1
      write(85, *) slicecontrol(theta, l, 1), " ", slicecontrol(theta, l, 2), " ", slicecontrol(theta, l, 3)
      theta = 1+nbtheta/2
      write(86, *) slicecontrol(theta, l, 1), " ", slicecontrol(theta, l, 2), " ", slicecontrol(theta, l, 3)
   enddo
   l = 120
   do theta = 1, nbtheta
      write(87, *) slicecontrol(theta, l, 2), " ", slicecontrol(theta, l, 3)

   enddo
   l = 120 
   write(89, *) points_courbe_equal(l, 2), " ", points_courbe_equal(l, 3)
   do theta = 1, nbtheta
      write(89, *) slice(theta, l, 2), " ", slice(theta, l, 3)
   enddo
   l = 120 
   write(88, *) points_courbe_equal(l, 2), " ", points_courbe_equal(l, 3)
   do theta = 1, nbthetaLarge
      write(88, *) slicetmp(theta, l, 2), " ", slicetmp(theta, l, 3)
   enddo
   theta = 34
   do l = itail, ihead
      write(90, *) l, " ", valTh(l-itail+1, theta)
   enddo
   do l=1, size(valTheta, 1)
      if (l<=size(points_courbe_equal, 1)) then
         write(91, *) points_courbe_equal(l, 1), " ", valTheta(l, theta)
      else
         write(91, *) points_courbe_equal(Ns, 1)+(l-Ns)*dsf, " ", valTheta(l, theta)
      endif
   enddo
   close(90)
   close(91)
   close(89)
   close(88)
   close(87)
   close(79)    
   close(85)
   close(86)
   close(81)
   close(82)
   close(83)
   close(84)
   deallocate(midline, midlinebis)
   deallocate(points_control)
   deallocate(slicecontrolm, slicecontrolf, slicecontroli, slicecontroltmpf, slicecontroltmpi)
   deallocate(slicecontrolLarge)
   deallocate(slicecontrol, slicecontroltmp)
   deallocate(valTh, valThLarge, valThtmp)

   deallocate(rhoSlices2, rhoSlices, gradPhi, xx, yy, zz, Nseed, skel, skel2, skel3, longTh)
   
   !! FILM
   N =  200 
   nx = 200 
   ny = 200 
   nz = 200

   dx = 0.001_pr*0.0256
   dy = 0.001_pr*0.0256
   dz = 0.001_pr*0.0256

   allocate(rhoSlices2(nx, ny), gradPhi(nx, ny))
   allocate(xx(nx), yy(ny), zz(nz))
   allocate(Nseed(nx, ny), skel(nx, ny), skel2(nx, ny), skel3(nx, ny), longTh(nx))

   do i=1, nx
      xx(i) = x0 + (float(i)-1)*dx
   enddo
   do j=1, ny
      yy(j) = x0 + (float(j)-1)*dy
   enddo
   do k=1, nz
      zz(k) = x0 + (float(k)-1)*dz
   enddo

   open(unit=80, file=trim(target_folder)//'/skeleton.txt', status='unknown')
   write(80, *) kt, "    ", nl, " ", long3, " ", long2, " ", long

   !**
   !! Loop over each experimental image to construct the deformed zebrafish shape (Lagrangian markers) based on the pre-built shape
   !**
   t = 0
   iter = 0
   idisplay = 10
   iUpdateDist = 5
   do kt = 1, n_frames
      write(*, *)"Image iteration num : ", kt
      !! FILM
      open(unit=78, file="D:\Users\Théo\GitHub\Repositories\deformation3d\IMAGES4\Image_"&
      //str(kt+picNum-1)//".dat", &
            status='unknown')
   !**
   !! Computation of the level-set 
   !**
      do k=1, nx
         do j=1, ny
            read(78, *) rhoSlices2(k, j)
            !read(78, *) rhoSlices2(nx-k+1, j)
         enddo
      enddo
      close(78)
      rhoSlices2 = rhoSlices2 - 0.5_pr*(maxval(rhoSlices2) + minval(rhoSlices2))
      !rhoSlices2 = (2*rhoSlices2-minval(rhoSlices2)-maxval(rhoSlices2))/(maxval(rhoSlices2)-minval(rhoSlices2))
      !dx = 0.001_pr*0.0301*1E6
      !dy = 0.001_pr*0.0301*1E6
      !dz = 0.001_pr*0.0301*1E6
      dx = 0.001_pr*0.0256*1E6
      dy = 0.001_pr*0.0256*1E6
      dz = 0.001_pr*0.0256*1E6
      call updateDistance(rhoSlices2, gradPhi)
      !rhoSlices2 = (2*rhoSlices2-minval(rhoSlices2)-maxval(rhoSlices2))/(maxval(rhoSlices2)-minval(rhoSlices2))
      !rhoSlices2 = rhoSlices2/maxval(rhoSlices2)
      !dx = 0.001_pr*0.0301
      !dy = 0.001_pr*0.0301
      !dz = 0.001_pr*0.0301
      dx = 0.001_pr*0.0256
      dy = 0.001_pr*0.0256
      dz = 0.001_pr*0.0256
      do i=1, nx
         xx(i) = x0 + (float(i)-1)*dx
      enddo
      do j=1, ny
         yy(j) = x0 + (float(j)-1)*dy
      enddo
      do k=1, nz
         zz(k) = x0 + (float(k)-1)*dz
      enddo

   !**
   !! Computation of the midline (skel: skeleton based one the level-set, skel2: gradient of the level-set)
   !**
      dir1 = 0
      dir2 = 0
      dir3 = 0
      dir4 = 0
      Nseed = 0
      skel = 0
      skel2 = 0
      do i=2, nx-1
         do j=2, ny-1

            !my_training_film
            if ((gradPhi(i, j)<0.62).and.(rhoSlices2(i, j)>0._pr)) skel(i, j) = 1
            if ((gradPhi(i, j)<0.94).and.(rhoSlices2(i, j)>0._pr)) skel2(i, j) = 1
            if ((kt==60).and.(gradPhi(i, j)<0.98).and.(rhoSlices2(i, j)>0._pr)) skel2(i, j) = 1
            if ((kt==365).and.(gradPhi(i, j)<0.98).and.(rhoSlices2(i, j)>0._pr)) skel2(i, j) = 1
            if ((kt==501).and.(gradPhi(i, j)<0.98).and.(rhoSlices2(i, j)>0._pr)) skel2(i, j) = 1      
            if ((kt==162).and.(gradPhi(i, j)<0.98).and.(rhoSlices2(i, j)>0._pr)) skel2(i, j) = 1      
            if ((kt==561).and.(gradPhi(i, j)<0.98).and.(rhoSlices2(i, j)>0._pr)) skel2(i, j) = 1  
            if ((kt==247).and.(gradPhi(i, j)<0.98).and.(rhoSlices2(i, j)>0._pr)) skel2(i, j) = 1
            if ((kt==303).and.(gradPhi(i, j)<0.98).and.(rhoSlices2(i, j)>0._pr)) skel2(i, j) = 1
            if ((kt==304).and.(gradPhi(i, j)<0.98).and.(rhoSlices2(i, j)>0._pr)) skel2(i, j) = 1
         !**
         !! First: Thresholds for extracting the skeleton and its gradient
         !**
         enddo
      enddo

      !**
      !! Second step: we reduce the gradient map
      !**
      skel3 = 0
      do i=2, nx-1
         do j=2, ny-1
            Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
            skel2(i, j-1) + skel2(i+1, j-1)
            !if ((skel2(i, j)==1).and.(skel2(i+1, j)==1).and.(skel2(i, j+1)==1).and.(skel2(i-1, j)==1).and.(skel2(i, j-1)==1)) skel22(i, j) = 1
            if ((skel2(i, j)==1).and.(Nseed(i, j)>5)) skel3(i, j) = 1
            !if ((skel2(i, j)==1).and.(Nseed(i, j)==6)) skel3(i, j) = 1
         enddo
      enddo
      skel2 = skel3

   !     l = 1
   !     pix1 = 0
   !     pix2 = 0
   !     skel2 = skel
   !     do i=2, N-1
   !        do j=2, ny-1
   !           Nseed(i, j) = skel(i+1, j) + skel(i+1, j+1) + skel(i, j+1) + skel(i-1, j+1) + skel(i-1, j) + skel(i-1, j-1) + skel(i, j-1) +&
   !                skel(i+1, j-1)
   !           !if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i-1, j)==1).and.(skel(i+1, j+1)==1)) skel2(i, j) = 1
   !           !if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i-1, j)==1).and.(skel(i+1, j-1)==1)) skel2(i, j) = 1
   !           !if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i, j-1)==1).and.(skel(i+1, j+1)==1)) skel2(i, j) = 1
   !           !if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i, j-1)==1).and.(skel(i-1, j+1)==1)) skel2(i, j) = 1
   !           !if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i, j-1)==1).and.(skel(i, j+1)==1)) skel2(i, j) = 1
   !           !if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i-1, j)==1).and.(skel(i+1, j)==1)) skel2(i, j) = 1
   !           !if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i-1, j-1)==1).and.(skel(i+1, j+1)==1)) skel2(i, j) = 1
   !           !if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i-1, j+1)==1).and.(skel(i+1, j-1)==1)) skel2(i, j) = 1
   !           if ((skel(i, j)==0).and.(skel(i-1, j)==1).and.(skel(i+1, j+1)==1)) skel2(i, j) = 1
   !           if ((skel(i, j)==0).and.(skel(i-1, j)==1).and.(skel(i+1, j-1)==1)) skel2(i, j) = 1
   !           if ((skel(i, j)==0).and.(skel(i, j-1)==1).and.(skel(i+1, j+1)==1)) skel2(i, j) = 1
   !           if ((skel(i, j)==0).and.(skel(i, j-1)==1).and.(skel(i-1, j+1)==1)) skel2(i, j) = 1
   !           if ((skel(i, j)==0).and.(skel(i, j-1)==1).and.(skel(i, j+1)==1)) skel2(i, j) = 1
   !           if ((skel(i, j)==0).and.(skel(i-1, j)==1).and.(skel(i+1, j)==1)) skel2(i, j) = 1
   !           if ((skel(i, j)==0).and.(skel(i-1, j-1)==1).and.(skel(i+1, j+1)==1)) skel2(i, j) = 1
   !           if ((skel(i, j)==0).and.(skel(i-1, j+1)==1).and.(skel(i+1, j-1)==1)) skel2(i, j) = 1
   !        enddo
   !     enddo
   !     skel = skel2
   !     skel2 = skel
   !     do i=2, N-1
   !        do j=2, ny-1
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i-1, j)==1).and.(skel2(i-1, j+1)==1)) skel2(i-1, j) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i-1, j+1)==1).and.(skel2(i, j+1)==1)) skel2(i, j+1) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i, j+1)==1).and.(skel2(i+1, j+1)==1)) skel2(i, j+1) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i+1, j+1)==1).and.(skel2(i+1, j)==1)) skel2(i+1, j) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i+1, j)==1).and.(skel2(i+1, j-1)==1)) skel2(i+1, j) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i+1, j-1)==1).and.(skel2(i, j-1)==1)) skel2(i, j-1) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i, j-1)==1).and.(skel2(i-1, j-1)==1)) skel2(i, j-1) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i-1, j-1)==1).and.(skel2(i-1, j)==1)) skel2(i-1, j) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)==0)) skel2(i, j) = 0
   !        enddo
   !     enddo
   !     skel = skel2
   !!     do while ((pix1==1).or.(pix2==1).or.(l==1))
   !!        pix1 = 0
   !!        pix2 = 0
   !!        l = l + 1
   !!        skel2 = skel
   !!        do i=2, N-1
   !!           do j=2, ny-1
   !!              maxv = 0
   !!              if ((skel(i+1, j)==0).and.(skel(i+1, j+1)==1)) maxv = maxv + 1
   !!              if ((skel(i+1, j+1)==0).and.(skel(i, j+1)==1)) maxv = maxv + 1
   !!              if ((skel(i, j+1)==0).and.(skel(i-1, j+1)==1)) maxv = maxv + 1
   !!              if ((skel(i-1, j+1)==0).and.(skel(i-1, j)==1)) maxv = maxv + 1
   !!              if ((skel(i-1, j)==0).and.(skel(i-1, j-1)==1)) maxv = maxv + 1
   !!              if ((skel(i-1, j-1)==0).and.(skel(i, j-1)==1)) maxv = maxv + 1
   !!              if ((skel(i, j-1)==0).and.(skel(i+1, j-1)==1)) maxv = maxv + 1
   !!              if ((skel(i+1, j-1)==0).and.(skel(i+1, j)==1)) maxv = maxv + 1
   !!              Nseed(i, j) = skel(i+1, j) + skel(i+1, j+1) + skel(i, j+1) + skel(i-1, j+1) + skel(i-1, j) + skel(i-1, j-1) + skel(i, j-1) +&
   !!                   skel(i+1, j-1)
   !!              if ((skel(i, j).eq.1).and.(maxv.eq.1).and.(Nseed(i, j)>=2).and.(Nseed(i, j)<=6).and.&
   !!                   ((skel(i+1, j)==0).or.(skel(i-1, j)==0).or.(skel(i, j-1)==0)).and.((skel(i, j-1)==0).or.(skel(i, j+1)==0)&
   !!                   !.or.(skel(i, j+1)==0))) then
   !!                   .or.(skel(i-1, j)==0))) then
   !!                 skel2(i, j) = 0
   !!                 pix1 = 1
   !!              endif
   !!              if ((skel(i, j).eq.1).and.(maxv.eq.1).and.(Nseed(i, j)>=2).and.(Nseed(i, j)<=6).and.&
   !!                   ((skel(i, j+1)==0).or.(skel(i-1, j)==0).or.(skel(i, j-1)==0)).and.((skel(i+1, j)==0).or.(skel(i, j+1)==0)&
   !!                   .or.(skel(i-1, j)==0))) then
   !!                 skel2(i, j) = 0
   !!                 pix2 = 1
   !!              endif
   !!           enddo
   !!        enddo
   !!        skel = skel2
   !!        skel2 = skel
   !!        do i=2, N-1
   !!           do j=2, ny-1
   !!              maxv = 0
   !!              if ((skel(i+1, j)==0).and.(skel(i+1, j+1)==1)) maxv = maxv + 1
   !!              if ((skel(i+1, j+1)==0).and.(skel(i, j+1)==1)) maxv = maxv + 1
   !!              if ((skel(i, j+1)==0).and.(skel(i-1, j+1)==1)) maxv = maxv + 1
   !!              if ((skel(i-1, j+1)==0).and.(skel(i-1, j)==1)) maxv = maxv + 1
   !!              if ((skel(i-1, j)==0).and.(skel(i-1, j-1)==1)) maxv = maxv + 1
   !!              if ((skel(i-1, j-1)==0).and.(skel(i, j-1)==1)) maxv = maxv + 1
   !!              if ((skel(i, j-1)==0).and.(skel(i+1, j-1)==1)) maxv = maxv + 1
   !!              if ((skel(i+1, j-1)==0).and.(skel(i+1, j)==1)) maxv = maxv + 1
   !!              Nseed(i, j) = skel(i+1, j) + skel(i+1, j+1) + skel(i, j+1) + skel(i-1, j+1) + skel(i-1, j) + skel(i-1, j-1) + skel(i, j-1) +&
   !!                   skel(i+1, j-1)
   !!              if ((skel(i, j).eq.1).and.(maxv.eq.1).and.(Nseed(i, j)>=2).and.(Nseed(i, j)<=6).and.&
   !!                   ((skel(i, j+1)==0).or.(skel(i+1, j)==0).or.(skel(i, j-1)==0)).and.((skel(i+1, j)==0).or.(skel(i, j-1)==0)&
   !!                   .or.(skel(i-1, j)==0))) then
   !!                 skel2(i, j) = 0
   !!                 pix2 = 1
   !!              endif
   !!           enddo
   !!        enddo
   !!        skel = skel2
   !!     enddo
   !     skel2 = skel
   !     do i=2, N-1
   !        do j=2, ny-1
   !           Nseed(i, j) = skel(i+1, j) + skel(i+1, j+1) + skel(i, j+1) + skel(i-1, j+1) + skel(i-1, j) + skel(i-1, j-1) + skel(i, j-1) +&
   !                skel(i+1, j-1)
   !           if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i-1, j)==1).and.(skel(i+1, j+1)==1)) skel2(i, j) = 1
   !           if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i-1, j)==1).and.(skel(i+1, j-1)==1)) skel2(i, j) = 1
   !           if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i, j-1)==1).and.(skel(i+1, j+1)==1)) skel2(i, j) = 1
   !           if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i, j-1)==1).and.(skel(i-1, j+1)==1)) skel2(i, j) = 1
   !           if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i, j-1)==1).and.(skel(i, j+1)==1)) skel2(i, j) = 1
   !           if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i-1, j)==1).and.(skel(i+1, j)==1)) skel2(i, j) = 1
   !           if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i-1, j-1)==1).and.(skel(i+1, j+1)==1)) skel2(i, j) = 1
   !           if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i-1, j+1)==1).and.(skel(i+1, j-1)==1)) skel2(i, j) = 1
   !        enddo
   !     enddo
   !     skel = skel2
   !     skel2 = skel
   !     do i=2, N-1
   !        do j=2, ny-1
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i-1, j)==1).and.(skel2(i-1, j+1)==1)) skel2(i-1, j) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i-1, j+1)==1).and.(skel2(i, j+1)==1)) skel2(i, j+1) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i, j+1)==1).and.(skel2(i+1, j+1)==1)) skel2(i, j+1) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i+1, j+1)==1).and.(skel2(i+1, j)==1)) skel2(i+1, j) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i+1, j)==1).and.(skel2(i+1, j-1)==1)) skel2(i+1, j) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i+1, j-1)==1).and.(skel2(i, j-1)==1)) skel2(i, j-1) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i, j-1)==1).and.(skel2(i-1, j-1)==1)) skel2(i, j-1) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)>2).and.(skel2(i-1, j-1)==1).and.(skel2(i-1, j)==1)) skel2(i-1, j) = 0
   !           Nseed(i, j) = skel2(i+1, j) + skel2(i+1, j+1) + skel2(i, j+1) + skel2(i-1, j+1) + skel2(i-1, j) + skel2(i-1, j-1) +&
   !                skel2(i, j-1) + skel2(i+1, j-1)
   !           if ((skel2(i, j)==1).and.(Nseed(i, j)==0)) skel2(i, j) = 0
   !        enddo
   !     enddo
   !     skel = skel2

      !**
      !! Third step: we complete the level-set map
      !**
   !     call ZST(skel)
      skel3 = skel
      do i=2, nx-1
         do j=2, ny-1
            Nseed(i, j) = skel(i+1, j) + skel(i+1, j+1) + skel(i, j+1) + skel(i-1, j+1) + skel(i-1, j) + &
            skel(i-1, j-1) + skel(i, j-1) +skel(i+1, j-1)
            if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i-1, j)==1).and.(skel(i+1, j+1)==1)) skel3(i, j) = 1
            if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i-1, j)==1).and.(skel(i+1, j-1)==1)) skel3(i, j) = 1
            if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i, j-1)==1).and.(skel(i+1, j+1)==1)) skel3(i, j) = 1
            if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i, j-1)==1).and.(skel(i-1, j+1)==1)) skel3(i, j) = 1
            if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i, j-1)==1).and.(skel(i, j+1)==1)) skel3(i, j) = 1
            if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i-1, j)==1).and.(skel(i+1, j)==1)) skel3(i, j) = 1
            if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i-1, j-1)==1).and.(skel(i+1, j+1)==1)) skel3(i, j) = 1
            if ((skel(i, j)==0).and.(Nseed(i, j)==2).and.(skel(i-1, j+1)==1).and.(skel(i+1, j-1)==1)) skel3(i, j) = 1
         enddo
      enddo

      !**
      !! Fourth step: MANUAL CORRECTIONS of level-set and gradient maps
      !**
         !! HERE for midline adjustments
      if(kt==370) skel3(158, 95) = 0
      if(kt==392) skel3(152, 92) = 0
      if(kt==395) then
         skel3(154, 93) = 0
         skel3(153, 92) = 0
         skel3(152, 93) = 0
      endif
      if(kt==440) skel3(160, 109) = 0
      if(kt==523) skel3(157, 95) = 0
      if(kt==247) skel3(152, 95) = 0
      if(kt==303) skel3(150, 112) = 0
      if(kt==70) skel3(159, 107) = 0
      if(kt==162) then
         skel3(148, 89) = 0
         skel3(147, 88) = 0
      endif

      if(kt == 306) skel3(150, 112) = 0
      if((kt >= 309).and.(kt <= 311)) skel3(150, 113) = 0
      if(kt == 315) then
         skel3(149, 113) = 0
         skel3(150, 114) = 0
      endif
      if(kt == 304) then
         skel3(148, 110) = 0
         skel3(149, 111) =0
      endif

      skel = skel3

      !**
      !! Fifth step: construction of a 1-pixel wide skeleton
      !**
      call filterskel(skel, rhoSlices2)
      tmpbool = 0

      sizeSkel = sum(skel)
      write(*, *) "Longlong 00  : ", long3, " ", long2, " ", long

      !**
      !! Sixth step: construction of the midline from the 2D-map (skeleton)
      !**
      allocate(midline(sizeSkel, 3))
      write(*, *) "okcut"
      l=1
      do i=2, nx-1
         do j=2, ny-1
            if ((skel(i, j)==1).and.(l==1)) then
               Nseed(i, j) = skel(i+1, j) + skel(i+1, j+1) + skel(i, j+1) + skel(i-1, j+1) + skel(i-1, j) + &
               skel(i-1, j-1) + skel(i, j-1) + skel(i+1, j-1)
               !!if ((Nseed(i, j)==1).and.((abs(xskelL-i)<100).and.(abs(yskelL-j)<100))) then 
               !              if ((((kt<121).or.(kt>140)).and.(Nseed(i, j)==1).and.((abs(xskelL-i)<100).and.(abs(yskelL-j)<100)))&
               !                   .or.((Nseed(i, j)==1).and.(kt==1))) then 
               !                 midline(l, 1) = i
               !                 midline(l, 2) = j
               !                 l = l+1
               !              else if (((Nseed(i, j)==1).and.((abs(xskelL-i)<20).and.(abs(yskelL-j)<20))).or.((Nseed(i, j)==1).and.(kt==1))) then 
               !if (((Nseed(i, j)==1).and.((abs(xskelL-i)<100).and.(abs(yskelL-j)<100))).or.((Nseed(i, j)==1).and.(kt==1))) then 
               if (((Nseed(i, j)==1).and.((abs(xskelL-i)<50).and.(abs(yskelL-j)<50))).or.((Nseed(i, j)==1).and.(kt==1))) then 
               !if (((Nseed(i, j)==1).and.((abs(xskelL-i)<10).and.(abs(yskelL-j)<10))).or.((Nseed(i, j)==1).and.(kt==1))) then 
                  midline(l, 1) = i
                  midline(l, 2) = j
                  l = l+1
               endif
            endif
         enddo
      enddo
      write(*, *) "okcut"
      !     if (kt==130) then
      !        skel(137, 78)=0
      !        midline(1, 1)=136
      !        midline(1, 2)=77
      !     endif
      midline(1, 3) = nint(zslice)
      write(*, *) "CHECK midline : ", midline(1, 1), " ", midline(1, 2), " ", midline(size(midline, 1), 1), " ", &
            midline(size(midline, 1), 2)
      xskelL = midline(1, 1)
      yskelL = midline(1, 2)
      boolskel=0
      nl = sum(skel)
      do l=2, sum(skel)
         do i=2, nx-1
            do j=2, ny-1
               if ((skel(i, j)==1).and.(boolskel==0)) then
                  Nseed(i, j) = skel(i+1, j) + skel(i+1, j+1) + skel(i, j+1) + skel(i-1, j+1) + skel(i-1, j) + &
                        skel(i-1, j-1) + skel(i, j-1) + skel(i+1, j-1)
                  if ((Nseed(i, j)>1).and.((midline(l-1, 1)==i).or.(midline(l-1, 1)==i-1).or.(midline(l-1, 1)==i+1)).and.&
                        ((midline(l-1, 2)==j).or.(midline(l-1, 2)==j-1).or.(midline(l-1, 2)==j+1)).and.&
                        (appartient(midline, i, j, l).neqv..true.)) then 
                     midline(l, 1) = i
                     midline(l, 2) = j
                  endif
                  !                 if ((Nseed(i, j)==1).and.(l==sum(skel))) then
                  if ((Nseed(i, j)==1).and.((midline(l-1, 1)==i).or.(midline(l-1, 1)==i-1).or.(midline(l-1, 1)==i+1)).and.&
                        ((midline(l-1, 2)==j).or.(midline(l-1, 2)==j-1).or.(midline(l-1, 2)==j+1)).and.&
                        (appartient(midline, i, j, l).neqv..true.)) then 
                     !                 if (Nseed(i, j)==1) then
                     midline(l, 1) = i
                     boolskel = 1
                     nl=l
                     midline(l, 2) = j
                  endif
               endif
            enddo
         enddo
         midline(l, 3) = nint(zslice)
      enddo
      write(*, *) "CHECK midline : ", midline(1, 1), " ", midline(1, 2), " ", midline(size(midline, 1), 1), " ", & 
            midline(size(midline, 1), 2) , " ", size(midline, 1), " ", sum(skel), " ", nl
      write(*, *) "CHECK midline BIS : ", midline(1, 1), " ", midline(1, 2), " ", midline(nl, 1), " ", midline(nl, 2), " ", nl

      if (nl<10) then
               write(*, *) "error midline queue fourche"

               do l=1, nl
                  skel(midline(l, 1), midline(l, 2)) = 0
               enddo

               l=1
               do i=2, nx-1
                  do j=2, ny-1
                     if ((skel(i, j)==1).and.(l==1)) then
                        Nseed(i, j) = skel(i+1, j) + skel(i+1, j+1) + skel(i, j+1) + skel(i-1, j+1) + &
                        skel(i-1, j) + skel(i-1, j-1) + skel(i, j-1) + skel(i+1, j-1)
                        !if (((Nseed(i, j)==1).and.((abs(xskelL-i)<100).and.(abs(yskelL-j)<100))).or.((Nseed(i, j)==1).and.(kt==1))) then 
                        if (((Nseed(i, j)==1).and.((abs(xskelL-i)<50).and.(abs(yskelL-j)<50))).or.&
                        ((Nseed(i, j)==1).and.(kt==1))) then 
                        !if (((Nseed(i, j)==1).and.((abs(xskelL-i)<10).and.(abs(yskelL-j)<10))).or.((Nseed(i, j)==1).and.(kt==1))) then 
                           midline(l, 1) = i
                           midline(l, 2) = j
                           l = l+1
                        endif
                     endif
                  enddo
               enddo
               midline(1, 3) = nint(zslice)
               write(*, *) "CHECK midline : ", midline(1, 1), " ", midline(1, 2), " ", midline(size(midline, 1), 1), " "&
               , midline(size(midline, 1), 2)
               xskelL = midline(1, 1)
               yskelL = midline(1, 2)
               boolskel=0
               do l=2, sum(skel)
                  do i=2, nx-1
                     do j=2, ny-1
                        if ((skel(i, j)==1).and.(boolskel==0)) then
                           Nseed(i, j) = skel(i+1, j) + skel(i+1, j+1) + skel(i, j+1) + skel(i-1, j+1) + skel(i-1, j) &
                           + skel(i-1, j-1) + skel(i, j-1) + skel(i+1, j-1)
                           if ((Nseed(i, j)>1).and.((midline(l-1, 1)==i).or.(midline(l-1, 1)==i-1).or.(midline(l-1, 1)==i+1)).and.&
                                 ((midline(l-1, 2)==j).or.(midline(l-1, 2)==j-1).or.(midline(l-1, 2)==j+1)).and.&
                                 (appartient(midline, i, j, l).neqv..true.)) then 
                              midline(l, 1) = i
                              midline(l, 2) = j
                           endif
                           if ((Nseed(i, j)==1).and.((midline(l-1, 1)==i).or.(midline(l-1, 1)==i-1).or.(midline(l-1, 1)==i+1)).and.&
                                 ((midline(l-1, 2)==j).or.(midline(l-1, 2)==j-1).or.(midline(l-1, 2)==j+1)).and.&
                                 (appartient(midline, i, j, l).neqv..true.)) then 
                              midline(l, 1) = i
                              boolskel = 1
                              nl=l
                              midline(l, 2) = j
                           endif
                        endif
                     enddo
                  enddo
                  midline(l, 3) = nint(zslice)
               enddo
      endif


















      !**
      !! Generating a 2D-map from the midline (for visualisation purpose)
      !**

   !     call filtermidline(midline, rhoSlices2)
   !     if (.not.((midline(1, 1)==midline(2, 1)).or.(midline(1, 2)==midline(2, 2)))) then
   !       do l=1, size(midline, 1)-1
   !         midline(l, 1) = midline(l+1, 1)
   !         midline(l, 2) = midline(l+1, 2)
   !       enddo
   !       nl=nl-1
   !     endif
   !     skel2 = skel
      skel = 0
      do l=1, nl
   !     do l=1, size(midline, 1)
         skel(midline(l, 1), midline(l, 2)) = 1
      enddo

      !**
      !! Seventh step: searching the endpoint of the midline (and cutting the rest) TRACKING STEP
      !**
      if (kt==1) then
               ic = -1
               jc = -1
      endif
      !if (kt==1) ic = midline(nl, 1)
      !if (kt==1) jc = midline(nl, 2)
      write(*, *) "okcut"
      !skel3 = skel2
      !call cutheadskel(midline, skel, skel2, nl, ic, jc)
      call cutheadskel(midline, skel, skel2, nl, ic, jc, kt)
      !skel2 = skel3
      write(*, *) "okcut"

      !**
      !! Heighth step: filling the midline
      !**
      xskelL = midline(1, 1)
      yskelL = midline(1, 2)
      boolskel=0
      do l=2, size(midline, 1)
         do i=2, nx-1
            do j=2, ny-1
               if ((skel(i, j)==1).and.(boolskel==0)) then
                  Nseed(i, j) = skel(i+1, j) + skel(i+1, j+1) + skel(i, j+1) + skel(i-1, j+1) + &
                        skel(i-1, j) + skel(i-1, j-1) + skel(i, j-1) + skel(i+1, j-1)
                  if ((Nseed(i, j)>1).and.((midline(l-1, 1)==i).or.(midline(l-1, 1)==i-1).or.(midline(l-1, 1)==i+1)).and.&
                        ((midline(l-1, 2)==j).or.(midline(l-1, 2)==j-1).or.(midline(l-1, 2)==j+1)).and.&
                        (appartient(midline, i, j, l).neqv..true.)) then 
                     midline(l, 1) = i
                     midline(l, 2) = j
                  endif
                  if ((Nseed(i, j)==1).and.((midline(l-1, 1)==i).or.(midline(l-1, 1)==i-1).or.(midline(l-1, 1)==i+1)).and.&
                        ((midline(l-1, 2)==j).or.(midline(l-1, 2)==j-1).or.(midline(l-1, 2)==j+1)).and.&
                        (appartient(midline, i, j, l).neqv..true.)) then 
                     midline(l, 1) = i
                     boolskel = 1
                     nl=l
                     midline(l, 2) = j
                  endif
               endif
            enddo
         enddo
      enddo

      !**
      !! Nineth step: smoothing considerably the midline
      !**
      !! FILM

      if (kt<2) then
               call filtermidline(midline, rhoSlices2)
      !call filtermidline(midline, rhoSlices2)

      !**
      !! Tenth step: deleting the first pixel if horizontal or vertical
      !**
      if (.not.((midline(1, 1)==midline(2, 1)).or.(midline(1, 2)==midline(2, 2)))) then
      do l=1, size(midline, 1)-1
      midline(l, 1) = midline(l+1, 1)
      midline(l, 2) = midline(l+1, 2)
      enddo
      nl=nl-1
      endif
      endif
      write(*, *) "CHECK midline : ", midline(1, 1), " ", midline(1, 2), " ", midline(size(midline, 1), 1), &
            " ", midline(size(midline, 1), 2), " ", size(midline, 1)
      write(*, *) "CHECK midline BIS : ", midline(1, 1), " ", midline(1, 2), " ", midline(nl, 1), " ", midline(nl, 2), " ", nl
      !call filtermidline(midline, rhoSlices2)
      !skel2 = skel

      !**
      !! Final 2D-map of the midline
      !**
      skel = 0
      do l=1, nl
         skel(midline(l, 1), midline(l, 2)) = 1
      enddo


      !!        rr = dx
      !!        x1 = xx(midline(nl-10, 1))
      !!        y1 = xx(midline(nl-10, 2))
      !!        x2 = xx(midline(nl, 1))
      !!        y2 = xx(midline(nl, 2))
      !!        sinPhi = (x2-x1)/(sqrt((x1-x2)**2+(y1-y2)**2))
      !!        cosPhi = (y2-y1)/(sqrt((x1-x2)**2+(y1-y2)**2)) 
   !!!        do while ((rhoSlices(nint(zslice), nint(x2 + rr*sinPhi), nint(y2 + rr*cosPhi))>0._pr).and.(rr<float(ny)))
      !!        do while ((rhoSlices2(nint(x2 + rr*sinPhi), nint(y2 + rr*cosPhi))>0._pr).and.(rr<float(ny)))
      !!           rr=rr+dx
      !!        enddo
      !!        xc = x2 + (rr-dx)*sinPhi +&
      !!sinPhi*rhoSlices2(nint(x2 + (rr-dx)*sinPhi), nint(y2 + (rr-dx)*cosPhi))*dy&
      !!/abs(rhoSlices2(nint(x2 + rr*sinPhi), nint(y2 + rr*cosPhi)) -&
      !!rhoSlices2(nint(x2 + (rr-dx)*sinPhi), nint(y2 + (rr-dx)*cosPhi)))
      !!        yc = y2 + (rr-dx)*cosPhi +&
      !!cosPhi*rhoSlices2(nint(x2 + (rr-dx)*sinPhi), nint(y2 + (rr-dx)*cosPhi))*dy&
      !!/abs(rhoSlices2(nint(x2 + rr*sinPhi), nint(y2 + rr*cosPhi)) -&
      !!rhoSlices2(nint(x2 + (rr-dx)*sinPhi), nint(y2 + (rr-dx)*cosPhi)))
      !!     write(*, *) "NEZ coord : ", xc, " ", yc, "    dist old ", disthead, " new ", dist(xc, yc, x2, y2), " ", sinPhi, " ", cosPhi, " ", x1, " ", y1&
      !!, " ", x2, " ", y2


      !**
      !! Definition of the control points of the midline
      !**
      !allocate(points_control(sum(skel), 2))
      !allocate(points_control(size(midline, 1), 3))
      allocate(points_control(nl, 3))
      !do l=1, sizeSkel
      do l=1, nl
         points_control(l, 1) = xx(midline(l, 1))
         points_control(l, 2) = yy(midline(l, 2))
         points_control(l, 3) = zz(nint(zslice))
      enddo
      !!TESTNEZ
      !!     points_control(nl, 1) = xc
      !!     points_control(nl, 2) = yc
      !!     points_control(nl, 3) = zslice

      !**
      !! Spline approximation of the midline
      !**
      points_courbe(1, 1) = points_control(1, 1)
      points_courbe(1, 2) = points_control(1, 2)
      points_courbe(1, 3) = points_control(1, 3)
      tb = dtb
      l = 1
      long = 0._pr
      do while ((tb<1._pr).and.(l+1<Ns))
         l = l+1
         !        call pointsBezierN(points_control, tb, px, py)
         call pointsBezierN3D(points_control, tb, px, py, pz)
         points_courbe(l, 1) = px
         points_courbe(l, 2) = py
         points_courbe(l, 3) = pz
         long = long + sqrt((points_courbe(l, 1)-points_courbe(l-1, 1))**2 + (points_courbe(l, 2)-points_courbe(l-1, 2))**2 +&
               (points_courbe(l, 3)-points_courbe(l-1, 3))**2)
         tb = tb+dtb
      enddo
      if (l==Ns-1) then
         l=l+1
         write(*, *) " ET VOILA"
         points_courbe(size(points_courbe_equal, 1), 1) = points_control(size(points_control, 1), 1)
         points_courbe(size(points_courbe_equal, 1), 2) = points_control(size(points_control, 1), 2)
         points_courbe(size(points_courbe_equal, 1), 3) = points_control(size(points_control, 1), 3)
         long = long + sqrt((points_courbe(l, 1)-points_courbe(l-1, 1))**2 + (points_courbe(l, 2)-points_courbe(l-1, 2))**2 +&
               (points_courbe(l, 3)-points_courbe(l-1, 3))**2)
      endif
      if (.not.(l==Ns)) write(*, *) "WARNING1  p ", l
      long2 = long
      write(*, *) "POINTS_COURBE_EQUAL  L==1  ", points_control(1, 1), " ", points_control(1, 2)
      write(*, *) "POINTS_COURBE_EQUAL  L==Ns  ", points_control(size(points_control, 1), 1), " ", &
            points_control(size(points_control, 1), 2)
   !     if (kt==1) long00 = long2

      !**
      !! Tail extrapolation 
      !**
      if (long00>long2) then
      !!x1 = points_courbe(1, 1)
      !!y1 = points_courbe(1, 2)
      !!x1 = points_courbe(nint(0.02*Ns), 1)
      !!y1 = points_courbe(nint(0.02*Ns), 2)
      !!x2 = points_courbe(nint(0.03*Ns), 1)
      !!y2 = points_courbe(nint(0.03*Ns), 2)
      !!x2 = points_courbe(nint(0.04*Ns), 1)
      !!y2 = points_courbe(nint(0.04*Ns), 2)
      !if (kt>80) then
      !        l=nint(0.1*nl)-1
      !else
      !        l=1
               l=5
      !endif
      x1=0.0
      y1=0.0
      do i=1, l!1!2!nint(0.01*nl)-1 !!0.1
         x1 = x1+points_control(i, 1)
         y1 = y1+points_control(i, 2)
      enddo
      x1=x1*1.0/l
      y1=y1*1.0/l
      do i=1, 1
      !if (kt>80) then
      !x2 = points_control(nint(0.1*nl), 1)
      !y2 = points_control(nint(0.1*nl), 2)
      !else
      !x2 = points_control(2, 1)!nint(0.01*nl), 1)
      !y2 = points_control(2, 2)!nint(0.01*nl), 2)
      x2 = points_control(6, 1)!nint(0.01*nl), 1)
      y2 = points_control(6, 2)!nint(0.01*nl), 2)
      !endif
      !rt = abs(long00-long2)
   !     rt = abs(long00-long2) + sqrt((points_control(i, 1)-points_control(1, 1))**2 + (points_control(i, 2)-points_control(1, 2))**2)
      rt = abs(long00-long2) + sqrt((x1-points_control(1, 1))**2 + (y1-points_control(1, 2))**2)
      !if (abs(x1-x2)<eepsilon*dx) then
      if (abs(x1-x2)/dx<eepsilon) then
         xt = x1
         !if (rt**2<(y2-y1-rt)**2) yt = y1+rt
         !if (rt**2<(y2-y1+rt)**2) yt = y1-rt
         yt = y1+rt

         if (((y2-y1)*(yt-y1))>0._pr) then
            write(*, *) "JEPASSEICIAUSSI"
            yt = y1-rt
         endif
      else
         xt = rt/(sqrt(1+((y2-y1)**2)/((x2-x1)**2))) + x1
         yt = (y2-y1)/(x2-x1)*(xt-x1) + y1
         !if (.not.(rt**2<(x2-xt)**2+(y2-yt)**2)) then
         !   xt = -rt/(sqrt(1+((y2-y1)**2)/((x2-x1)**2))) + x1
         !   yt = (y2-y1)/(x2-x1)*(xt-x1) + y1
         !endif

         if (((x2-x1)*(xt-x1)+(y2-y1)*(yt-y1))>0._pr) then
            write(*, *) "JEPASSEICI"
            xt = -rt/(sqrt(1+((y2-y1)**2)/((x2-x1)**2))) + x1
            yt = (y2-y1)/(x2-x1)*(xt-x1) + y1
         endif
      endif
      write(*, *) "CHECKCJECKi  ", i, "    ", x1, " ", y1, " ", xt, " ", yt, "  RT  ", rt, " ", long00, " ", long2&
            , "  PS  ", (x2-x1)*(xt-x1)+(y2-y1)*(yt-y1)
      if (i==1) points_control(1, 1) = 0._pr
      if (i==1) points_control(1, 2) = 0._pr
      points_control(1, 1) = points_control(1, 1)+xt
      points_control(1, 2) = points_control(1, 2)+yt
      enddo
      !if (kt>80) then
      !points_control(1, 1) = points_control(1, 1)/(nint(0.1*nl)-1._pr)
      !points_control(1, 2) = points_control(1, 2)/(nint(0.1*nl)-1._pr)
      !endif
      endif

      write(*, *) "POINTS_COURBE_EQUAL  L==1  ", points_control(1, 1), " ", points_control(1, 2)
      points_control(1, 3) = zz(nint(zslice))
      points_courbe(1, 1) = points_control(1, 1)
      points_courbe(1, 2) = points_control(1, 2)
      points_courbe(1, 3) = zz(nint(zslice))


      tb = dtb
      l = 1
      long = 0._pr
      do while ((tb<1._pr).and.(l+1<Ns+1))
         l = l+1
         call pointsBezierN3D(points_control, tb, px, py, pz)
         points_courbe(l, 1) = px
         points_courbe(l, 2) = py
         points_courbe(l, 3) = pz
         long = long + sqrt((points_courbe(l, 1)-points_courbe(l-1, 1))**2 + (points_courbe(l, 2)-points_courbe(l-1, 2))**2 +&
               (points_courbe(l, 3)-points_courbe(l-1, 3))**2)
         tb = tb+dtb
      enddo
      if (l==Ns-1) then
         write(*, *) " ET VOILA"
         l=l+1
         points_courbe(size(points_courbe_equal, 1), 1) = points_control(size(points_control, 1), 1)
         points_courbe(size(points_courbe_equal, 1), 2) = points_control(size(points_control, 1), 2)
         points_courbe(size(points_courbe_equal, 1), 3) = points_control(size(points_control, 1), 3)
         long = long + sqrt((points_courbe(l, 1)-points_courbe(l-1, 1))**2 + (points_courbe(l, 2)-points_courbe(l-1, 2))**2 +&
               (points_courbe(l, 3)-points_courbe(l-1, 3))**2)
      endif
      if (.not.(l==Ns)) write(*, *) "WARNING2  p ", l
      !    write(*, *) "IMPORTANT  ", x1, " ", y1, " ", x2, " ", y2, " ", xt, " ", yt, " ", dist(x1, y1, xt, yt), " ", dist(x2, y2, xt, yt), " ", delta, " ", rt&
      !, " ", rt**2
      write(*, *) "IMPORTANT  ", long00, " ", long2, " ", long, "  RT ", rt, " ", sqrt((x1-xt)**2+(y1-yt)**2)&
            , "  coord ", x1, " ", y1, " ", x2, " ", y2
      long2 = long

      !**
      !! Tail extrapolation 
      !**
      l=1
      do i=1, l!1!2!nint(0.01*nl)-1 !!0.1
      x1 = points_control(i, 1)
      y1 = points_control(i, 2)
      x2 = points_control(2, 1)!nint(0.01*nl), 1)
      y2 = points_control(2, 2)!nint(0.01*nl), 2)
      rt = abs(long00-long2)
      !if (abs(x1-x2)<eepsilon*dx) then
      if (abs(x1-x2)/dx<eepsilon) then
         xt = x1
         !if (rt**2<(y2-y1-rt)**2) yt = y1+rt
         !if (rt**2<(y2-y1+rt)**2) yt = y1-rt
         yt = y1+rt

         if (((y2-y1)*(yt-y1))>0._pr) then
            write(*, *) "JEPASSEICIAUSSI"
            yt = y1-rt
         endif
      else
         xt = rt/(sqrt(1+((y2-y1)**2)/((x2-x1)**2))) + x1
         yt = (y2-y1)/(x2-x1)*(xt-x1) + y1
         !if (.not.(rt**2<(x2-xt)**2+(y2-yt)**2)) then
         !   xt = -rt/(sqrt(1+((y2-y1)**2)/((x2-x1)**2))) + x1
         !   yt = (y2-y1)/(x2-x1)*(xt-x1) + y1
         !endif

         if (((x2-x1)*(xt-x1)+(y2-y1)*(yt-y1))>0._pr) then
            write(*, *) "JEPASSEICI"
            xt = -rt/(sqrt(1+((y2-y1)**2)/((x2-x1)**2))) + x1
            yt = (y2-y1)/(x2-x1)*(xt-x1) + y1
         endif
      endif
      write(*, *) "CHECKCJECKi  ", i, "    ", x1, " ", y1, " ", xt, " ", yt, "  RT  ", rt, " ", long00, " ", long2&
            , "  PS  ", (x2-x1)*(xt-x1)+(y2-y1)*(yt-y1)
      if (i==1) points_control(1, 1) = 0._pr
      if (i==1) points_control(1, 2) = 0._pr
      points_control(1, 1) = points_control(1, 1)+xt
      points_control(1, 2) = points_control(1, 2)+yt
      enddo
      points_control(1, 3) = zz(nint(zslice))
      points_courbe(1, 1) = points_control(1, 1)
      points_courbe(1, 2) = points_control(1, 2)
      points_courbe(1, 3) = zz(nint(zslice))

      tb = dtb
      l = 1
      long = 0._pr
      do while ((tb<1._pr).and.(l+1<Ns+1))
         l = l+1
         call pointsBezierN3D(points_control, tb, px, py, pz)
         points_courbe(l, 1) = px
         points_courbe(l, 2) = py
         points_courbe(l, 3) = pz
         long = long + sqrt((points_courbe(l, 1)-points_courbe(l-1, 1))**2 + (points_courbe(l, 2)-points_courbe(l-1, 2))**2 +&
               (points_courbe(l, 3)-points_courbe(l-1, 3))**2)
         tb = tb+dtb
      enddo
      if (l==Ns-1) then
         l=l+1
         points_courbe(size(points_courbe_equal, 1), 1) = points_control(size(points_control, 1), 1)
         points_courbe(size(points_courbe_equal, 1), 2) = points_control(size(points_control, 1), 2)
         points_courbe(size(points_courbe_equal, 1), 3) = points_control(size(points_control, 1), 3)
         long = long + sqrt((points_courbe(l, 1)-points_courbe(l-1, 1))**2 + (points_courbe(l, 2)-points_courbe(l-1, 2))**2 +&
               (points_courbe(l, 3)-points_courbe(l-1, 3))**2)
      endif
      if (.not.(l==Ns)) write(*, *) "WARNING2  p ", l
      write(*, *) "IMPORTANTBIS  ", long00, " ", long2, " ", long, "  RT ", rt, " ", sqrt((x1-xt)**2+(y1-yt)**2)&
            , "  coord ", x1, " ", y1, " ", x2, " ", y2


      !**
      !! Uniform spline approximation 
      !**
      long2 = long
      ds = long/(Ns-1)
      !ds = area*long2/(meshRatio*Ns-1)
      points_courbe_equal =0._pr
      points_courbe_equal(1, 1) = points_control(1, 1)
      points_courbe_equal(1, 2) = points_control(1, 2)
      points_courbe_equal(1, 3) = points_control(1, 3)
      l = 1
      long = 0._pr
      tb = 0._pr
      do while ((tb<1._pr).and.(l+1<Ns)) !+1))
      !do while ((tb<1._pr).and.(l+1<meshRatio*Ns+1)) !+1))
         l = l+1
         nt = 1
         s = 0._pr
         do while ((l-1)*ds-s>0._pr) 
            nt = nt+1
            s = s + sqrt((points_courbe(nt, 1)-points_courbe(nt-1, 1))**2 + (points_courbe(nt, 2)-points_courbe(nt-1, 2))**2 +&
                  (points_courbe(nt, 3)-points_courbe(nt-1, 3))**2)
         enddo
         tbm = (nt-2)*dtb
         tbp = (nt-1)*dtb
         tb = tbm
         s0 = s
         sinit = s - sqrt((points_courbe(nt, 1)-points_courbe(nt-1, 1))**2 + (points_courbe(nt, 2)-points_courbe(nt-1, 2))**2 +&
               (points_courbe(nt, 3)-points_courbe(nt-1, 3))**2)
         s = sinit
         bool = 0
         !do while ((abs((l-1)*ds-s)>eepsilon*dx).and.(bool==0)) !.and.((tbm+tbp)*0.5_pr<1._pr))
         do while ((abs((l-1)*ds-s)/dx>eepsilon).and.(bool==0)) !.and.((tbm+tbp)*0.5_pr<1._pr))
            tb = (tbm + tbp)*0.5_pr
            !           call pointsBezierN(points_control, tb, px, py)
            call pointsBezierN3D(points_control, tb, px, py, pz)
            s = sinit + sqrt((px-points_courbe(nt-1, 1))**2 + (py-points_courbe(nt-1, 2))**2 + (pz-points_courbe(nt-1, 3))**2)
            if ((l-1)*ds-s>0._pr) then
               tbm = tb
            else
               tbp = tb
            endif
            if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then
               bool = 1
            endif
         enddo
         !        call pointsBezierN(points_control, tb, px, py)
         call pointsBezierN3D(points_control, tb, px, py, pz)
         points_courbe_equal(l, 1) = px
         points_courbe_equal(l, 2) = py
         points_courbe_equal(l, 3) = pz
         long = long +&
               sqrt((points_courbe_equal(l, 1)-points_courbe_equal(l-1, 1))**2 +&
               (points_courbe_equal(l, 2)-points_courbe_equal(l-1, 2))**2 +&
               (points_courbe_equal(l, 3)-points_courbe_equal(l-1, 3))**2)
      enddo
   !     lp = l
   !     !ds = 0.9*long2/(0.8*Ns-1)
   !     !ds = 0.9*long2/(Ns-lp)
   !     ds = (1._pr-2._pr*area)*long2/((1._pr-meshRatio)*Ns-lp)
   !     tp = tb
   !     tb = 0._pr !dtb
   !     l = lp !1
   !     !do while ((tb<1._pr).and.(l+1<Ns+1))
   !     do while ((tb<1._pr).and.(l+1<(1._pr-meshRatio)*Ns+1))
   !        l = l+1
   !        nt = 1
   !        s = 0._pr
   !        !do while ((lp-1)*0.1*long2/(0.2*Ns-1)+(l-lp-1)*ds-s>0._pr) 
   !        do while ((lp-1)*area*long2/(meshRatio*Ns-1)+(l-lp)*ds-s>0._pr) 
   !           nt = nt+1
   !           s = s + sqrt((points_courbe(nt, 1)-points_courbe(nt-1, 1))**2 + (points_courbe(nt, 2)&
   !                -points_courbe(nt-1, 2))**2 + (points_courbe(nt-1, 3)-points_courbe(nt-1, 3))**2)
   !        enddo
   !        tbm = (nt-2)*dtb!-tp+dtb
   !        tbp = (nt-1)*dtb!-tp+dtb
   !        tb = tbm !(tbm + tbp)*0.5_pr
   !        s0 = s
   !        sinit = s - sqrt((points_courbe(nt, 1)-points_courbe(nt-1, 1))**2 + (points_courbe(nt, 2)&
   !             -points_courbe(nt-1, 2))**2 + (points_courbe(nt, 3)-points_courbe(nt-1, 3))**2)
   !        s = sinit
   !        bool = 0
   !        !do while ((abs((lp-1)*0.1*long2/(0.2*Ns-1)+(l-lp-1)*ds-s)>eepsilon).and.(bool==0)) !.and.((tbm+tbp)*0.5_pr<1._pr))
   !        do while ((abs((lp-1)*area*long2/(meshRatio*Ns-1)+(l-lp)*ds-s)>eepsilon).and.(bool==0)) !.and.((tbm+tbp)*0.5_pr<1._pr))
   !           tb = (tbm + tbp)*0.5_pr
   !           call pointsBezierN3D(points_control, tb, px, py, pz)
   !           s = sinit + sqrt((px-points_courbe(nt-1, 1))**2 + (py-points_courbe(nt-1, 2))**2&
   !                + (pz-points_courbe(nt-1, 3))**2)
   !           !if ((lp-1)*0.1*long2/(0.2*Ns-1)+(l-lp-1)*ds-s>0._pr) then
   !           if ((lp-1)*area*long2/(meshRatio*Ns-1)+(l-lp)*ds-s>0._pr) then
   !              tbm = tb
   !           else
   !              tbp = tb
   !           endif
   !           if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then
   !              bool = 1
   !           endif
   !        enddo
   !        call pointsBezierN3D(points_control, tb, px, py, pz)
   !        points_courbe_equal(l, 1) = px
   !        points_courbe_equal(l, 2) = py
   !        points_courbe_equal(l, 3) = pz
   !        long = long +&
   !             sqrt((points_courbe_equal(l, 1)-points_courbe_equal(l-1, 1))**2 +&
   !             (points_courbe_equal(l, 2)-points_courbe_equal(l-1, 2))**2 +&
   !             (points_courbe_equal(l, 3)-points_courbe_equal(l-1, 3))**2)
   !     enddo
   !     lpt = l
   !     ds = area*long2/(Ns-lpt)
   !     tp = tb
   !     tb = 0._pr !dtb
   !     l = lpt !1
   !     do while ((tb<1._pr).and.(l+1<Ns))
   !        l = l+1
   !        nt = 1
   !        s = 0._pr
   !        do while ((lp-1)*area*long2/(meshRatio*Ns-1)+(lpt-lp)*(1._pr-2._pr*area)*long2/((1._pr-meshRatio)*Ns-lp)+(l-lpt)*ds-s>0._pr) 
   !           nt = nt+1
   !           s = s + sqrt((points_courbe(nt, 1)-points_courbe(nt-1, 1))**2 + (points_courbe(nt, 2)&
   !                -points_courbe(nt-1, 2))**2 + (points_courbe(nt-1, 3)-points_courbe(nt-1, 3))**2)
   !        enddo
   !        tbm = (nt-2)*dtb!-tp+dtb
   !        tbp = (nt-1)*dtb!-tp+dtb
   !        tb = tbm !(tbm + tbp)*0.5_pr
   !        s0 = s
   !        sinit = s - sqrt((points_courbe(nt, 1)-points_courbe(nt-1, 1))**2 + (points_courbe(nt, 2)&
   !             -points_courbe(nt-1, 2))**2 + (points_courbe(nt, 3)-points_courbe(nt-1, 3))**2)
   !        s = sinit
   !        bool = 0
   !        do while ((abs((lp-1)*area*long2/(meshRatio*Ns-1)+(lpt-lp)*(1._pr-2._pr*area)*long2/((1._pr-meshRatio)*Ns-lp)&
   !             +(l-lpt)*ds-s)>eepsilon).and.(bool==0)) !.and.((tbm+tbp)*0.5_pr<1._pr))
   !           tb = (tbm + tbp)*0.5_pr
   !           call pointsBezierN3D(points_control, tb, px, py, pz)
   !           s = sinit + sqrt((px-points_courbe(nt-1, 1))**2 + (py-points_courbe(nt-1, 2))**2&
   !                + (pz-points_courbe(nt-1, 3))**2)
   !           if ((lp-1)*area*long2/(meshRatio*Ns-1)+(lpt-lp)*(1._pr-2._pr*area)*long2/((1._pr-meshRatio)*Ns-lp)+(l-lpt)*ds-s>0._pr)&
   !                then
   !              tbm = tb
   !           else
   !              tbp = tb
   !           endif
   !           if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then
   !              bool = 1
   !           endif
   !        enddo
   !        call pointsBezierN3D(points_control, tb, px, py, pz)
   !        points_courbe_equal(l, 1) = px
   !        points_courbe_equal(l, 2) = py
   !        points_courbe_equal(l, 3) = pz
   !        long = long +&
   !             sqrt((points_courbe_equal(l, 1)-points_courbe_equal(l-1, 1))**2 +&
   !             (points_courbe_equal(l, 2)-points_courbe_equal(l-1, 2))**2 +&
   !             (points_courbe_equal(l, 3)-points_courbe_equal(l-1, 3))**2)
   !     enddo
      if (l==Ns-1) then
         l = l+1
         write(*, *) " ET VOILA j y passe  "
         points_courbe_equal(l, 1) = points_control(size(points_control, 1), 1)
         points_courbe_equal(l, 2) = points_control(size(points_control, 1), 2)
         points_courbe_equal(l, 3) = points_control(size(points_control, 1), 3)
         long = long +&
               sqrt((points_courbe_equal(l, 1)-points_courbe_equal(l-1, 1))**2 +&
               (points_courbe_equal(l, 2)-points_courbe_equal(l-1, 2))**2 +&
               (points_courbe_equal(l, 3)-points_courbe_equal(l-1, 3))**2)
      endif
      if (.not.(l==Ns)) write(*, *) "WARNING  s ", l
      write(*, *) "CHECKPOINT 1 ", long2, " ", long, " ", l, " ", tb, " POINTS  ", points_courbe_equal(Ns, 1), " ", &
            points_courbe_equal(Ns, 2), " ", points_courbe(Ns, 1), " ", points_courbe(Ns, 2), " ", &
            points_control(size(points_control, 1), 1), " ", points_control(size(points_control, 1), 2)
      do l=nint(0.04*Ns), 1, -1
         x2 = points_courbe_equal(Ns-l-1, 1)
         y2 = points_courbe_equal(Ns-l-1, 2)
         x1 = points_courbe_equal(Ns-l, 1)
         y1 = points_courbe_equal(Ns-l, 2)
         !! RTMODIF
         rt = dist(x1, y1, x2, y2)
         !if (abs(x1-x2)<eepsilon*dx) then
         if (abs(x1-x2)/dx<eepsilon) then
            xt = x1
            !if (rt**2<(y2-y1-rt)**2) yt = y1+rt
            !if (rt**2<(y2-y1+rt)**2) yt = y1-rt
            yt = y1+rt

            if (((y2-y1)*(yt-y1))>0._pr) then
               yt = y1-rt
            endif
         else
            xt = rt/(sqrt(1+((y2-y1)**2)/((x2-x1)**2))) + x1
            yt = (y2-y1)/(x2-x1)*(xt-x1) + y1
            !if (.not.(rt**2<(x2-xt)**2+(y2-yt)**2)) then
            !   xt = -rt/(sqrt(1+((y2-y1)**2)/((x2-x1)**2))) + x1
            !   yt = (y2-y1)/(x2-x1)*(xt-x1) + y1
            !endif

            if (((x2-x1)*(xt-x1)+(y2-y1)*(yt-y1))>0._pr) then
               xt = -rt/(sqrt(1+((y2-y1)**2)/((x2-x1)**2))) + x1
               yt = (y2-y1)/(x2-x1)*(xt-x1) + y1
            endif
         endif
         points_courbe_equal(size(points_courbe_equal, 1)-l+1, 1) = xt
         points_courbe_equal(size(points_courbe_equal, 1)-l+1, 2) = yt
         points_courbe_equal(size(points_courbe_equal, 1)-l+1, 3) = zz(nint(zslice))
      enddo

      !call body_masscenter(points_courbe_equal, xg, yg)
      !call body_masscentering(points_courbe_equal, xg, yg, xgref, ygref)
      !call body_masscenter(points_courbe_equal, xg, yg)
      !call body_rotationdefPS(points_courbe_equal_ref, points_courbe_equal, xg, yg, alphadef, dx)

      !**
      !! Head and Tail linear approximations 
      !**
      dti = dsi/ds !distslice(1+Ni, 1)/(Ni-1)/dsi
      do l=2, Ni
         tb = - (l-1)*dti
         !        call courbeBezierN(points_courbe_equal, tb, Px, Py)
         call courbeBezierN3D(points_courbe_equal, tb, Px, Py, Pz)
         tail_courbe(l, 1) = Px
         tail_courbe(l, 2) = Py
         tail_courbe(l, 3) = Pz
      enddo
      tail_courbe(1, 1) = points_courbe_equal(1, 1) !points_control(1, 1)
      tail_courbe(1, 2) = points_courbe_equal(1, 2) !points_control(1, 2)
      tail_courbe(1, 3) = points_courbe_equal(1, 3) !points_control(1, 2)
      dtf = dsf/ds !*90) !dist(points_control(size(points_control, 1), 1), points_control(size(points_control, 1), 2), points_courbe(size(points_courbe, 1)-1, 1), points_courbe(size(points_courbe, 1)-1, 2)) !ds!/(100)!*ds)!*100 !distslice(Ns+Ni, 2)/(Nf-1)/dsf !(10*ds)
      do l=2, Nf
         tb = 1._pr + (l-1)*dtf
         !        call courbeBezierN(points_courbe_equal, tb, Px, Py)
         call courbeBezierN3D(points_courbe_equal, tb, Px, Py, Pz)
         head_courbe(l, 1) = Px
         head_courbe(l, 2) = Py
         head_courbe(l, 3) = Pz
      enddo
      head_courbe(1, 1) = points_courbe_equal(Ns, 1) !points_control(size(points_control, 1), 1) !points_courbe_equal(Ns, 1) !points_control(Ns, 1)
      head_courbe(1, 2) = points_courbe_equal(Ns, 2) !points_control(size(points_control, 1), 2) !points_courbe_equal(Ns, 2) !points_control(Ns, 2)
      head_courbe(1, 3) = points_courbe_equal(Ns, 3) !points_control(size(points_control, 1), 2) !points_courbe_equal(Ns, 2) !points_control(Ns, 2)
      long3 = 0._pr
      !do l=1, sum(skel)-1
      do l=1, nl-1
         long3 = long3 +sqrt((xx(midline(l+1, 1))-xx(midline(l, 1)))**2 + (yy(midline(l+1, 2))-yy(midline(l, 2)))**2)
      enddo

      !write(*, *) "Longlong : ", long, " ", long3, " ", sum(skel), " ", sizeSkel
      write(*, *) "Longlong : ", long, " ", long3, " ", nl, " ", sizeSkel, "          "&
      , dsi, " ", dsf, " ", ds, " ", long, " ", (Nf-1)*dsf, " ", disthead, "     ", Px-head_courbe(1, 1), " "&
      , points_courbe_equal(Ns, 1)-points_courbe_equal(Ns-1, 1)




      !**
      !! Shape reconstruction (Lagrangian markers) based on the deformed midline: each transverse slice remains orthogonal to the midline
      !**

      l=1
      indextab = 0
      !     do theta=1, nbtheta
      !        thetatab(theta) = 1
   !!!        indextab(theta, thetatab(theta)) = 1
      !        indextab(theta, 1) = 1
      !        slice(theta, 1, 1) = tail_courbe(Ni, 1)
      !        slice(theta, 1, 2) = tail_courbe(Ni, 2)
      !        slice(theta, 1, 3) = tail_courbe(Ni, 3)
      !        slice2(theta, 1, 1) = slice(theta, 1, 1)
      !        slice2(theta, 1, 2) = slice(theta, 1, 2)
      !        slice2(theta, 1, 3) = slice(theta, 1, 3)
      !     enddo
      do l=1, Ni-1
         if (l==1) then
            !call intersection(tail_courbe(Ni-l+1, 1), tail_courbe(Ni-l+1-1, 1), tail_courbe(Ni-l+1, 2), tail_courbe(Ni-l+1-1, 2)&
            !     , tail_courbe(Ni-l+1, 1), tail_courbe(Ni-l+1, 2), distslice(l, 1), distslice(l, 2)&
            !     , distslice(Ns+Ni+Nf+l, 1), distslice(Ns+Ni+Nf+l, 2)&
            !     , xTheta(1), yTheta(1), xTheta(2), yTheta(2), errorl, errorr, deltal, deltar) 
            call intersection(tail_courbe(Ni, 1), tail_courbe(Ni-2, 1), tail_courbe(Ni, 2), tail_courbe(Ni-2, 2)&
                  , tail_courbe(Ni-1, 1), tail_courbe(Ni-1, 2), distslice(2, 1), distslice(2, 2)&
                  , distslice(Ns+Ni+Nf+2, 1), distslice(Ns+Ni+Nf+2, 2)&
                  , xTheta(1), yTheta(1), xTheta(2), yTheta(2), errorl, errorr, deltal, deltar)
         else
            call intersection(tail_courbe(Ni-l+1+1, 1), tail_courbe(Ni-l+1-1, 1), tail_courbe(Ni-l+1+1, 2), &
                  tail_courbe(Ni-l+1-1, 2), tail_courbe(Ni-l+1, 1), tail_courbe(Ni-l+1, 2), &
                  distslice(l, 1), distslice(l, 2), distslice(Ns+Ni+Nf+l, 1), distslice(Ns+Ni+Nf+l, 2)&
                  , xTheta(1), yTheta(1), xTheta(2), yTheta(2), errorl, errorr, deltal, deltar) 
         endif
         if (l==1) then
            sinPhi = (xTheta(1)-tail_courbe(Ni-1, 1))/(&
                  sqrt((xTheta(1)-tail_courbe(Ni-1, 1))**2+(yTheta(1)-tail_courbe(Ni-1, 2))**2))
            cosPhi = (yTheta(1)-tail_courbe(Ni-1, 2))/(&
                  sqrt((xTheta(1)-tail_courbe(Ni-1, 1))**2+(yTheta(1)-tail_courbe(Ni-1, 2))**2))
         else
            sinPhi = (xTheta(1)-tail_courbe(Ni-l+1, 1))/(&
                  sqrt((xTheta(1)-tail_courbe(Ni-l+1, 1))**2+(yTheta(1)-tail_courbe(Ni-l+1, 2))**2))
            cosPhi = (yTheta(1)-tail_courbe(Ni-l+1, 2))/(&
                  sqrt((xTheta(1)-tail_courbe(Ni-l+1, 1))**2+(yTheta(1)-tail_courbe(Ni-l+1, 2))**2)) 
         endif
         do theta=1, nbtheta
            !           slice(theta, l, 1) = tail_courbe(Ni-l+1, 1) + valDist(l, theta)*cos(valTheta(l, theta))*sinPhi
            !           slice(theta, l, 2) = tail_courbe(Ni-l+1, 2) + valDist(l, theta)*cos(valTheta(l, theta))*cosPhi
            !           slice(theta, l, 3) = zslice + valDist(l, theta)*sin(valTheta(l, theta))
            slice(theta, l, 1) = tail_courbe(Ni-l+1, 1) + valDist(l, theta)*cosTheta_tab(l, theta)*sinPhi
            slice(theta, l, 2) = tail_courbe(Ni-l+1, 2) + valDist(l, theta)*cosTheta_tab(l, theta)*cosPhi
            slice(theta, l, 3) = zz(nint(zslice)) + valDist(l, theta)*sinTheta_tab(l, theta)
         enddo

         do theta=1, nbtheta
            if (l==1) then
               thetatab(theta) = 1
            else
               thetatab(theta) = thetatab(theta) + 1
            endif
            !!           indextab(theta, thetatab(theta)) = l
            indextab(theta, l) = 1
            !!           indextab(theta, l) = indextab(theta, l) +1
            slice2(theta, thetatab(theta), 1) = slice(theta, l, 1)
            slice2(theta, thetatab(theta), 2) = slice(theta, l, 2)
            slice2(theta, thetatab(theta), 3) = slice(theta, l, 3)
         enddo
   !!!        slice(1, l, 1) = xTheta(1)
   !!!        slice(1, l, 2) = yTheta(1)
   !!!        if (errorl==1) then
   !!!           slice(1, l, 1) = slice(1, l-1, 1)
   !!!           slice(1, l, 2) = slice(1, l-1, 2)
   !!!        endif
   !!!        if ((abs(slice(1, l, 1)-slice(1, l-1, 1))<1.e-2).and.(abs(slice(1, l, 2)-slice(1, l-1, 2))<1.e-2)) then
   !!!           slice(2, l, 1) = slice(2, l-1, 1)
   !!!           slice(2, l, 2) = slice(2, l-1, 2)
   !!!        else
   !!!           slice(2, l, 1) = -slice(1, l, 1) + 2*tail_courbe(Ni-l+1, 1)
   !!!           slice(2, l, 2) = -slice(1, l, 2) + 2*tail_courbe(Ni-l+1, 2)
   !!!        endif
   !!!        do theta=1, nbtheta
   !!!            slice(theta, l, 1) = slice(theta, l-1, 1)
   !!!            slice(theta, l, 2) = slice(theta, l-1, 2)
   !!!            slice(theta, l, 3) = slice(theta, l-1, 3)
   !!!        enddo
         !        if (det(tail_courbe(Ni-l+1+1, 1), tail_courbe(Ni-l+1+1, 2), slice2(lpli, 1, 1), slice2(lpli, 1, 2), tail_courbe(Ni-l+1, 1), tail_courbe(Ni-l+1, 2))*det(tail_courbe(Ni-l+1+1, 1), tail_courbe(Ni-l+1+1, 2), slice2(lpli, 1, 1), slice2(lpli, 1, 2), xLeft, yLeft)>= 0._pr) then
   !!!        do theta=1, nbtheta
   !!!           thetatab(theta) = thetatab(theta) + 1
   !!!        enddo
   !!!        slice2(1, thetatab(1), 1) = xTheta(1)
   !!!        slice2(1, thetatab(1), 2) = yTheta(1)
   !!!        slice2(2, thetatab(2), 1) = slice(2, l, 1)
   !!!        slice2(2, thetatab(2), 2) = slice(2, l, 2)
   !!!        do theta=2, nbtheta
   !!!           slice2(theta, thetatab(2), 1) = slice(theta, l, 1)
   !!!           slice2(theta, thetatab(2), 2) = slice(theta, l, 2)
   !!!           slice2(theta, thetatab(2), 3) = slice(theta, l, 3)
   !!!        enddo
         !        endif
   !!!        if (errorl==1) then
   !!!           write(*, *) "blabla 1", l
   !!!           slice2(1, thetatab(1), 1) = slice2(1, thetatab(1)-1, 1)
   !!!           slice2(1, thetatab(1), 2) = slice2(1, thetatab(1)-1, 2)
   !!!           slice2(1, thetatab(1), 3) = slice2(1, thetatab(1)-1, 3)
   !!!        endif
   !!!        if (errorr==1) then
   !!!           write(*, *) "blabla 2", l
   !!!           slice2(2, thetatab(2), 1) = slice2(2, thetatab(2)-1, 1)
   !!!           slice2(2, thetatab(2), 2) = slice2(2, thetatab(2)-1, 2)
   !!!           slice2(2, thetatab(2), 3) = slice2(2, thetatab(2)-1, 3)
   !!!        endif
         if ((l>1).and.((errorl==1).or.(errorr==1))) then
            do theta=1, nbtheta
               slice(theta, l, 1) = slice(theta, l-1, 1)
               slice(theta, l, 2) = slice(theta, l-1, 2)
               slice(theta, l, 3) = slice(theta, l-1, 3)
               slice2(theta, thetatab(theta), 1) = slice2(theta, thetatab(theta)-1, 1)
               slice2(theta, thetatab(theta), 2) = slice2(theta, thetatab(theta)-1, 2)
               slice2(theta, thetatab(theta), 3) = slice2(theta, thetatab(theta)-1, 3)
            enddo
            write(*, *) "error initialsegment : ", theta, " ", l, " ", errorl, " ", errorr
         endif


         !if (errorl==1) write(*, *) "WARNING : ", errorl, " ", l
      enddo
      write(*, *) "VOILA VOILA"
      call intersection(tail_courbe(2, 1), points_courbe_equal(2, 1), tail_courbe(2, 2), points_courbe_equal(2, 2), &
            points_courbe_equal(1, 1), points_courbe_equal(1, 2), distslice(Ni, 1), distslice(Ni+1, 2), &
            distslice(Ns+Ni+Nf+Ni+1, 1), distslice(Ns+Ni+Nf+Ni+1, 2), xTheta(1), yTheta(1), xTheta(2), &
            yTheta(2), errorl, errorr, deltal, deltar) 
      l = Ni
      sinPhi = (xTheta(1)-points_courbe_equal(1, 1))/(&
            sqrt((xTheta(1)-points_courbe_equal(1, 1))**2+(yTheta(1)-points_courbe_equal(1, 2))**2))
      cosPhi = (yTheta(1)-points_courbe_equal(1, 2))/(&
            sqrt((xTheta(1)-points_courbe_equal(1, 1))**2+(yTheta(1)-points_courbe_equal(1, 2))**2)) 
   !!!     slice(1, l, 1) = points_courbe_equal(1, 1) + valDist(l, 1)*cos(valTheta(l, 1))*sinPhi
   !!!     slice(1, l, 2) = points_courbe_equal(1, 2) + valDist(l, 1)*cos(valTheta(l, 1))*cosPhi
   !!!     slice(2, l, 1) = points_courbe_equal(1, 1) + valDist(l, 2)*cos(valTheta(l, 2))*sinPhi 
   !!!     slice(2, l, 2) = points_courbe_equal(1, 2) + valDist(l, 2)*cos(valTheta(l, 2))*cosPhi 
      do theta=1, nbtheta
         !        slice(theta, l, 1) = points_courbe_equal(1, 1) + valDist(l, theta)*cos(valTheta(l, theta))*sinPhi
         !        slice(theta, l, 2) = points_courbe_equal(1, 2) + valDist(l, theta)*cos(valTheta(l, theta))*cosPhi
         !        slice(theta, l, 3) = zslice + valDist(l, theta)*sin(valTheta(l, theta))
         slice(theta, l, 1) = points_courbe_equal(1, 1) + valDist(l, theta)*cosTheta_tab(l, theta)*sinPhi
         slice(theta, l, 2) = points_courbe_equal(1, 2) + valDist(l, theta)*cosTheta_tab(l, theta)*cosPhi
         slice(theta, l, 3) = zz(nint(zslice)) + valDist(l, theta)*sinTheta_tab(l, theta)
      enddo

   !!!     if ((det(tail_courbe(2, 1), tail_courbe(2, 2), slice2(1, thetatab(1), 1), slice2(1, thetatab(1), 2), points_courbe_equal(1, 1), &
   !!!points_courbe_equal(1, 2))*det(tail_courbe(2, 1), tail_courbe(2, 2), slice2(1, thetatab(1), 1), slice2(1, thetatab(1), 2), slice(1, l, 1), &
   !!!slice(1, l, 2))>= 0._pr).or.(thetatab(1)<Ni)) then
   !!!        thetatab(1) = thetatab(1) + 1
   !!!        slice2(1, thetatab(1), 1) = slice(1, l, 1)
   !!!        slice2(1, thetatab(1), 2) = slice(1, l, 2)
   !!!     endif
      do theta=1, nbtheta
         if ((det(tail_courbe(2, 1), tail_courbe(2, 2), slice2(theta, thetatab(theta), 1), slice2(theta, thetatab(theta), 2), &
               points_courbe_equal(1, 1), points_courbe_equal(1, 2))*det(tail_courbe(2, 1), tail_courbe(2, 2)&
               , slice2(theta, thetatab(theta), 1), &
               slice2(theta, thetatab(theta), 2), slice(theta, l, 1), slice(theta, l, 2))>= 0._pr).or.(thetatab(theta)<Ni)) then
            thetatab(theta) = thetatab(theta) + 1
            !!           indextab(theta, thetatab(theta)) = l
            indextab(theta, l) = 1
            slice2(theta, thetatab(theta), 1) = slice(theta, l, 1)
            slice2(theta, thetatab(theta), 2) = slice(theta, l, 2)
            slice2(theta, thetatab(theta), 3) = slice(theta, l, 3)
         endif
      enddo
      if (errorl==1) write(*, *) "WARNING : ", errorl, " ", l
      do l=2, Ns-1
         sinPhir = sinPhi
         cosPhir = cosPhi
         call intersection(points_courbe_equal(l-1, 1), points_courbe_equal(l+1, 1), &
               points_courbe_equal(l-1, 2), points_courbe_equal(l+1, 2), points_courbe_equal(l, 1), &
               points_courbe_equal(l, 2), distslice(l+Ni, 1), distslice(l+Ni, 2), &
               distslice(Ns+Ni+Nf+l+Ni, 1), distslice(Ns+Ni+Nf+l+Ni, 2), xTheta(1), yTheta(1)&
               , xTheta(2), yTheta(2), errorl, errorr, deltal, deltar) 
         sinPhi = (xTheta(1)-points_courbe_equal(l, 1))/(&
               sqrt((xTheta(1)-points_courbe_equal(l, 1))**2+(yTheta(1)-points_courbe_equal(l, 2))**2))
         cosPhi = (yTheta(1)-points_courbe_equal(l, 2))/(&
               sqrt((xTheta(1)-points_courbe_equal(l, 1))**2+(yTheta(1)-points_courbe_equal(l, 2))**2)) 
         sinAlpha = sinPhi
         cosAlpha = cosPhi
         sinPhil = sinPhi
         cosPhil = cosPhi
         !!        slice(1, l+Ni-1, 1) = points_courbe_equal(l, 1) + valDist(l+Ni-1, 1)*cos(valTheta(l+Ni-1, 1))*sinPhi
         !!        slice(1, l+Ni-1, 2) = points_courbe_equal(l, 2) + valDist(l+Ni-1, 1)*cos(valTheta(l+Ni-1, 1))*cosPhi
         !!        slice(1, l+Ni-1, 3) = zslice + valDist(l+Ni-1, 1)*sin(valTheta(l+Ni-1, 1))
         !!        slice(2, l+Ni-1, 1) = points_courbe_equal(l, 1) + valDist(l+Ni-1, 2)*cos(valTheta(l+Ni-1, 2))*sinPhi 
         !!        slice(2, l+Ni-1, 2) = points_courbe_equal(l, 2) + valDist(l+Ni-1, 2)*cos(valTheta(l+Ni-1, 2))*cosPhi 
         do theta=1, nbtheta
            !           slice(theta, l+Ni-1, 1) = points_courbe_equal(l, 1) + valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*sinPhi
            !           slice(theta, l+Ni-1, 2) = points_courbe_equal(l, 2) + valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*cosPhi
            !           slice(theta, l+Ni-1, 3) = zslice + valDist(l+Ni-1, theta)*sin(valTheta(l+Ni-1, theta))
            slice(theta, l+Ni-1, 1) = points_courbe_equal(l, 1) + valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*sinPhi
            slice(theta, l+Ni-1, 2) = points_courbe_equal(l, 2) + valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*cosPhi
            slice(theta, l+Ni-1, 3) = zz(nint(zslice)) + valDist(l+Ni-1, theta)*sinTheta_tab(l+Ni-1, theta)
         enddo
         !!        if ((det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2), slice2(1, thetatab(1), 1), slice2(1, thetatab(1), 2)&
         !!, points_courbe_equal(l, 1), points_courbe_equal(l, 2))*det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2)&
         !!, slice2(1, thetatab(1), 1), slice2(1, thetatab(1), 2), slice(1, l+Ni-1, 1), slice(1, l+Ni-1, 2))>= 0._pr).or.(thetatab(1)<Ni)) then
         !!           thetatab(1) = thetatab(1) + 1
         !!           slice2(1, thetatab(1), 1) = slice(1, l+Ni-1, 1)
         !!           slice2(1, thetatab(1), 2) = slice(1, l+Ni-1, 2)
         !!        endif
         bool = 0
         do theta=1, nbtheta
            !           if (((theta==1).or.(theta==1+nbtheta/2)).and.(.not.((det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2)&
            !, slice(theta, l+Ni-1-1, 1)&
            !, slice(theta, l+Ni-1-1, 2), points_courbe_equal(l, 1), points_courbe_equal(l, 2))*det(points_courbe_equal(l-1, 1)&
            !, points_courbe_equal(l-1, 2), slice(theta, l+Ni-1-1, 1), slice(theta, l+Ni-1-1, 2), slice(theta, l+Ni-1, 1)&
            !, slice(theta, l+Ni-1, 2))>= 0._pr).or.(thetatab(theta)<Ni)))) then
            if (.not.(det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2)&
                                 !if (((theta==1).or.(theta==1+nbtheta/2)).and.(.not.(det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2)&
                                 !, points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*sinPhir&
                                 !, points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*cosPhir&
                                 !, points_courbe_equal(l, 1), points_courbe_equal(l, 2))*det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2)&
                                 !, points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*sinPhir&
                                 !, points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*cosPhir&
                                 !, points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*sinPhil&
                                 !, points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*cosPhil)>= 0._pr)) then
                  , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*sinPhir&
                  , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*cosPhir&
                  , points_courbe_equal(l, 1), points_courbe_equal(l, 2))*det(points_courbe_equal(l-1, 1)&
                  , points_courbe_equal(l-1, 2)&
                  , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*sinPhir&
                  , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*cosPhir&
                  , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*sinPhil&
                  , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*cosPhil)>= 0._pr)) then

               !              sinTheta = -det(points_courbe_equal(l, 1), points_courbe_equal(l, 2), slice(theta, l+Ni-1, 1), slice(theta, l+Ni-1, 2)&
               !, slice(theta, l+Ni-1-1, 1), slice(theta, l+Ni-1-1, 2))/(norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
               !, slice(theta, l+Ni-1, 1), slice(theta, l+Ni-1, 2))*norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
               !, slice(theta, l+Ni-1-1, 1), slice(theta, l+Ni-1-1, 2)))
               !              cosTheta = dotProd(points_courbe_equal(l, 1), points_courbe_equal(l, 2), slice(theta, l+Ni-1, 1), slice(theta, l+Ni-1, 2)&
               !, slice(theta, l+Ni-1-1, 1), slice(theta, l+Ni-1-1, 2))/(norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
               !, slice(theta, l+Ni-1, 1), slice(theta, l+Ni-1, 2))*norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
               !, slice(theta, l+Ni-1-1, 1), slice(theta, l+Ni-1-1, 2)))
               !sinTheta = -det(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
               !     , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*sinPhil&
               !     , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*cosPhil&
               !     , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*sinPhir&
               !     , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*cosPhir)&
               !     /(norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
               !     , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*sinPhil&
               !     , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*cosPhil)&
               !     *norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
               !     , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*sinPhir&
               !     , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*cosPhir))
               !cosTheta = dotProd(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
               !     , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*sinPhil&
               !     , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*cosPhil&
               !     , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*sinPhir&
               !     , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*cosPhir)&
               !     /(norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
               !     , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*sinPhil&
               !     , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*cosPhil)&
               !     *norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
               !     , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*sinPhir&
               !     , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*cosPhir))
               sinTheta = -det(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
                     , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*sinPhil&
                     , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*cosPhil&
                     , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*sinPhir&
                     , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*cosPhir)&
                     /(norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
                     , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*sinPhil&
                     , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*cosPhil)&
                     *norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
                     , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*sinPhir&
                     , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*cosPhir))
               cosTheta = dotProd(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
                     , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*sinPhil&
                     , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*cosPhil&
                     , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*sinPhir&
                     , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*cosPhir)&
                     /(norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
                     , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*sinPhil&
                     , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*cosPhil)&
                     *norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
                     , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*sinPhir&
                     , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*cosPhir))
               if (cosTheta*oldC<0._pr) write(*, *) "ATTENTION COS NEGATIF ", l, " ", l+Ni-1, " ", theta, "     kt ", kt
               !               if (bool==0) then
               !               endif
               if ((abs(sinTheta/cosTheta)<abs(oldS/oldC)).and.(bool==1)) then
                  sinTheta = oldS
                  cosTheta = oldC
               endif
               bool = 1
               oldS = sinTheta
               oldC = cosTheta
               write(*, *) "cccTESTEST  ", theta, " ", l, " ", l+Ni-1, " ", sinTheta, " ", cosTheta, " ", sinTheta/cosTheta
            endif
         enddo
         do theta=1, nbtheta
            !           if ((det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2), slice2(theta, thetatab(theta), 1)&
            !, slice2(theta, thetatab(theta), 2), points_courbe_equal(l, 1), points_courbe_equal(l, 2))*det(points_courbe_equal(l-1, 1)&
            !, points_courbe_equal(l-1, 2), slice2(theta, thetatab(theta), 1), slice2(theta, thetatab(theta), 2), slice(theta, l+Ni-1, 1)&
            !, slice(theta, l+Ni-1, 2))>= 0._pr).or.(thetatab(theta)<Ni)) then
            !              xr = slice2(theta, thetatab(theta), 1)
            !              yr = slice2(theta, thetatab(theta), 2)
            xr = slice(theta, l+Ni-1, 1)
            yr = slice(theta, l+Ni-1, 2)
            thetatab(theta) = thetatab(theta) + 1
            !!              indextab(theta, thetatab(theta)) = l+Ni-1
            indextab(theta, l+Ni-1) = 1
            slice2(theta, thetatab(theta), 1) = slice(theta, l+Ni-1, 1)
            slice2(theta, thetatab(theta), 2) = slice(theta, l+Ni-1, 2)
            slice2(theta, thetatab(theta), 3) = slice(theta, l+Ni-1, 3)
            if (bool==1) then
               slice(theta, l+Ni-1, 1) = cosTheta*(xr-points_courbe_equal(l, 1)) + sinTheta*(yr-points_courbe_equal(l, 2)) +&
                     points_courbe_equal(l, 1)
               slice(theta, l+Ni-1, 2) = -sinTheta*(xr-points_courbe_equal(l, 1)) + cosTheta*(yr-points_courbe_equal(l, 2)) +&
                     points_courbe_equal(l, 2)
               sinPhi = sinAlpha*cosTheta + sinTheta*cosAlpha
               cosPhi = cosAlpha*cosTheta - sinAlpha*sinTheta
            endif
            !           else
            !
            !           if (bool==1) then
            !              xr = slice2(theta, thetatab(theta), 1)
            !              yr = slice2(theta, thetatab(theta), 2)
            !              thetatab(theta) = thetatab(theta) + 1
            !              indextab(theta, l+Ni-1) = 1
            !              slice2(theta, thetatab(theta), 1) = cosTheta*xr + sinTheta*yr
            !              slice2(theta, thetatab(theta), 2) = -sinTheta*xr + cosTheta*yr
            !            endif
            !           endif
         enddo

         if (errorl==1) write(*, *) "WARNING11 : ", errorl, " ", l
      enddo
      l = Ns
      sinPhir = sinPhi
      cosPhir = cosPhi
      call intersection(points_courbe_equal(l-1, 1), head_courbe(1+1, 1), points_courbe_equal(Ns-1, 2), head_courbe(1+1, 2)&
                                 !, head_courbe(1, 1), head_courbe(1, 2), distslice(Ni+Ns, 1), distslice(1+Ni+Ns, 2), distslice(Ns+Ni+Nf+Ni+Ns, 1)&
            , head_courbe(1, 1), head_courbe(1, 2), distslice(Ni+Ns, 1), distslice(Ni+Ns+1, 2), distslice(Ns+Ni+Nf+Ni+Ns, 1)&
            , distslice(Ns+Ni+Nf+1+Ni+Ns, 2), xTheta(1), yTheta(1), xTheta(2), yTheta(2), errorl, errorr, deltal, deltar) 
      sinPhi = (xTheta(1)-head_courbe(1, 1))/(sqrt((xTheta(1)-head_courbe(1, 1))**2+(yTheta(1)-head_courbe(1, 2))**2))
      cosPhi = (yTheta(1)-head_courbe(1, 2))/(sqrt((xTheta(1)-head_courbe(1, 1))**2+(yTheta(1)-head_courbe(1, 2))**2)) 
      sinAlpha = sinPhi
      cosAlpha = cosPhi
      sinPhil = sinPhi
      cosPhil = cosPhi
   !!!     slice(1, l+Ni-1, 1) = head_courbe(1, 1) + valDist(l+Ni-1, 1)*cos(valTheta(l+Ni-1, 1))*sinPhi
   !!!     slice(1, l+Ni-1, 2) = head_courbe(1, 2) + valDist(l+Ni-1, 1)*cos(valTheta(l+Ni-1, 1))*cosPhi
   !!!     slice(2, l+Ni-1, 1) = head_courbe(1, 1) + valDist(l+Ni-1, 2)*cos(valTheta(l+Ni-1, 2))*sinPhi
   !!!     slice(2, l+Ni-1, 2) = head_courbe(1, 2) + valDist(l+Ni-1, 2)*cos(valTheta(l+Ni-1, 2))*cosPhi
      do theta=1, nbtheta
         !        slice(theta, l+Ni-1, 1) = head_courbe(1, 1) + valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*sinPhi
         !        slice(theta, l+Ni-1, 2) = head_courbe(1, 2) + valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*cosPhi
         !        slice(theta, l+Ni-1, 3) = zslice + valDist(l+Ni-1, theta)*sin(valTheta(l+Ni-1, theta))
         slice(theta, l+Ni-1, 1) = head_courbe(1, 1) + valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*sinPhi
         slice(theta, l+Ni-1, 2) = head_courbe(1, 2) + valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*cosPhi
         slice(theta, l+Ni-1, 3) = zz(nint(zslice)) + valDist(l+Ni-1, theta)*sinTheta_tab(l+Ni-1, theta)
      enddo

   !!!     if (det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2), slice2(1, thetatab(1), 1), slice2(1, thetatab(1), 2)&
   !!!, points_courbe_equal(l, 1), points_courbe_equal(l, 2))*det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2)&
   !!!, slice2(1, thetatab(1), 1), slice2(1, thetatab(1), 2), slice(1, l+Ni-1, 1), slice(1, l+Ni-1, 2))>= 0._pr) then
   !!!        thetatab(1) = thetatab(1) + 1
   !!!        slice2(1, thetatab(1), 1) = slice(1, l+Ni-1, 1)
   !!!        slice2(1, thetatab(1), 2) = slice(1, l+Ni-1, 2)
   !!!     endif
      bool = 0
      do theta=1, nbtheta
         !        if (((theta==1).or.(theta==1+nbtheta/2)).and.(.not.(det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2)&
         !, slice(theta, l+Ni-1-1, 1)&
         !, slice(theta, l+Ni-1-1, 2), points_courbe_equal(l, 1), points_courbe_equal(l, 2))*det(points_courbe_equal(l-1, 1)&
         !, points_courbe_equal(l-1, 2), slice(theta, l+Ni-1-1, 1), slice(theta, l+Ni-1-1, 2), slice(theta, l+Ni-1, 1)&
         !, slice(theta, l+Ni-1, 2))>= 0._pr))) then
         if (.not.(det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2)&
                                 !, points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*sinPhir&
                                 !, points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*cosPhir&
                                 !, points_courbe_equal(l, 1), points_courbe_equal(l, 2))*det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2)&
                                 !, points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*sinPhir&
                                 !, points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*cosPhir&
                                 !, points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*sinPhil&
                                 !, points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*cosPhil)>= 0._pr)) then
               , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*sinPhir&
               , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*cosPhir&
               , points_courbe_equal(l, 1), points_courbe_equal(l, 2))*det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2)&
               , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*sinPhir&
               , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*cosPhir&
               , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*sinPhil&
               , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*cosPhil)>= 0._pr)) then

            !              sinTheta = -det(points_courbe_equal(l, 1), points_courbe_equal(l, 2), slice(theta, l+Ni-1, 1), slice(theta, l+Ni-1, 2)&
            !, slice(theta, l+Ni-1-1, 1), slice(theta, l+Ni-1-1, 2))/(norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
            !, slice(theta, l+Ni-1, 1), slice(theta, l+Ni-1, 2))*norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
            !, slice(theta, l+Ni-1-1, 1), slice(theta, l+Ni-1-1, 2)))
            !              cosTheta = dotProd(points_courbe_equal(l, 1), points_courbe_equal(l, 2), slice(theta, l+Ni-1, 1), slice(theta, l+Ni-1, 2)&
            !, slice(theta, l+Ni-1-1, 1), slice(theta, l+Ni-1-1, 2))/(norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
            !, slice(theta, l+Ni-1, 1), slice(theta, l+Ni-1, 2))*norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
            !, slice(theta, l+Ni-1-1, 1), slice(theta, l+Ni-1-1, 2)))
            !sinTheta = -det(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
            !     , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*sinPhil&
            !     , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*cosPhil&
            !     , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*sinPhir&
            !     , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*cosPhir)&
            !     /(norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
            !     , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*sinPhil&
            !     , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*cosPhil)&
            !     *norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
            !     , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*sinPhir&
            !     , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*cosPhir))
            !cosTheta = dotProd(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
            !     , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*sinPhil&
            !     , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*cosPhil&
            !     , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*sinPhir&
            !     , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*cosPhir)&
            !     /(norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
            !     , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*sinPhil&
            !     , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*cosPhil)&
            !     *norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
            !     , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*sinPhir&
            !     , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cos(valTheta(l+Ni-1-1, theta))*cosPhir))
            sinTheta = -det(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
                  , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*sinPhil&
                  , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*cosPhil&
                  , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*sinPhir&
                  , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*cosPhir)&
                  /(norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
                  , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*sinPhil&
                  , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*cosPhil)&
                  *norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
                  , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*sinPhir&
                  , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*cosPhir))
            cosTheta = dotProd(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
                  , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*sinPhil&
                  , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*cosPhil&
                  , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*sinPhir&
                  , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*cosPhir)&
                  /(norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
                  , points_courbe_equal(l, 1)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*sinPhil&
                  , points_courbe_equal(l, 2)+valDist(l+Ni-1, theta)*cosTheta_tab(l+Ni-1, theta)*cosPhil)&
                  *norme(points_courbe_equal(l, 1), points_courbe_equal(l, 2)&
                  , points_courbe_equal(l-1, 1)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*sinPhir&
                  , points_courbe_equal(l-1, 2)+valDist(l+Ni-1-1, theta)*cosTheta_tab(l+Ni-1-1, theta)*cosPhir))
            if (cosTheta*oldC<0._pr) write(*, *) "ATTENTION COS NEGATIF ", l, " ", l+Ni-1, " ", theta, "     kt ", kt
            !               if (bool==0) then
            !               endif
            if ((abs(sinTheta/cosTheta)<abs(oldS/oldC)).and.(bool==1)) then
               sinTheta = oldS
               cosTheta = oldC
            endif
            bool = 1
            oldS = sinTheta
            oldC = cosTheta
            write(*, *) "bbbTESTEST  ", theta, " ", l, " ", l+Ni-1, " ", sinTheta, " ", cosTheta, " ", sinTheta/cosTheta
         endif
      enddo
      !boolPhi = 0
      do theta=1, nbtheta
         !        if (det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2), slice2(theta, thetatab(theta), 1)&
         !, slice2(theta, thetatab(theta), 2), points_courbe_equal(l, 1), points_courbe_equal(l, 2))*det(points_courbe_equal(l-1, 1)&
         !, points_courbe_equal(l-1, 2), slice2(theta, thetatab(theta), 1), slice2(theta, thetatab(theta), 2), slice(theta, l+Ni-1, 1)&
         !, slice(theta, l+Ni-1, 2))>= 0._pr) then
         !              xr = slice2(theta, thetatab(theta), 1)
         !              yr = slice2(theta, thetatab(theta), 2)
         xr = slice(theta, l+Ni-1, 1)
         yr = slice(theta, l+Ni-1, 2)
         thetatab(theta) = thetatab(theta) + 1
         !!           indextab(theta, thetatab(theta)) = l+Ni-1
         indextab(theta, l+Ni-1) = 1
         slice2(theta, thetatab(theta), 1) = slice(theta, l+Ni-1, 1)
         slice2(theta, thetatab(theta), 2) = slice(theta, l+Ni-1, 2)
         slice2(theta, thetatab(theta), 3) = slice(theta, l+Ni-1, 3)
         if (bool==1) then
            slice(theta, l+Ni-1, 1) = cosTheta*(xr-points_courbe_equal(l, 1)) + sinTheta*(yr-points_courbe_equal(l, 2)) +&
                  points_courbe_equal(l, 1)
            slice(theta, l+Ni-1, 2) = -sinTheta*(xr-points_courbe_equal(l, 1)) + cosTheta*(yr-points_courbe_equal(l, 2)) +&
                  points_courbe_equal(l, 2)
            !              if (boolPhi==0) then
            sinPhi = sinAlpha*cosTheta + sinTheta*cosAlpha
            cosPhi = cosAlpha*cosTheta - sinAlpha*sinTheta
            !              sinPhi = sinPhi*cosTheta + sinTheta*cosPhi
            !              cosPhi = cosPhi*cosTheta - sinPhi*sinTheta
            !                boolPhi = 1
            !                endif
         endif
         !           else
         !
         !           if (bool==1) then
         !              xr = slice2(theta, thetatab(theta), 1)
         !              yr = slice2(theta, thetatab(theta), 2)
         !              thetatab(theta) = thetatab(theta) + 1
         !              indextab(theta, l+Ni-1) = 1
         !              slice2(theta, thetatab(theta), 1) = cosTheta*xr + sinTheta*yr
         !              slice2(theta, thetatab(theta), 2) = -sinTheta*xr + cosTheta*yr
         !            endif
         !           endif
      enddo
      if (errorl==1) write(*, *) "WARNING22 : ", errorl, " ", l
      do l=2, Nf
         sinPhir = sinPhi
         cosPhir = cosPhi
         if (l==Nf) then
            !call intersection(head_courbe(l-1, 1), head_courbe(l, 1), head_courbe(l-1, 2), head_courbe(l, 2), head_courbe(l, 1)&
            !     , head_courbe(l, 2), distslice(l+Ni+Ns, 1), distslice(l+Ni+Ns, 2)&
            !     , distslice(Ns+Ni+Nf+l+Ni+Ns, 1), distslice(Ns+Ni+Nf+l+Ni+Ns, 2), xTheta(1)&
            !     , yTheta(1), xTheta(2), yTheta(2), errorl, errorr, deltal, deltar) 
            call intersection(head_courbe(Nf-2, 1), head_courbe(Nf-1, 1), head_courbe(Nf-2, 2), head_courbe(Nf-1, 2)&
                  , head_courbe(Nf-1, 1), head_courbe(Nf-1, 2), distslice(Nf-1+Ni+Ns, 1), distslice(Nf-1+Ni+Ns, 2)&
                  , distslice(Ns+Ni+Nf+Nf-1+Ni+Ns, 1), distslice(Ns+Ni+Nf+Nf-1+Ni+Ns, 2), xTheta(1)&
                  , yTheta(1), xTheta(2), yTheta(2), errorl, errorr, deltal, deltar)
         else
            call intersection(head_courbe(l-1, 1), head_courbe(l+1, 1), head_courbe(l-1, 2), head_courbe(l+1, 2), head_courbe(l, 1)&
                  , head_courbe(l, 2), distslice(l+Ni+Ns, 1), distslice(l+Ni+Ns, 2)&
                  , distslice(Ns+Ni+Nf+l+Ni+Ns, 1), distslice(Ns+Ni+Nf+l+Ni+Ns, 2), xTheta(1)&
                  , yTheta(1), xTheta(2), yTheta(2), errorl, errorr, deltal, deltar) 
         endif
         if (l==Nf) then
            sinPhi = (xTheta(1)-head_courbe(Nf-1, 1))/(sqrt((xTheta(1)-head_courbe(Nf-1, 1))**2+(yTheta(1)-&
            head_courbe(Nf-1, 2))**2))
            cosPhi = (yTheta(1)-head_courbe(Nf-1, 2))/(sqrt((xTheta(1)-head_courbe(Nf-1, 1))**2+(yTheta(1)-&
            head_courbe(Nf-1, 2))**2))
         else
            sinPhi = (xTheta(1)-head_courbe(l, 1))/(sqrt((xTheta(1)-head_courbe(l, 1))**2+(yTheta(1)-&
            head_courbe(l, 2))**2))
            cosPhi = (yTheta(1)-head_courbe(l, 2))/(sqrt((xTheta(1)-head_courbe(l, 1))**2+(yTheta(1)-&
            head_courbe(l, 2))**2)) 
         endif
         sinAlpha = sinPhi
         cosAlpha = cosPhi
         sinPhil = sinPhi
         cosPhil = cosPhi
   !!!        slice(1, l+Ni-1+Ns-1, 1) = head_courbe(l, 1) + valDist(l+Ni-1+Ns-1, 1)*cos(valTheta(l+Ni-1+Ns-1, 1))*sinPhi
   !!!        slice(1, l+Ni-1+Ns-1, 2) = head_courbe(l, 2) + valDist(l+Ni-1+Ns-1, 1)*cos(valTheta(l+Ni-1+Ns-1, 1))*cosPhi
   !!!        slice(2, l+Ni-1+Ns-1, 1) = head_courbe(l, 1) + valDist(l+Ni-1+Ns-1, 2)*cos(valTheta(l+Ni-1+Ns-1, 2))*sinPhi
   !!!        slice(2, l+Ni-1+Ns-1, 2) = head_courbe(l, 2) + valDist(l+Ni-1+Ns-1, 2)*cos(valTheta(l+Ni-1+Ns-1, 2))*cosPhi
         do theta=1, nbtheta
            !           slice(theta, l+Ni-1+Ns-1, 1) = head_courbe(l, 1) + valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*sinPhi
            !           slice(theta, l+Ni-1+Ns-1, 2) = head_courbe(l, 2) + valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*cosPhi
            !           slice(theta, l+Ni-1+Ns-1, 3) = zslice + valDist(l+Ni-1+Ns-1, theta)*sin(valTheta(l+Ni-1+Ns-1, theta))
            slice(theta, l+Ni-1+Ns-1, 1) = head_courbe(l, 1) + valDist(l+Ni-1+Ns-1, theta)*cosTheta_tab(l+Ni-1+Ns-1, theta)*sinPhi
            slice(theta, l+Ni-1+Ns-1, 2) = head_courbe(l, 2) + valDist(l+Ni-1+Ns-1, theta)*cosTheta_tab(l+Ni-1+Ns-1, theta)*cosPhi
            slice(theta, l+Ni-1+Ns-1, 3) = zz(nint(zslice)) + valDist(l+Ni-1+Ns-1, theta)*sinTheta_tab(l+Ni-1+Ns-1, theta)
         enddo
   !!!        if (det(head_courbe(l-1, 1), head_courbe(l-1, 2), slice2(1, thetatab(1), 1), slice2(1, thetatab(1), 2), head_courbe(l, 1), &
   !!!head_courbe(l, 2))*det(head_courbe(l-1, 1), head_courbe(l-1, 2), slice2(1, thetatab(1), 1), slice2(1, thetatab(1), 2), slice(1, l+Ni-1+Ns-1, 1), &
   !!!slice(1, l+Ni-1+Ns-1, 2))>= 0._pr) then
   !!!           thetatab(1) = thetatab(1) + 1
   !!!           slice2(1, thetatab(1), 1) = slice(1, l+ni-1+Ns-1, 1)
   !!!           slice2(1, thetatab(1), 2) = slice(1, l+Ni-1+Ns-1, 2)
   !!!        endif
         bool = 0
         do theta=1, nbtheta
            !           if (.not.(det(head_courbe(l-1, 1), head_courbe(l-1, 2), slice2(theta, thetatab(theta), 1), slice2(theta, thetatab(theta), 2), &
            !head_courbe(l, 1), head_courbe(l, 2))*det(head_courbe(l-1, 1), head_courbe(l-1, 2), slice2(theta, thetatab(theta), 1), &
            !slice2(theta, thetatab(theta), 2), slice(theta, l+Ni-1+Ns-1, 1), slice(theta, l+Ni-1+Ns-1, 2))>= 0._pr)) then
            !              bool = 1
            !              sinTheta = det(head_courbe(l, 1), head_courbe(l, 2), slice(theta, l+Ni-1+Ns-1, 1), slice(theta, l+Ni-1+Ns-1, 2)&
            !, slice2(theta, thetatab(theta), 1), slice2(theta, thetatab(theta), 2))/(norme(head_courbe(l, 1), head_courbe(l, 2)&
            !, slice(theta, l+Ni-1+Ns-1, 1), slice(theta, l+Ni-1+Ns-1, 2))*norme(head_courbe(l, 1), head_courbe(l, 2)&
            !, slice2(theta, thetatab(theta), 1), slice2(theta, thetatab(theta), 2)))
            !cosTheta = dotProd(head_courbe(l, 1), head_courbe(l, 2), slice(theta, l+Ni-1+Ns-1, 1), slice(theta, l+Ni-1+Ns-1, 2)&
            !, slice2(theta, thetatab(theta), 1), slice2(theta, thetatab(theta), 2))/(norme(head_courbe(l, 1), head_courbe(l, 2)&
            !, slice(theta, l+Ni-1+Ns-1, 1), slice(theta, l+Ni-1+Ns-1, 2))*norme(head_courbe(l, 1), head_courbe(l, 2)&
            !, slice2(theta, thetatab(theta), 1), slice2(theta, thetatab(theta), 2)))

            !xr = dist3D(head_courbe(l-1, 1), head_courbe(l-1, 2), head_courbe(l-1, 3), slice(theta, l+Ni-1+Ns-1-1, 1), slice(theta, l+Ni-1+Ns-1-1, 2)&
            !, slice(theta, l+Ni-1+Ns-1-1, 3))
            !xl = dist3D(head_courbe(l, 1), head_courbe(l, 2), head_courbe(l, 3), slice(theta, l+Ni-1+Ns-1, 1), slice(theta, l+Ni-1+Ns-1, 2)&
            !, slice(theta, l+Ni-1+Ns-1, 3))
            !           if (theta==1+nbtheta/2) write(*, *) "mytheta     ", l, "    ", slice(theta, l+Ni-1+Ns-1, 1), " ", slice(theta, l+Ni-1+Ns-1, 2)&
            !, " ", head_courbe(l, 1)+valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*sinPhil, " ", head_courbe(l, 2)+&
            !valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*cosPhil, " ", head_courbe(l, 1)+xl*cos(valTheta(l+Ni-1+Ns-1, theta))&
            !*sinPhil, " ", head_courbe(l, 2)+xl*cos(valTheta(l+Ni-1+Ns-1, theta))*cosPhil, "  dist  ", valDist(l+Ni-1+Ns-1, theta), " ", xl, " "&
            !, valTheta(l+Ni-1+Ns-1, theta)
            !           if (theta==1+nbtheta/2) write(*, *) "mytheta     ", l, "    ", slice(theta, l+Ni-1+Ns-1-1, 1), " ", slice(theta, l+Ni-1+Ns-1-1, 2)&
            !, " ", head_courbe(l-1, 1)+valDist(l+Ni-1+Ns-1-1, theta)*cos(valTheta(l+Ni-1+Ns-1-1, theta))*sinPhir, " ", head_courbe(l-1, 2)+&
            !valDist(l+Ni-1+Ns-1-1, theta)*cos(valTheta(l+Ni-1+Ns-1-1, theta))*cosPhir, " ", head_courbe(l-1, 1)+&
            !xr*cos(valTheta(l+Ni-1+Ns-1-1, theta))*sinPhir, " ", head_courbe(l-1, 2)+xr*cos(valTheta(l+Ni-1+Ns-1-1, theta))*cosPhir, "  dist  "&
            !, valDist(l+Ni-1+Ns-1-1, theta), " ", xr, " ", valTheta(l+Ni-1+Ns-1-1, theta)

            !           if (((theta==1).or.(theta==1+nbtheta/2)).and.(.not.(det(head_courbe(l-1, 1), head_courbe(l-1, 2)&
            if (.not.(det(head_courbe(l-1, 1), head_courbe(l-1, 2)&
                                 !, head_courbe(l-1, 1)+valDist(l+Ni-1+Ns-1-1, theta)*cos(valTheta(l+Ni-1+Ns-1-1, theta))*sinPhir&
                                 !, head_courbe(l-1, 2)+valDist(l+Ni-1+Ns-1-1, theta)*cos(valTheta(l+Ni-1+Ns-1-1, theta))*cosPhir&
                                 !, head_courbe(l, 1), head_courbe(l, 2))*det(head_courbe(l-1, 1), head_courbe(l-1, 2)&
                                 !, head_courbe(l-1, 1)+valDist(l+Ni-1+Ns-1-1, theta)*cos(valTheta(l+Ni-1+Ns-1-1, theta))*sinPhir&
                                 !, head_courbe(l-1, 2)+valDist(l+Ni-1+Ns-1-1, theta)*cos(valTheta(l+Ni-1+Ns-1-1, theta))*cosPhir&
                                 !, head_courbe(l, 1)+valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*sinPhil&
                                 !, head_courbe(l, 2)+valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*cosPhil)&
                  , head_courbe(l-1, 1)+valDist(l+Ni-1+Ns-1-1, theta)*cosTheta_tab(l+Ni-1+Ns-1-1, theta)*sinPhir&
                  , head_courbe(l-1, 2)+valDist(l+Ni-1+Ns-1-1, theta)*cosTheta_tab(l+Ni-1+Ns-1-1, theta)*cosPhir&
                  , head_courbe(l, 1), head_courbe(l, 2))*det(head_courbe(l-1, 1), head_courbe(l-1, 2)&
                  , head_courbe(l-1, 1)+valDist(l+Ni-1+Ns-1-1, theta)*cosTheta_tab(l+Ni-1+Ns-1-1, theta)*sinPhir&
                  , head_courbe(l-1, 2)+valDist(l+Ni-1+Ns-1-1, theta)*cosTheta_tab(l+Ni-1+Ns-1-1, theta)*cosPhir&
                  , head_courbe(l, 1)+valDist(l+Ni-1+Ns-1, theta)*cosTheta_tab(l+Ni-1+Ns-1, theta)*sinPhil&
                  , head_courbe(l, 2)+valDist(l+Ni-1+Ns-1, theta)*cosTheta_tab(l+Ni-1+Ns-1, theta)*cosPhil)&
                  >= 0._pr)) then
               !           if (((theta==1).or.(theta==1+nbtheta/2)).and.(.not.(det(head_courbe(l-1, 1), head_courbe(l-1, 2)&
               !, head_courbe(l-1, 1)+xr*cos(valTheta(l+Ni-1+Ns-1-1, theta))*sinPhir, head_courbe(l-1, 2)+xr*cos(valTheta(l+Ni-1+Ns-1-1, theta))*cosPhir, &
               !head_courbe(l, 1), head_courbe(l, 2))*det(head_courbe(l-1, 1), head_courbe(l-1, 2), head_courbe(l-1, 1)+&
               !xr*cos(valTheta(l+Ni-1+Ns-1-1, theta))*sinPhir&
               !, head_courbe(l-1, 2)+xr*cos(valTheta(l+Ni-1+Ns-1-1, theta))*cosPhir, head_courbe(l, 1)+xl*cos(valTheta(l+Ni-1+Ns-1, theta))*sinPhil&
               !, head_courbe(l, 2)+xl*cos(valTheta(l+Ni-1+Ns-1, theta))*cosPhil)&
               !>= 0._pr))) then
               !           if (((theta==1).or.(theta==1+nbtheta/2)).and.(.not.(det(head_courbe(l-1, 1), head_courbe(l-1, 2)&
               !, slice(theta, l+Ni-1+Ns-1-1, 1), slice(theta, l+Ni-1+Ns-1-1, 2), &
               !head_courbe(l, 1), head_courbe(l, 2))*det(head_courbe(l-1, 1), head_courbe(l-1, 2), slice(theta, l+Ni-1+Ns-1-1, 1), &
               !slice(theta, l+Ni-1+Ns-1-1, 2), slice(theta, l+Ni-1+Ns-1, 1), slice(theta, l+Ni-1+Ns-1, 2))>= 0._pr))) then
               !sinTheta = -det(head_courbe(l, 1), head_courbe(l, 2)&
               !     , head_courbe(l, 1)+valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*sinPhil&
               !     , head_courbe(l, 2)+valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*cosPhil&
               !     , head_courbe(l-1, 1)+valDist(l+Ni-1+Ns-1-1, theta)*cos(valTheta(l+Ni-1+Ns-1-1, theta))*sinPhir&
               !     , head_courbe(l-1, 2)+valDist(l+Ni-1+Ns-1-1, theta)*cos(valTheta(l+Ni-1+Ns-1-1, theta))*cosPhir)&
               !     /(norme(head_courbe(l, 1), head_courbe(l, 2)&
               !     , head_courbe(l, 1)+valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*sinPhil&
               !     , head_courbe(l, 2)+valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*cosPhil)&
               !     *norme(head_courbe(l, 1), head_courbe(l, 2)&
               !     , head_courbe(l-1, 1)+valDist(l+Ni-1+Ns-1-1, theta)*cos(valTheta(l+Ni-1+Ns-1-1, theta))*sinPhir&
               !     , head_courbe(l-1, 2)+valDist(l+Ni-1+Ns-1-1, theta)*cos(valTheta(l+Ni-1+Ns-1-1, theta))*cosPhir))
               !cosTheta = dotProd(head_courbe(l, 1), head_courbe(l, 2)&
               !     , head_courbe(l, 1)+valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*sinPhil&
               !     , head_courbe(l, 2)+valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*cosPhil&
               !     , head_courbe(l-1, 1)+valDist(l+Ni-1+Ns-1-1, theta)*cos(valTheta(l+Ni-1+Ns-1-1, theta))*sinPhir&
               !     , head_courbe(l-1, 2)+valDist(l+Ni-1+Ns-1-1, theta)*cos(valTheta(l+Ni-1+Ns-1-1, theta))*cosPhir)&
               !     /(norme(head_courbe(l, 1), head_courbe(l, 2)&
               !     , head_courbe(l, 1)+valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*sinPhil&
               !     , head_courbe(l, 2)+valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*cosPhil)&
               !     *norme(head_courbe(l, 1), head_courbe(l, 2)&
               !     , head_courbe(l-1, 1)+valDist(l+Ni-1+Ns-1-1, theta)*cos(valTheta(l+Ni-1+Ns-1-1, theta))*sinPhir&
               !     , head_courbe(l-1, 2)+valDist(l+Ni-1+Ns-1-1, theta)*cos(valTheta(l+Ni-1+Ns-1-1, theta))*cosPhir))
               sinTheta = -det(head_courbe(l, 1), head_courbe(l, 2)&
                     , head_courbe(l, 1)+valDist(l+Ni-1+Ns-1, theta)*cosTheta_tab(l+Ni-1+Ns-1, theta)*sinPhil&
                     , head_courbe(l, 2)+valDist(l+Ni-1+Ns-1, theta)*cosTheta_tab(l+Ni-1+Ns-1, theta)*cosPhil&
                     , head_courbe(l-1, 1)+valDist(l+Ni-1+Ns-1-1, theta)*cosTheta_tab(l+Ni-1+Ns-1-1, theta)*sinPhir&
                     , head_courbe(l-1, 2)+valDist(l+Ni-1+Ns-1-1, theta)*cosTheta_tab(l+Ni-1+Ns-1-1, theta)*cosPhir)&
                     /(norme(head_courbe(l, 1), head_courbe(l, 2)&
                     , head_courbe(l, 1)+valDist(l+Ni-1+Ns-1, theta)*cosTheta_tab(l+Ni-1+Ns-1, theta)*sinPhil&
                     , head_courbe(l, 2)+valDist(l+Ni-1+Ns-1, theta)*cosTheta_tab(l+Ni-1+Ns-1, theta)*cosPhil)&
                     *norme(head_courbe(l, 1), head_courbe(l, 2)&
                     , head_courbe(l-1, 1)+valDist(l+Ni-1+Ns-1-1, theta)*cosTheta_tab(l+Ni-1+Ns-1-1, theta)*sinPhir&
                     , head_courbe(l-1, 2)+valDist(l+Ni-1+Ns-1-1, theta)*cosTheta_tab(l+Ni-1+Ns-1-1, theta)*cosPhir))
               cosTheta = dotProd(head_courbe(l, 1), head_courbe(l, 2)&
                     , head_courbe(l, 1)+valDist(l+Ni-1+Ns-1, theta)*cosTheta_tab(l+Ni-1+Ns-1, theta)*sinPhil&
                     , head_courbe(l, 2)+valDist(l+Ni-1+Ns-1, theta)*cosTheta_tab(l+Ni-1+Ns-1, theta)*cosPhil&
                     , head_courbe(l-1, 1)+valDist(l+Ni-1+Ns-1-1, theta)*cosTheta_tab(l+Ni-1+Ns-1-1, theta)*sinPhir&
                     , head_courbe(l-1, 2)+valDist(l+Ni-1+Ns-1-1, theta)*cosTheta_tab(l+Ni-1+Ns-1-1, theta)*cosPhir)&
                     /(norme(head_courbe(l, 1), head_courbe(l, 2)&
                     , head_courbe(l, 1)+valDist(l+Ni-1+Ns-1, theta)*cosTheta_tab(l+Ni-1+Ns-1, theta)*sinPhil&
                     , head_courbe(l, 2)+valDist(l+Ni-1+Ns-1, theta)*cosTheta_tab(l+Ni-1+Ns-1, theta)*cosPhil)&
                     *norme(head_courbe(l, 1), head_courbe(l, 2)&
                     , head_courbe(l-1, 1)+valDist(l+Ni-1+Ns-1-1, theta)*cosTheta_tab(l+Ni-1+Ns-1-1, theta)*sinPhir&
                     , head_courbe(l-1, 2)+valDist(l+Ni-1+Ns-1-1, theta)*cosTheta_tab(l+Ni-1+Ns-1-1, theta)*cosPhir))
               if (cosTheta*oldC<0._pr) write(*, *) "ATTENTION COS NEGATIF ", l, " ", l+Ni-1+Ns-1, " ", theta, "    kt ", kt
               !               if (bool==0) then
               !               endif
               if ((abs(sinTheta/cosTheta)<abs(oldS/oldC)).and.(bool==1)) then
                  sinTheta = oldS
                  cosTheta = oldC
               endif
               bool = 1
               oldS = sinTheta
               oldC = cosTheta

               !              sinTheta = -det(head_courbe(l, 1), head_courbe(l, 2), slice(theta, l+Ni-1+Ns-1, 1), slice(theta, l+Ni-1+Ns-1, 2)&
               !, slice(theta, l+Ni-1+Ns-1-1, 1), slice(theta, l+Ni-1+Ns-1-1, 2))/(norme(head_courbe(l, 1), head_courbe(l, 2)&
               !, slice(theta, l+Ni-1+Ns-1, 1), slice(theta, l+Ni-1+Ns-1, 2))*norme(head_courbe(l, 1), head_courbe(l, 2)&
               !, slice(theta, l+Ni-1+Ns-1-1, 1), slice(theta, l+Ni-1+Ns-1-1, 2)))
               !              cosTheta = dotProd(head_courbe(l, 1), head_courbe(l, 2), slice(theta, l+Ni-1+Ns-1, 1), slice(theta, l+Ni-1+Ns-1, 2)&
               !, slice(theta, l+Ni-1+Ns-1-1, 1), slice(theta, l+Ni-1+Ns-1-1, 2))/(norme(head_courbe(l, 1), head_courbe(l, 2)&
               !, slice(theta, l+Ni-1+Ns-1, 1), slice(theta, l+Ni-1+Ns-1, 2))*norme(head_courbe(l, 1), head_courbe(l, 2)&
               !, slice(theta, l+Ni-1+Ns-1-1, 1), slice(theta, l+Ni-1+Ns-1-1, 2)))
               write(*, *) "aaaTESTEST  ", theta, " ", l, " ", l+Ni-1+Ns-1, " ", sinTheta, " ", cosTheta, " ", &
                     sinPhir, " ", cosPhir, " ", sinPhil, " ", cosPhil, " ", sinTheta/cosTheta
               !              write(*, *) "aaaTESTEST  ", theta, " ", l, " ", l+Ni-1+Ns-1, " ", sinTh(theta), " ", cosTh(theta), " ", sinPhir, " ", cosPhir&
               !, " ", sinPhil, " ", cosPhil
            endif
         enddo
         !boolPhi = 0
         do theta=1, nbtheta
            !           if (det(head_courbe(l-1, 1), head_courbe(l-1, 2), slice2(theta, thetatab(theta), 1), slice2(theta, thetatab(theta), 2), &
            !head_courbe(l, 1), head_courbe(l, 2))*det(head_courbe(l-1, 1), head_courbe(l-1, 2), slice2(theta, thetatab(theta), 1), &
            !slice2(theta, thetatab(theta), 2), slice(theta, l+Ni-1+Ns-1, 1), slice(theta, l+Ni-1+Ns-1, 2))>= 0._pr) then
            !              xr = slice2(theta, thetatab(theta), 1)
            !              yr = slice2(theta, thetatab(theta), 2)
            xr = slice(theta, l+Ni-1+Ns-1, 1)
            yr = slice(theta, l+Ni-1+Ns-1, 2)
            !           if (theta==101) write(*, *) "CHECK ouch appartient  ", theta, " ", l, " ", l+Ni-1+Ns-1, " ", slice(theta, l+Ni-1+Ns-1, 1)&
            !, " ", slice(theta, l+Ni-1+Ns-1, 2), " ", slice(theta, l+Ni-1+Ns-1, 3), " ", cos(valTheta(l+Ni-1+Ns-1, theta))&
            !, " ", sin(valTheta(l+Ni-1+Ns-1, theta)), " ", indextheta(l)
            thetatab(theta) = thetatab(theta) + 1
            !!              indextab(theta, thetatab(theta)) = l+Ni-1+Ns-1
            indextab(theta, l+Ni-1+Ns-1) = 1
            slice2(theta, thetatab(theta), 1) = slice(theta, l+Ni-1+Ns-1, 1)
            slice2(theta, thetatab(theta), 2) = slice(theta, l+Ni-1+Ns-1, 2)
            slice2(theta, thetatab(theta), 3) = slice(theta, l+Ni-1+Ns-1, 3)
            !              xr = slice2(theta, thetatab(theta), 1)
            !              yr = slice2(theta, thetatab(theta), 2)
            if (bool==1) then
               !              write(*, *) "oint : ", slice2(theta, thetatab(theta), 1), " ", slice2(theta, thetatab(theta), 2)
               !              slice2(theta, thetatab(theta), 1) = cosTheta*xr + sinTheta*yr
               !              slice2(theta, thetatab(theta), 2) = -sinTheta*xr + cosTheta*yr
               !write(*, *) "AH oui ", theta, " ", l, " ", l+Ni-1+Ns-1, " ", slice(theta, l+Ni-1+Ns-1, 1), " ", slice(theta, l+Ni-1+Ns-1, 2)
               slice(theta, l+Ni-1+Ns-1, 1) = cosTheta*(xr-head_courbe(l, 1)) + sinTheta*(yr-head_courbe(l, 2)) + head_courbe(l, 1)
               slice(theta, l+Ni-1+Ns-1, 2) = -sinTheta*(xr-head_courbe(l, 1)) + cosTheta*(yr-head_courbe(l, 2)) + head_courbe(l, 2)
               !write(*, *) "ah oui ", theta, " ", l, " ", l+Ni-1+Ns-1, " ", slice(theta, l+Ni-1+Ns-1, 1), " ", slice(theta, l+Ni-1+Ns-1, 2)
               !            slice(theta, l+Ni-1+Ns-1, 1) = cosTh(theta)*(xr-head_courbe(l, 1)) + sinTh(theta)*(yr-head_courbe(l, 2)) + head_courbe(l, 1)
               !            slice(theta, l+Ni-1+Ns-1, 2) = -sinTh(theta)*(xr-head_courbe(l, 1)) + cosTh(theta)*(yr-head_courbe(l, 2)) + head_courbe(l, 2)
               !              write(*, *) "oint : ", slice2(theta, thetatab(theta), 1), " ", slice2(theta, thetatab(theta), 2)
               !              slice2(theta, thetatab(theta), 1) = cosTheta*xr - sinTheta*yr
               !              slice2(theta, thetatab(theta), 2) = sinTheta*xr + cosTheta*yr
               !              write(*, *) "oint : ", slice2(theta, thetatab(theta), 1), " ", slice2(theta, thetatab(theta), 2)
               !if (boolPhi==0) then
               sinPhi = sinAlpha*cosTheta + sinTheta*cosAlpha
               cosPhi = cosAlpha*cosTheta - sinAlpha*sinTheta
               !              sinPhi = sinPhi*cosTheta + sinTheta*cosPhi
               !              cosPhi = cosPhi*cosTheta - sinPhi*sinTheta
               !              sinPhi = sinAlpha*cosTh(theta) + sinTh(theta)*cosAlpha
               !              cosPhi = cosAlpha*cosTh(theta) - sinAlpha*sinTh(theta)
               !  boolPhi = 1
               !  endif
            endif
            !           else
            !
            !           if (bool==1) then
            !              xr = slice2(theta, thetatab(theta), 1)
            !              yr = slice2(theta, thetatab(theta), 2)
            !!              xr = slice(theta, l+Ni-1+Ns-1, 1)
            !!              yr = slice(theta, +Ni-1+Ns-1, 2)
            !              thetatab(theta) = thetatab(theta) + 1
            !              indextab(theta, l+Ni-1+Ns-1) = 1
            !              slice2(theta, thetatab(theta), 1) = cosTheta*xr + sinTheta*yr
            !              slice2(theta, thetatab(theta), 2) = -sinTheta*xr + cosTheta*yr
            !            endif
            !           endif
         enddo
         if (errorl==1) write(*, *) "WARNING33 : ", errorl, " ", l
      enddo
      !     l = Nf
      !     do theta=1, nbtheta
      !        thetatab(theta) = thetatab(theta) + 1
   !!!        indextab(theta, thetatab(theta)) = Nf+Ni-1+Ns-1
      !        indextab(theta, l+Ni-1+Ns-1) = 1
      !     slice(theta, l+Ni-1+Ns-1, 1) = head_courbe(l, 1)
      !     slice(theta, l+Ni-1+Ns-1, 2) = head_courbe(l, 2)
      !     slice(theta, l+Ni-1+Ns-1, 3) = head_courbe(l, 3)
      !        slice2(theta, thetatab(theta), 1) = head_courbe(l, 1)
      !        slice2(theta, thetatab(theta), 2) = head_courbe(l, 2)
      !        slice2(theta, thetatab(theta), 3) = head_courbe(l, 3)
      !     enddo
      !!   do theta=1, nbtheta
      !!        indextheta(1) = 1
   !!!!!        slice(theta, 1, 1) = tail_courbe(Ni, 1)
   !!!!!        slice(theta, 1, 2) = tail_courbe(Ni, 2)
   !!!!!        slice(theta, 1, 3) = tail_courbe(Ni, 3)
      !!        vect(theta, 1, 1) = slice(theta, 1, 1)
      !!        vect(theta, 1, 2) = slice(theta, 1, 2)
      !!        vect(theta, 1, 3) = slice(theta, 1, 3)
      !!     do l=2, Ni-1
   !!!!!        call intersection(tail_courbe(Ni-l+1+1, 1), tail_courbe(Ni-l+1-1, 1), tail_courbe(Ni-l+1+1, 2), tail_courbe(Ni-l+1-1, 2)&
   !!!!!, tail_courbe(Ni-l+1, 1), tail_courbe(Ni-l+1, 2), distslice(l, 1), distslice(l, 2), distslice(Ns+Ni+Nf+l, 1), distslice(Ns+Ni+Nf+l, 2)&
   !!!!!, xTheta(1), yTheta(1), xTheta(2), yTheta(2), errorl, errorr, deltal, deltar) 
   !!!!!        sinPhi = (xTheta(1)-tail_courbe(Ni-l+1, 1))/(&
   !!!!!sqrt((xTheta(1)-tail_courbe(Ni-l+1, 1))**2+(yTheta(1)-tail_courbe(Ni-l+1, 2))**2))
   !!!!!        cosPhi = (yTheta(1)-tail_courbe(Ni-l+1, 2))/(&
   !!!!!sqrt((xTheta(1)-tail_courbe(Ni-l+1, 1))**2+(yTheta(1)-tail_courbe(Ni-l+1, 2))**2)) 
   !!!!!           slice(theta, l, 1) = tail_courbe(Ni-l+1, 1) + valDist(l, theta)*cos(valTheta(l, theta))*sinPhi
   !!!!!           slice(theta, l, 2) = tail_courbe(Ni-l+1, 2) + valDist(l, theta)*cos(valTheta(l, theta))*cosPhi
   !!!!!           slice(theta, l, 3) = zslice + valDist(l, theta)*sin(valTheta(l, theta))
      !!
      !!           indextheta(l) = indextheta(l) + 1
      !!           vect(indextheta(l), l, 1) = slice(theta, l, 1)
      !!           vect(indextheta(l), l, 2) = slice(theta, l, 2)
      !!           vect(indextheta(l), l, 3) = slice(theta, l, 3)
      !!        if ((errorl==1).or.(errorr==1)) then
      !!           vect(indextheta(l), l, 1) = vect(indextheta(l), l-1, 1)
      !!           vect(indextheta(l), l, 2) = vect(indextheta(l), l-1, 2)
      !!           vect(indextheta(l), l, 3) = vect(indextheta(l), l-1, 3)
      !!   endif
      !!     enddo
      !!     l = Ni
   !!!!!     call intersection(tail_courbe(2, 1), points_courbe_equal(2, 1), tail_courbe(2, 2), points_courbe_equal(2, 2), points_courbe_equal(1, 1)&
   !!!!!, points_courbe_equal(1, 2), distslice(Ni, 1), distslice(Ni+1, 2), distslice(Ns+Ni+Nf+Ni, 1), distslice(Ns+Ni+Nf+Ni+1, 2), xTheta(1), yTheta(1)&
   !!!!!, xTheta(2), yTheta(2), errorl, errorr, deltal, deltar) 
   !!!!!     sinPhi = (xTheta(1)-points_courbe_equal(1, 1))/(&
   !!!!!sqrt((xTheta(1)-points_courbe_equal(1, 1))**2+(yTheta(1)-points_courbe_equal(1, 2))**2))
   !!!!!     cosPhi = (yTheta(1)-points_courbe_equal(1, 2))/(&
   !!!!!sqrt((xTheta(1)-points_courbe_equal(1, 1))**2+(yTheta(1)-points_courbe_equal(1, 2))**2)) 
   !!!!!        slice(theta, l, 1) = points_courbe_equal(1, 1) + valDist(l, theta)*cos(valTheta(l, theta))*sinPhi
   !!!!!        slice(theta, l, 2) = points_courbe_equal(1, 2) + valDist(l, theta)*cos(valTheta(l, theta))*cosPhi
   !!!!!        slice(theta, l, 3) = zslice + valDist(l, theta)*sin(valTheta(l, theta))
      !!        if ((det(tail_courbe(2, 1), tail_courbe(2, 2), vect(indextheta(l), l, 1), vect(indextheta(l), l, 2), &
      !!points_courbe_equal(1, 1), points_courbe_equal(1, 2))*det(tail_courbe(2, 1), tail_courbe(2, 2), vect(indextheta(l), l, 1), &
      !!vect(indextheta(l), l, 2), slice(theta, l, 1), slice(theta, l, 2))>= 0._pr).or.(l<Ni)) then
      !!           indextheta(l) = indextheta(l) + 1
      !!           vect(indextheta(l), l, 1) = slice(theta, l, 1)
      !!           vect(indextheta(l), l, 2) = slice(theta, l, 2)
      !!           vect(indextheta(l), l, 3) = slice(theta, l, 3)
      !!        endif
      !!     do l=2, Ns-1
   !!!!!        call intersection(points_courbe_equal(l-1, 1), points_courbe_equal(l+1, 1), points_courbe_equal(l-1, 2), &
   !!!!!points_courbe_equal(l+1, 2), points_courbe_equal(l, 1), points_courbe_equal(l, 2), distslice(l+Ni, 1), distslice(l+Ni, 2), &
   !!!!!distslice(Ns+Ni+Nf+l+Ni, 1), distslice(Ns+Ni+Nf+l+Ni, 2), xTheta(1), yTheta(1), xTheta(2), yTheta(2), errorl, errorr, deltal, deltar) 
   !!!!!        sinPhi = (xTheta(1)-points_courbe_equal(l, 1))/(&
   !!!!!sqrt((xTheta(1)-points_courbe_equal(l, 1))**2+(yTheta(1)-points_courbe_equal(l, 2))**2))
   !!!!!        cosPhi = (yTheta(1)-points_courbe_equal(l, 2))/(&
   !!!!!sqrt((xTheta(1)-points_courbe_equal(l, 1))**2+(yTheta(1)-points_courbe_equal(l, 2))**2)) 
   !!!!!           slice(theta, l+Ni-1, 1) = points_courbe_equal(l, 1) + valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*sinPhi
   !!!!!           slice(theta, l+Ni-1, 2) = points_courbe_equal(l, 2) + valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*cosPhi
   !!!!!           slice(theta, l+Ni-1, 3) = zslice + valDist(l+Ni-1, theta)*sin(valTheta(l+Ni-1, theta))
      !!           if ((det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2), vect(indextheta(l+Ni-1), l+Ni-1, 1)&
      !!, vect(indextheta(l+Ni-1), l+Ni-1, 2), points_courbe_equal(l, 1), points_courbe_equal(l, 2))*det(points_courbe_equal(l-1, 1)&
      !!, points_courbe_equal(l-1, 2), vect(indextheta(l+Ni-1), l+Ni-1, 1), vect(indextheta(l+Ni-1), l+Ni-1, 2), slice(theta, l+Ni-1, 1)&
      !!, slice(theta, l+Ni-1, 2))>= 0._pr).or.(l<Ni)) then
      !!              indextheta(l+Ni-1) = indextheta(l+Ni-1) + 1
      !!              vect(indextheta(l+Ni-1), l+Ni-1, 1) = slice(theta, l+Ni-1, 1)
      !!              vect(indextheta(l+Ni-1), l+Ni-1, 2) = slice(theta, l+Ni-1, 2)
      !!              vect(indextheta(l+Ni-1), l+Ni-1, 3) = slice(theta, l+Ni-1, 3)
      !!           endif
      !!     enddo
      !!     l = Ns
   !!!!!     call intersection(points_courbe_equal(l-1, 1), head_courbe(1+1, 1), points_courbe_equal(Ns-1, 2), head_courbe(1+1, 2)&
   !!!!!, head_courbe(1, 1), head_courbe(1, 2), distslice(Ni+Ns, 1), distslice(1+Ni+Ns, 2), distslice(Ns+Ni+Nf+Ni+Ns, 1)&
   !!!!!, distslice(Ns+Ni+Nf+1+Ni+Ns, 2), xTheta(1), yTheta(1), xTheta(2), yTheta(2), errorl, errorr, deltal, deltar) 
   !!!!!     sinPhi = (xTheta(1)-head_courbe(1, 1))/(sqrt((xTheta(1)-head_courbe(1, 1))**2+(yTheta(1)-head_courbe(1, 2))**2))
   !!!!!     cosPhi = (yTheta(1)-head_courbe(1, 2))/(sqrt((xTheta(1)-head_courbe(1, 1))**2+(yTheta(1)-head_courbe(1, 2))**2)) 
   !!!!!        slice(theta, l+Ni-1, 1) = head_courbe(1, 1) + valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*sinPhi
   !!!!!        slice(theta, l+Ni-1, 2) = head_courbe(1, 2) + valDist(l+Ni-1, theta)*cos(valTheta(l+Ni-1, theta))*cosPhi
   !!!!!        slice(theta, l+Ni-1, 3) = zslice + valDist(l+Ni-1, theta)*sin(valTheta(l+Ni-1, theta))
      !!        if (det(points_courbe_equal(l-1, 1), points_courbe_equal(l-1, 2), vect(indextheta(l+Ni-1), l+Ni-1, 1)&
      !!, vect(indextheta(l+Ni-1), l+Ni-1, 2), points_courbe_equal(l, 1), points_courbe_equal(l, 2))*det(points_courbe_equal(l-1, 1)&
      !!, points_courbe_equal(l-1, 2), vect(indextheta(l+Ni-1), l+Ni-1, 1), vect(indextheta(l+Ni-1), l+Ni-1, 2), slice(theta, l+Ni-1, 1)&
      !!, slice(theta, l+Ni-1, 2))>= 0._pr) then
      !!           indextheta(l+Ni-1) = indextheta(l+Ni-1) + 1
      !!           vect(indextheta(l+Ni-1), l+Ni-1, 1) = slice(theta, l+Ni-1, 1)
      !!           vect(indextheta(l+Ni-1), l+Ni-1, 2) = slice(theta, l+Ni-1, 2)
      !!           vect(indextheta(l+Ni-1), l+Ni-1, 3) = slice(theta, l+Ni-1, 3)
      !!        endif
      !!     do l=2, Nf-1
   !!!!!        call intersection(head_courbe(l-1, 1), head_courbe(l+1, 1), head_courbe(l-1, 2), head_courbe(l+1, 2), head_courbe(l, 1)&
   !!!!!, head_courbe(l, 2), distslice(l+Ni+Ns, 1), distslice(l+Ni+Ns, 2), distslice(Ns+Ni+Nf+l+Ni+Ns, 1), distslice(Ns+Ni+Nf+l+Ni+Ns, 2), xTheta(1)&
   !!!!!, yTheta(1), xTheta(2), yTheta(2), errorl, errorr, deltal, deltar) 
   !!!!!        sinPhi = (xTheta(1)-head_courbe(l, 1))/(sqrt((xTheta(1)-head_courbe(l, 1))**2+(yTheta(1)-head_courbe(l, 2))**2))
   !!!!!        cosPhi = (yTheta(1)-head_courbe(l, 2))/(sqrt((xTheta(1)-head_courbe(l, 1))**2+(yTheta(1)-head_courbe(l, 2))**2)) 
   !!!!!           slice(theta, l+Ni-1+Ns-1, 1) = head_courbe(l, 1) + valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*sinPhi
   !!!!!           slice(theta, l+Ni-1+Ns-1, 2) = head_courbe(l, 2) + valDist(l+Ni-1+Ns-1, theta)*cos(valTheta(l+Ni-1+Ns-1, theta))*cosPhi
   !!!!!           slice(theta, l+Ni-1+Ns-1, 3) = zslice + valDist(l+Ni-1+Ns-1, theta)*sin(valTheta(l+Ni-1+Ns-1, theta))
      !!           if (det(head_courbe(l-1, 1), head_courbe(l-1, 2), vect(indextheta(l+Ni-1+Ns-1), l+Ni-1+Ns-1, 1), vect(indextheta(l+Ni-1+Ns-1)&
      !!, l+Ni-1+Ns-1, 2), head_courbe(l, 1), head_courbe(l, 2))*det(head_courbe(l-1, 1), head_courbe(l-1, 2), vect(indextheta(l+Ni-1+Ns-1)&
      !!, l+Ni-1+Ns-1, 1), vect(indextheta(l+NI-1+Ns-1), l+Ni-1+Ns-1, 2), slice(theta, l+Ni-1+Ns-1, 1), slice(theta, l+Ni-1+Ns-1, 2))>= 0._pr) then
      !!              indextheta(l+Ni-1+Ns-1) = indextheta(l+Ni-1+Ns-1) + 1
      !!              vect(indextheta(l+Ni-1+Ns-1), l+Ni-1+Ns-1, 1) = slice(theta, l+Ni-1+Ns-1, 1)
      !!              vect(indextheta(l+Ni-1+Ns-1), l+Ni-1+Ns-1, 2) = slice(theta, l+Ni-1+Ns-1, 2)
      !!              vect(indextheta(l+Ni-1+Ns-1), l+Ni-1+Ns-1, 3) = slice(theta, l+Ni-1+Ns-1, 3)
      !!           endif
      !!     enddo
      !!     l = Nf
      !!        indextheta(l+Ni-1+Ns-1) = indextheta(l+Ni-1+Ns-1) + 1
      !!        vect(indextheta(l+Ni-1+Ns-1), l+Ni-1+Ns-1, 1) = head_courbe(l, 1)
      !!        vect(indextheta(l+Ni-1+Ns-1), l+Ni-1+Ns-1, 2) = head_courbe(l, 2)
      !!        vect(indextheta(l+Ni-1+Ns-1), l+Ni-1+Ns-1, 3) = head_courbe(l, 3)
      !!   enddo
      do l=1, Ni+Ns+Nf-2
         indextheta(l) = 0
         do theta = 1, nbtheta
   !!!           if (l==indextab(theta, l)) then
            !           if (theta==101) write(*, *) "CHECK ERROR101 ", theta, " ", l, " ", l+Ni-1+Ns-1, " ", slice(theta, l+Ni-1+Ns-1, 1)&
            !, " ", slice(theta, l+Ni-1+Ns-1, 2), " ", slice(theta, l+Ni-1+Ns-1, 3), " ", cos(valTheta(l+Ni-1+Ns-1, theta))&
            !, " ", sin(valTheta(l+Ni-1+Ns-1, theta)), " ", indextheta(l)
            !!           if (appartientElmt(indextab, theta, l).eqv..true.) then
            if (indextab(theta, l)==1) then
               indextheta(l) = indextheta(l) + 1
   !!!              vect(l, indextheta(l)) = theta
   !!!              vect(indextheta(l), l, 1) = slice2(theta, thetatab(theta), 1)
   !!!              vect(indextheta(l), l, 2) = slice2(theta, thetatab(theta), 2)
   !!!              vect(indextheta(l), l, 3) = slice2(theta, thetatab(theta), 3)
               !           if (theta==101) write(*, *) "CHECK ERROR101 appartient  ", theta, " ", l, " ", l+Ni-1+Ns-1, " ", slice(theta, l+Ni-1+Ns-1, 1)&
               !, " ", slice(theta, l+Ni-1+Ns-1, 2), " ", slice(theta, l+Ni-1+Ns-1, 3), " ", cos(valTheta(l+Ni-1+Ns-1, theta))&
               !, " ", sin(valTheta(l+Ni-1+Ns-1, theta)), " ", indextheta(l)
               vect(indextheta(l), l, 1) = slice(theta, l, 1)
               vect(indextheta(l), l, 2) = slice(theta, l, 2)
               !              if (vect(101, l, 2)>110.0) write(*, *) "WHHHYY  ", theta, " ", l, " ", indextheta(l)
               vect(indextheta(l), l, 3) = slice(theta, l, 3)
            endif
         enddo
      enddo

      !     write(*, *) "lpl lpr finaux : ", thetatab(1), " ", thetatab(2)

      deallocate(midline)
      deallocate(points_control)

      !**
      !! OUTPUTS 
      !**

      !! FILM

      open(unit=85, file=trim(target_folder)//'/skelhead'//str(kt+picNum-1)//'.txt', status='unknown')
      open(unit=84, file=trim(target_folder)//'/right'//str(kt+picNum-1)//'.txt', status='unknown')
      open(unit=83, file=trim(target_folder)//'/left'//str(kt+picNum-1)//'.txt', status='unknown')
      open(unit=82, file=trim(target_folder)//'/skeletteq'//str(kt+picNum-1)//'.txt', status='unknown')
      open(unit=81, file=trim(target_folder)//'/skelett'//str(kt+picNum-1)//'.txt', status='unknown')
      open(unit=79, file=trim(target_folder)//'/skelett'//str(kt+picnum-1)//'.vtk', status='unknown')
      file_id = 79
      n_dim = 2
      call write_vtk_header(file_id, nx, ny, nz, dx, dy, dz, n_dim)

      do k=1, ny
         do j=1, nx
            write(79, *) rhoSlices2(j, k)*(1-skel(j, k)) + minval(rhoSlices2)*skel(j, k)

         enddo
      enddo
      write(*, *) "CHECK size nbtheta : ", sum(thetatab), " ", sum(indextheta), "  THETA  ", thetatab(nbtheta/2+1)
      
      !! FILM

      open(unit=78, file=trim(target_folder)//'/surf'//str(kt+picNum-1)//'.vts', status='unknown')
      write(78, '(a)') "<?xml version=""1.0""?>"
      write(78, '(a)') "<VTKFile type=""StructuredGrid"" version=""0.1"" byte_order=""LittleEndian"" &
      &     compressor=""vtkZLibDataCompressor"">"
      write(78, '(a, I3, a, I3, a)') "<StructuredGrid WholeExtent=""0 ", nbtheta, " 0 ", Ni-1+Ns+Nf-1-1, " 0 0"">"
      write(78, '(a, I3, a, I3, a)') "<Piece Extent=""0 ", nbtheta, " 0 ", Ni-1+Ns+Nf-1-1, " 0 0"">"
      write(78, '(a)') "<PointData >"
      write(78, '(a)') "</PointData>"
      write(78, '(a)') "<CellData>"
      write(78, '(a)') "</CellData>"
      write(78, '(a)') "<Points>"
      write(78, '(a)') "<DataArray NumberOfComponents=""3"" type=""Float64"" format=""ascii"" >"  
      write(78, '(a, I4, I4, I4)') 'DIMENSIONS', nx, ny, nz
      write(78, '(a, E23.15, E23.15, E23.15)') 'ORIGIN', 1., 1., 1.
      write(78, '(a, E23.15, E23.15, E23.15)') 'SPACING', dx, dy, dz
      write(78, '(a, I9)') 'POINT_DATA' , nx*ny*nz
      write(78, '(a)') 'SCALARS values double'
      write(78, '(a)') 'LOOKUP_TABLE default'
!!        do theta=1, nbtheta
!!     do l=1, thetatab(theta)
!!           write(78, *) slice2(theta, l, 1), " ", slice2(theta, l, 2), " ", slice2(theta, l, 3)
!!     enddo
!!        enddo
!!        do theta=1, nbtheta
!!        do l=1, thetatab(theta)
!!        write(*, *) "OHOH th ", 
      open(unit=92, file=trim(target_folder)//'/surf/surf'//str(kt+picNum-1)//'.dat', status='unknown')

      do l=1, Ni+Ns+Nf-2
         do theta=1, nbtheta 
            if (indextab(theta, l)==1) write(78, *) slice(theta, l, 1), " ", slice(theta, l, 2), " ", slice(theta, l, 3)
            if (indextab(theta, l)==1) write(92, *) slice(theta, l, 1), " ", slice(theta, l, 2), " ", slice(theta, l, 3)
         enddo
            write(92, *) slice(1, l, 1), " ", slice(1, l, 2), " ", slice(1, l, 3)
      enddo

      write(78, '(a)') "</DataArray>"
      write(78, '(a)') "</Points>"
      write(78, '(a)') "</Piece>"
      write(78, '(a)') "</StructuredGrid>"
      write(78, '(a)') "</VTKFile>"
      close(78)

            do l=1, Ni-1
               write(82, *) tail_courbe(Ni-l+1, 1), " ", tail_courbe(Ni-l+1, 2), " ", tail_courbe(Ni-l+1, 3)
            enddo

            do l=2, Nf
               write(85, *) head_courbe(l, 1), " ", head_courbe(l, 2), " ", head_courbe(l, 3)
            enddo

            write(80, *) kt, "    ", nl, " ", long3, " ", long2, " ", long
            do l=1, Ns
               write(81, *) points_courbe(l, 1), " ", points_courbe(l, 2), " ", points_courbe(l, 3)
               write(82, *) points_courbe_equal(l, 1), " ", points_courbe_equal(l, 2), " ", points_courbe_equal(l, 3)
            enddo

            do l=1, Ni+Ns+Nf-2 
               if (indextab(1+nbtheta/2, l)==1) write(84, *) slice(1+nbtheta/2, l, 1), " ", slice(1+nbtheta/2, l, 2)&
                     , " ", slice(1+nbtheta/2, l, 3)
               if (indextab(1, l)==1) write(83, *) slice(1, l, 1), " ", slice(1, l, 2), " ", slice(1, l, 3)
            enddo
      close(79)     
      close(92)    
      close(81)
      close(82)
      close(83)
      close(84)
      close(85)
   !     write(*, *) "okcut"
   !     call cutheadskel(midline, skel, skel2, nl, ic, jc)
   !     write(*, *) "okcut"
   !
   !     deallocate(midline)
   enddo
   close(80)
   write(*, *) "************************************************" 
   deallocate(rhoSlices2)
   !deallocate(rhoSlices)
   !  deallocate(rho0) 
   deallocate(xx, yy)
   !  deallocate(u, v, rho, rhou, rhov) 
   !  deallocate(rhoup, rhovp, rhop)
   !  deallocate(uu, vv, rhouu, rhovv, rhouup, rhovvp)
   !  deallocate(u0, v0)
   deallocate(tmp1, tmp2, tmpbool)
   deallocate(dir1, dir2, dir3, dir4, Nseed, skel, skel2, skel3)
   !deallocate(un, zero)
   deallocate(distslice)
   deallocate(points_courbe, points_courbe_equal, points_courbe_equal_ref)
   deallocate(slice_courbe, slice_courbe_equal, slice_courbemidLarge)
   deallocate(tail_courbe, head_courbe)
   deallocate(slice, slice2, thetatab, xTheta, yTheta, valTheta, valDist, slicetmp, longTh, longTheta, sth, dsth)
   deallocate(longslice, longslicei, longslicef, sslice, dsslice, slicemid)
   deallocate(indextab, indextheta, vect)
   deallocate(gradPhi) !, valTh)
   deallocate(cosTheta_tab, sinTheta_tab)
end program def3D
