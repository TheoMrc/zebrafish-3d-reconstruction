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
   integer :: i, j, k, kt, idisplay, iUpdateDist, iter, picNum, l, ll, pix1, pix2, l0, pos, sizeSkel, nl, lp, bool, booltmp, &
   &    boolPhi, loopbool, xskelL, yskelL, xskelR, yskelR, ii, lph, lpt, lpf, hh, ic, jc, ic0, jc0, niter, file_id, n_dim
   real(pr) :: Tps, threshold, surf, maxv, minimL, minimL_m, minimL_p, minimR, long, tb, dtb, long2, long3, dti, dtf, dsi, dsf, &
   &    xr, yr, xl, yl, long00, tp, longexp, longratio, xg, yg, xgref, ygref
   real(pr), dimension(:), allocatable :: longslice, longslicei, longslicef, longTh, longTheta
   real(pr) :: t, tPic, xi, yj, distW, distWb, xLeft, yLeft, xRight, yRight
   real(pr), dimension(:,:), allocatable :: fish_shape, rhoSlices, gradPhi
   real(pr), dimension(:,:,:), allocatable :: rhoSlices3D,rhoSlices3D2
   integer, dimension(:,:), allocatable :: midline,midlinebis
   real(pr), dimension(:,:), allocatable :: tmp, tmp1, tmp2, un, zero,distslice,voisin
   integer, dimension(:,:), allocatable :: dir1, dir2, dir3, dir4, Nseed, skel, skel2, tmpbool, skel3, skeltmp
   real(pr) :: px,py,pz,s0,sinit,tbm,tbp,s,ds,rr,pxx,pyy,pzz
   real(pr),dimension(:),allocatable :: sslice,dsslice,stheta,dstheta,sth,dsth
   real(pr),dimension(:,:),allocatable :: points_control, points_courbe, points_courbe_equal, tail_courbe, head_courbe&
   & , points_courbe_equal_ref
   real(pr),dimension(:,:,:),allocatable :: slicecontrolm, slicecontrolf, slicecontrol,slicemid,slicecontroli,slicecontroltmp
   real(pr),dimension(:,:,:),allocatable :: slice_courbe, slice_courbe_equal
   integer :: Ns,Ni,Nf,nt,errorl,errorr,itail,ihead
   real(pr) :: deltal,deltar,disttail,disthead,rhead,rtail,disttaillY,disttailrY,distheadlY,distheadrY
   real(pr),dimension(:,:,:),allocatable :: slice,slice2,vect,slicetmp
   integer,dimension(:),allocatable :: thetatab, indextheta
   integer, dimension(:,:), allocatable :: indextab
   real(pr),dimension(:),allocatable :: xTheta,yTheta
   real(pr),dimension(:,:),allocatable :: valDist,valTheta,valTh,valThtmp
   integer :: th, theta, boolskel
   real(pr) :: PI, sigma, sigma2 
   real(pr) :: cosPhi,sinPhi,cosTheta,sinTheta,oldS,oldC,cosAlpha,sinAlpha,cosPhil,cosPhir,sinPhil,sinPhir
   real(pr), dimension(:,:), allocatable :: cosTheta_tab,sinTheta_tab
   real(pr) :: zslice,area,meshRatio,seuil,seuil0,inc
   real(pr) :: x1,x2,y1,y2,z1,z2,xc,yc,zc,delta,rt,xt,yt,xp,yp,x0
   real(pr) :: LS1,LS2,LS3,LS4,LSp,LSrr,rp,alpha,alphadef,calpha,salpha, start_time, end_time
   integer :: ip1,ip2,jp1,jp2
   integer :: index, sum_grad_thresh, reduced_grad_thresh, n_frames
   integer, dimension(3) :: file_ids
   character(len=256) :: target_folder, subfolder
 
   !**
   !! Main Initialisation
   !**

   call cpu_time(start_time)

   target_folder = "D:/data_modelisation/results/training_video"
   subfolder = "full_10th_corr"      ! CHANGE DESTINATION SUBFOLDER HERE

   target_folder = trim(trim(adjustl(target_folder))//"/"//trim(adjustl(subfolder)))

   call create_directory(target_folder//"/surf")

   n_frames = 580
   meshRatio = 0.2_pr
   area = 0.2_pr
   PI = acos(-1.0_pr)
   sigma = 1._pr
   sigma2 = 1._pr
   nx = 1602 
   ny = 300 
   nz = 300
   x0 = 1._pr
   eepsilon = 1.e-6_pr
   dt = 1._pr
   threshold = 0.001

   n_dim = 2

   dx = 2.4 
   dy = 2.25171 
   dz = 2.25171 
   Tps = 1._pr

   Ns = 263 
   Ni = 4 
   Nf = 35 
   dtb = 1._pr/(Ns-1)
 

   !**
   !! Main Allocation of arrays
   !**
   allocate(cosTheta_tab(Ns+Ni+Nf-2,2),sinTheta_tab(Ns+Ni+Nf-2,2))
   allocate(slicemid(2,Ns,2))
   allocate(longslice(2),longslicei(2),longslicef(2),sslice(2),dsslice(2),stheta(Ns+Ni+Nf-2)&
        ,dstheta(Ns+Ni+Nf-2))
   allocate(thetatab(2),indextheta(Ns+Ni+Nf-2),indextab(2,Ns+Ni+Nf-2))
   allocate(xTheta(2))
   allocate(yTheta(2))
   allocate(valTheta(Ns+Ni+Nf-2,2),longTheta(Ns+Ni+Nf-2))
   allocate(valDist(Ns+Ni+Nf-2,2))
   allocate(head_courbe(Nf,2),tail_courbe(Ni,2))
   allocate(points_courbe(Ns,2))
   allocate(points_courbe_equal(Ns,2))
   allocate(points_courbe_equal_ref(Ns,2))
   allocate(rhoSlices(nx,ny), gradPhi(nx,ny), fish_shape(nx, ny))
   allocate(rhoSlices3D(nz,nx,ny))
   allocate(rhoSlices3D2(nz,ny,nx))

   allocate(tmp(nx,ny),tmp1(nx,ny),tmp2(nx,ny),tmpbool(nx,ny))
   allocate(xx(nx))
   allocate(yy(ny)) 
   allocate(dir1(nx,ny),dir2(nx,ny),dir3(nx,ny),dir4(nx,ny),Nseed(nx,ny),skel(nx,ny),skel2(nx,ny),skel3(nx,ny),skeltmp(nx,ny))
   allocate(distslice(2*(Ns+Ni+Nf),2))
   allocate(slice_courbe(2,Ns+Ni+Nf-2,2))
   allocate(slice_courbe_equal(2,Ns+Ni+Nf-2,2))
   allocate(slice(2,Ns+Ni+Nf-2,2),slicetmp(2,Ns+Ni+Nf-2,2))
   allocate(slice2(2,Ns+Ni+Nf-2,2),vect(2,Ns+Ni+Nf-2,2))
 
 
   x0 = 0._pr
   do i=1,nx
      xx(i) = x0 + (float(i)-1)*dx
   enddo
   do j=1,ny
      yy(j) = x0 + (float(j)-1)*dy
   enddo
 

   !! FILM
   picNum = 1

   !**
   !! Flatten the 3D shape into 2D shape (time consuming)
   !*

   ! shape_file_path = 'fish_shapes/rhoSlices_v2.dat'
   ! call flatten3Dshape(shape_file_path, rhoSlices, rhoSlices3D, gradPhi)
   ! deallocate(rhoSlices3D)
   
   !**
   !! Loading the 2D silhouette
   !**

   open(unit=79,file='fish_shapes/2Dshape.dat',status='unknown')
   do j=1,ny
      do k=1,nx
        read(79,*) fish_shape(k,j)
      enddo
   enddo
   close(79)
   
   call cpu_time(end_time)
   write(*,*) "******Initialisation and shape loading time = ", end_time-start_time, " sec******"  ! negligible (0,35 sec)

   
   !**
   !! Now, the level-set is computed in fish_shape
   !**
   call cpu_time(start_time)
   call updateDistanceINI(fish_shape, gradPhi)

   dx = 2.4*0.000001_pr 
   dy = 2.25171*0.000001_pr 
   dz = 2.25171*0.000001_pr

   do i=1,nx
      xx(i) = x0 + (float(i)-1)*dx
   enddo
   do j=1,ny
      yy(j) = x0 + (float(j)-1)*dy
   enddo

   open(unit=79,file= trim(target_folder)//'/sol2D00.txt',status='unknown')
   open(unit=80,file= trim(target_folder)//'/sol2D00.vtk',status='unknown')
   file_id = 80
   call write_vtk_header(file_id, nx, ny, nz, dx, dy, dz, n_dim)
   do k=1,nx
      do j=1,ny
         if (fish_shape(k,j)>0._pr) write(79,*) xx(k)," ",yy(j)
         write(80,*) fish_shape(k,j)
     enddo
   enddo
   close(79)
   close(80)

   call cpu_time(end_time)
   write(*,*) "Shape processing in Rho_slices time = ", end_time-start_time, " sec" ! Around 45 sec of exec

!    !**
!    !! Definition of the zebrafish length
!    !**

   call cpu_time(start_time)

!    !! FILM
   long00 = 3.8424001526832577E-003 - 0._pr  
   longexp = 164*0.001_pr*0.0256
   longratio = longexp/long00
   write(*,*) "longueurs initiales ",long00," ",longexp," ",longratio
  
   !**
   !! Construction of the initial midline
   !**
 
   dir1 = 0
   dir2 = 0
   dir3 = 0
   dir4 = 0
   Nseed = 0
   skel = 0
   do i=2,ny-1
      do j=2,nx-1
         if ((gradPhi(j,i)<0.74).and.(fish_shape(j,i)>0._pr)) skel(j,i) = 1
         enddo
      enddo
 
   tmpbool = 0
   nl = sum(skel)
   sizeSkel = sum(skel)
   allocate(midlinebis(sizeSkel,2))
   l=1
   midlinebis = 0
   do j=1,nx
      do k=1,ny
         if ((skel(j,k)==1).and.(l==1)) then
            midlinebis(l,1) = j
            midlinebis(l,2) = k
            if (j==nx/2) l0 = l
            l = l+1
         endif
      enddo
   enddo
   !**
   !! first point of the midline (tail)
   !**
   xskelL = midlinebis(1,1)
   yskelL = midlinebis(1,2)
   lp = 1
   boolskel=0
   do i=2,nx-1

     !**
     !! horizontal midline
     !**
      if (lp<=1462) then !nombre de points sur la première midline
         lp=lp+1
         midlinebis(lp,1) = i
         midlinebis(lp,2) = midlinebis(1,2)
      endif
   enddo
   nl = lp
   
   !**
   !! Now we have the midline coordinates
   !**
   allocate(midline(nl,2))
   midline(1:nl,1) = midlinebis(1:nl,1)
   midline(1:nl,2) = midlinebis(1:nl,2)
  
   !**
   !! Definition of the control points of the midline
   !**
   allocate(points_control(nl,2))!,leftcontrolm(nl,2),rightcontrolm(nl,2))
   do l=1,nl
      points_control(l,1) = xx(midline(l,1))
      points_control(l,2) = yy(midline(1,2))!yy(midline(l,2))
      if (l==1) write(*,*) "POINT_CONTROL L==1  ",points_control(l,1)," ",points_control(l,2)&
      ,"-- ",midline(1,1)," ",midline(1,2)
      if (l==nl) write(*,*) "POINT_CONTROL L==NL  ",points_control(l,1)," ",points_control(l,2)&
      ,"  ",nl,"-- ",midline(nl,1)," ",midline(1,2)
   enddo
   !points_control(1,1) = points_control(1,1) - 1
   !! TESTNEZ
   !!  points_control(nl,1) = xc
   !!  points_control(nl,2) = yc
   !!  points_control(nl,3) = zslice
 
   !**
   !! Spline approximation of the midline
   !**
   points_courbe(1,1) = points_control(1,1)
   points_courbe(1,2) = points_control(1,2)
 
   tb = dtb
   l = 1
   long = 0._pr
   do while ((tb<1._pr).and.(l+1<Ns+1))
      l = l+1
      call pointsBezierN(points_control,tb,px,py)
      points_courbe(l,1) = px
      points_courbe(l,2) = py
      long = long + sqrt((points_courbe(l,1)-points_courbe(l-1,1))**2 + (points_courbe(l,2)-points_courbe(l-1,2))**2)
      tb = tb+dtb
   enddo
   if (l==Ns-1) then
      l = l+1
      write(*,*) "ET VOILA"
      points_courbe(Ns,1) = points_control(size(points_control,1),1)
      points_courbe(Ns,2) = points_control(size(points_control,1),2)
      long = long + sqrt((points_courbe(l,1)-points_courbe(l-1,1))**2 + (points_courbe(l,2)-points_courbe(l-1,2))**2)
   endif
   if (.not.(l==Ns)) write(*,*) "WARNIINNGG ",l
   write(*,*) "longueur ==  ",(points_control(size(points_control,1),1)-points_control(1,1))," ",long," ",(points_courbe(Ns,1)-&
        points_courbe(1,1))," ",points_courbe(1,1)," ",points_courbe(Ns,1),"       ",nl
   long2 = long
   ds = long/(Ns-1)
 
   !**
   !! Uniform spline approximation of the midline
   !**
   points_courbe_equal(1,1) = points_control(1,1)
   points_courbe_equal(1,2) = points_control(1,2)
   l = 1
   long = 0._pr
   dtb = 1._pr/(Ns-1)
   tb = 0._pr
   do while ((tb<1._pr).and.(l+1<Ns))
      l = l+1
      nt = 1
      s = 0._pr
      do while ((l-1)*ds-s>0._pr) 
         nt = nt+1
         s = s + sqrt((points_courbe(nt,1)-points_courbe(nt-1,1))**2 + (points_courbe(nt,2)-points_courbe(nt-1,2))**2)
      enddo
      tbm = (nt-2)*dtb
      tbp = (nt-1)*dtb
      tb = tbm 
      s0 = s
      sinit = s - sqrt((points_courbe(nt,1)-points_courbe(nt-1,1))**2 + (points_courbe(nt,2)-points_courbe(nt-1,2))**2)
      s = sinit
      bool = 0
      do while ((abs((l-1)*ds-s)/dx>eepsilon).and.(bool==0))
         tb = (tbm + tbp)*0.5_pr
         call pointsBezierN(points_control,tb,px,py)
         s = sinit + sqrt((px-points_courbe(nt-1,1))**2 + (py-points_courbe(nt-1,2))**2)
         if ((l-1)*ds-s>0._pr) then
            tbm = tb
         else
            tbp = tb
            if (tbp>1._pr) tbp = 1._pr
         endif
         if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then 
            bool = 1
         endif
      enddo
      write(*,*) "TBpoints  ",l," ",tb," ",ds," ",tbm," ",tbp
      call pointsBezierN(points_control,tb,px,py)
      points_courbe_equal(l,1) = px
      points_courbe_equal(l,2) = py
      long = long +&
           sqrt((points_courbe_equal(l,1)-points_courbe_equal(l-1,1))**2 + (points_courbe_equal(l,2)-points_courbe_equal(l-1,2))**2)
      write(*,*) "LLfirst  ",l,"       ",px," ",py," ",ds
   enddo
   lp = l
 
   if (l==Ns-1) then
      l = l+1
      write(*,*) "ET VOILA1"
      points_courbe_equal(Ns,1) = points_control(size(points_control,1),1)
      points_courbe_equal(Ns,2) = points_control(size(points_control,1),2)
      long = long + sqrt((points_courbe(l,1)-points_courbe(l-1,1))**2 + (points_courbe(l,2)-points_courbe(l-1,2))**2)
   endif
   if (.not.(l==Ns)) write(*,*) "WARNIINNGG2 ",l
   !if (.not.(l+lp-1==Ns)) write(*,*) "WARNIINNGG2 ",l
   write(*,*) "longueur =  ",(points_control(size(points_control,1),1)-points_control(1,1))," ",long," ",& 
   & (points_courbe_equal(Ns,1)-points_courbe_equal(1,1)),"  LL  ",l
   long3 = 0._pr
   do l=1,nl-1
      long3 = long3 +sqrt((xx(midline(l+1,1))-xx(midline(l,1)))**2 + (yy(midline(l+1,2))-yy(midline(l,2)))**2)
   enddo
   long3 = long3+1
   write(*,*) "longueur =  ",long3
 
   !**
   !! Interpolation of endpoints of the real midline 
   !**
   rtail = xx(1)
   itail = floor((rtail-x0)/dx+1)
   rhead = xx(nx)
   ihead = ceiling((rhead-x0)/dx+1)
   write(*,*) fish_shape(1,midline(1,2))
   write(*,*) fish_shape(2,midline(1,2))
   do i=2,nx-1
      if ((fish_shape(i,midline(1,2))>0._pr).and.(fish_shape(i-1,midline(1,2))<0._pr)) then 
         itail = i
         rtail = xx(i-1) - fish_shape(i,midline(1,2))*dx/(fish_shape(i-1,midline(1,2)) - fish_shape(i,midline(1,2)))
         itail = floor((rtail-x0)/dx+1)
      endif
      if ((fish_shape(i,midline(1,2))>0._pr).and.(fish_shape(i+1,midline(1,2))<0._pr)) then
         ihead = i
         rhead = xx(i) - fish_shape(i,midline(1,2))*dx/(fish_shape(i+1,midline(1,2)) - fish_shape(i,midline(1,2)))
         ihead = ceiling((rhead-x0)/dx+1)
      endif
   enddo
   disttail = abs(rtail - points_courbe_equal(1,1))
   disthead = abs(rhead-points_courbe_equal(size(points_courbe_equal,1),1))
   write(*,*) "longueur =  ",rhead-rtail," ",rhead-rtail-disttail," ",rhead-rtail-disthead," ",rhead-rtail-disttail-disthead
   long00 = rhead-rtail-disttail-disthead
   write(*,*) "longueur ",long3," ",long2," ",long," ",long00," ",longexp
   write(*,*) "DIST : TAIL ",disttail," HEAD ",disthead," NINT ",nint(disttail/dx)," ",nint(disthead/dx)," ",ceiling(disthead/dx)&
        ," ",floor(disthead/dx)," ",rtail," ",rhead&
        ," ",itail," ",ihead," ",ceiling(points_courbe_equal(1,1)/dx)," ",floor(points_courbe_equal(Ns,1)/dx)!," 2 ",2
 
   !**
   !! Allocation of control arrays 
   !**
   allocate(slicecontrolm(2,nint(points_courbe_equal(Ns,1)/dx)-nint(points_courbe_equal(1,1)/dx)+1,2))
   allocate(slicecontrolf(2,nint(disthead/dx)+1,2))
   allocate(slicecontroli(2,nint(disttail/dx)+1,2))
   allocate(slicecontrol(2,ihead-itail+1,2),slicecontroltmp(2,ihead-itail+1,2)) 
   allocate(valTh(ihead-itail+1,2))
   allocate(valThtmp(2,ihead-itail+1))
   write(*,*) "SIZE CONTROLLARGE  ",points_courbe_equal(Ns,1)+disthead-points_courbe_equal(1,1)&
        -disttail," ",rhead," ",rtail
   do l=1,size(valTh,1)
      do theta=1,2
         valTh(l,theta) = (theta-1)*2*PI/2
         valThtmp(theta,l) = (theta-1)*2*PI/2
      enddo
   enddo
 
   !**
   !! Filling control arrays (based on level-set zero searching - dichotomy) 
   !**
   lp = 0
   lph = 0
   lpf = 0
   lpt = 0
   do l=itail,ihead
      do theta=1,2
         rr = 0._pr
         bool = 0
         do while (((LSrr<0._pr).and.(rr<float(ny))).or.(bool==0))
            rp = rr
            if (bool==1) then
               LSp = LSrr
            else
               if (cos(valTh(l-itail+1,theta))<0._pr) then
                 xp = -rp*cos(valTh(l-itail+1,theta))+yy(1)
               else
                 xp = -rp*cos(valTh(l-itail+1,theta))+yy(ny)
               endif
               ip1 = int((xp-x0)/dy)+1
               if (cos(valTh(l-itail+1,theta))<0._pr) then
                 ip2 = ip1+1
               else
                 ip2 = ip1-1
               endif
               LS1 = fish_shape(l,ip1)
               LS2 = fish_shape(l,ip2)
               x1 = yy(ip1)
               if (cos(valTh(l-itail+1,theta))<0._pr) then
                 x2 = x1 + dy
               else
                 x2 = x1 - dy
               endif
 
               LSp = interpLS2D(LS1,x1,LS2,x2,xp) 
            endif
 
            bool = 1
            rr=rr+dy
            if (cos(valTh(l-itail+1,theta))<0._pr) then
              xp = -rr*cos(valTh(l-itail+1,theta))+yy(1)
            else
              xp = -rr*cos(valTh(l-itail+1,theta))+yy(ny)
            endif
            ip1 = int((xp-x0)/dy)+1

            if (cos(valTh(l-itail+1,theta))<0._pr) then
              ip2 = ip1+1
            else
              ip2 = ip1-1
            endif
 
            LS1 = fish_shape(l,ip1)
            LS2 = fish_shape(l,ip2)

            x1 = yy(ip1)

            if (cos(valTh(l-itail+1,theta))<0._pr) then
              x2 = x1 + dy
            else
              x2 = x1 - dy
            endif
 
            LSrr = interpLS2D(LS1,x1,LS2,x2,xp) 

         enddo
         rr = rp
         LSrr = LSp
         do while ((LSrr<0._pr).and.(rr<float(ny)))
            rp = rr
            LSp = LSrr
 
            rr=rr+0.1*dy
            if (cos(valTh(l-itail+1,theta))<0._pr) then
              xp = -rr*cos(valTh(l-itail+1,theta))+yy(1)
            else
              xp = -rr*cos(valTh(l-itail+1,theta))+yy(ny)
            endif
            ip1 = int((xp-x0)/dy)+1
            if (cos(valTh(l-itail+1,theta))<0._pr) then
              ip2 = ip1+1
            else
              ip2 = ip1-1
            endif
            LS1 = fish_shape(l,ip1)
            LS2 = fish_shape(l,ip2)
            x1 = yy(ip1)
            if (cos(valTh(l-itail+1,theta))<0._pr) then
              x2 = x1 + dy
            else
              x2 = x1 - dy
            endif
 
            LSrr = interpLS2D(LS1,x1,LS2,x2,xp) 
         enddo
         px = xx(l)

         if (cos(valTh(l-itail+1,theta))<0._pr) then
           py = -rp*cos(valTh(l-itail+1,theta))+yy(1) +&
              -cos(valTh(l-itail+1,theta))*LSP*(rr-rp)/abs(LSrr - LSp) !rr+1
         else
           py = -rp*cos(valTh(l-itail+1,theta))+yy(ny) +&
              -cos(valTh(l-itail+1,theta))*LSP*(rr-rp)/abs(LSrr - LSp) !rr+1
         endif
 
         if (l<=nint((points_courbe_equal(1,1)-x0)/dx+1)) then
            if (lpt==0) lpt = l
            write(*,*) "LPT ",lpt,l-lpt+1
            slicecontroli(theta,l-lpt+1,1) = px
            slicecontroli(theta,l-lpt+1,2) = py

         endif

         if ((l>=nint((points_courbe_equal(1,1)-x0)/dx+1)).and.(l<=nint((points_courbe_equal(Ns,1)-x0)/dx+1))) then
            if (lp==0) lpf = l
            if (lp==0) lp = l
            slicecontrolm(theta,l-lp+1,1) = px
            slicecontrolm(theta,l-lp+1,2) = py

         endif

         if (l>=nint((points_courbe_equal(Ns,1)-x0)/dx+1)) then
            if (lph==0) lph = l
            slicecontrolf(theta,l-lph+1,1) = px
            slicecontrolf(theta,l-lph+1,2) = py

         endif
         slicecontrol(theta,l-itail+1,1) = px
         slicecontrol(theta,l-itail+1,2) = py
         if (l==itail) write(*,*) "SLICECONTROL L==itail  ",px," ",py," ",l," ",xx(l)
         if (l==itail+1) write(*,*) "SLICECONTROL L==itail+1  ",px," ",py," ",l," ",xx(l)
         if (l==itail+2) write(*,*) "SLICECONTROL L==itail+2  ",px," ",py," ",l," ",xx(l)
      enddo
 
      if (lp==0) then
      else
         lpf = lpf+1
      endif
   enddo
 
   !**
   !! Filling control arrays (based on level-set zero searching - dichotomy) for TAIL/HEAD slices
   !**
   do theta=1,2
      l=itail+1
      rr = 0._pr
      bool = 0
      do while (((LSrr<0._pr).and.(rr<float(ny))).or.(bool==0))
         rp = rr
         if (bool==1) then
            LSp = LSrr
         else

            if (cos(valTh(l-itail+1,theta))<0._pr) then
              xp = -rp*cos(valTh(l-itail+1,theta))+yy(1)
            else
              xp = -rp*cos(valTh(l-itail+1,theta))+yy(ny)
            endif
            ip1 = int((xp-x0)/dy)+1
            if (cos(valTh(l-itail+1,theta))<0._pr) then
              ip2 = ip1+1
            else
              ip2 = ip1-1
            endif
            LS1 = fish_shape(l,ip1)
            LS2 = fish_shape(l,ip2)
            !x1 = xx(ip1)
            !x2 = x1 + dx
            x1 = yy(ip1)
            !x2 = x1 + dy
            if (cos(valTh(l-itail+1,theta))<0._pr) then
              x2 = x1 + dy
            else
              x2 = x1 - dy
            endif
            LSp = interpLS2D(LS1,x1,LS2,x2,xp) 
            if (theta==1) write(*,*) "checknantail  ",x1," ",LSrr," ",LSp
         endif
 
         bool = 1
         rr=rr+dy

         if (cos(valTh(l-itail+1,theta))<0._pr) then
           xp = -rr*cos(valTh(l-itail+1,theta))+yy(1)
         else
           xp = -rr*cos(valTh(l-itail+1,theta))+yy(ny)
         endif
         ip1 = int((xp-x0)/dy)+1
         if (cos(valTh(l-itail+1,theta))<0._pr) then
           ip2 = ip1+1
         else
           ip2 = ip1-1
         endif
         LS1 = fish_shape(l,ip1)
         LS2 = fish_shape(l,ip2)
         x1 = yy(ip1)
         if (cos(valTh(l-itail+1,theta))<0._pr) then
           x2 = x1 + dy
         else
           x2 = x1 - dy
         endif
         LSrr = interpLS2D(LS1,x1,LS2,x2,xp) 
      enddo
      rr = rp
      LSrr = LSp
      do while ((LSrr<0._pr).and.(rr<float(ny)))
         rp = rr
         LSp = LSrr
         rr=rr+0.1*dy
         if (cos(valTh(l-itail+1,theta))<0._pr) then
           xp = -rr*cos(valTh(l-itail+1,theta))+yy(1)
         else
           xp = -rr*cos(valTh(l-itail+1,theta))+yy(ny)
         endif
         ip1 = int((xp-x0)/dy)+1
         if (cos(valTh(l-itail+1,theta))<0._pr) then
           ip2 = ip1+1
         else
           ip2 = ip1-1
         endif
         LS1 = fish_shape(l,ip1)
         LS2 = fish_shape(l,ip2)
         x1 = yy(ip1)
         if (cos(valTh(l-itail+1,theta))<0._pr) then
           x2 = x1 + dy
         else
           x2 = x1 - dy
         endif
         LSrr = interpLS2D(LS1,x1,LS2,x2,xp) 
      enddo
      px = rtail
      py = 0.01_pr*cos(valTh(l-itail+1,theta))*abs(slicecontrol(theta,2,2)-points_courbe_equal(1,2)) +&
      points_courbe_equal(1,2)
      slicecontroli(theta,1,1) = px 
      slicecontroli(theta,1,2) = py 
      slicecontrol(theta,1,1) = px 
      slicecontrol(theta,1,2) = py 
      if (theta==1) write(*,*) "checknantail  ",slicecontrol(theta,1,2)," ",LSrr," ",LSp
      l=ihead-1!+1
      rr = 0._pr
      bool = 0
      do while (((LSrr>0._pr).and.(rr<float(ny))).or.(bool==0))
         rp = rr
         if (bool==1) then
            LSp = LSrr
         else
            if (cos(valTh(l-itail+1,theta))<0._pr) then
              xp = -rp*cos(valTh(l-itail+1,theta))+yy(1)
            else
              xp = -rp*cos(valTh(l-itail+1,theta))+yy(ny)
            endif
            ip1 = int((xp-x0)/dy)+1
            if (cos(valTh(l-itail+1,theta))<0._pr) then
              ip2 = ip1+1
            else
              ip2 = ip1-1
            endif
            LS1 = fish_shape(l,ip1)
            LS2 = fish_shape(l,ip2)

            x1 = yy(ip1)
            if (cos(valTh(l-itail+1,theta))<0._pr) then
              x2 = x1 + dy
            else
              x2 = x1 - dy
            endif
            LSp = interpLS2D(LS1,x1,LS2,x2,xp) 
         endif
 
         bool = 1
         rr=rr+dy

         if (cos(valTh(l-itail+1,theta))<0._pr) then
           xp = -rr*cos(valTh(l-itail+1,theta))+yy(1)
         else
           xp = -rr*cos(valTh(l-itail+1,theta))+yy(ny)
         endif
         ip1 = int((xp-x0)/dy)+1
         if (cos(valTh(l-itail+1,theta))<0._pr) then
           ip2 = ip1+1
         else
           ip2 = ip1-1
         endif
         LS1 = fish_shape(l,ip1)
         LS2 = fish_shape(l,ip2)

         x1 = yy(ip1)
         if (cos(valTh(l-itail+1,theta))<0._pr) then
           x2 = x1 + dy
         else
           x2 = x1 - dy
         endif
         LSrr = interpLS2D(LS1,x1,LS2,x2,xp) 
      enddo
      rr = rp
      LSrr = LSp
      do while ((LSrr<0._pr).and.(rr<float(ny)))
         rp = rr
         LSp = LSrr
         rr=rr+0.1*dy
         if (cos(valTh(l-itail+1,theta))<0._pr) then
           xp = -rr*cos(valTh(l-itail+1,theta))+yy(1)
         else
           xp = -rr*cos(valTh(l-itail+1,theta))+yy(ny)
         endif
         ip1 = int((xp-x0)/dy)+1
         if (cos(valTh(l-itail+1,theta))<0._pr) then
           ip2 = ip1+1
         else
           ip2 = ip1-1
         endif
         LS1 = fish_shape(l,ip1)
         LS2 = fish_shape(l,ip2)
         x1 = yy(ip1)
         if (cos(valTh(l-itail+1,theta))<0._pr) then
           x2 = x1 + dy
         else
           x2 = x1 - dy
         endif
         LSrr = interpLS2D(LS1,x1,LS2,x2,xp) 
      enddo
      px = rhead
      py = 0.01_pr*cos(valTh(l-itail+1,theta))*abs(slicecontrol(theta,size(slicecontrol,2)-1,2)-points_courbe_equal(1,2)) +&
      points_courbe_equal(1,2)
      slicecontrolf(theta,size(slicecontrolf,2),1) = px 
      slicecontrolf(theta,size(slicecontrolf,2),2) = py 
      slicecontrol(theta,size(slicecontrol,2),1) = px 
      slicecontrol(theta,size(slicecontrol,2),2) = py 
   enddo
 
   !**
   !! Spline approximation of contour/surface control points
   !**
   dtb = 1._pr/(Ni-1)
   tb = dtb
   do theta=1,2
      slice_courbe(theta,1,1) = slicecontroli(theta,1,1)
      slice_courbe(theta,1,2) = slicecontroli(theta,1,2)
      valTheta(1,theta) = valTh(1,theta)
   enddo
   dsi = disttail/(Ni-1)
   dsf = disthead/(Nf-1)
   do theta=1,2
      tb = 0._pr
      do l=1,Ni-1
         nt = 1
         tbm = 0._pr
         tbp = 1._pr
         tb = tbm 
         bool = 0
         do while (((abs(sslice(theta))/(dx*dx)>eepsilon).and.(bool==0)).or.(nt==1))
            tb = (tbm + tbp)*0.5_pr
            nt = nt+1
            call pointsBezierN(slicecontroli(theta,:,:),tb,px,py)
            sslice(theta) = dotProd(points_courbe_equal(1,1)-(Ni-l)*dsi,points_courbe_equal(1,2)&
                 ,points_courbe_equal(1,1)-(Ni-l-1)*dsi,points_courbe_equal(1,2),px,py)
 
            if (sslice(theta)<0._pr) then
               tbm = tb
            else
               tbp = tb
            endif
            if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then!*dx) then
               write(*,*) "TBslice0  ",l," ",tb
               bool = 1
            endif
         enddo
         if (theta==1) write(*,*) "TBsliceI  ",l," ",tb," ",sslice(theta)," ",px," ",py," ",&
         points_courbe_equal(1,1)-(Ni-l-1)*dsi," ",points_courbe_equal(1,1)-(Ni-l)*dsi
         if (l<Ni) then
            call pointsBezierN(slicecontrol(theta,2:3,:),0.1_pr,pxx,pyy)
            tp = (slicecontrol(theta,2,1)-slicecontrol(theta,1,1))/(pxx-slicecontrol(theta,1,1))
            call poly2(slicecontrol(theta,1,1),slicecontrol(theta,2,1),pxx,tb*tp,tp,px) 
            call poly2(slicecontrol(theta,1,2),slicecontrol(theta,2,2),pyy,tb*tp,tp,py) 
            if (theta==1) write(*,*) "TBsliceII  ",l," ",tb*tp," ",tp," ",px," ",py," ",slicecontrol(theta,1,1)&
                 ," ",slicecontrol(theta,2,1)," ",pxx
         else
            call pointsBezierN(slicecontroli(theta,:,:),tb,px,py)
         endif
         slice_courbe(theta,l,1) = px
         slice_courbe(theta,l,2) = py
      enddo
      nt = 1
      tbm = 0._pr
      tbp = 1._pr
      tb = tbm 
      bool = 0
      do while (((abs(sslice(theta))/(dx*dx)>eepsilon).and.(bool==0)).or.(nt==1))
         tb = (tbm + tbp)*0.5_pr
         nt = nt+1
         call pointsBezierN(slicecontroli(theta,:,:),tb,px,py)
         sslice(theta) = dotProd(points_courbe_equal(1,1),points_courbe_equal(1,2)&
              ,points_courbe_equal(2,1),points_courbe_equal(2,2),px,py)
 
         if (sslice(theta)<0._pr) then
            tbm = tb
         else
            tbp = tb
         endif
         if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then
            bool = 1
         endif
      enddo
      if (theta==1) write(*,*) "TBsliceI  ",Ni," ",tb
      call pointsBezierN(slicecontroli(theta,:,:),tb,px,py)
      slice_courbe(theta,Ni,1) = px
      slice_courbe(theta,Ni,2) = py
   enddo
 
 
   dtb = 1._pr/(Ns-1)
   tb = dtb
   longslice = 0._pr
   do theta=1,2
      tb = 0._pr
      do l=1,Ns-1
         nt = 1
         tbm = 0._pr
         tbp = 1._pr
         tb = tbm 
         bool = 0
         do while (((abs(sslice(theta))/(dx*dx)>eepsilon).and.(bool==0)).or.(nt==1))
         if (theta==1) write(*,*) "TBslice  ",Ns," ",tb," ",nt," ",px," ",py
            tb = (tbm + tbp)*0.5_pr
            nt = nt+1
            call pointsBezierN(slicecontrolm(theta,:,:),tb,px,py)
            sslice(theta) = dotProd(points_courbe_equal(l,1),points_courbe_equal(l,2)&
                 ,points_courbe_equal(l+1,1),points_courbe_equal(l+1,2),px,py)
 
            if (sslice(theta)<0._pr) then
               tbm = tb
            else
               tbp = tb
            endif
            if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then
               bool = 1
            endif
         if (theta==1) write(*,*) "TBslice  ",Ns," ",tb," ",nt," ",px," ",py
         enddo
         if (theta==1) write(*,*) "TBslicem  ",l," ",tb
         call pointsBezierN(slicecontrolm(theta,:,:),tb,px,py)
         slice_courbe(theta,l+Ni-1,1) = px
         slice_courbe(theta,l+Ni-1,2) = py
      enddo
      nt = 1
      tbm = 0._pr
      tbp = 1._pr
      tb = tbm 
      bool = 0
      do while (((abs(sslice(theta))/(dx*dx)>eepsilon).and.(bool==0)).or.(nt==1))
         tb = (tbm + tbp)*0.5_pr
         nt = nt+1
         call pointsBezierN(slicecontrolm(theta,:,:),tb,px,py)
         sslice(theta) = dotProd(points_courbe_equal(Ns,1),points_courbe_equal(Ns,2)&
              ,points_courbe_equal(Ns,1)+dsf,points_courbe_equal(1,2),px,py)
 
         if (sslice(theta)<0._pr) then
            tbm = tb
         else
            tbp = tb
         endif
         if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then
            bool = 1
         endif
      enddo
      if (theta==1) write(*,*) "TBslice  ",Ns," ",tb
      call pointsBezierN(slicecontrolm(theta,:,:),tb,px,py)
      slice_courbe(theta,Ns+Ni-1,1) = px
      slice_courbe(theta,Ns+Ni-1,2) = py
   enddo
   dtb = 1._pr/(Nf-1)
   tb = dtb
   longslicef = 0._pr
   do theta=1,2
      slice_courbe(theta,size(slice_courbe,2),1) = slicecontrolf(theta,size(slicecontrolf,2),1)
      slice_courbe(theta,size(slice_courbe,2),2) = slicecontrolf(theta,size(slicecontrolf,2),2)
   enddo
    do theta=1,2
      tb = 0._pr
      do l=1,Nf-1
         nt = 1
         tbm = 0._pr
         tbp = 1._pr
         tb = tbm 
         bool = 0
         do while (((abs(sslice(theta))/(dx*dx)>eepsilon).and.(bool==0)).or.(nt==1))
            tb = (tbm + tbp)*0.5_pr
            nt = nt+1
            call pointsBezierN(slicecontrolf(theta,:,:),tb,px,py)
            sslice(theta) = dotProd(points_courbe_equal(Ns,1)+(l-1)*dsf,points_courbe_equal(1,2)&
                 ,points_courbe_equal(Ns,1)+l*dsf,points_courbe_equal(1,2),px,py)
 
            if (sslice(theta)<0._pr) then
               tbm = tb
            else
               tbp = tb
            endif
            if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then
               bool = 1
            endif
         enddo
         if (theta==1) write(*,*) "TBsliceF  ",l," ",tb
         call pointsBezierN(slicecontrolf(theta,:,:),tb,px,py)
         slice_courbe(theta,l+Ni-1+Ns-1,1) = px
         slice_courbe(theta,l+Ni-1+Ns-1,2) = py
      enddo
   enddo
   l = Nf-1+Ns-1+Ni-1
   if (l==Ni-1+Ns-1+Nf-1) then
      l = l+1
      do theta=1,2
         slice_courbe(theta,l,1) = slicecontrolf(theta,size(slicecontrolf,2),1)
         slice_courbe(theta,l,2) = slicecontrolf(theta,size(slicecontrolf,2),2)
         call pointsBezierN1D(valThtmp(theta,lph-itail+1:ihead-itail+1)&
              ,tb,alpha)
         valTheta(l,theta) = alpha
         longslicef(theta) = longslicef(theta) +&
              sqrt((slice_courbe(theta,l,1)-slice_courbe(theta,l-1,1))**2 + (slice_courbe(theta,l,2)-slice_courbe(theta,l-1,2))**2)
      enddo
      write(*,*) "ET VOILA F "!,l," ",slice_courbe(1,l,1)," ",slice_courbe(1,l,2)," ",slice_courbe(1,l,3)
   endif
   if (.not.(l==Ni-1+Ns-1+Nf)) write(*,*) "WARNIINNGG  F ",l
 
   !**
   !! Scaling the contour/surface Lagrangian markers
   !**
   do l=1,Ni-1+Ns-1+Nf
     do theta=1,2
       slice(theta,l,1) = (slice_courbe(theta,l,1)-x0)*longratio+x0
       slice(theta,l,2) = (slice_courbe(theta,l,2)-x0)*longratio+x0
     enddo
   enddo
   do l=1,size(points_courbe_equal,1)
     points_courbe_equal(l,1) = (points_courbe_equal(l,1)-x0)*longratio+x0
     points_courbe_equal(l,2) = (points_courbe_equal(l,2)-x0)*longratio+x0
   enddo

   do theta=1,2
     slice(theta,1,1) = (rtail-x0)*longratio+x0
     slice(theta,Ni-1+Ns-1+Nf,1) = (rhead-x0)*longratio+x0
   enddo
   long00 = long00*longratio
   disthead = disthead*longratio
   disttail = disttail*longratio
   dsi = disttail/(Ni-1)
   dsf = disthead/(Nf-1)
   dtb = 1._pr/(Ns-1)

   do l =1,Ni+Ns+Nf-2
      do theta=1,2
         if (theta==2/2) write(*,*) "FINAL THETA  ",l," ",theta," ",valTheta(l,theta)
      enddo
      !     tb = tb + dtb
   enddo
   slicetmp = slice
   
   call cpu_time(end_time)
   write(*,*) "******Frame 1 processing part 1 time = ", end_time-start_time, " sec******" ! Around 25 sec
 
   !**
   !! Computing the distance and angles between each Lagrangian markers and associated midline points
   !**
   call cpu_time(start_time)

   distslice(1+Ni,1) = dist(slice(1,1+Ni-1,1),slice(1,1+Ni-1,2)&
        ,points_courbe_equal(1,1),points_courbe_equal(1,2))
   distslice(1+Ni,2) = dist(slice(1,1+Ni-1,1),slice(1,1+Ni-1,2)&
        ,points_courbe_equal(2,1),points_courbe_equal(2,2))
   do theta=1,2
      valDist(1+Ni-1,theta) = dist(slice(theta,1+Ni-1,1),slice(theta,1+Ni-1,2)&
           ,points_courbe_equal(1,1),points_courbe_equal(1,2))
      cosTheta_tab(1+Ni-1,theta) = (slice(theta,1+Ni-1,2)-points_courbe_equal(1,2))/valDist(1+Ni-1,theta)  
   enddo
   do l=2,Ns-1
      distslice(l+Ni,1) = dist(slice(1,l+Ni-1,1),slice(1,l+Ni-1,2)&
           ,points_courbe_equal(l-1,1),points_courbe_equal(l-1,2))
      distslice(l+Ni,2) = dist(slice(1,l+Ni-1,1),slice(1,l+Ni-1,2)&
           ,points_courbe_equal(l+1,1),points_courbe_equal(l+1,2))
      write(*,*) "DOTPROD  ",l," ", dotProd(points_courbe_equal(l,1),points_courbe_equal(l,2),slice(1,l+Ni-1,1)&
           ,slice(1,l+Ni-1,2),points_courbe_equal(l+1,1),points_courbe_equal(l+1,2))
      do theta=1,2
         valDist(l+Ni-1,theta) = dist(slice(theta,l+Ni-1,1),slice(theta,l+Ni-1,2)&
              ,points_courbe_equal(l,1),points_courbe_equal(l,2))
         cosTheta_tab(l+Ni-1,theta) = (slice(theta,l+Ni-1,2)-points_courbe_equal(l,2))/valDist(l+Ni-1,theta) !(sqrt((slice(theta,l+Ni-1,2)-points_courbe_equal(l,2))**2+(slice(theta,l+Ni-1,3)-points_courbe_equal(l,3))**2)) 
      enddo
   enddo
   distslice(Ns+Ni,1) = dist(slice(1,Ns+Ni-1,1),slice(1,Ns+Ni-1,2)&
        ,points_courbe_equal(Ns-1,1),points_courbe_equal(Ns-1,2))
   distslice(Ns+Ni,2) = dist(slice(1,Ns+Ni-1,1),slice(1,Ns+Ni-1,2)&
        ,points_courbe_equal(Ns,1),points_courbe_equal(Ns,2))
   do theta=1,2
      valDist(Ns+Ni-1,theta) =&
           dist(slice(theta,Ns+Ni-1,1),slice(theta,Ns+Ni-1,2)&
           ,points_courbe_equal(Ns,1),points_courbe_equal(Ns,2))
      cosTheta_tab(Ns+Ni-1,theta) = (slice(theta,Ns+Ni-1,2)-points_courbe_equal(Ns,2))/valDist(Ns+Ni-1,theta) !(sqrt((slice(theta,Ns+Ni-1,2)-points_courbe_equal(Ns,2))**2+(slice(theta,Ns+Ni-1,3)-points_courbe_equal(Ns,3))**2)) 
   enddo
 
   dsi = disttail/(Ni-1)
   distslice(1,1) = dist(slice(1,1,1),slice(1,1,2)&
        ,points_courbe_equal(1,1)-(Ni-1)*dsi,points_courbe_equal(1,2))
   distslice(1,2) = dist(slice(1,1,1),slice(1,1,2)&
        ,points_courbe_equal(1,1)-(Ni-2)*dsi,points_courbe_equal(1,2))
   do theta=1,2 
      valDist(1,theta) = dist(slice(theta,1,1),slice(theta,1,2)&
           ,points_courbe_equal(1,1)-(Ni-1)*dsi,points_courbe_equal(1,2))
      write(*,*) "VALDIST  L==1  ",theta," ",valDist(1,theta)," ",slice(theta,1,1)," ",slice(theta,1,2)
      cosTheta_tab(1,theta) = abs(slice(theta,1,2)-points_courbe_equal(1,2))/valDist(1,theta) !(sqrt((slice(theta,l+Ni-1,2)-points_courbe_equal(l,2))**2+(slice(theta,l+Ni-1,3)-points_courbe_equal(l,3))**2)) 
      if (theta==2) cosTheta_tab(1,2) = -1*cosTheta_tab(1,1)
   enddo
   do l=2,Ni-1
      distslice(l,1) = dist(slice(1,l,1),slice(1,l,2)&
           ,points_courbe_equal(1,1)-(Ni-l+1)*dsi,points_courbe_equal(1,2))
      distslice(l,2) = dist(slice(1,l,1),slice(1,l,2)&
           ,points_courbe_equal(1,1)-(Ni-l-1)*dsi,points_courbe_equal(1,2))
      write(*,*) "DOTPRODi  ",l," ", dotProd(points_courbe_equal(1,1)-(Ni-l)*dsi,points_courbe_equal(1,2),slice(1,l,1)&
           ,slice(1,l,2),points_courbe_equal(1,1)-(Ni-l-1)*dsi,points_courbe_equal(1,2))
      do theta=1,2
         valDist(l,theta) = dist(slice(theta,l,1),slice(theta,l,2)&
              ,points_courbe_equal(1,1)-(Ni-l)*dsi,points_courbe_equal(1,2))
         cosTheta_tab(l,theta) = (slice(theta,l,2)-points_courbe_equal(1,2))/valDist(l,theta)  
      enddo
   enddo
   distslice(Ni,1) = dist(slice(1,Ni,1),slice(1,Ni,2)&
        ,points_courbe_equal(1,1)-dsi,points_courbe_equal(1,2))
   distslice(Ni,2) = dist(slice(1,Ni,1),slice(1,Ni,2)&
        ,points_courbe_equal(1,1),points_courbe_equal(1,2))
   do theta=1,2
      valDist(Ni,theta) = dist(slice(theta,Ni,1),slice(theta,Ni,2)&
           ,points_courbe_equal(1,1),points_courbe_equal(1,2))
      write(*,*) "VALDIST  L==2  ",theta," ",valDist(Ni,theta)," ",slice(theta,Ni,1)," ",slice(theta,Ni,2)
      cosTheta_tab(Ni,theta) = (slice(theta,Ni,2)-points_courbe_equal(1,2))/valDist(Ni,theta)  
   enddo
 
   dsf = disthead/(Nf-1)
   distslice(1+Ni+Ns,1) = dist(slice(1,1+Ni-1+Ns-1,1),slice(1,1+Ni-1+Ns-1,2)&
        ,points_courbe_equal(Ns,1),points_courbe_equal(1,2))
   distslice(1+Ni+Ns,2) = dist(slice(1,1+Ni-1+Ns-1,1),slice(1,1+Ni-1+Ns-1,2)&
        ,points_courbe_equal(Ns,1)+dsf,points_courbe_equal(1,2))
   do theta=1,2
      valDist(1+Ni-1+Ns-1,theta) =&
           dist(slice(theta,1+Ni-1+Ns-1,1),slice(theta,1+Ni-1+Ns-1,2)&
           ,points_courbe_equal(Ns,1),points_courbe_equal(1,2))
      cosTheta_tab(1+Ni-1+Ns-1,theta) = (slice(theta,1+Ni-1+Ns-1,2)-points_courbe_equal(1,2))/valDist(1+Ni-1+Ns-1,theta)  
   enddo
   do l=2,Nf-1
      distslice(l+Ni+Ns,1) =&
           dist(slice(1,l+Ni-1+Ns-1,1),slice(1,l+Ni-1+Ns-1,2)&
           ,points_courbe_equal(Ns,1)+(l-2)*dsf,points_courbe_equal(1,2))
      distslice(l+Ni+Ns,2) =&
           dist(slice(1,l+Ni-1+Ns-1,1),slice(1,l+Ni-1+Ns-1,2)&
           ,points_courbe_equal(Ns,1)+l*dsf,points_courbe_equal(1,2))
      write(*,*) "DOTPRODf  ",l," ", dotProd(points_courbe_equal(Ns,1)+(l-1)*dsf,points_courbe_equal(1,2),slice(1,l+Ni-1+Ns-1,1)&
           ,slice(1,l+Ni-1+Ns-1,2),points_courbe_equal(Ns,1)+l*dsf,points_courbe_equal(1,2))
      do theta=1,2
         valDist(l+Ni-1+Ns-1,theta) =&
              dist(slice(theta,l+Ni-1+Ns-1,1),slice(theta,l+Ni-1+Ns-1,2)&
              ,points_courbe_equal(Ns,1)+(l-1)*dsf,points_courbe_equal(1,2))
         cosTheta_tab(l+Ni-1+Ns-1,theta) = (slice(theta,l+Ni-1+Ns-1,2)-points_courbe_equal(1,2))/valDist(l+Ni-1+Ns-1,theta)  
      enddo
   enddo
   distslice(Nf+Ni+Ns,1) =&
        dist(slice(1,Nf+Ni-1+Ns-1,1),slice(1,Nf+Ni-1+Ns-1,2)&
        ,points_courbe_equal(Ns,1)+(Nf-2)*dsf,points_courbe_equal(1,2))
   distslice(Nf+Ni+Ns,2) =&
        dist(slice(1,Nf+Ni-1+Ns-1,1),slice(1,Nf+Ni-1+Ns-1,2)&
        ,points_courbe_equal(Ns,1)+(Nf-1)*dsf,points_courbe_equal(1,2))
   do theta=1,2
      valDist(Nf+Ni-1+Ns-1,theta) =&
           dist(slice(theta,Nf+Ni-1+Ns-1,1),slice(theta,Nf+Ni-1+Ns-1,2)&
           ,points_courbe_equal(Ns,1)+(Nf-1)*dsf,points_courbe_equal(1,2))
      cosTheta_tab(Nf+Ni-1+Ns-1,theta) = (slice(theta,Nf+Ni-1+Ns-1,2)-points_courbe_equal(1,2))/valDist(Nf+Ni-1+Ns-1,theta)  
   enddo
 
   call cpu_time(end_time)
   write(*,*) "Frame 1 Langrangian markers processing time = ", end_time-start_time, " sec" ! Negligible (0 sec measured???)

   !**
   !! Now, we can write in output the initial shape of the zebrafish (Lagrangian markers / midline)
   !**
 
   ! call cpu_time(start_time)

   ! open(unit=91,file= trim(target_folder)//'/theta2D00.txt',status='unknown')
   ! open(unit=84,file= trim(target_folder)//'/right2D00.txt',status='unknown')
   ! open(unit=83,file= trim(target_folder)//'/left2D00.txt',status='unknown')
   ! open(unit=86,file= trim(target_folder)//'/controlright2D00.txt',status='unknown')
   ! open(unit=85,file= trim(target_folder)//'/controlleft2D00.txt',status='unknown')
   ! open(unit=79,file= trim(target_folder)//'/skelet2D00.vtk',status='unknown')
   ! open(unit=81,file= trim(target_folder)//'/skelett2D00.txt',status='unknown')
   ! open(unit=82,file= trim(target_folder)//'/skeletteq2D00.txt',status='unknown')
   ! file_id = 79
   ! call write_vtk_header(file_id, nx, ny, nz, dx, dy, dz, n_dim)

   ! open(unit=78,file= trim(target_folder)//'/surf2D00.vts',status='unknown')
   ! open(unit=92,file= trim(target_folder)//'/surf2D00.dat',status='unknown')
   ! write(78,'(a)') "<?xml version=""1.0""?>"
   ! write(78,'(a)') "<VTKFile type=""StructuredGrid"" version=""0.1"" byte_order=""LittleEndian"" compressor=&
   ! &     ""vtkZLibDataCompressor"">"
   ! write(78,'(a,I3,a,I3,a)') "<StructuredGrid WholeExtent=""0 ",1," 0 ",Ni-1+Ns+Nf-1-1," 0 0"">"
   ! write(78,'(a,I3,a,I3,a)') "<Piece Extent=""0 ",1," 0 ",Ni-1+Ns+Nf-1-1," 0 0"">"
   ! write(78,'(a)') "<PointData >"
   ! write(78,'(a)') "</PointData>"
   ! write(78,'(a)') "<CellData>"
   ! write(78,'(a)') "</CellData>"
   ! write(78,'(a)') "<Points>"
   ! write(78,'(a)') "<DataArray NumberOfComponents=""3"" type=""Float64"" format=""ascii"" >"  
   ! do l=1,Ni+Ns+Nf-2
   !    do theta=1,2 !indextheta(l)
   !       write(78,*) slice(theta,l,1)," ",slice(theta,l,2)," ",1
   !       write(92,*) slice(theta,l,1)," ",slice(theta,l,2)
   !       if (l==Ni+1+10) write(*,*) "cmpslice10  ",theta," ", dist(slice(theta,l,1),slice(theta,l,2)&
   !            ,slicetmp(theta,l,1),slicetmp(theta,l,2))
   !       if (l==Ni+1+5) write(*,*) "cmpslice5  ",theta," ", dist(slice(theta,l,1),slice(theta,l,2)&
   !            ,slicetmp(theta,l,1),slicetmp(theta,l,2))
   !       if (l==Ni+1+2) write(*,*) "cmpslice2  ",theta," ", dist(slice(theta,l,1),slice(theta,l,2)&
   !            ,slicetmp(theta,l,1),slicetmp(theta,l,2))
   !       if (l==Ni+1) write(*,*) "cmpsliceNI  ",theta," ", dist(slice(theta,l,1),slice(theta,l,2)&
   !            ,slicetmp(theta,l,1),slicetmp(theta,l,2))
   !    enddo
   ! enddo
   ! write(78,'(a)') "</DataArray>"
   ! write(78,'(a)') "</Points>"
   ! write(78,'(a)') "</Piece>"
   ! write(78,'(a)') "</StructuredGrid>"
   ! write(78,'(a)') "</VTKFile>"
   ! close(78)    
   ! close(92)    

   ! call cpu_time(end_time) !negligible time
   ! write(*,*) "Save surf2D00.dat and .vts time = ", end_time-start_time, " sec from savings start" ! Negligible (0.016 sec)
 
   ! do k=1,ny
   !    do j=1,nx
   !       write(79,*) fish_shape(j,k)*(1-skel(j,k)) + minval(fish_shape)*skel(j,k)
   !    enddo
   ! enddo
 
   ! do l=1,Ns
   !    write(81,*) points_courbe(l,1)," ",points_courbe(l,2)
   !    write(82,*) points_courbe_equal(l,1)," ",points_courbe_equal(l,2)
   ! enddo
 
   ! do l=1,size(slice,2) 
   !    theta = 1
   !    write(83,*) slice(theta,l,1)," ",slice(theta,l,2)
   !    theta = 1+2/2
   !    write(84,*) slice(theta,l,1)," ",slice(theta,l,2)
   ! enddo
  
   ! do l=1,size(slicecontrol,2) 
   !    theta = 1
   !    write(85,*) slicecontrol(theta,l,1)," ",slicecontrol(theta,l,2)
   !    theta = 1+2/2
   !    write(86,*) slicecontrol(theta,l,1)," ",slicecontrol(theta,l,2)
   ! enddo
   ! do l=1,size(valTheta,1)
   !    if (l<=size(points_courbe_equal,1)) then
   !       write(91,*) points_courbe_equal(l,1)," ",valTheta(l,theta)
   !    else
   !       write(91,*) points_courbe_equal(Ns,1)+(l-Ns)*dsf," ",valTheta(l,theta)
   !    endif
   ! enddo
   ! close(91)
   ! close(79)    
   ! close(85)
   ! close(86)
   ! close(81)
   ! close(82)
   ! close(83)
   ! close(84)
   deallocate(midline,midlinebis)
   deallocate(points_control)
   deallocate(slicecontrolm,slicecontrolf,slicecontroli)
   deallocate(slicecontrol,slicecontroltmp)
   deallocate(valTh,valThtmp)

   ! call cpu_time(end_time)
   ! write(*,*) "Save first frame files time = ", end_time-start_time, " sec" ! Really slow : Around 4 minutes to iterate and save files 

   !!FILM

   N =  200 
   nx = 200 
   ny = 200 
   nz = 200
   dx = 0.001_pr*0.0256
   dy = 0.001_pr*0.0256
   dz = 0.001_pr*0.0256
 
   deallocate(fish_shape, rhoSlices, gradPhi,xx,yy,Nseed,skel,skel2,skel3,skeltmp)
   allocate(rhoSlices(nx,ny), gradPhi(nx,ny))
   allocate(xx(nx))
   allocate(yy(ny)) 
   allocate(Nseed(nx,ny),skel(nx,ny),skel2(nx,ny),skel3(nx,ny),skeltmp(nx,ny))
 
   do i=1,nx
      xx(i) = x0 + (float(i)-1)*dx
   enddo
   do j=1,ny
      yy(j) = x0 + (float(j)-1)*dy
   enddo
 
   open(unit=80,file= trim(target_folder)//'/skeleton.txt',status='unknown')
   write(80,*) kt,"    ",nl," ",long3," ",long2," ",long
 
   !**
   !! Loop over each experimental image to construct the deformed zebrafish shape (Lagrangian markers) based on the pre-built shape
   !**

   call cpu_time(start_time)

   t = 0
   iter = 0
   idisplay = 10
   iUpdateDist = 5
   !! FILM
   do kt = 1, n_frames 
      write(*,*) ""
      write(*,*) "Image iteration num: ", kt
      !! FILM
      open(unit=78,file="IMAGES4/Image_"&
      //str(kt+picNum-1)//".dat",&
           status='unknown')
 
   !**
   !! Computation of the level-set 
   !**
      do k=1,nx
         do j=1,ny
            read(78,*) rhoSlices(k,j)
         enddo
      enddo
      close(78)
      rhoSlices = rhoSlices - 0.5_pr*(maxval(rhoSlices) + minval(rhoSlices))
      !! FILM
      dx = 0.001_pr*0.0256*1E6
      dy = 0.001_pr*0.0256*1E6
      call updateDistance(rhoSlices,gradPhi)
      dx = 0.001_pr*0.0256
      dy = 0.001_pr*0.0256
      do i=1,nx
         xx(i) = x0 + (float(i)-1)*dx
      enddo
      do j=1,ny
         yy(j) = x0 + (float(j)-1)*dy
      enddo
 
   !**
   !! Computation of the midline (skel: skeleton based one the level-set, skel2: gradient of the level-set)
   !**
      Nseed = 0
      skel = 0
      skel2 = 0
      sum_grad_thresh = 0
      do i=2,ny-1
         do j=2,nx-1
            !! FILM
 
            !serie5_18_angle2
            if ((gradPhi(j,i)<0.62).and.(rhoSlices(j,i)>0._pr)) skel(j,i) = 1
            if ((gradPhi(j,i)<0.94).and.(rhoSlices(j,i)>0._pr)) then
               skel2(j,i) = 1
               sum_grad_thresh = sum_grad_thresh + 1
            endif
            if ((kt==60).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1
            if ((kt==365).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1
            if ((kt==501).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1      
            if ((kt==162).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1      
            if ((kt==561).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1  
            if ((kt==247).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1
            if ((kt==303).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1
            if ((kt==304).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1


            ! if ((kt==60).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1
            ! if ((kt==365).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1
            ! if ((kt==561).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1
            ! if ((kt==162).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1                                                                                         
            ! if ((kt==247).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1                                                                                         
            ! if ((kt==303).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1                                                                                         
            ! if ((kt==304).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1                                                                                         
            ! if ((kt==501).and.(gradPhi(j,i)<0.98).and.(rhoSlices(j,i)>0._pr)) skel2(j,i) = 1      
      
 
         !**
         !! First: Thresholds for extracting the skeleton and its gradient
         !**
         enddo
      enddo
      write(*,*) "Len of sum_grad_thresh at frame ", kt, " = ", sum_grad_thresh

      !!FILM (threshold no more used)

      seuil0 = 0.98
      !**
      !! Second step: we reduce the gradient map
      !**

      reduced_grad_thresh = 0

      skel3 = 0
      do i=2,nx-1
         do j=2,ny-1
            Nseed(i,j) = skel2(i+1,j) + skel2(i+1,j+1) + skel2(i,j+1) + skel2(i-1,j+1) + skel2(i-1,j) + skel2(i-1,j-1) +&
            skel2(i,j-1) + skel2(i+1,j-1)
            if ((skel2(i,j)==1).and.(Nseed(i,j)>5)) then
               skel3(i,j) = 1
               reduced_grad_thresh  = reduced_grad_thresh + 1
            endif
         enddo
      enddo
      write(*,*) "Len of reduced_grad_thresh at frame ", kt, " = ", reduced_grad_thresh 

      skel2 = skel3
 
      l = 1
      pix1 = 0
      pix2 = 0

      !**
      !! Third step: we complete the level-set map
      !**
   !!FILM
      skel3 = skel
      do i=2,nx-1
         do j=2,ny-1
            Nseed(i,j) = skel(i+1,j) + skel(i+1,j+1) + skel(i,j+1) + skel(i-1,j+1) + skel(i-1,j) + skel(i-1,j-1) + skel(i,j-1) +&
                 skel(i+1,j-1)
            if ((skel(i,j)==0).and.(Nseed(i,j)==2).and.(skel(i-1,j)==1).and.(skel(i+1,j+1)==1)) skel3(i,j) = 1
            if ((skel(i,j)==0).and.(Nseed(i,j)==2).and.(skel(i-1,j)==1).and.(skel(i+1,j-1)==1)) skel3(i,j) = 1
            if ((skel(i,j)==0).and.(Nseed(i,j)==2).and.(skel(i,j-1)==1).and.(skel(i+1,j+1)==1)) skel3(i,j) = 1
            if ((skel(i,j)==0).and.(Nseed(i,j)==2).and.(skel(i,j-1)==1).and.(skel(i-1,j+1)==1)) skel3(i,j) = 1
            if ((skel(i,j)==0).and.(Nseed(i,j)==2).and.(skel(i,j-1)==1).and.(skel(i,j+1)==1)) skel3(i,j) = 1
            if ((skel(i,j)==0).and.(Nseed(i,j)==2).and.(skel(i-1,j)==1).and.(skel(i+1,j)==1)) skel3(i,j) = 1
            if ((skel(i,j)==0).and.(Nseed(i,j)==2).and.(skel(i-1,j-1)==1).and.(skel(i+1,j+1)==1)) skel3(i,j) = 1
            if ((skel(i,j)==0).and.(Nseed(i,j)==2).and.(skel(i-1,j+1)==1).and.(skel(i+1,j-1)==1)) skel3(i,j) = 1
         enddo
      enddo

      !**
      !! Fourth step: MANUAL CORRECTIONS of level-set and gradient maps
      !**
      !! FILM
      !!!! filme5_18_angle2
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


      ! if(kt==79)skel3(155,109)=0
      ! if(kt==219)skel3(136,85)=0
      ! if(kt==370)skel3(158,95)=0 ok
      ! if(kt==395)skel3(152,93)=0 ok
      ! if(kt==395)skel3(154,93)=0 ok
      ! if(kt==395)skel3(153,92)=0 ok
      ! if(kt==392)skel3(152,92)=0 ok
      ! if(kt==440)skel3(160,109)=0 ok
      ! if(kt==480)skel3(157,107)=0
      ! if(kt==523)skel3(157,95)=0
      ! if(kt==103)skel3(151,110)=0
 
      skel = skel3
      !skeltmp = skel
      skel3 = skel
      skeltmp = skel
      tmpbool = skel
      
      !**
      !! Fifth step: construction of a 1-pixel wide skeleton
      !**
      call filterskel(skel,rhoSlices)
      !!FILM

      tmpbool = skel
      !skeltmp = skel

      sizeSkel = sum(skel)
      write(*,*) "Longlong 00  : ",long3," ",long2," ",long
 
      !**
      !! Sixth step: construction of the midline from the 2D-map (skeleton)
      !**
      allocate(midline(sizeSkel,2))
      write(*,*) "okcut"
      l=1
      do i=2,nx-1
         do j=2,ny-1
            if ((skel(i,j)==1).and.(l==1)) then
               Nseed(i,j) = skel(i+1,j) + skel(i+1,j+1) + skel(i,j+1) + skel(i-1,j+1) + skel(i-1,j) + skel(i-1,j-1) + skel(i,j-1) +&
                    skel(i+1,j-1)

               if (((Nseed(i,j)==1).and.((abs(xskelL-i)<50).and.(abs(yskelL-j)<50))).or.((Nseed(i,j)==1).and.(kt==1))) then 
                  midline(l,1) = i
                  midline(l,2) = j
                  l = l+1
               endif
            endif
         enddo
      enddo
      !FILM

      write(*,*) "okcut"

      write(*,*) "CHECK midline : ",midline(1,1)," ",midline(1,2)," ",midline(size(midline,1),1)," ",midline(size(midline,1),2)
      xskelL = midline(1,1)
      yskelL = midline(1,2)
      boolskel=0
      nl = sum(skel)
      do l=2,sum(skel)
         do i=2,nx-1
            do j=2,ny-1
               if ((skel(i,j)==1).and.(boolskel==0)) then
                  Nseed(i,j) = skel(i+1,j) + skel(i+1,j+1) + skel(i,j+1) + skel(i-1,j+1) + skel(i-1,j) + skel(i-1,j-1) + &
                  & skel(i,j-1) + skel(i+1,j-1)
                  if ((Nseed(i,j)>1).and.((midline(l-1,1)==i).or.(midline(l-1,1)==i-1).or.(midline(l-1,1)==i+1)).and.&
                       ((midline(l-1,2)==j).or.(midline(l-1,2)==j-1).or.(midline(l-1,2)==j+1)).and.(appartient(midline,i,j,l)&
                       .neqv..true.)) then 
                     midline(l,1) = i
                     midline(l,2) = j
                  endif
                  if ((Nseed(i,j)==1).and.((midline(l-1,1)==i).or.(midline(l-1,1)==i-1).or.(midline(l-1,1)==i+1)).and.&
                       ((midline(l-1,2)==j).or.(midline(l-1,2)==j-1).or.(midline(l-1,2)==j+1)).and.(appartient(midline,i,j,l)&
                       .neqv..true.)) then 
                     midline(l,1) = i
                     boolskel = 1
                     nl=l
                     midline(l,2) = j
                  if (Nseed(i,j)==1) write(*,*) "extrem  ",i," ",j," ",l
                  endif
               endif
            enddo
         enddo
      enddo
      write(*,*) "CHECK midline : ",midline(1,1)," ",midline(1,2)," ",midline(size(midline,1),1)," ",midline(size(midline,1),2)&
           ," ",size(midline,1)," ",sum(skel)," ",nl
      write(*,*) "CHECK midline BIS : ",midline(1,1)," ",midline(1,2)," ",midline(nl,1)," ",midline(nl,2)," ",nl
 
      if (nl<10) then
              write(*,*) "error midline queue fourche"
 
              do l=1,nl
                skel(midline(l,1),midline(l,2)) = 0
              enddo
 
              l=1
              do i=2,N-1
                 do j=2,ny-1
                    if ((skel(i,j)==1).and.(l==1)) then
                       Nseed(i,j) = skel(i+1,j) + skel(i+1,j+1) + skel(i,j+1) + skel(i-1,j+1) + skel(i-1,j) + skel(i-1,j-1) +&
                       skel(i,j-1) + skel(i+1,j-1)
                       if (((Nseed(i,j)==1).and.((abs(xskelL-i)<50).and.(abs(yskelL-j)<50))).or.((Nseed(i,j)==1).and.(kt==1))) then 
                          midline(l,1) = i
                          midline(l,2) = j
                          l = l+1
                       endif
                    endif
                 enddo
              enddo
              write(*,*) "CHECK midline : ",midline(1,1)," ",midline(1,2)," ",midline(size(midline,1),1)," ",&
              midline(size(midline,1),2)
              xskelL = midline(1,1)
              yskelL = midline(1,2)
              boolskel=0
              do l=2,sum(skel)
                 do i=2,N-1
                    do j=2,ny-1
                       if ((skel(i,j)==1).and.(boolskel==0)) then
                          Nseed(i,j) = skel(i+1,j) + skel(i+1,j+1) + skel(i,j+1) + skel(i-1,j+1) + skel(i-1,j) + skel(i-1,j-1) +&
                          skel(i,j-1) + skel(i+1,j-1)
                          if ((Nseed(i,j)>1).and.((midline(l-1,1)==i).or.(midline(l-1,1)==i-1).or.(midline(l-1,1)==i+1)).and.&
                               ((midline(l-1,2)==j).or.(midline(l-1,2)==j-1).or.(midline(l-1,2)==j+1)).and.&
                               (appartient(midline,i,j,l).neqv..true.)) then 
                             midline(l,1) = i
                             midline(l,2) = j
                          endif
                          if ((Nseed(i,j)==1).and.((midline(l-1,1)==i).or.(midline(l-1,1)==i-1).or.(midline(l-1,1)==i+1)).and.&
                               ((midline(l-1,2)==j).or.(midline(l-1,2)==j-1).or.(midline(l-1,2)==j+1)).and.&
                               (appartient(midline,i,j,l).neqv..true.)) then 
                             midline(l,1) = i
                             boolskel = 1
                             nl=l
                             midline(l,2) = j
                          endif
                       endif
                    enddo
                 enddo
              enddo
      endif
 
 
 
 
 
 
      !**
      !! Generating a 2D-map from the midline (for visualisation purpose)
      !**
 
      skel = 0
      do l=1,nl
        skel(midline(l,1),midline(l,2)) = 1
      enddo
 
      !**
      !! Seventh step: searching the endpoint of the midline (and cutting the rest) TRACKING STEP
      !**
      if (kt==1) then
              ic = -1
              jc = -1
      endif
      if (kt==1) ic = midline(nl,1)
      if (kt==1) jc = midline(nl,2)
      write(*,*) "okcut"
      inc = 0.01 !incremental step no more used)
      ic0 = ic
      jc0 = jc
      !skeltmp = skel
      skeltmp = skel2
      call cutheadskel(midline,skel,skel2,nl,ic,jc,kt) !! Note all variables are modified herein !
 !     !! CONTACT ATTENTION
 !     if (kt==384) ic = 134
 !     if (kt==384) jc = 77
      !if (kt>1) call cutheadskel(midline,skel,skel2,nl,ic,jc,kt)
 
      !**
      !! Generating 2D-maps for visualisation purpose
      !**
      skel3 = 0
      do l=1,nl
        skel3(midline(l,1),midline(l,2)) = skel2(midline(l,1),midline(l,2))
      enddo
      skel2 = skel3
      !skeltmp = skel2
      niter=1
      seuil = seuil0 !no more used (deprecated)
      !!FILM
 
      if (sum(skel2)>11) write(*,*) "ATTENTIONSKEL20  ",kt,"     ",sum(skel2)
      if (sum(skel2)<4) write(*,*) "ATTENTIONSKEL02  ",kt,"     ",sum(skel2)
      !skel2 = skel3
      write(*,*) "okcut"
 
      !**
      !! Heighth step: filling the midline
      !**
      xskelL = midline(1,1)
      yskelL = midline(1,2)
      boolskel=0
      do l=2,size(midline,1)
         do i=2,N-1
            do j=2,ny-1
               if ((skel(i,j)==1).and.(boolskel==0)) then
                  Nseed(i,j) = skel(i+1,j) + skel(i+1,j+1) + skel(i,j+1) + skel(i-1,j+1) + skel(i-1,j) + skel(i-1,j-1) + &
                  &     skel(i,j-1) + skel(i+1,j-1)
                  if ((Nseed(i,j)>1).and.((midline(l-1,1)==i).or.(midline(l-1,1)==i-1).or.(midline(l-1,1)==i+1)).and.&
                       ((midline(l-1,2)==j).or.(midline(l-1,2)==j-1).or.(midline(l-1,2)==j+1)).and.(appartient(midline,i,j,l)&
                       .neqv..true.)) then 
                     midline(l,1) = i
                     midline(l,2) = j
                  endif
                  if ((Nseed(i,j)==1).and.((midline(l-1,1)==i).or.(midline(l-1,1)==i-1).or.(midline(l-1,1)==i+1)).and.&
                       ((midline(l-1,2)==j).or.(midline(l-1,2)==j-1).or.(midline(l-1,2)==j+1)).and.(appartient(midline,i,j,l)&
                       .neqv..true.)) then 
                     midline(l,1) = i
                     boolskel = 1
                     nl=l
                     midline(l,2) = j
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
         call filtermidline(midline,rhoSlices)
      pix1 = midline(nl,1)
      pix2 = midline(nl,2)
 
      !**
      !! Tenth step: deleting the first pixel if horizontal or vertical
      !**
      if (.not.((midline(1,1)==midline(2,1)).or.(midline(1,2)==midline(2,2)))) then
      do l=1,size(midline,1)-1
      midline(l,1) = midline(l+1,1)
      midline(l,2) = midline(l+1,2)
      enddo
      nl=nl-1
      endif
      endif
      write(*,*) "CHECK midline : ",midline(1,1)," ",midline(1,2)," ",midline(size(midline,1),1)," ",midline(size(midline,1),2)&
           ," ",size(midline,1)
      write(*,*) "CHECK midline BIS : ",midline(1,1)," ",midline(1,2)," ",midline(nl,1)," ",midline(nl,2)," ",nl
 
      !**
      !! Final 2D-map of the midline
      !**
      !skel2 = skel
      skel = 0
      do l=1,nl
        skel(midline(l,1),midline(l,2)) = 1
      enddo

!! From HERE
!       allocate(points_control(nl,2))
!       do l=1,nl
!          points_control(l,1) = xx(midline(l,1))
!          points_control(l,2) = yy(midline(l,2))
!       enddo
!       !!TESTNEZ
!       !!     points_control(nl,1) = xc
!       !!     points_control(nl,2) = yc
 
!       !**
!       !! Spline approximation of the midline
!       !**
!       points_courbe(1,1) = points_control(1,1)
!       points_courbe(1,2) = points_control(1,2)
!       tb = dtb
!       l = 1
!       long = 0._pr
!       do while ((tb<1._pr).and.(l+1<Ns+1))
!          l = l+1
!          call pointsBezierN(points_control,tb,px,py)
!          points_courbe(l,1) = px
!          points_courbe(l,2) = py
!          long = long + sqrt((points_courbe(l,1)-points_courbe(l-1,1))**2 + (points_courbe(l,2)-points_courbe(l-1,2))**2)
!          tb = tb+dtb
!       enddo
!       if (l==Ns-1) then
!          l=l+1
!          write(*,*) " ET VOILA"
!          points_courbe(size(points_courbe_equal,1),1) = points_control(size(points_control,1),1)
!          points_courbe(size(points_courbe_equal,1),2) = points_control(size(points_control,1),2)
!          long = long + sqrt((points_courbe(l,1)-points_courbe(l-1,1))**2 + (points_courbe(l,2)-points_courbe(l-1,2))**2)
!       endif
!       if (.not.(l==Ns)) write(*,*) "WARNING1  p ",l
!       long2 = long
!       write(*,*) "POINTS_COURBE_EQUAL  L==1  ",points_control(1,1)," ",points_control(1,2)
 
!       !**
!       !! Tail extrapolation 
!       !**
!       !! ICI ON RALLONGE LA QUEUE POUR NORMALISER LA LONGUEUR : ON MOYENNE LA
!       !DIRECTION SUR LES X DERNIERS POINTS

!       if (long00>long2) then
!       !!x1 = points_courbe(1,1)
!       !!y1 = points_courbe(1,2)
!       !!x1 = points_courbe(nint(0.02*Ns),1)
!       !!y1 = points_courbe(nint(0.02*Ns),2)
!       !!x2 = points_courbe(nint(0.03*Ns),1)
!       !!y2 = points_courbe(nint(0.03*Ns),2)
!       !!x2 = points_courbe(nint(0.04*Ns),1)
!       !!y2 = points_courbe(nint(0.04*Ns),2)
!          l=5
!       !endif
!       x1=0.0
!       y1=0.0
!       do i=1,l
!          x1 = x1+points_control(i,1)
!          y1 = y1+points_control(i,2)
!       enddo
!       x1=x1*1.0/l
!       y1=y1*1.0/l
!       write(*,*) "X1Y1  ",x1," ",y1
!       do i=1,1

!       x2 = points_control(6,1)
!       y2 = points_control(6,2)

!       rt = abs(long00-long2) + sqrt((x1-points_control(1,1))**2 + (y1-points_control(1,2))**2)
!       if (abs(x1-x2)/dx<eepsilon) then
!          xt = x1
!          yt = y1+rt
!          if (((y2-y1)*(yt-y1))>0._pr) then
!            write(*,*) "JEPASSEICIAUSSI"
!            yt = y1-rt
!          endif
!       else
!          xt = rt/(sqrt(1+((y2-y1)**2)/((x2-x1)**2))) + x1
!          yt = (y2-y1)/(x2-x1)*(xt-x1) + y1
  
!          if (((x2-x1)*(xt-x1)+(y2-y1)*(yt-y1))>0._pr) then
!            write(*,*) "JEPASSEICI"
!            xt = -rt/(sqrt(1+((y2-y1)**2)/((x2-x1)**2))) + x1
!            yt = (y2-y1)/(x2-x1)*(xt-x1) + y1
 
!          endif
!       endif
!       write(*,*) "CHECKCJECKi  ",i,"    ",x1," ",y1," ",xt," ",yt,"  RT  ",rt," ",long00," ",long2&
!       ,"  PS  ",(x2-x1)*(xt-x1)+(y2-y1)*(yt-y1)
!       if (i==1) points_control(1,1) = 0._pr
!       if (i==1) points_control(1,2) = 0._pr
!       points_control(1,1) = points_control(1,1)+xt
!       points_control(1,2) = points_control(1,2)+yt
!       enddo
 
!       endif
!       write(*,*) "POINTS_COURBE_EQUAL  L==1  ",points_control(1,1)," ",points_control(1,2)
!       points_courbe(1,1) = points_control(1,1)
!       points_courbe(1,2) = points_control(1,2)
 
 
!       tb = dtb
!       l = 1
!       long = 0._pr
!       do while ((tb<1._pr).and.(l+1<Ns+1))
!          l = l+1
!          call pointsBezierN(points_control,tb,px,py)
!          points_courbe(l,1) = px
!          points_courbe(l,2) = py
!          long = long + sqrt((points_courbe(l,1)-points_courbe(l-1,1))**2 + (points_courbe(l,2)-points_courbe(l-1,2))**2)
!          tb = tb+dtb
!       enddo
!       if (l==Ns-1) then
!          write(*,*) " ET VOILA"
!          l=l+1
!          points_courbe(size(points_courbe_equal,1),1) = points_control(size(points_control,1),1)
!          points_courbe(size(points_courbe_equal,1),2) = points_control(size(points_control,1),2)
!          long = long + sqrt((points_courbe(l,1)-points_courbe(l-1,1))**2 + (points_courbe(l,2)-points_courbe(l-1,2))**2)
!       endif
!       if (.not.(l==Ns)) write(*,*) "WARNING2  p ",l
   
!       write(*,*) "IMPORTANT  ",long00," ",long2," ",long,"  RT ",rt," ",sqrt((x1-xt)**2+(y1-yt)**2)&
!            ,"  coord ",x1," ",y1," ",x2," ",y2
!       long2 = long
 
!       !**
!       !! Tail extrapolation 
!       !**
!       !! ICI ON RALLONGE LA QUEUE POUR NORMALISER LA LONGUEUR : ON MOYENNE LA
!       !DIRECTION SUR LES X DERNIERS POINTS
!       l=1
!       do i=1,l
!       x1 = points_control(i,1)
!       y1 = points_control(i,2)
!       x2 = points_control(2,1)
!       y2 = points_control(2,2)
!       rt = abs(long00-long2)
!       if (abs(x1-x2)/dx<eepsilon) then
!          xt = x1
!          yt = y1+rt
 
!          if (((y2-y1)*(yt-y1))>0._pr) then
!            write(*,*) "JEPASSEICIAUSSI"
!            yt = y1-rt
!          endif
!       else
!          xt = rt/(sqrt(1+((y2-y1)**2)/((x2-x1)**2))) + x1
!          yt = (y2-y1)/(x2-x1)*(xt-x1) + y1
 
!          if (((x2-x1)*(xt-x1)+(y2-y1)*(yt-y1))>0._pr) then
!            write(*,*) "JEPASSEICI"
!            xt = -rt/(sqrt(1+((y2-y1)**2)/((x2-x1)**2))) + x1
!            yt = (y2-y1)/(x2-x1)*(xt-x1) + y1
!          endif
!       endif
!       write(*,*) "CHECKCJECKi  ",i,"    ",x1," ",y1," ",xt," ",yt,"  RT  ",rt," ",long00," ",long2&
!       ,"  PS  ",(x2-x1)*(xt-x1)+(y2-y1)*(yt-y1)
!       if (i==1) points_control(1,1) = 0._pr
!       if (i==1) points_control(1,2) = 0._pr
!       points_control(1,1) = points_control(1,1)+xt
!       points_control(1,2) = points_control(1,2)+yt
!       enddo
!       points_courbe(1,1) = points_control(1,1)
!       points_courbe(1,2) = points_control(1,2)
 
!       tb = dtb
!       l = 1
!       long = 0._pr
!       do while ((tb<1._pr).and.(l+1<Ns+1))
!          l = l+1
!          call pointsBezierN(points_control,tb,px,py)
!          points_courbe(l,1) = px
!          points_courbe(l,2) = py
!          long = long + sqrt((points_courbe(l,1)-points_courbe(l-1,1))**2 + (points_courbe(l,2)-points_courbe(l-1,2))**2)
!          tb = tb+dtb
!       enddo
!       if (l==Ns-1) then
!          l=l+1
!          points_courbe(size(points_courbe_equal,1),1) = points_control(size(points_control,1),1)
!          points_courbe(size(points_courbe_equal,1),2) = points_control(size(points_control,1),2)
!          long = long + sqrt((points_courbe(l,1)-points_courbe(l-1,1))**2 + (points_courbe(l,2)-points_courbe(l-1,2))**2)
!       endif
!       if (.not.(l==Ns)) write(*,*) "WARNING2  p ",l
!       write(*,*) "IMPORTANTBIS  ",long00," ",long2," ",long,"  RT ",rt," ",sqrt((x1-xt)**2+(y1-yt)**2)&
!            ,"  coord ",x1," ",y1," ",x2," ",y2
 
 
!       !**
!       !! Uniform spline approximation 
!       !**
!       long2 = long
!       ds = long/(Ns-1)
!       points_courbe_equal =0._pr
!       points_courbe_equal(1,1) = points_control(1,1)
!       points_courbe_equal(1,2) = points_control(1,2)
!       write(*,*) "POINTS_COURBE_EQUAL  L==1  ",points_courbe_equal(1,1)," ",points_courbe_equal(1,2)
!       l = 1
!       long = 0._pr
!       tb = 0._pr
!       do while ((tb<1._pr).and.(l+1<Ns))
!          l = l+1
!          nt = 1
!          s = 0._pr
!          do while ((l-1)*ds-s>0._pr) 
!             nt = nt+1
!             s = s + sqrt((points_courbe(nt,1)-points_courbe(nt-1,1))**2 + (points_courbe(nt,2)-points_courbe(nt-1,2))**2)
!          enddo
!          tbm = (nt-2)*dtb
!          tbp = (nt-1)*dtb
!          tb = tbm
!          s0 = s
!          sinit = s - sqrt((points_courbe(nt,1)-points_courbe(nt-1,1))**2 + (points_courbe(nt,2)-points_courbe(nt-1,2))**2)
!          s = sinit
!          bool = 0
!          do while ((abs((l-1)*ds-s)/dx>eepsilon).and.(bool==0))
!             tb = (tbm + tbp)*0.5_pr
!             call pointsBezierN(points_control,tb,px,py)
!             s = sinit + sqrt((px-points_courbe(nt-1,1))**2 + (py-points_courbe(nt-1,2))**2)
!             if ((l-1)*ds-s>0._pr) then
!                tbm = tb
!             else
!                tbp = tb
!             endif
!             if (abs(tb-(tbm+tbp)*0.5_pr)<eepsilon) then
!                bool = 1
!             endif
!          enddo
!          call pointsBezierN(points_control,tb,px,py)
!          points_courbe_equal(l,1) = px
!          points_courbe_equal(l,2) = py
!          if (l==2) write(*,*) "POINTS_COURBE_EQUAL  L==2  ",points_courbe_equal(l,1)," ",points_courbe_equal(l,2)
!          long = long +&
!               sqrt((points_courbe_equal(l,1)-points_courbe_equal(l-1,1))**2 +&
!               (points_courbe_equal(l,2)-points_courbe_equal(l-1,2))**2)
!       enddo
!       lp = l
!       if (l==Ns-1) then
!          l = l+1
!          write(*,*) " ET VOILA j y passe  "
!          points_courbe_equal(l,1) = points_control(size(points_control,1),1)
!          points_courbe_equal(l,2) = points_control(size(points_control,1),2)
!          long = long +&
!               sqrt((points_courbe_equal(l,1)-points_courbe_equal(l-1,1))**2 +&
!               (points_courbe_equal(l,2)-points_courbe_equal(l-1,2))**2) 
!       endif
!       if (.not.(l==Ns)) write(*,*) "WARNING  s ",l
!       write(*,*) "CHECKPOINT 1 ",long2," ",long," ",l," ",tb," POINTS  ",points_courbe_equal(Ns,1)," ",points_courbe_equal(Ns,2)&
!            ," ",points_courbe(Ns,1)," ",points_courbe(Ns,2)," ",points_control(size(points_control,1),1)&
!            ," ",points_control(size(points_control,1),2)
!       do l=nint(0.04*Ns),1,-1
!          x2 = points_courbe_equal(Ns-l-1,1)
!          y2 = points_courbe_equal(Ns-l-1,2)
!          x1 = points_courbe_equal(Ns-l,1)
!          y1 = points_courbe_equal(Ns-l,2)
!          rt = dist(x1,y1,x2,y2)
!          if (abs(x1-x2)/dx<eepsilon) then
!             xt = x1
!             yt = y1+rt
 
!             if (((y2-y1)*(yt-y1))>0._pr) then
!               yt = y1-rt
!             endif
!          else
!             xt = rt/(sqrt(1+((y2-y1)**2)/((x2-x1)**2))) + x1
!             yt = (y2-y1)/(x2-x1)*(xt-x1) + y1
 
!             if (((x2-x1)*(xt-x1)+(y2-y1)*(yt-y1))>0._pr) then
!               xt = -rt/(sqrt(1+((y2-y1)**2)/((x2-x1)**2))) + x1
!               yt = (y2-y1)/(x2-x1)*(xt-x1) + y1
!             endif
!          endif
!          points_courbe_equal(size(points_courbe_equal,1)-l+1,1) = xt
!          points_courbe_equal(size(points_courbe_equal,1)-l+1,2) = yt
!       enddo
 
!       !**
!       !! Compute the bend amplitude 
!       !**
!       alphadef = 0._pr
!       call compute_thetadef(points_courbe_equal,alphadef,calpha,salpha,kt)
 
!       !**
!       !! Head and Tail linear approximations 
!       !**
!       dti = dsi/ds 
!       do l=2,Ni
!          tb = - (l-1)*dti
!          call courbeBezierN(points_courbe_equal,tb,Px,Py)
!          tail_courbe(l,1) = Px
!          tail_courbe(l,2) = Py
!       enddo
!       write(*,*) "POINTS_COURBE_EQUAL  L==1  ",points_courbe_equal(1,1)," ",points_courbe_equal(1,2)
!       tail_courbe(1,1) = points_courbe_equal(1,1) 
!       tail_courbe(1,2) = points_courbe_equal(1,2) 
!       dtf = dsf/ds 
!       do l=2,Nf
!          tb = 1._pr + (l-1)*dtf
!          call courbeBezierN(points_courbe_equal,tb,Px,Py)
!          head_courbe(l,1) = Px
!          head_courbe(l,2) = Py
!       enddo
!       head_courbe(1,1) = points_courbe_equal(Ns,1) 
!       head_courbe(1,2) = points_courbe_equal(Ns,2) 
!       long3 = 0._pr
!       do l=1,nl-1
!          long3 = long3 +sqrt((xx(midline(l+1,1))-xx(midline(l,1)))**2 + (yy(midline(l+1,2))-yy(midline(l,2)))**2)
!       enddo
 
!       write(*,*) "Longlong : ",long," ",long3," ",nl," ",sizeSkel
!       !**
!       !! Compute the bend amplitude 
!       !**
!       call compute_thetadef(head_courbe,alphadef,calpha,salpha,kt)
!       write(*,*) "ALPHA  ",alphadef*180/PI," ",cos(alphadef)," ",sin(alphadef)
 
 
!      l=1
!      indextab = 0
!      !     do theta=1,2
!      !        thetatab(theta) = 1
!  !!!        indextab(theta,thetatab(theta)) = 1
!      !        indextab(theta,1) = 1
!      !        slice(theta,1,1) = tail_courbe(Ni,1)
!      !        slice(theta,1,2) = tail_courbe(Ni,2)
!      !        slice(theta,1,3) = tail_courbe(Ni,3)
!      !        slice2(theta,1,1) = slice(theta,1,1)
!      !        slice2(theta,1,2) = slice(theta,1,2)
!      !        slice2(theta,1,3) = slice(theta,1,3)
!      !     enddo
!      do l=1,Ni-1
!         write(*,*) "intersection L : ",l
!         if (l==1) then
!  !           call intersection(tail_courbe(Ni-l+1,1),tail_courbe(Ni-l+1-1,1),tail_courbe(Ni-l+1,2),tail_courbe(Ni-l+1-1,2)&
!  !                ,tail_courbe(Ni-l+1,1),tail_courbe(Ni-l+1,2),distslice(l,1),distslice(l,2)&
!  !                ,distslice(Ns+Ni+Nf+l,1),distslice(Ns+Ni+Nf+l,2)&
!  !                ,xTheta(1),yTheta(1),xTheta(2),yTheta(2),errorl,errorr,deltal,deltar) 
!  !           write(*,*) "INITIALslice:  ",tail_courbe(Ni-l+1,1)," ",tail_courbe(Ni-l+1-1,1)&
!  !           ," ",distslice(l,1)," ",distslice(l,2)," ",yTheta(1)," ",deltal
!            call intersection(tail_courbe(Ni,1),tail_courbe(Ni-2,1),tail_courbe(Ni,2),tail_courbe(Ni-2,2)&
!                 ,tail_courbe(Ni-1,1),tail_courbe(Ni-1,2),distslice(2,1),distslice(2,2)&
!                 ,distslice(Ns+Ni+Nf+2,1),distslice(Ns+Ni+Nf+2,2)&
!                 ,xTheta(1),yTheta(1),xTheta(2),yTheta(2),errorl,errorr,deltal,deltar) 
!         else
!            call intersection(tail_courbe(Ni-l+1+1,1),tail_courbe(Ni-l+1-1,1),tail_courbe(Ni-l+1+1,2),tail_courbe(Ni-l+1-1,2)&
!                 ,tail_courbe(Ni-l+1,1),tail_courbe(Ni-l+1,2),distslice(l,1),distslice(l,2)&
!                 ,distslice(Ns+Ni+Nf+l,1),distslice(Ns+Ni+Nf+l,2)&
!                 ,xTheta(1),yTheta(1),xTheta(2),yTheta(2),errorl,errorr,deltal,deltar) 
!         endif
!         if (l==1) then
!           sinPhi = (xTheta(1)-tail_courbe(Ni-1,1))/(&
!                sqrt((xTheta(1)-tail_courbe(Ni-1,1))**2+(yTheta(1)-tail_courbe(Ni-1,2))**2))
!           cosPhi = (yTheta(1)-tail_courbe(Ni-1,2))/(&
!                sqrt((xTheta(1)-tail_courbe(Ni-1,1))**2+(yTheta(1)-tail_courbe(Ni-1,2))**2)) 
!         else
!           sinPhi = (xTheta(1)-tail_courbe(Ni-l+1,1))/(&
!                sqrt((xTheta(1)-tail_courbe(Ni-l+1,1))**2+(yTheta(1)-tail_courbe(Ni-l+1,2))**2))
!           cosPhi = (yTheta(1)-tail_courbe(Ni-l+1,2))/(&
!                sqrt((xTheta(1)-tail_courbe(Ni-l+1,1))**2+(yTheta(1)-tail_courbe(Ni-l+1,2))**2)) 
!         endif
!         do theta=1,2
!            !           slice(theta,l,1) = tail_courbe(Ni-l+1,1) + valDist(l,theta)*cos(valTheta(l,theta))*sinPhi
!            !           slice(theta,l,2) = tail_courbe(Ni-l+1,2) + valDist(l,theta)*cos(valTheta(l,theta))*cosPhi
!            slice(theta,l,1) = tail_courbe(Ni-l+1,1) + valDist(l,theta)*cosTheta_tab(l,theta)*sinPhi
!            slice(theta,l,2) = tail_courbe(Ni-l+1,2) + valDist(l,theta)*cosTheta_tab(l,theta)*cosPhi
!         enddo
 
!         do theta=1,2
!            if (l==1) then
!               thetatab(theta) = 1
!            else
!               thetatab(theta) = thetatab(theta) + 1
!            endif
!            !!           indextab(theta,thetatab(theta)) = l
!            indextab(theta,l) = 1
!            !!           indextab(theta,l) = indextab(theta,l) +1
!            slice2(theta,thetatab(theta),1) = slice(theta,l,1)
!            slice2(theta,thetatab(theta),2) = slice(theta,l,2)
 
!  !           if (l==1) write(*,*) "CHECKslice L==1 ",theta," ",slice(theta,l,1)," ",slice(theta,l,2)," ",cosPhi&
!  !           ," ",tail_courbe(Ni-l+1,2)," ",yTheta(1)," ",cosTheta_tab(l,theta)
!  !           if (l==2) write(*,*) "CHECKslice L==2 ",theta," ",slice(theta,l,1)," ",slice(theta,l,2)
!         enddo
!  !!!        slice(1,l,1) = xTheta(1)
!  !!!        slice(1,l,2) = yTheta(1)
!  !!!        if (errorl==1) then
!  !!!           slice(1,l,1) = slice(1,l-1,1)
!  !!!           slice(1,l,2) = slice(1,l-1,2)
!  !!!        endif
!  !!!        if ((abs(slice(1,l,1)-slice(1,l-1,1))<1.e-2).and.(abs(slice(1,l,2)-slice(1,l-1,2))<1.e-2)) then
!  !!!           slice(2,l,1) = slice(2,l-1,1)
!  !!!           slice(2,l,2) = slice(2,l-1,2)
!  !!!        else
!  !!!           slice(2,l,1) = -slice(1,l,1) + 2*tail_courbe(Ni-l+1,1)
!  !!!           slice(2,l,2) = -slice(1,l,2) + 2*tail_courbe(Ni-l+1,2)
!  !!!        endif
!  !!!        do theta=1,2
!  !!!            slice(theta,l,1) = slice(theta,l-1,1)
!  !!!            slice(theta,l,2) = slice(theta,l-1,2)
!  !!!            slice(theta,l,3) = slice(theta,l-1,3)
!  !!!        enddo
!         !        if (det(tail_courbe(Ni-l+1+1,1),tail_courbe(Ni-l+1+1,2),slice2(lpli,1,1),slice2(lpli,1,2),tail_courbe(Ni-l+1,1),tail_courbe(Ni-l+1,2))*det(tail_courbe(Ni-l+1+1,1),tail_courbe(Ni-l+1+1,2),slice2(lpli,1,1),slice2(lpli,1,2),xLeft,yLeft)>= 0._pr) then
!  !!!        do theta=1,2
!  !!!           thetatab(theta) = thetatab(theta) + 1
!  !!!        enddo
!  !!!        slice2(1,thetatab(1),1) = xTheta(1)
!  !!!        slice2(1,thetatab(1),2) = yTheta(1)
!  !!!        slice2(2,thetatab(2),1) = slice(2,l,1)
!  !!!        slice2(2,thetatab(2),2) = slice(2,l,2)
!  !!!        do theta=2,2
!  !!!           slice2(theta,thetatab(2),1) = slice(theta,l,1)
!  !!!           slice2(theta,thetatab(2),2) = slice(theta,l,2)
!  !!!           slice2(theta,thetatab(2),3) = slice(theta,l,3)
!  !!!        enddo
!         !        endif
!  !!!        if (errorl==1) then
!  !!!           write(*,*) "blabla 1",l
!  !!!           slice2(1,thetatab(1),1) = slice2(1,thetatab(1)-1,1)
!  !!!           slice2(1,thetatab(1),2) = slice2(1,thetatab(1)-1,2)
!  !!!           slice2(1,thetatab(1),3) = slice2(1,thetatab(1)-1,3)
!  !!!        endif
!  !!!        if (errorr==1) then
!  !!!           write(*,*) "blabla 2",l
!  !!!           slice2(2,thetatab(2),1) = slice2(2,thetatab(2)-1,1)
!  !!!           slice2(2,thetatab(2),2) = slice2(2,thetatab(2)-1,2)
!  !!!           slice2(2,thetatab(2),3) = slice2(2,thetatab(2)-1,3)
!  !!!        endif
!         if ((l>1).and.((errorl==1).or.(errorr==1))) then
!            do theta=1,2
!               slice(theta,l,1) = slice(theta,l-1,1)
!               slice(theta,l,2) = slice(theta,l-1,2)
!               slice2(theta,thetatab(theta),1) = slice2(theta,thetatab(theta)-1,1)
!               slice2(theta,thetatab(theta),2) = slice2(theta,thetatab(theta)-1,2)
!            enddo
!            write(*,*) "error initialsegment : ",theta," ",l," ",errorl," ",errorr
!         endif
 
 
!         !if (errorl==1) write(*,*) "WARNING : ",errorl," ",l
!      enddo
!      write(*,*) "VOILA VOILA"
!      call intersection(tail_courbe(2,1),points_courbe_equal(2,1),tail_courbe(2,2),points_courbe_equal(2,2),points_courbe_equal(1,1)&
!           !,points_courbe_equal(1,2),distslice(Ni+1,1),distslice(Ni+1,2),distslice(Ns+Ni+Nf+Ni+1,1),distslice(Ns+Ni+Nf+Ni+1,2)&
!           ,points_courbe_equal(1,2),distslice(Ni,1),distslice(Ni+1,2),distslice(Ns+Ni+Nf+Ni+1,1),distslice(Ns+Ni+Nf+Ni+1,2)&
!           ,xTheta(1),yTheta(1),xTheta(2),yTheta(2),errorl,errorr,deltal,deltar) 
!      l = Ni
!      sinPhi = (xTheta(1)-points_courbe_equal(1,1))/(&
!           sqrt((xTheta(1)-points_courbe_equal(1,1))**2+(yTheta(1)-points_courbe_equal(1,2))**2))
!      cosPhi = (yTheta(1)-points_courbe_equal(1,2))/(&
!           sqrt((xTheta(1)-points_courbe_equal(1,1))**2+(yTheta(1)-points_courbe_equal(1,2))**2)) 
!  !     if (yTheta(1)<1.e-6_pr) cosPhi = 1._pr
!  !!!     slice(1,l,1) = points_courbe_equal(1,1) + valDist(l,1)*cos(valTheta(l,1))*sinPhi
!  !!!     slice(1,l,2) = points_courbe_equal(1,2) + valDist(l,1)*cos(valTheta(l,1))*cosPhi
!  !!!     slice(2,l,1) = points_courbe_equal(1,1) + valDist(l,2)*cos(valTheta(l,2))*sinPhi 
!  !!!     slice(2,l,2) = points_courbe_equal(1,2) + valDist(l,2)*cos(valTheta(l,2))*cosPhi 
!      do theta=1,2
!         !        slice(theta,l,1) = points_courbe_equal(1,1) + valDist(l,theta)*cos(valTheta(l,theta))*sinPhi
!         !        slice(theta,l,2) = points_courbe_equal(1,2) + valDist(l,theta)*cos(valTheta(l,theta))*cosPhi
!         !        slice(theta,l,3) = zslice + valDist(l,theta)*sin(valTheta(l,theta))
!         slice(theta,l,1) = points_courbe_equal(1,1) + valDist(l,theta)*cosTheta_tab(l,theta)*sinPhi
!         slice(theta,l,2) = points_courbe_equal(1,2) + valDist(l,theta)*cosTheta_tab(l,theta)*cosPhi
!         !if (l==2) write(*,*) "CHECKslice L==2 ",theta," ",slice(theta,l,1)," ",slice(theta,l,2)," ",cosPhi&
!         write(*,*) "CHECKslice L==2 ",theta," ",slice(theta,l,1)," ",slice(theta,l,2)," ",cosPhi&
!         ," ",yTheta(1)," ",distslice(Ni+1,1)," ",distslice(Ni+1,2)&
!         ," ",distslice(Ns+Ni+Nf+Ni+1,1)," ",distslice(Ns+Ni+Nf+Ni+1,2)&
!         ," ",tail_courbe(2,1)," ",points_courbe_equal(2,1)&
!         ," ",tail_courbe(2,2)," ",points_courbe_equal(2,2)&
!         ," ",points_courbe_equal(1,1)," ",points_courbe_equal(1,2)
!      enddo
 
!  !!!     if ((det(tail_courbe(2,1),tail_courbe(2,2),slice2(1,thetatab(1),1),slice2(1,thetatab(1),2),points_courbe_equal(1,1),&
!  !!!points_courbe_equal(1,2))*det(tail_courbe(2,1),tail_courbe(2,2),slice2(1,thetatab(1),1),slice2(1,thetatab(1),2),slice(1,l,1),&
!  !!!slice(1,l,2))>= 0._pr).or.(thetatab(1)<Ni)) then
!  !!!        thetatab(1) = thetatab(1) + 1
!  !!!        slice2(1,thetatab(1),1) = slice(1,l,1)
!  !!!        slice2(1,thetatab(1),2) = slice(1,l,2)
!  !!!     endif
!      do theta=1,2
!         if ((det(tail_courbe(2,1),tail_courbe(2,2),slice2(theta,thetatab(theta),1),slice2(theta,thetatab(theta),2),&
!              points_courbe_equal(1,1),points_courbe_equal(1,2))*det(tail_courbe(2,1),tail_courbe(2,2)&
!              ,slice2(theta,thetatab(theta),1),&
!              slice2(theta,thetatab(theta),2),slice(theta,l,1),slice(theta,l,2))>= 0._pr).or.(thetatab(theta)<Ni)) then
!            thetatab(theta) = thetatab(theta) + 1
!            !!           indextab(theta,thetatab(theta)) = l
!            indextab(theta,l) = 1
!            slice2(theta,thetatab(theta),1) = slice(theta,l,1)
!            slice2(theta,thetatab(theta),2) = slice(theta,l,2)
!         endif
!      enddo
!      if (errorl==1) write(*,*) "WARNING : ",errorl," ",l
!      do l=2,Ns-1
!         sinPhir = sinPhi
!         cosPhir = cosPhi
!         call intersection(points_courbe_equal(l-1,1),points_courbe_equal(l+1,1),points_courbe_equal(l-1,2),&
!              points_courbe_equal(l+1,2),points_courbe_equal(l,1),points_courbe_equal(l,2),distslice(l+Ni,1),distslice(l+Ni,2),&
!              distslice(Ns+Ni+Nf+l+Ni,1),distslice(Ns+Ni+Nf+l+Ni,2),xTheta(1),yTheta(1)&
!              ,xTheta(2),yTheta(2),errorl,errorr,deltal,deltar) 
!         sinPhi = (xTheta(1)-points_courbe_equal(l,1))/(&
!              sqrt((xTheta(1)-points_courbe_equal(l,1))**2+(yTheta(1)-points_courbe_equal(l,2))**2))
!         cosPhi = (yTheta(1)-points_courbe_equal(l,2))/(&
!              sqrt((xTheta(1)-points_courbe_equal(l,1))**2+(yTheta(1)-points_courbe_equal(l,2))**2)) 
!         sinAlpha = sinPhi
!         cosAlpha = cosPhi
!         sinPhil = sinPhi
!         cosPhil = cosPhi
!         !!        slice(1,l+Ni-1,1) = points_courbe_equal(l,1) + valDist(l+Ni-1,1)*cos(valTheta(l+Ni-1,1))*sinPhi
!         !!        slice(1,l+Ni-1,2) = points_courbe_equal(l,2) + valDist(l+Ni-1,1)*cos(valTheta(l+Ni-1,1))*cosPhi
!         !!        slice(1,l+Ni-1,3) = zslice + valDist(l+Ni-1,1)*sin(valTheta(l+Ni-1,1))
!         !!        slice(2,l+Ni-1,1) = points_courbe_equal(l,1) + valDist(l+Ni-1,2)*cos(valTheta(l+Ni-1,2))*sinPhi 
!         !!        slice(2,l+Ni-1,2) = points_courbe_equal(l,2) + valDist(l+Ni-1,2)*cos(valTheta(l+Ni-1,2))*cosPhi 
!         do theta=1,2
!            !           slice(theta,l+Ni-1,1) = points_courbe_equal(l,1) + valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*sinPhi
!            !           slice(theta,l+Ni-1,2) = points_courbe_equal(l,2) + valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*cosPhi
!            slice(theta,l+Ni-1,1) = points_courbe_equal(l,1) + valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*sinPhi
!            slice(theta,l+Ni-1,2) = points_courbe_equal(l,2) + valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*cosPhi
!  !           if (l==2) write(*,*) "CHECKslice L==2 ",theta," ",slice(theta,l,1)," ",slice(theta,l,2)
!         enddo
!         !!        if ((det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2),slice2(1,thetatab(1),1),slice2(1,thetatab(1),2)&
!         !!,points_courbe_equal(l,1),points_courbe_equal(l,2))*det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2)&
!         !!,slice2(1,thetatab(1),1),slice2(1,thetatab(1),2),slice(1,l+Ni-1,1),slice(1,l+Ni-1,2))>= 0._pr).or.(thetatab(1)<Ni)) then
!         !!           thetatab(1) = thetatab(1) + 1
!         !!           slice2(1,thetatab(1),1) = slice(1,l+Ni-1,1)
!         !!           slice2(1,thetatab(1),2) = slice(1,l+Ni-1,2)
!         !!        endif
!         bool = 0
!         do theta=1,2
!            !           if (((theta==1).or.(theta==1+2/2)).and.(.not.((det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2)&
!            !,slice(theta,l+Ni-1-1,1)&
!            !,slice(theta,l+Ni-1-1,2),points_courbe_equal(l,1),points_courbe_equal(l,2))*det(points_courbe_equal(l-1,1)&
!            !,points_courbe_equal(l-1,2),slice(theta,l+Ni-1-1,1),slice(theta,l+Ni-1-1,2),slice(theta,l+Ni-1,1)&
!            !,slice(theta,l+Ni-1,2))>= 0._pr).or.(thetatab(theta)<Ni)))) then
!            if (.not.(det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2)&
!                                 !if (((theta==1).or.(theta==1+2/2)).and.(.not.(det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2)&
!                                 !,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*sinPhir&
!                                 !,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*cosPhir&
!                                 !,points_courbe_equal(l,1),points_courbe_equal(l,2))*det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2)&
!                                 !,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*sinPhir&
!                                 !,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*cosPhir&
!                                 !,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*sinPhil&
!                                 !,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*cosPhil)>= 0._pr)) then
!                 ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*sinPhir&
!                 ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*cosPhir&
!                 ,points_courbe_equal(l,1),points_courbe_equal(l,2))*det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2)&
!                 ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*sinPhir&
!                 ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*cosPhir&
!                 ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*sinPhil&
!                 ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*cosPhil)>= 0._pr)) then
 
!  !              !              sinTheta = -det(points_courbe_equal(l,1),points_courbe_equal(l,2),slice(theta,l+Ni-1,1),slice(theta,l+Ni-1,2)&
!  !              !,slice(theta,l+Ni-1-1,1),slice(theta,l+Ni-1-1,2))/(norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !              !,slice(theta,l+Ni-1,1),slice(theta,l+Ni-1,2))*norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !              !,slice(theta,l+Ni-1-1,1),slice(theta,l+Ni-1-1,2)))
!  !              !              cosTheta = dotProd(points_courbe_equal(l,1),points_courbe_equal(l,2),slice(theta,l+Ni-1,1),slice(theta,l+Ni-1,2)&
!  !              !,slice(theta,l+Ni-1-1,1),slice(theta,l+Ni-1-1,2))/(norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !              !,slice(theta,l+Ni-1,1),slice(theta,l+Ni-1,2))*norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !              !,slice(theta,l+Ni-1-1,1),slice(theta,l+Ni-1-1,2)))
!  !              !sinTheta = -det(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !              !     ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*sinPhil&
!  !              !     ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*cosPhil&
!  !              !     ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*sinPhir&
!  !              !     ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*cosPhir)&
!  !              !     /(norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !              !     ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*sinPhil&
!  !              !     ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*cosPhil)&
!  !              !     *norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !              !     ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*sinPhir&
!  !              !     ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*cosPhir))
!  !              !cosTheta = dotProd(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !              !     ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*sinPhil&
!  !              !     ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*cosPhil&
!  !              !     ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*sinPhir&
!  !              !     ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*cosPhir)&
!  !              !     /(norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !              !     ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*sinPhil&
!  !              !     ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*cosPhil)&
!  !              !     *norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !              !     ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*sinPhir&
!  !              !     ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*cosPhir))
!  !              sinTheta = -det(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !                   ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*sinPhil&
!  !                   ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*cosPhil&
!  !                   ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*sinPhir&
!  !                   ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*cosPhir)&
!  !                   /(norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !                   ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*sinPhil&
!  !                   ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*cosPhil)&
!  !                   *norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !                   ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*sinPhir&
!  !                   ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*cosPhir))
!  !              cosTheta = dotProd(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !                   ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*sinPhil&
!  !                   ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*cosPhil&
!  !                   ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*sinPhir&
!  !                   ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*cosPhir)&
!  !                   /(norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !                   ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*sinPhil&
!  !                   ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*cosPhil)&
!  !                   *norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!  !                   ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*sinPhir&
!  !                   ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*cosPhir))
!  !              if (cosTheta*oldC<0._pr) write(*,*) "ATTENTION COS NEGATIF ",l," ",l+Ni-1," ",theta,"     kt ",kt
!  !              !               if (bool==0) then
!  !              !               endif
!  !              if ((abs(sinTheta/cosTheta)<abs(oldS/oldC)).and.(bool==1)) then
!  !                 sinTheta = oldS
!  !                 cosTheta = oldC
!  !              endif
!  !              bool = 1
!  !              oldS = sinTheta
!  !              oldC = cosTheta
!               write(*,*) "cccTESTEST  ",theta," ",l," ",l+Ni-1," ",sinTheta," ",cosTheta," ",sinTheta/cosTheta&
!               ," ",errorr," ",errorl," ",deltar," ",deltal
!            endif
!         enddo
!         do theta=1,2
!            !           if ((det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2),slice2(theta,thetatab(theta),1)&
!            !,slice2(theta,thetatab(theta),2),points_courbe_equal(l,1),points_courbe_equal(l,2))*det(points_courbe_equal(l-1,1)&
!            !,points_courbe_equal(l-1,2),slice2(theta,thetatab(theta),1),slice2(theta,thetatab(theta),2),slice(theta,l+Ni-1,1)&
!            !,slice(theta,l+Ni-1,2))>= 0._pr).or.(thetatab(theta)<Ni)) then
!            !              xr = slice2(theta,thetatab(theta),1)
!            !              yr = slice2(theta,thetatab(theta),2)
!            xr = slice(theta,l+Ni-1,1)
!            yr = slice(theta,l+Ni-1,2)
!            thetatab(theta) = thetatab(theta) + 1
!            !!              indextab(theta,thetatab(theta)) = l+Ni-1
!            indextab(theta,l+Ni-1) = 1
!            slice2(theta,thetatab(theta),1) = slice(theta,l+Ni-1,1)
!            slice2(theta,thetatab(theta),2) = slice(theta,l+Ni-1,2)
!            if (bool==1) then
!               slice(theta,l+Ni-1,1) = cosTheta*(xr-points_courbe_equal(l,1)) + sinTheta*(yr-points_courbe_equal(l,2)) +&
!                    points_courbe_equal(l,1)
!               slice(theta,l+Ni-1,2) = -sinTheta*(xr-points_courbe_equal(l,1)) + cosTheta*(yr-points_courbe_equal(l,2)) +&
!                    points_courbe_equal(l,2)
!               sinPhi = sinAlpha*cosTheta + sinTheta*cosAlpha
!               cosPhi = cosAlpha*cosTheta - sinAlpha*sinTheta
!            endif
!            !           else
!            !
!            !           if (bool==1) then
!            !              xr = slice2(theta,thetatab(theta),1)
!            !              yr = slice2(theta,thetatab(theta),2)
!            !              thetatab(theta) = thetatab(theta) + 1
!            !              indextab(theta,l+Ni-1) = 1
!            !              slice2(theta,thetatab(theta),1) = cosTheta*xr + sinTheta*yr
!            !              slice2(theta,thetatab(theta),2) = -sinTheta*xr + cosTheta*yr
!            !            endif
!            !           endif
!         enddo
 
!         if (errorl==1) write(*,*) "WARNING11 : ",errorl," ",l
!      enddo
!      l = Ns
!      sinPhir = sinPhi
!      cosPhir = cosPhi
!      call intersection(points_courbe_equal(l-1,1),head_courbe(1+1,1),points_courbe_equal(Ns-1,2),head_courbe(1+1,2)&
!                                 !,head_courbe(1,1),head_courbe(1,2),distslice(Ni+Ns,1),distslice(1+Ni+Ns,2),distslice(Ns+Ni+Nf+Ni+Ns,1)&
!           !,head_courbe(1,1),head_courbe(1,2),distslice(Ni+Ns,1),distslice(Ni+Ns,2),distslice(Ns+Ni+Nf+Ni+Ns,1)&
!           ,head_courbe(1,1),head_courbe(1,2),distslice(Ni+Ns,1),distslice(Ni+Ns+1,2),distslice(Ns+Ni+Nf+Ni+Ns,1)&
!           ,distslice(Ns+Ni+Nf+1+Ni+Ns,2),xTheta(1),yTheta(1),xTheta(2),yTheta(2),errorl,errorr,deltal,deltar) 
!      sinPhi = (xTheta(1)-head_courbe(1,1))/(sqrt((xTheta(1)-head_courbe(1,1))**2+(yTheta(1)-head_courbe(1,2))**2))
!      cosPhi = (yTheta(1)-head_courbe(1,2))/(sqrt((xTheta(1)-head_courbe(1,1))**2+(yTheta(1)-head_courbe(1,2))**2)) 
!      sinAlpha = sinPhi
!      cosAlpha = cosPhi
!      sinPhil = sinPhi
!      cosPhil = cosPhi
!  !!!     slice(1,l+Ni-1,1) = head_courbe(1,1) + valDist(l+Ni-1,1)*cos(valTheta(l+Ni-1,1))*sinPhi
!  !!!     slice(1,l+Ni-1,2) = head_courbe(1,2) + valDist(l+Ni-1,1)*cos(valTheta(l+Ni-1,1))*cosPhi
!  !!!     slice(2,l+Ni-1,1) = head_courbe(1,1) + valDist(l+Ni-1,2)*cos(valTheta(l+Ni-1,2))*sinPhi
!  !!!     slice(2,l+Ni-1,2) = head_courbe(1,2) + valDist(l+Ni-1,2)*cos(valTheta(l+Ni-1,2))*cosPhi
!      do theta=1,2
!         !        slice(theta,l+Ni-1,1) = head_courbe(1,1) + valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*sinPhi
!         !        slice(theta,l+Ni-1,2) = head_courbe(1,2) + valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*cosPhi
!         slice(theta,l+Ni-1,1) = head_courbe(1,1) + valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*sinPhi
!         slice(theta,l+Ni-1,2) = head_courbe(1,2) + valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*cosPhi
!      enddo
 
!  !!!     if (det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2),slice2(1,thetatab(1),1),slice2(1,thetatab(1),2)&
!  !!!,points_courbe_equal(l,1),points_courbe_equal(l,2))*det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2)&
!  !!!,slice2(1,thetatab(1),1),slice2(1,thetatab(1),2),slice(1,l+Ni-1,1),slice(1,l+Ni-1,2))>= 0._pr) then
!  !!!        thetatab(1) = thetatab(1) + 1
!  !!!        slice2(1,thetatab(1),1) = slice(1,l+Ni-1,1)
!  !!!        slice2(1,thetatab(1),2) = slice(1,l+Ni-1,2)
!  !!!     endif
!      bool = 0
!      do theta=1,2
!         !        if (((theta==1).or.(theta==1+2/2)).and.(.not.(det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2)&
!         !,slice(theta,l+Ni-1-1,1)&
!         !,slice(theta,l+Ni-1-1,2),points_courbe_equal(l,1),points_courbe_equal(l,2))*det(points_courbe_equal(l-1,1)&
!         !,points_courbe_equal(l-1,2),slice(theta,l+Ni-1-1,1),slice(theta,l+Ni-1-1,2),slice(theta,l+Ni-1,1)&
!         !,slice(theta,l+Ni-1,2))>= 0._pr))) then
!         if (.not.(det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2)&
!                                 !,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*sinPhir&
!                                 !,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*cosPhir&
!                                 !,points_courbe_equal(l,1),points_courbe_equal(l,2))*det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2)&
!                                 !,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*sinPhir&
!                                 !,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*cosPhir&
!                                 !,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*sinPhil&
!                                 !,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*cosPhil)>= 0._pr)) then
!              ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*sinPhir&
!              ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*cosPhir&
!              ,points_courbe_equal(l,1),points_courbe_equal(l,2))*det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2)&
!              ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*sinPhir&
!              ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*cosPhir&
!              ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*sinPhil&
!              ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*cosPhil)>= 0._pr)) then
 
!            !              sinTheta = -det(points_courbe_equal(l,1),points_courbe_equal(l,2),slice(theta,l+Ni-1,1),slice(theta,l+Ni-1,2)&
!            !,slice(theta,l+Ni-1-1,1),slice(theta,l+Ni-1-1,2))/(norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!            !,slice(theta,l+Ni-1,1),slice(theta,l+Ni-1,2))*norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!            !,slice(theta,l+Ni-1-1,1),slice(theta,l+Ni-1-1,2)))
!            !              cosTheta = dotProd(points_courbe_equal(l,1),points_courbe_equal(l,2),slice(theta,l+Ni-1,1),slice(theta,l+Ni-1,2)&
!            !,slice(theta,l+Ni-1-1,1),slice(theta,l+Ni-1-1,2))/(norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!            !,slice(theta,l+Ni-1,1),slice(theta,l+Ni-1,2))*norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!            !,slice(theta,l+Ni-1-1,1),slice(theta,l+Ni-1-1,2)))
!            !sinTheta = -det(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!            !     ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*sinPhil&
!            !     ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*cosPhil&
!            !     ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*sinPhir&
!            !     ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*cosPhir)&
!            !     /(norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!            !     ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*sinPhil&
!            !     ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*cosPhil)&
!            !     *norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!            !     ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*sinPhir&
!            !     ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*cosPhir))
!            !cosTheta = dotProd(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!            !     ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*sinPhil&
!            !     ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*cosPhil&
!            !     ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*sinPhir&
!            !     ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*cosPhir)&
!            !     /(norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!            !     ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*sinPhil&
!            !     ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*cosPhil)&
!            !     *norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!            !     ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*sinPhir&
!            !     ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cos(valTheta(l+Ni-1-1,theta))*cosPhir))
!            sinTheta = -det(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!                 ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*sinPhil&
!                 ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*cosPhil&
!                 ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*sinPhir&
!                 ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*cosPhir)&
!                 /(norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!                 ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*sinPhil&
!                 ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*cosPhil)&
!                 *norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!                 ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*sinPhir&
!                 ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*cosPhir))
!            cosTheta = dotProd(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!                 ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*sinPhil&
!                 ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*cosPhil&
!                 ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*sinPhir&
!                 ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*cosPhir)&
!                 /(norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!                 ,points_courbe_equal(l,1)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*sinPhil&
!                 ,points_courbe_equal(l,2)+valDist(l+Ni-1,theta)*cosTheta_tab(l+Ni-1,theta)*cosPhil)&
!                 *norme(points_courbe_equal(l,1),points_courbe_equal(l,2)&
!                 ,points_courbe_equal(l-1,1)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*sinPhir&
!                 ,points_courbe_equal(l-1,2)+valDist(l+Ni-1-1,theta)*cosTheta_tab(l+Ni-1-1,theta)*cosPhir))
!            if (cosTheta*oldC<0._pr) write(*,*) "ATTENTION COS NEGATIF ",l," ",l+Ni-1," ",theta,"     kt ",kt
!            !               if (bool==0) then
!            !               endif
!            if ((abs(sinTheta/cosTheta)<abs(oldS/oldC)).and.(bool==1)) then
!               sinTheta = oldS
!               cosTheta = oldC
!            endif
!            bool = 1
!            oldS = sinTheta
!            oldC = cosTheta
!            write(*,*) "bbbTESTEST  ",theta," ",l," ",l+Ni-1," ",sinTheta," ",cosTheta," ",sinTheta/cosTheta
!         endif
!      enddo
!      !boolPhi = 0
!      do theta=1,2
!         !        if (det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2),slice2(theta,thetatab(theta),1)&
!         !,slice2(theta,thetatab(theta),2),points_courbe_equal(l,1),points_courbe_equal(l,2))*det(points_courbe_equal(l-1,1)&
!         !,points_courbe_equal(l-1,2),slice2(theta,thetatab(theta),1),slice2(theta,thetatab(theta),2),slice(theta,l+Ni-1,1)&
!         !,slice(theta,l+Ni-1,2))>= 0._pr) then
!         !              xr = slice2(theta,thetatab(theta),1)
!         !              yr = slice2(theta,thetatab(theta),2)
!         xr = slice(theta,l+Ni-1,1)
!         yr = slice(theta,l+Ni-1,2)
!         thetatab(theta) = thetatab(theta) + 1
!         !!           indextab(theta,thetatab(theta)) = l+Ni-1
!         indextab(theta,l+Ni-1) = 1
!         slice2(theta,thetatab(theta),1) = slice(theta,l+Ni-1,1)
!         slice2(theta,thetatab(theta),2) = slice(theta,l+Ni-1,2)
!         if (bool==1) then
!            slice(theta,l+Ni-1,1) = cosTheta*(xr-points_courbe_equal(l,1)) + sinTheta*(yr-points_courbe_equal(l,2)) +&
!                 points_courbe_equal(l,1)
!            slice(theta,l+Ni-1,2) = -sinTheta*(xr-points_courbe_equal(l,1)) + cosTheta*(yr-points_courbe_equal(l,2)) +&
!                 points_courbe_equal(l,2)
!            !              if (boolPhi==0) then
!            sinPhi = sinAlpha*cosTheta + sinTheta*cosAlpha
!            cosPhi = cosAlpha*cosTheta - sinAlpha*sinTheta
!            !              sinPhi = sinPhi*cosTheta + sinTheta*cosPhi
!            !              cosPhi = cosPhi*cosTheta - sinPhi*sinTheta
!            !                boolPhi = 1
!            !                endif
!         endif
!         !           else
!         !
!         !           if (bool==1) then
!         !              xr = slice2(theta,thetatab(theta),1)
!         !              yr = slice2(theta,thetatab(theta),2)
!         !              thetatab(theta) = thetatab(theta) + 1
!         !              indextab(theta,l+Ni-1) = 1
!         !              slice2(theta,thetatab(theta),1) = cosTheta*xr + sinTheta*yr
!         !              slice2(theta,thetatab(theta),2) = -sinTheta*xr + cosTheta*yr
!         !            endif
!         !           endif
!      enddo
!      if (errorl==1) write(*,*) "WARNING22 : ",errorl," ",l
!      do l=2,Nf
!         sinPhir = sinPhi
!         cosPhir = cosPhi
!         if (l==Nf) then
!  !           call intersection(head_courbe(l-1,1),head_courbe(l,1),head_courbe(l-1,2),head_courbe(l,2),head_courbe(l,1)&
!  !                ,head_courbe(l,2),distslice(l+Ni+Ns,1),distslice(l+Ni+Ns,2)&
!  !                ,distslice(Ns+Ni+Nf+l+Ni+Ns,1),distslice(Ns+Ni+Nf+l+Ni+Ns,2),xTheta(1)&
!  !                ,yTheta(1),xTheta(2),yTheta(2),errorl,errorr,deltal,deltar) 
!            call intersection(head_courbe(Nf-2,1),head_courbe(Nf-1,1),head_courbe(Nf-2,2),head_courbe(Nf-1,2),head_courbe(Nf-1,1)&
!                 ,head_courbe(Nf-1,2),distslice(Nf-1+Ni+Ns,1),distslice(Nf-1+Ni+Ns,2)&
!                 ,distslice(Ns+Ni+Nf+Nf-1+Ni+Ns,1),distslice(Ns+Ni+Nf+Nf-1+Ni+Ns,2),xTheta(1)&
!                 ,yTheta(1),xTheta(2),yTheta(2),errorl,errorr,deltal,deltar) 
!         else
!            call intersection(head_courbe(l-1,1),head_courbe(l+1,1),head_courbe(l-1,2),head_courbe(l+1,2),head_courbe(l,1)&
!                 ,head_courbe(l,2),distslice(l+Ni+Ns,1),distslice(l+Ni+Ns,2)&
!                 ,distslice(Ns+Ni+Nf+l+Ni+Ns,1),distslice(Ns+Ni+Nf+l+Ni+Ns,2),xTheta(1)&
!                 ,yTheta(1),xTheta(2),yTheta(2),errorl,errorr,deltal,deltar) 
!         endif
!         if (l==Nf) then
!           sinPhi = (xTheta(1)-head_courbe(Nf-1,1))/(sqrt((xTheta(1)-head_courbe(Nf-1,1))**2+(yTheta(1)-head_courbe(Nf-1,2))**2))
!           cosPhi = (yTheta(1)-head_courbe(Nf-1,2))/(sqrt((xTheta(1)-head_courbe(Nf-1,1))**2+(yTheta(1)-head_courbe(Nf-1,2))**2)) 
!         else
!           sinPhi = (xTheta(1)-head_courbe(l,1))/(sqrt((xTheta(1)-head_courbe(l,1))**2+(yTheta(1)-head_courbe(l,2))**2))
!           cosPhi = (yTheta(1)-head_courbe(l,2))/(sqrt((xTheta(1)-head_courbe(l,1))**2+(yTheta(1)-head_courbe(l,2))**2)) 
!         endif
!         sinAlpha = sinPhi
!         cosAlpha = cosPhi
!         sinPhil = sinPhi
!         cosPhil = cosPhi
!  !!!        slice(1,l+Ni-1+Ns-1,1) = head_courbe(l,1) + valDist(l+Ni-1+Ns-1,1)*cos(valTheta(l+Ni-1+Ns-1,1))*sinPhi
!  !!!        slice(1,l+Ni-1+Ns-1,2) = head_courbe(l,2) + valDist(l+Ni-1+Ns-1,1)*cos(valTheta(l+Ni-1+Ns-1,1))*cosPhi
!  !!!        slice(2,l+Ni-1+Ns-1,1) = head_courbe(l,1) + valDist(l+Ni-1+Ns-1,2)*cos(valTheta(l+Ni-1+Ns-1,2))*sinPhi
!  !!!        slice(2,l+Ni-1+Ns-1,2) = head_courbe(l,2) + valDist(l+Ni-1+Ns-1,2)*cos(valTheta(l+Ni-1+Ns-1,2))*cosPhi
!         do theta=1,2
!            !           slice(theta,l+Ni-1+Ns-1,1) = head_courbe(l,1) + valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*sinPhi
!            !           slice(theta,l+Ni-1+Ns-1,2) = head_courbe(l,2) + valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*cosPhi
!            slice(theta,l+Ni-1+Ns-1,1) = head_courbe(l,1) + valDist(l+Ni-1+Ns-1,theta)*cosTheta_tab(l+Ni-1+Ns-1,theta)*sinPhi
!            slice(theta,l+Ni-1+Ns-1,2) = head_courbe(l,2) + valDist(l+Ni-1+Ns-1,theta)*cosTheta_tab(l+Ni-1+Ns-1,theta)*cosPhi
!         enddo
!  !!!        if (det(head_courbe(l-1,1),head_courbe(l-1,2),slice2(1,thetatab(1),1),slice2(1,thetatab(1),2),head_courbe(l,1),&
!  !!!head_courbe(l,2))*det(head_courbe(l-1,1),head_courbe(l-1,2),slice2(1,thetatab(1),1),slice2(1,thetatab(1),2),slice(1,l+Ni-1+Ns-1,1),&
!  !!!slice(1,l+Ni-1+Ns-1,2))>= 0._pr) then
!  !!!           thetatab(1) = thetatab(1) + 1
!  !!!           slice2(1,thetatab(1),1) = slice(1,l+ni-1+Ns-1,1)
!  !!!           slice2(1,thetatab(1),2) = slice(1,l+Ni-1+Ns-1,2)
!  !!!        endif
!         bool = 0
!         do theta=1,2
!            !           if (.not.(det(head_courbe(l-1,1),head_courbe(l-1,2),slice2(theta,thetatab(theta),1),slice2(theta,thetatab(theta),2),&
!            !head_courbe(l,1),head_courbe(l,2))*det(head_courbe(l-1,1),head_courbe(l-1,2),slice2(theta,thetatab(theta),1),&
!            !slice2(theta,thetatab(theta),2),slice(theta,l+Ni-1+Ns-1,1),slice(theta,l+Ni-1+Ns-1,2))>= 0._pr)) then
!            !              bool = 1
!            !              sinTheta = det(head_courbe(l,1),head_courbe(l,2),slice(theta,l+Ni-1+Ns-1,1),slice(theta,l+Ni-1+Ns-1,2)&
!            !,slice2(theta,thetatab(theta),1),slice2(theta,thetatab(theta),2))/(norme(head_courbe(l,1),head_courbe(l,2)&
!            !,slice(theta,l+Ni-1+Ns-1,1),slice(theta,l+Ni-1+Ns-1,2))*norme(head_courbe(l,1),head_courbe(l,2)&
!            !,slice2(theta,thetatab(theta),1),slice2(theta,thetatab(theta),2)))
!            !cosTheta = dotProd(head_courbe(l,1),head_courbe(l,2),slice(theta,l+Ni-1+Ns-1,1),slice(theta,l+Ni-1+Ns-1,2)&
!            !,slice2(theta,thetatab(theta),1),slice2(theta,thetatab(theta),2))/(norme(head_courbe(l,1),head_courbe(l,2)&
!            !,slice(theta,l+Ni-1+Ns-1,1),slice(theta,l+Ni-1+Ns-1,2))*norme(head_courbe(l,1),head_courbe(l,2)&
!            !,slice2(theta,thetatab(theta),1),slice2(theta,thetatab(theta),2)))
 
!            !xr = dist3D(head_courbe(l-1,1),head_courbe(l-1,2),head_courbe(l-1,3),slice(theta,l+Ni-1+Ns-1-1,1),slice(theta,l+Ni-1+Ns-1-1,2)&
!            !,slice(theta,l+Ni-1+Ns-1-1,3))
!            !xl = dist3D(head_courbe(l,1),head_courbe(l,2),head_courbe(l,3),slice(theta,l+Ni-1+Ns-1,1),slice(theta,l+Ni-1+Ns-1,2)&
!            !,slice(theta,l+Ni-1+Ns-1,3))
!            !           if (theta==1+2/2) write(*,*) "mytheta     ",l,"    ",slice(theta,l+Ni-1+Ns-1,1)," ",slice(theta,l+Ni-1+Ns-1,2)&
!            !," ",head_courbe(l,1)+valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*sinPhil," ",head_courbe(l,2)+&
!            !valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*cosPhil," ",head_courbe(l,1)+xl*cos(valTheta(l+Ni-1+Ns-1,theta))&
!            !*sinPhil," ",head_courbe(l,2)+xl*cos(valTheta(l+Ni-1+Ns-1,theta))*cosPhil,"  dist  ",valDist(l+Ni-1+Ns-1,theta)," ",xl," "&
!            !,valTheta(l+Ni-1+Ns-1,theta)
!            !           if (theta==1+2/2) write(*,*) "mytheta     ",l,"    ",slice(theta,l+Ni-1+Ns-1-1,1)," ",slice(theta,l+Ni-1+Ns-1-1,2)&
!            !," ",head_courbe(l-1,1)+valDist(l+Ni-1+Ns-1-1,theta)*cos(valTheta(l+Ni-1+Ns-1-1,theta))*sinPhir," ",head_courbe(l-1,2)+&
!            !valDist(l+Ni-1+Ns-1-1,theta)*cos(valTheta(l+Ni-1+Ns-1-1,theta))*cosPhir," ",head_courbe(l-1,1)+&
!            !xr*cos(valTheta(l+Ni-1+Ns-1-1,theta))*sinPhir," ",head_courbe(l-1,2)+xr*cos(valTheta(l+Ni-1+Ns-1-1,theta))*cosPhir,"  dist  "&
!            !,valDist(l+Ni-1+Ns-1-1,theta)," ",xr," ",valTheta(l+Ni-1+Ns-1-1,theta)
 
!            !           if (((theta==1).or.(theta==1+2/2)).and.(.not.(det(head_courbe(l-1,1),head_courbe(l-1,2)&
!            if (.not.(det(head_courbe(l-1,1),head_courbe(l-1,2)&
!                                 !,head_courbe(l-1,1)+valDist(l+Ni-1+Ns-1-1,theta)*cos(valTheta(l+Ni-1+Ns-1-1,theta))*sinPhir&
!                                 !,head_courbe(l-1,2)+valDist(l+Ni-1+Ns-1-1,theta)*cos(valTheta(l+Ni-1+Ns-1-1,theta))*cosPhir&
!                                 !,head_courbe(l,1),head_courbe(l,2))*det(head_courbe(l-1,1),head_courbe(l-1,2)&
!                                 !,head_courbe(l-1,1)+valDist(l+Ni-1+Ns-1-1,theta)*cos(valTheta(l+Ni-1+Ns-1-1,theta))*sinPhir&
!                                 !,head_courbe(l-1,2)+valDist(l+Ni-1+Ns-1-1,theta)*cos(valTheta(l+Ni-1+Ns-1-1,theta))*cosPhir&
!                                 !,head_courbe(l,1)+valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*sinPhil&
!                                 !,head_courbe(l,2)+valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*cosPhil)&
!                 ,head_courbe(l-1,1)+valDist(l+Ni-1+Ns-1-1,theta)*cosTheta_tab(l+Ni-1+Ns-1-1,theta)*sinPhir&
!                 ,head_courbe(l-1,2)+valDist(l+Ni-1+Ns-1-1,theta)*cosTheta_tab(l+Ni-1+Ns-1-1,theta)*cosPhir&
!                 ,head_courbe(l,1),head_courbe(l,2))*det(head_courbe(l-1,1),head_courbe(l-1,2)&
!                 ,head_courbe(l-1,1)+valDist(l+Ni-1+Ns-1-1,theta)*cosTheta_tab(l+Ni-1+Ns-1-1,theta)*sinPhir&
!                 ,head_courbe(l-1,2)+valDist(l+Ni-1+Ns-1-1,theta)*cosTheta_tab(l+Ni-1+Ns-1-1,theta)*cosPhir&
!                 ,head_courbe(l,1)+valDist(l+Ni-1+Ns-1,theta)*cosTheta_tab(l+Ni-1+Ns-1,theta)*sinPhil&
!                 ,head_courbe(l,2)+valDist(l+Ni-1+Ns-1,theta)*cosTheta_tab(l+Ni-1+Ns-1,theta)*cosPhil)&
!                 >= 0._pr)) then
!               !           if (((theta==1).or.(theta==1+2/2)).and.(.not.(det(head_courbe(l-1,1),head_courbe(l-1,2)&
!               !,head_courbe(l-1,1)+xr*cos(valTheta(l+Ni-1+Ns-1-1,theta))*sinPhir,head_courbe(l-1,2)+xr*cos(valTheta(l+Ni-1+Ns-1-1,theta))*cosPhir,&
!               !head_courbe(l,1),head_courbe(l,2))*det(head_courbe(l-1,1),head_courbe(l-1,2),head_courbe(l-1,1)+&
!               !xr*cos(valTheta(l+Ni-1+Ns-1-1,theta))*sinPhir&
!               !,head_courbe(l-1,2)+xr*cos(valTheta(l+Ni-1+Ns-1-1,theta))*cosPhir,head_courbe(l,1)+xl*cos(valTheta(l+Ni-1+Ns-1,theta))*sinPhil&
!               !,head_courbe(l,2)+xl*cos(valTheta(l+Ni-1+Ns-1,theta))*cosPhil)&
!               !>= 0._pr))) then
!               !           if (((theta==1).or.(theta==1+2/2)).and.(.not.(det(head_courbe(l-1,1),head_courbe(l-1,2)&
!               !,slice(theta,l+Ni-1+Ns-1-1,1),slice(theta,l+Ni-1+Ns-1-1,2),&
!               !head_courbe(l,1),head_courbe(l,2))*det(head_courbe(l-1,1),head_courbe(l-1,2),slice(theta,l+Ni-1+Ns-1-1,1),&
!               !slice(theta,l+Ni-1+Ns-1-1,2),slice(theta,l+Ni-1+Ns-1,1),slice(theta,l+Ni-1+Ns-1,2))>= 0._pr))) then
!               !sinTheta = -det(head_courbe(l,1),head_courbe(l,2)&
!               !     ,head_courbe(l,1)+valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*sinPhil&
!               !     ,head_courbe(l,2)+valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*cosPhil&
!               !     ,head_courbe(l-1,1)+valDist(l+Ni-1+Ns-1-1,theta)*cos(valTheta(l+Ni-1+Ns-1-1,theta))*sinPhir&
!               !     ,head_courbe(l-1,2)+valDist(l+Ni-1+Ns-1-1,theta)*cos(valTheta(l+Ni-1+Ns-1-1,theta))*cosPhir)&
!               !     /(norme(head_courbe(l,1),head_courbe(l,2)&
!               !     ,head_courbe(l,1)+valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*sinPhil&
!               !     ,head_courbe(l,2)+valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*cosPhil)&
!               !     *norme(head_courbe(l,1),head_courbe(l,2)&
!               !     ,head_courbe(l-1,1)+valDist(l+Ni-1+Ns-1-1,theta)*cos(valTheta(l+Ni-1+Ns-1-1,theta))*sinPhir&
!               !     ,head_courbe(l-1,2)+valDist(l+Ni-1+Ns-1-1,theta)*cos(valTheta(l+Ni-1+Ns-1-1,theta))*cosPhir))
!               !cosTheta = dotProd(head_courbe(l,1),head_courbe(l,2)&
!               !     ,head_courbe(l,1)+valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*sinPhil&
!               !     ,head_courbe(l,2)+valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*cosPhil&
!               !     ,head_courbe(l-1,1)+valDist(l+Ni-1+Ns-1-1,theta)*cos(valTheta(l+Ni-1+Ns-1-1,theta))*sinPhir&
!               !     ,head_courbe(l-1,2)+valDist(l+Ni-1+Ns-1-1,theta)*cos(valTheta(l+Ni-1+Ns-1-1,theta))*cosPhir)&
!               !     /(norme(head_courbe(l,1),head_courbe(l,2)&
!               !     ,head_courbe(l,1)+valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*sinPhil&
!               !     ,head_courbe(l,2)+valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*cosPhil)&
!               !     *norme(head_courbe(l,1),head_courbe(l,2)&
!               !     ,head_courbe(l-1,1)+valDist(l+Ni-1+Ns-1-1,theta)*cos(valTheta(l+Ni-1+Ns-1-1,theta))*sinPhir&
!               !     ,head_courbe(l-1,2)+valDist(l+Ni-1+Ns-1-1,theta)*cos(valTheta(l+Ni-1+Ns-1-1,theta))*cosPhir))
!               sinTheta = -det(head_courbe(l,1),head_courbe(l,2)&
!                    ,head_courbe(l,1)+valDist(l+Ni-1+Ns-1,theta)*cosTheta_tab(l+Ni-1+Ns-1,theta)*sinPhil&
!                    ,head_courbe(l,2)+valDist(l+Ni-1+Ns-1,theta)*cosTheta_tab(l+Ni-1+Ns-1,theta)*cosPhil&
!                    ,head_courbe(l-1,1)+valDist(l+Ni-1+Ns-1-1,theta)*cosTheta_tab(l+Ni-1+Ns-1-1,theta)*sinPhir&
!                    ,head_courbe(l-1,2)+valDist(l+Ni-1+Ns-1-1,theta)*cosTheta_tab(l+Ni-1+Ns-1-1,theta)*cosPhir)&
!                    /(norme(head_courbe(l,1),head_courbe(l,2)&
!                    ,head_courbe(l,1)+valDist(l+Ni-1+Ns-1,theta)*cosTheta_tab(l+Ni-1+Ns-1,theta)*sinPhil&
!                    ,head_courbe(l,2)+valDist(l+Ni-1+Ns-1,theta)*cosTheta_tab(l+Ni-1+Ns-1,theta)*cosPhil)&
!                    *norme(head_courbe(l,1),head_courbe(l,2)&
!                    ,head_courbe(l-1,1)+valDist(l+Ni-1+Ns-1-1,theta)*cosTheta_tab(l+Ni-1+Ns-1-1,theta)*sinPhir&
!                    ,head_courbe(l-1,2)+valDist(l+Ni-1+Ns-1-1,theta)*cosTheta_tab(l+Ni-1+Ns-1-1,theta)*cosPhir))
!               cosTheta = dotProd(head_courbe(l,1),head_courbe(l,2)&
!                    ,head_courbe(l,1)+valDist(l+Ni-1+Ns-1,theta)*cosTheta_tab(l+Ni-1+Ns-1,theta)*sinPhil&
!                    ,head_courbe(l,2)+valDist(l+Ni-1+Ns-1,theta)*cosTheta_tab(l+Ni-1+Ns-1,theta)*cosPhil&
!                    ,head_courbe(l-1,1)+valDist(l+Ni-1+Ns-1-1,theta)*cosTheta_tab(l+Ni-1+Ns-1-1,theta)*sinPhir&
!                    ,head_courbe(l-1,2)+valDist(l+Ni-1+Ns-1-1,theta)*cosTheta_tab(l+Ni-1+Ns-1-1,theta)*cosPhir)&
!                    /(norme(head_courbe(l,1),head_courbe(l,2)&
!                    ,head_courbe(l,1)+valDist(l+Ni-1+Ns-1,theta)*cosTheta_tab(l+Ni-1+Ns-1,theta)*sinPhil&
!                    ,head_courbe(l,2)+valDist(l+Ni-1+Ns-1,theta)*cosTheta_tab(l+Ni-1+Ns-1,theta)*cosPhil)&
!                    *norme(head_courbe(l,1),head_courbe(l,2)&
!                    ,head_courbe(l-1,1)+valDist(l+Ni-1+Ns-1-1,theta)*cosTheta_tab(l+Ni-1+Ns-1-1,theta)*sinPhir&
!                    ,head_courbe(l-1,2)+valDist(l+Ni-1+Ns-1-1,theta)*cosTheta_tab(l+Ni-1+Ns-1-1,theta)*cosPhir))
!               if (cosTheta*oldC<0._pr) write(*,*) "ATTENTION COS NEGATIF ",l," ",l+Ni-1+Ns-1," ",theta,"    kt ",kt
!               !               if (bool==0) then
!               !               endif
!               if ((abs(sinTheta/cosTheta)<abs(oldS/oldC)).and.(bool==1)) then
!                  sinTheta = oldS
!                  cosTheta = oldC
!               endif
!               bool = 1
!               oldS = sinTheta
!               oldC = cosTheta
 
!               !              sinTheta = -det(head_courbe(l,1),head_courbe(l,2),slice(theta,l+Ni-1+Ns-1,1),slice(theta,l+Ni-1+Ns-1,2)&
!               !,slice(theta,l+Ni-1+Ns-1-1,1),slice(theta,l+Ni-1+Ns-1-1,2))/(norme(head_courbe(l,1),head_courbe(l,2)&
!               !,slice(theta,l+Ni-1+Ns-1,1),slice(theta,l+Ni-1+Ns-1,2))*norme(head_courbe(l,1),head_courbe(l,2)&
!               !,slice(theta,l+Ni-1+Ns-1-1,1),slice(theta,l+Ni-1+Ns-1-1,2)))
!               !              cosTheta = dotProd(head_courbe(l,1),head_courbe(l,2),slice(theta,l+Ni-1+Ns-1,1),slice(theta,l+Ni-1+Ns-1,2)&
!               !,slice(theta,l+Ni-1+Ns-1-1,1),slice(theta,l+Ni-1+Ns-1-1,2))/(norme(head_courbe(l,1),head_courbe(l,2)&
!               !,slice(theta,l+Ni-1+Ns-1,1),slice(theta,l+Ni-1+Ns-1,2))*norme(head_courbe(l,1),head_courbe(l,2)&
!               !,slice(theta,l+Ni-1+Ns-1-1,1),slice(theta,l+Ni-1+Ns-1-1,2)))
!               write(*,*) "aaaTESTEST  ",theta," ",l," ",l+Ni-1+Ns-1," ",sinTheta," ",cosTheta," ",sinPhir," ",cosPhir&
!                    ," ",sinPhil," ",cosPhil," ",sinTheta/cosTheta
!               !              write(*,*) "aaaTESTEST  ",theta," ",l," ",l+Ni-1+Ns-1," ",sinTh(theta)," ",cosTh(theta)," ",sinPhir," ",cosPhir&
!               !," ",sinPhil," ",cosPhil
!            endif
!         enddo
!         !boolPhi = 0
!         do theta=1,2
!            !           if (det(head_courbe(l-1,1),head_courbe(l-1,2),slice2(theta,thetatab(theta),1),slice2(theta,thetatab(theta),2),&
!            !head_courbe(l,1),head_courbe(l,2))*det(head_courbe(l-1,1),head_courbe(l-1,2),slice2(theta,thetatab(theta),1),&
!            !slice2(theta,thetatab(theta),2),slice(theta,l+Ni-1+Ns-1,1),slice(theta,l+Ni-1+Ns-1,2))>= 0._pr) then
!            !              xr = slice2(theta,thetatab(theta),1)
!            !              yr = slice2(theta,thetatab(theta),2)
!            xr = slice(theta,l+Ni-1+Ns-1,1)
!            yr = slice(theta,l+Ni-1+Ns-1,2)
!            !           if (theta==101) write(*,*) "CHECK ouch appartient  ",theta," ",l," ",l+Ni-1+Ns-1," ",slice(theta,l+Ni-1+Ns-1,1)&
!            !," ",slice(theta,l+Ni-1+Ns-1,2)," ",slice(theta,l+Ni-1+Ns-1,3)," ",cos(valTheta(l+Ni-1+Ns-1,theta))&
!            !," ",sin(valTheta(l+Ni-1+Ns-1,theta))," ",indextheta(l)
!            thetatab(theta) = thetatab(theta) + 1
!            !!              indextab(theta,thetatab(theta)) = l+Ni-1+Ns-1
!            indextab(theta,l+Ni-1+Ns-1) = 1
!            slice2(theta,thetatab(theta),1) = slice(theta,l+Ni-1+Ns-1,1)
!            slice2(theta,thetatab(theta),2) = slice(theta,l+Ni-1+Ns-1,2)
!            !              xr = slice2(theta,thetatab(theta),1)
!            !              yr = slice2(theta,thetatab(theta),2)
!            if (bool==1) then
!               !              write(*,*) "oint : ",slice2(theta,thetatab(theta),1)," ",slice2(theta,thetatab(theta),2)
!               !              slice2(theta,thetatab(theta),1) = cosTheta*xr + sinTheta*yr
!               !              slice2(theta,thetatab(theta),2) = -sinTheta*xr + cosTheta*yr
!               !write(*,*) "AH oui ",theta," ",l," ",l+Ni-1+Ns-1," ",slice(theta,l+Ni-1+Ns-1,1)," ",slice(theta,l+Ni-1+Ns-1,2)
!               slice(theta,l+Ni-1+Ns-1,1) = cosTheta*(xr-head_courbe(l,1)) + sinTheta*(yr-head_courbe(l,2)) + head_courbe(l,1)
!               slice(theta,l+Ni-1+Ns-1,2) = -sinTheta*(xr-head_courbe(l,1)) + cosTheta*(yr-head_courbe(l,2)) + head_courbe(l,2)
!               !write(*,*) "ah oui ",theta," ",l," ",l+Ni-1+Ns-1," ",slice(theta,l+Ni-1+Ns-1,1)," ",slice(theta,l+Ni-1+Ns-1,2)
!               !            slice(theta,l+Ni-1+Ns-1,1) = cosTh(theta)*(xr-head_courbe(l,1)) + sinTh(theta)*(yr-head_courbe(l,2)) + head_courbe(l,1)
!               !            slice(theta,l+Ni-1+Ns-1,2) = -sinTh(theta)*(xr-head_courbe(l,1)) + cosTh(theta)*(yr-head_courbe(l,2)) + head_courbe(l,2)
!               !              write(*,*) "oint : ",slice2(theta,thetatab(theta),1)," ",slice2(theta,thetatab(theta),2)
!               !              slice2(theta,thetatab(theta),1) = cosTheta*xr - sinTheta*yr
!               !              slice2(theta,thetatab(theta),2) = sinTheta*xr + cosTheta*yr
!               !              write(*,*) "oint : ",slice2(theta,thetatab(theta),1)," ",slice2(theta,thetatab(theta),2)
!               !if (boolPhi==0) then
!               sinPhi = sinAlpha*cosTheta + sinTheta*cosAlpha
!               cosPhi = cosAlpha*cosTheta - sinAlpha*sinTheta
!               !              sinPhi = sinPhi*cosTheta + sinTheta*cosPhi
!               !              cosPhi = cosPhi*cosTheta - sinPhi*sinTheta
!               !              sinPhi = sinAlpha*cosTh(theta) + sinTh(theta)*cosAlpha
!               !              cosPhi = cosAlpha*cosTh(theta) - sinAlpha*sinTh(theta)
!               !  boolPhi = 1
!               !  endif
!            endif
!            !           else
!            !
!            !           if (bool==1) then
!            !              xr = slice2(theta,thetatab(theta),1)
!            !              yr = slice2(theta,thetatab(theta),2)
!            !!              xr = slice(theta,l+Ni-1+Ns-1,1)
!            !!              yr = slice(theta,+Ni-1+Ns-1,2)
!            !              thetatab(theta) = thetatab(theta) + 1
!            !              indextab(theta,l+Ni-1+Ns-1) = 1
!            !              slice2(theta,thetatab(theta),1) = cosTheta*xr + sinTheta*yr
!            !              slice2(theta,thetatab(theta),2) = -sinTheta*xr + cosTheta*yr
!            !            endif
!            !           endif
!         enddo
!         if (errorl==1) write(*,*) "WARNING33 : ",errorl," ",l
!      enddo
!      !     l = Nf
!      !     do theta=1,2
!      !        thetatab(theta) = thetatab(theta) + 1
!  !!!        indextab(theta,thetatab(theta)) = Nf+Ni-1+Ns-1
!      !        indextab(theta,l+Ni-1+Ns-1) = 1
!      !     slice(theta,l+Ni-1+Ns-1,1) = head_courbe(l,1)
!      !     slice(theta,l+Ni-1+Ns-1,2) = head_courbe(l,2)
!      !     slice(theta,l+Ni-1+Ns-1,3) = head_courbe(l,3)
!      !        slice2(theta,thetatab(theta),1) = head_courbe(l,1)
!      !        slice2(theta,thetatab(theta),2) = head_courbe(l,2)
!      !        slice2(theta,thetatab(theta),3) = head_courbe(l,3)
!      !     enddo
!      !!   do theta=1,2
!      !!        indextheta(1) = 1
!  !!!!!        slice(theta,1,1) = tail_courbe(Ni,1)
!  !!!!!        slice(theta,1,2) = tail_courbe(Ni,2)
!  !!!!!        slice(theta,1,3) = tail_courbe(Ni,3)
!      !!        vect(theta,1,1) = slice(theta,1,1)
!      !!        vect(theta,1,2) = slice(theta,1,2)
!      !!        vect(theta,1,3) = slice(theta,1,3)
!      !!     do l=2,Ni-1
!  !!!!!        call intersection(tail_courbe(Ni-l+1+1,1),tail_courbe(Ni-l+1-1,1),tail_courbe(Ni-l+1+1,2),tail_courbe(Ni-l+1-1,2)&
!  !!!!!,tail_courbe(Ni-l+1,1),tail_courbe(Ni-l+1,2),distslice(l,1),distslice(l,2),distslice(Ns+Ni+Nf+l,1),distslice(Ns+Ni+Nf+l,2)&
!  !!!!!,xTheta(1),yTheta(1),xTheta(2),yTheta(2),errorl,errorr,deltal,deltar) 
!  !!!!!        sinPhi = (xTheta(1)-tail_courbe(Ni-l+1,1))/(&
!  !!!!!sqrt((xTheta(1)-tail_courbe(Ni-l+1,1))**2+(yTheta(1)-tail_courbe(Ni-l+1,2))**2))
!  !!!!!        cosPhi = (yTheta(1)-tail_courbe(Ni-l+1,2))/(&
!  !!!!!sqrt((xTheta(1)-tail_courbe(Ni-l+1,1))**2+(yTheta(1)-tail_courbe(Ni-l+1,2))**2)) 
!  !!!!!           slice(theta,l,1) = tail_courbe(Ni-l+1,1) + valDist(l,theta)*cos(valTheta(l,theta))*sinPhi
!  !!!!!           slice(theta,l,2) = tail_courbe(Ni-l+1,2) + valDist(l,theta)*cos(valTheta(l,theta))*cosPhi
!  !!!!!           slice(theta,l,3) = zslice + valDist(l,theta)*sin(valTheta(l,theta))
!      !!
!      !!           indextheta(l) = indextheta(l) + 1
!      !!           vect(indextheta(l),l,1) = slice(theta,l,1)
!      !!           vect(indextheta(l),l,2) = slice(theta,l,2)
!      !!           vect(indextheta(l),l,3) = slice(theta,l,3)
!      !!        if ((errorl==1).or.(errorr==1)) then
!      !!           vect(indextheta(l),l,1) = vect(indextheta(l),l-1,1)
!      !!           vect(indextheta(l),l,2) = vect(indextheta(l),l-1,2)
!      !!           vect(indextheta(l),l,3) = vect(indextheta(l),l-1,3)
!      !!   endif
!      !!     enddo
!      !!     l = Ni
!  !!!!!     call intersection(tail_courbe(2,1),points_courbe_equal(2,1),tail_courbe(2,2),points_courbe_equal(2,2),points_courbe_equal(1,1)&
!  !!!!!,points_courbe_equal(1,2),distslice(Ni,1),distslice(Ni+1,2),distslice(Ns+Ni+Nf+Ni,1),distslice(Ns+Ni+Nf+Ni+1,2),xTheta(1),yTheta(1)&
!  !!!!!,xTheta(2),yTheta(2),errorl,errorr,deltal,deltar) 
!  !!!!!     sinPhi = (xTheta(1)-points_courbe_equal(1,1))/(&
!  !!!!!sqrt((xTheta(1)-points_courbe_equal(1,1))**2+(yTheta(1)-points_courbe_equal(1,2))**2))
!  !!!!!     cosPhi = (yTheta(1)-points_courbe_equal(1,2))/(&
!  !!!!!sqrt((xTheta(1)-points_courbe_equal(1,1))**2+(yTheta(1)-points_courbe_equal(1,2))**2)) 
!  !!!!!        slice(theta,l,1) = points_courbe_equal(1,1) + valDist(l,theta)*cos(valTheta(l,theta))*sinPhi
!  !!!!!        slice(theta,l,2) = points_courbe_equal(1,2) + valDist(l,theta)*cos(valTheta(l,theta))*cosPhi
!  !!!!!        slice(theta,l,3) = zslice + valDist(l,theta)*sin(valTheta(l,theta))
!      !!        if ((det(tail_courbe(2,1),tail_courbe(2,2),vect(indextheta(l),l,1),vect(indextheta(l),l,2),&
!      !!points_courbe_equal(1,1),points_courbe_equal(1,2))*det(tail_courbe(2,1),tail_courbe(2,2),vect(indextheta(l),l,1),&
!      !!vect(indextheta(l),l,2),slice(theta,l,1),slice(theta,l,2))>= 0._pr).or.(l<Ni)) then
!      !!           indextheta(l) = indextheta(l) + 1
!      !!           vect(indextheta(l),l,1) = slice(theta,l,1)
!      !!           vect(indextheta(l),l,2) = slice(theta,l,2)
!      !!           vect(indextheta(l),l,3) = slice(theta,l,3)
!      !!        endif
!      !!     do l=2,Ns-1
!  !!!!!        call intersection(points_courbe_equal(l-1,1),points_courbe_equal(l+1,1),points_courbe_equal(l-1,2),&
!  !!!!!points_courbe_equal(l+1,2),points_courbe_equal(l,1),points_courbe_equal(l,2),distslice(l+Ni,1),distslice(l+Ni,2),&
!  !!!!!distslice(Ns+Ni+Nf+l+Ni,1),distslice(Ns+Ni+Nf+l+Ni,2),xTheta(1),yTheta(1),xTheta(2),yTheta(2),errorl,errorr,deltal,deltar) 
!  !!!!!        sinPhi = (xTheta(1)-points_courbe_equal(l,1))/(&
!  !!!!!sqrt((xTheta(1)-points_courbe_equal(l,1))**2+(yTheta(1)-points_courbe_equal(l,2))**2))
!  !!!!!        cosPhi = (yTheta(1)-points_courbe_equal(l,2))/(&
!  !!!!!sqrt((xTheta(1)-points_courbe_equal(l,1))**2+(yTheta(1)-points_courbe_equal(l,2))**2)) 
!  !!!!!           slice(theta,l+Ni-1,1) = points_courbe_equal(l,1) + valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*sinPhi
!  !!!!!           slice(theta,l+Ni-1,2) = points_courbe_equal(l,2) + valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*cosPhi
!  !!!!!           slice(theta,l+Ni-1,3) = zslice + valDist(l+Ni-1,theta)*sin(valTheta(l+Ni-1,theta))
!      !!           if ((det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2),vect(indextheta(l+Ni-1),l+Ni-1,1)&
!      !!,vect(indextheta(l+Ni-1),l+Ni-1,2),points_courbe_equal(l,1),points_courbe_equal(l,2))*det(points_courbe_equal(l-1,1)&
!      !!,points_courbe_equal(l-1,2),vect(indextheta(l+Ni-1),l+Ni-1,1),vect(indextheta(l+Ni-1),l+Ni-1,2),slice(theta,l+Ni-1,1)&
!      !!,slice(theta,l+Ni-1,2))>= 0._pr).or.(l<Ni)) then
!      !!              indextheta(l+Ni-1) = indextheta(l+Ni-1) + 1
!      !!              vect(indextheta(l+Ni-1),l+Ni-1,1) = slice(theta,l+Ni-1,1)
!      !!              vect(indextheta(l+Ni-1),l+Ni-1,2) = slice(theta,l+Ni-1,2)
!      !!              vect(indextheta(l+Ni-1),l+Ni-1,3) = slice(theta,l+Ni-1,3)
!      !!           endif
!      !!     enddo
!      !!     l = Ns
!  !!!!!     call intersection(points_courbe_equal(l-1,1),head_courbe(1+1,1),points_courbe_equal(Ns-1,2),head_courbe(1+1,2)&
!  !!!!!,head_courbe(1,1),head_courbe(1,2),distslice(Ni+Ns,1),distslice(1+Ni+Ns,2),distslice(Ns+Ni+Nf+Ni+Ns,1)&
!  !!!!!,distslice(Ns+Ni+Nf+1+Ni+Ns,2),xTheta(1),yTheta(1),xTheta(2),yTheta(2),errorl,errorr,deltal,deltar) 
!  !!!!!     sinPhi = (xTheta(1)-head_courbe(1,1))/(sqrt((xTheta(1)-head_courbe(1,1))**2+(yTheta(1)-head_courbe(1,2))**2))
!  !!!!!     cosPhi = (yTheta(1)-head_courbe(1,2))/(sqrt((xTheta(1)-head_courbe(1,1))**2+(yTheta(1)-head_courbe(1,2))**2)) 
!  !!!!!        slice(theta,l+Ni-1,1) = head_courbe(1,1) + valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*sinPhi
!  !!!!!        slice(theta,l+Ni-1,2) = head_courbe(1,2) + valDist(l+Ni-1,theta)*cos(valTheta(l+Ni-1,theta))*cosPhi
!  !!!!!        slice(theta,l+Ni-1,3) = zslice + valDist(l+Ni-1,theta)*sin(valTheta(l+Ni-1,theta))
!      !!        if (det(points_courbe_equal(l-1,1),points_courbe_equal(l-1,2),vect(indextheta(l+Ni-1),l+Ni-1,1)&
!      !!,vect(indextheta(l+Ni-1),l+Ni-1,2),points_courbe_equal(l,1),points_courbe_equal(l,2))*det(points_courbe_equal(l-1,1)&
!      !!,points_courbe_equal(l-1,2),vect(indextheta(l+Ni-1),l+Ni-1,1),vect(indextheta(l+Ni-1),l+Ni-1,2),slice(theta,l+Ni-1,1)&
!      !!,slice(theta,l+Ni-1,2))>= 0._pr) then
!      !!           indextheta(l+Ni-1) = indextheta(l+Ni-1) + 1
!      !!           vect(indextheta(l+Ni-1),l+Ni-1,1) = slice(theta,l+Ni-1,1)
!      !!           vect(indextheta(l+Ni-1),l+Ni-1,2) = slice(theta,l+Ni-1,2)
!      !!           vect(indextheta(l+Ni-1),l+Ni-1,3) = slice(theta,l+Ni-1,3)
!      !!        endif
!      !!     do l=2,Nf-1
!  !!!!!        call intersection(head_courbe(l-1,1),head_courbe(l+1,1),head_courbe(l-1,2),head_courbe(l+1,2),head_courbe(l,1)&
!  !!!!!,head_courbe(l,2),distslice(l+Ni+Ns,1),distslice(l+Ni+Ns,2),distslice(Ns+Ni+Nf+l+Ni+Ns,1),distslice(Ns+Ni+Nf+l+Ni+Ns,2),xTheta(1)&
!  !!!!!,yTheta(1),xTheta(2),yTheta(2),errorl,errorr,deltal,deltar) 
!  !!!!!        sinPhi = (xTheta(1)-head_courbe(l,1))/(sqrt((xTheta(1)-head_courbe(l,1))**2+(yTheta(1)-head_courbe(l,2))**2))
!  !!!!!        cosPhi = (yTheta(1)-head_courbe(l,2))/(sqrt((xTheta(1)-head_courbe(l,1))**2+(yTheta(1)-head_courbe(l,2))**2)) 
!  !!!!!           slice(theta,l+Ni-1+Ns-1,1) = head_courbe(l,1) + valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*sinPhi
!  !!!!!           slice(theta,l+Ni-1+Ns-1,2) = head_courbe(l,2) + valDist(l+Ni-1+Ns-1,theta)*cos(valTheta(l+Ni-1+Ns-1,theta))*cosPhi
!  !!!!!           slice(theta,l+Ni-1+Ns-1,3) = zslice + valDist(l+Ni-1+Ns-1,theta)*sin(valTheta(l+Ni-1+Ns-1,theta))
!      !!           if (det(head_courbe(l-1,1),head_courbe(l-1,2),vect(indextheta(l+Ni-1+Ns-1),l+Ni-1+Ns-1,1),vect(indextheta(l+Ni-1+Ns-1)&
!      !!,l+Ni-1+Ns-1,2),head_courbe(l,1),head_courbe(l,2))*det(head_courbe(l-1,1),head_courbe(l-1,2),vect(indextheta(l+Ni-1+Ns-1)&
!      !!,l+Ni-1+Ns-1,1),vect(indextheta(l+NI-1+Ns-1),l+Ni-1+Ns-1,2),slice(theta,l+Ni-1+Ns-1,1),slice(theta,l+Ni-1+Ns-1,2))>= 0._pr) then
!      !!              indextheta(l+Ni-1+Ns-1) = indextheta(l+Ni-1+Ns-1) + 1
!      !!              vect(indextheta(l+Ni-1+Ns-1),l+Ni-1+Ns-1,1) = slice(theta,l+Ni-1+Ns-1,1)
!      !!              vect(indextheta(l+Ni-1+Ns-1),l+Ni-1+Ns-1,2) = slice(theta,l+Ni-1+Ns-1,2)
!      !!              vect(indextheta(l+Ni-1+Ns-1),l+Ni-1+Ns-1,3) = slice(theta,l+Ni-1+Ns-1,3)
!      !!           endif
!      !!     enddo
!      !!     l = Nf
!      !!        indextheta(l+Ni-1+Ns-1) = indextheta(l+Ni-1+Ns-1) + 1
!      !!        vect(indextheta(l+Ni-1+Ns-1),l+Ni-1+Ns-1,1) = head_courbe(l,1)
!      !!        vect(indextheta(l+Ni-1+Ns-1),l+Ni-1+Ns-1,2) = head_courbe(l,2)
!      !!        vect(indextheta(l+Ni-1+Ns-1),l+Ni-1+Ns-1,3) = head_courbe(l,3)
!      !!   enddo
!      do l=1,Ni+Ns+Nf-2
!         indextheta(l) = 0
!         do theta = 1,2
!  !!!           if (l==indextab(theta,l)) then
!            !           if (theta==101) write(*,*) "CHECK ERROR101 ",theta," ",l," ",l+Ni-1+Ns-1," ",slice(theta,l+Ni-1+Ns-1,1)&
!            !," ",slice(theta,l+Ni-1+Ns-1,2)," ",slice(theta,l+Ni-1+Ns-1,3)," ",cos(valTheta(l+Ni-1+Ns-1,theta))&
!            !," ",sin(valTheta(l+Ni-1+Ns-1,theta))," ",indextheta(l)
!            !!           if (appartientElmt(indextab,theta,l).eqv..true.) then
!            if (indextab(theta,l)==1) then
!               indextheta(l) = indextheta(l) + 1
!  !!!              vect(l,indextheta(l)) = theta
!  !!!              vect(indextheta(l),l,1) = slice2(theta,thetatab(theta),1)
!  !!!              vect(indextheta(l),l,2) = slice2(theta,thetatab(theta),2)
!  !!!              vect(indextheta(l),l,3) = slice2(theta,thetatab(theta),3)
!               !           if (theta==101) write(*,*) "CHECK ERROR101 appartient  ",theta," ",l," ",l+Ni-1+Ns-1," ",slice(theta,l+Ni-1+Ns-1,1)&
!               !," ",slice(theta,l+Ni-1+Ns-1,2)," ",slice(theta,l+Ni-1+Ns-1,3)," ",cos(valTheta(l+Ni-1+Ns-1,theta))&
!               !," ",sin(valTheta(l+Ni-1+Ns-1,theta))," ",indextheta(l)
!               vect(indextheta(l),l,1) = slice(theta,l,1)
!               vect(indextheta(l),l,2) = slice(theta,l,2)
!               !              if (vect(101,l,2)>110.0) write(*,*) "WHHHYY  ",theta," ",l," ",indextheta(l)
!            endif
!         enddo
!      enddo
 
!      !     write(*,*) "lpl lpr finaux : ",thetatab(1)," ",thetatab(2)
!  !     do l=2,Ni+Ns+Nf-3
!  !        do theta = 1,2
!  !           slicetmp(theta,l,1) = 0._pr
!  !           slicetmp(theta,l,2) = 0._pr
!  !           if (((slice(theta,l,1)>slice(theta,l-1,1)).and.(slice(theta,l,1)<=slice(theta,l+1,1))).or.&
!  !           ((slice(theta,l,1)>slice(theta,l+1,1)).and.(slice(theta,l,1)<=slice(theta,l-1,1))))&
!  !           slicetmp(theta,l,1)=slice(theta,l,1)
!  !           if (((slice(theta,l+1,1)>slice(theta,l-1,1)).and.(slice(theta,l+1,1)<=slice(theta,l,1))).or.&
!  !           ((slice(theta,l+1,1)>slice(theta,l,1)).and.(slice(theta,l+1,1)<=slice(theta,l-1,1))))&
!  !           slicetmp(theta,l,1)=slice(theta,l+1,1)
!  !           if (((slice(theta,l-1,1)>slice(theta,l,1)).and.(slice(theta,l-1,1)<=slice(theta,l+1,1))).or.&
!  !           ((slice(theta,l-1,1)>slice(theta,l+1,1)).and.(slice(theta,l-1,1)<=slice(theta,l,1))))&
!  !           slicetmp(theta,l,1)=slice(theta,l-1,1)
!  !       
!  !           if (((slice(theta,l,2)>slice(theta,l-1,2)).and.(slice(theta,l,2)<=slice(theta,l+1,2))).or.&
!  !           ((slice(theta,l,2)>slice(theta,l+1,2)).and.(slice(theta,l,2)<=slice(theta,l-1,2))))&
!  !           slicetmp(theta,l,2)=slice(theta,l,2)
!  !           if (((slice(theta,l+1,2)>slice(theta,l-1,2)).and.(slice(theta,l+1,2)<=slice(theta,l,2))).or.&
!  !           ((slice(theta,l+1,2)>slice(theta,l,2)).and.(slice(theta,l+1,2)<=slice(theta,l-1,2))))&
!  !           slicetmp(theta,l,2)=slice(theta,l+1,2)
!  !           if (((slice(theta,l-1,2)>slice(theta,l,2)).and.(slice(theta,l-1,2)<=slice(theta,l+1,2))).or.&
!  !           ((slice(theta,l-1,2)>slice(theta,l+1,2)).and.(slice(theta,l-1,2)<=slice(theta,l,2))))&
!  !           slicetmp(theta,l,2)=slice(theta,l-1,2)
!  !        enddo
!  !     enddo
!  !     do l=2,Ni+Ns+Nf-3
!  !        do theta = 1,2
!  !          slice(theta,l,1)=slicetmp(theta,l,1)
!  !          slice(theta,l,2)=slicetmp(theta,l,2)
!  !        enddo
!  !     enddo
 
 
!        deallocate(points_control)
!! To HERE

       deallocate(midline)
 
      !**
      !! OUTPUTS 
      !**

      open(unit=85,file= trim(target_folder)//'/skelhead2D'//str(kt+picNum-1)//'.txt',status='unknown')

      open(unit=84,file= trim(target_folder)//'/right2D'//str(kt+picNum-1)//'.txt',status='unknown')

      open(unit=83,file= trim(target_folder)//'/left2D'//str(kt+picNum-1)//'.txt',status='unknown')

      open(unit=82,file= trim(target_folder)//'/skeletteq2D'//str(kt+picNum-1)//'.txt',status='unknown')

      open(unit=81,file= trim(target_folder)//'/skelett2D'//str(kt+picNum-1)//'.txt',status='unknown')

      open(unit=79,file= trim(target_folder)//'/skelett2D'//str(kt+picnum-1)//'.vtk',status='unknown')

      open(unit=86,file= trim(target_folder)//'/gradPhi2D'//str(kt+picnum-1)//'.vtk',status='unknown')

      open(unit=87,file= trim(target_folder)//'/gradPhiBis2D'//str(kt+picnum-1)//'.vtk',status='unknown')

      file_ids = (/ 79, 86, 87 /)
      do index = 1, size(file_ids)
         file_id = file_ids(index)
         call write_vtk_header(file_id, nx, ny, nz, dx, dy, dz, n_dim)
      enddo

      do k=1,ny
         do j=1,nx
            write(86,*) gradPhi(j,k)*(1-skel(j,k)) + 2*skel(j,k)
            write(79,*) gradPhi(j,k)*(1-skeltmp(j,k)) + 2*skeltmp(j,k)
            write(87,*) gradPhi(j,k)*(1-tmpbool(j,k)) + 2*tmpbool(j,k)
         enddo
      enddo
      write(*,*) "CHECK size 2 : ",sum(thetatab)," ",sum(indextheta),"  THETA  ",thetatab(2/2+1)

      open(unit=78,file= trim(target_folder)//'/surf2D'//str(kt+picNum-1)//'.vts',status='unknown')

           write(78,'(a)') "<?xml version=""1.0""?>"
           write(78,'(a)') "<VTKFile type=""StructuredGrid"" version=""0.1"" byte_order=""LittleEndian"" compressor=&
           &     ""vtkZLibDataCompressor"">"
           write(78,'(a,I3,a,I3,a)') "<StructuredGrid WholeExtent=""0 ",1," 0 ",Ni-1+Ns+Nf-1-1," 0 0"">"
           write(78,'(a,I3,a,I3,a)') "<Piece Extent=""0 ",1," 0 ",Ni-1+Ns+Nf-1-1," 0 0"">"
           write(78,'(a)') "<PointData >"
           write(78,'(a)') "</PointData>"
           write(78,'(a)') "<CellData>"
           write(78,'(a)') "</CellData>"
           write(78,'(a)') "<Points>"
           write(78,'(a)') "<DataArray NumberOfComponents=""3"" type=""Float64"" format=""ascii"" >"  

      do l=1,Ni+Ns+Nf-2
         do theta=1,2

         enddo
      enddo

      theta=1
      do l=1,Ni+Ns+Nf-2
         if (indextab(theta,l)==1) write(78,*) slice(theta,l,1)," ",slice(theta,l,2)," ",1.
      enddo
      theta=2
      do l=Ni+Ns+Nf-2,1,-1
         if (indextab(theta,l)==1) write(78,*) slice(theta,l,1)," ",slice(theta,l,2)," ",1.
      enddo
           write(78,'(a)') "</DataArray>"
           write(78,'(a)') "</Points>"
           write(78,'(a)') "</Piece>"
           write(78,'(a)') "</StructuredGrid>"
           write(78,'(a)') "</VTKFile>" 
           do l=1,Ni-1
              write(82,*) tail_courbe(Ni-l+1,1)," ",tail_courbe(Ni-l+1,2)
           enddo
           do l=2,Nf
              write(85,*) head_courbe(l,1)," ",head_courbe(l,2)
           enddo

           write(80,*) kt,"    ",nl," ",long3," ",long2," ",long
           do l=1,Ns
              write(81,*) points_courbe(l,1)," ",points_courbe(l,2)
              write(82,*) points_courbe_equal(l,1)," ",points_courbe_equal(l,2)
           enddo

           do l=1,Ni+Ns+Nf-2 
              if (indextab(1+2/2,l)==1) write(84,*) slice(1+2/2,l,1)," ",slice(1+2/2,l,2)
              if (indextab(1,l)==1) write(83,*) slice(1,l,1)," ",slice(1,l,2)

           enddo
      close(79)    
      close(86)    
      close(87)    
      close(78)    
      close(81)
      close(82)
      close(83)
      close(84)
      close(85)
   
   enddo
   close(80)
   write(*,*) "************************************************" 
   deallocate(rhoSlices)
   deallocate(xx,yy)

   deallocate(tmp1,tmp2,tmpbool)
   deallocate(dir1,dir2,dir3,dir4,Nseed,skel,skel2,skel3,skeltmp)
   deallocate(distslice)
   deallocate(points_courbe,points_courbe_equal,points_courbe_equal_ref)
   deallocate(slice_courbe,slice_courbe_equal)
   deallocate(tail_courbe,head_courbe)
   deallocate(slice,slice2,thetatab,xTheta,yTheta,valTheta,valDist,slicetmp)
   deallocate(longslice,longslicei,longslicef,sslice,dsslice,slicemid)
   deallocate(indextab,indextheta,vect)
   deallocate(gradPhi) 
   deallocate(cosTheta_tab,sinTheta_tab)

   call cpu_time(end_time)
   write(*,*) "Full movie iterations time = ", end_time-start_time, " sec"
   write(*,*) "Full movie iterations time = ", (end_time-start_time)/n_frames, " sec per iteration"

end program def3D
 