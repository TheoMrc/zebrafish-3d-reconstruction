module midline2D_modules
    use doubler
    use AdvectionProblem

    implicit none

    contains
        subroutine flatten3Dshape(shape_file_path, rhoSlices, rhoSlices3D, gradPhi)
            implicit none 
            character(len=*), intent(in) :: shape_file_path
            real(pr), dimension(:,:,:), intent(inout) :: rhoSlices3D
            real(pr), dimension(:,:), intent(inout) :: rhoSlices, gradPhi
            integer :: i, j, k, nx, ny, nz

            nx = size(rhoSlices3D, 1)
            ny = size(rhoSlices3D, 2)
            nz = size(rhoSlices3D, 3)
            open(unit=78,file=shape_file_path, status='unknown')
            
            do k=1,nx
                do j=1,ny
                    do i=1,nz
                        read(78,*) rhoSlices3D(i,k,j)
                    enddo
                enddo
            enddo
            close(78)
            do k=1,nx
                rhoSlices3D(:, k, :) = rhoSlices3D(:, k, :)/maxval(rhoSlices3D(:,k,:))
            enddo
            
            rhoSlices3D = 10*rhoSlices3D
            do i=1,nz
                do k=1,nx
                    do j=1,ny
                        if (rhoSlices3D(i,k,j)>1._pr) then
                            rhoSlices3D(i,k,j) = 1._pr
                        endif
                    enddo
                enddo
            enddo
            rhoSlices3D = rhoSlices3D - 0.5_pr*(maxval(rhoSlices3D)+minval(rhoSlices3D))
            dy = 2.4
            dx = 2.25171
            dz = 2.25171
            call updateDistance3D(rhoSlices3D,gradPhi,150._pr)
            where(rhoSlices3D<0._pr) rhoSlices3D = 0._pr
            rhoSlices(:,:) = sum(rhoSlices3D(:,:,:),1)
            where(rhoSlices>0._pr) rhoSlices = 1._pr            
            rhoSlices = rhoSlices - 0.5_pr*(maxval(rhoSlices) + minval(rhoSlices))

        end subroutine flatten3Dshape


        subroutine write_vtk_header(file_id, nx, ny, nz, dx, dy, dz, n_dim)
            implicit none
            integer :: nx, ny, nz,  n_dim, nz_bis
            real(pr) :: dx, dy, dz, dz_bis
            integer, intent(in) :: file_id
            
            if (n_dim == 2) then
                nz_bis = 1
                dz_bis = 1
            elseif (n_dim == 3) then
                nz_bis = nz
                dz_bis = dz
            endif

            write(file_id,'(1A26)') '# vtk DataFile Version 2.0'
            write(file_id,'(a)') 'rho'
            write(file_id,'(a)') 'ASCII'
            write(file_id,'(a)') 'DATASET STRUCTURED_POINTS'
            write(file_id,'(a,I4,I4,I4)') 'DIMENSIONS ', nx, ny, nz_bis
            write(file_id,'(a,E23.15,E23.15,E23.15)') 'ORIGIN', 1., 1., 1.
            ! write(file_id,'(a,E23.15,E23.15,E23.15)') 'SPACING', dx, dy, dz_bis
            write(file_id,'(a,E23.15,E23.15,E23.15)') 'SPACING', 1., 1., 1.

            write(file_id,'(a,I9)') 'POINT_DATA' , nx * ny * nz_bis
            write(file_id,'(a)') 'SCALARS values double'
            write(file_id,'(a)') 'LOOKUP_TABLE default'

        end subroutine write_vtk_header


        subroutine create_directory(newDirPath)       
            implicit none
        
            character(len=*), intent(in) :: newDirPath
            character(len=256)           :: mkdirCmd
            logical                      :: dirExists
        
            ! Check if the directory exists first
            inquire(file=trim(newDirPath)//'/.', exist=dirExists)        
        
            if (dirExists) then
                write (*,*) "Directory already exists: '"//trim(newDirPath)//"'"
            else
                mkdirCmd = 'mkdir -p '//trim(newDirPath)
                write(*,'(a)') "Creating new directory: '"//trim(mkdirCmd)//"'"
                call system(mkdirCmd)
            endif
        end subroutine create_directory


end module midline2D_modules