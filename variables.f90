module variables
  use doubler  
  implicit none
 
  integer :: N
  integer :: nx, ny, nz
  real(pr) :: dx, dy, dz, dt 
  real(pr) :: a, bbb !!,LL
  real(pr), dimension(4) :: beta, c, alp
  real(pr) :: eepsilon
  real(pr), dimension(:), allocatable :: xx, yy, zz 

end module variables
