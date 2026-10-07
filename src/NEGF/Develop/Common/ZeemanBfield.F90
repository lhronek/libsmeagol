module mBfield

  use mNegfOutput, only: negf_master, negf_out_unit
use mConstants
use mTypes

implicit none

private

real(kdp), public :: ZeemanBx, ZeemanBy, ZeemanBz

public AddZeemanBfield

contains

subroutine AddZeemanBfield(v,n,nspin)

integer, intent(in) :: n,nspin
real, intent(inout) :: v(n,nspin) ! the potential is a real valued quantity, not double precision

if(negf_master().and.(ZeemanBx .ne. 0.0_kdp.or.ZeemanBy .ne. 0.0_kdp.or.ZeemanBz .ne. 0.0_kdp)) write(negf_out_unit(),*)"ZeemanBxyz=",ZeemanBx,ZeemanBy,ZeemanBz

if(ZeemanBz .ne. 0.0_kdp)then
!    write(*,*)"ZeemanBz=",ZeemanBz
  v(:,1)=v(:,1)-0.5_kdp * ZeemanBz
  v(:,2)=v(:,2)+0.5_kdp * ZeemanBz
endif

if(nspin>2)then
  if(ZeemanBx .ne. 0.0_kdp)then
!    write(*,*)"ZeemanBx=",ZeemanBx
    v(:,3)=v(:,3)+0.5_kdp * ZeemanBx
  endif
  if(ZeemanBy .ne. 0.0_kdp)then
!    write(*,*)"ZeemanBx=",ZeemanBx
    v(:,4)=v(:,4)+0.5_kdp * ZeemanBy !xxx (check sign)
  endif
endif

end subroutine AddZeemanBfield

end module mBfield
