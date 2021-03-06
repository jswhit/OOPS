! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Multiply streamfunction by sqrt(C) - Adjoint

subroutine qg_3d_covar_sqrt_mult_ad(kx,ky,xincr,xctl,config)

use iso_c_binding
use qg_fields
use qg_3d_covar_configs
use kinds

implicit none
integer(c_int), intent(in)    :: kx            !< Zonal grid spacing
integer(c_int), intent(in)    :: ky            !< Meridional grid spacing
real(c_double), intent(inout) :: xctl(kx,ky,2) !< Result
type(qg_field), intent(in)    :: xincr        !< Streamfunction: psi
type(qg_3d_covar_config), intent(in) :: config !< covar config structure

integer :: ifax(13)
real(kind=kind_real), allocatable :: trigs(:), xext(:,:,:), work(:), xout_copy(:,:,:)
integer :: i, j, k, iri, m
real(kind=kind_real) :: zc, zero, one

! ------------------------------------------------------------------------------

!--- adjoint of multiplication by standard deviation

allocate(xout_copy(kx,ky,2))

do k=1,2
  do j=1,ky
    do i=1,kx
      xout_copy(i,j,k) = xincr%x(i,j,k) * config%sigma
    enddo
  enddo
enddo

!--- adjoint of multiplication by square-root of zonal correlation matrix

allocate (trigs(3*kx/2+1))
call set99 (trigs,ifax,kx)

allocate(xext(kx+2,ky,2))
xext(1:kx,:,:) = xout_copy(:,:,:)
xext(kx+1:kx+2,:,:) = 0.0_kind_real

allocate(work(2*ky*(kx+2)))

call fft991 (xext,work,trigs,ifax,1,kx+2,kx,2*ky,-1)

do k=1,2
  do j=1,ky
    do m=0,kx/2
      do iri=1,2
        xext(2*m+iri,j,k) = xext(2*m+iri,j,k) * config%sqrt_zonal(m)
      enddo
    enddo
  enddo
enddo

call fft991 (xext,work,trigs,ifax,1,kx+2,kx,2*ky,+1)
xout_copy(:,:,:) = xext(1:kx,:,:)

deallocate(work)
deallocate(xext)
deallocate(trigs)

!--- adjoint of multiplication by square-root of meridional correlation matrix

zero = 0.0_kind_real
one = 1.0_kind_real
allocate(work(ky))
work(:)=0.0_kind_real
do k=1,2
  do i=1,kx
    call DSYMV('L',ky,ONE,CONFIg%sqrt_merid,ky,xout_copy(i,1,k),kx,zero,work,1)
    do j=1,ky
      xout_copy(i,j,k) = work(j)
    enddo
  enddo
enddo
deallocate(work)

!--- adjoint multiplication by symmetric square-root of vert correlation matrix

zc = sqrt(1.0_kind_real-config%vert_corr*config%vert_corr)
do j=1,ky
  do i=1,kx
    xctl(i,j,1) = xctl(i,j,1) + xout_copy(i,j,1) &
              & +config%vert_corr * xout_copy(i,j,2)
    xctl(i,j,2) = xctl(i,j,2) + xout_copy(i,j,2) * zc
  enddo
enddo

deallocate(xout_copy)

end subroutine qg_3d_covar_sqrt_mult_ad
