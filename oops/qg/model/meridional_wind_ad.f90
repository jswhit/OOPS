! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Calculate meridional wind component - Adjoint

subroutine meridional_wind_ad (v,x,nx,ny,deltax)

use kinds

implicit none
integer, intent(in) :: nx         !< Zonal grid dimension
integer, intent(in) :: ny         !< Meridional grid dimension
real(kind=kind_real), intent(inout) :: v(nx,ny,2) !< Meridional wind adjoint variable
real(kind=kind_real), intent(inout) :: x(nx,ny,2) !< Streamfunction adjoint variable
real(kind=kind_real), intent(in) :: deltax        !< Zonal grid spacing (non-dimensional)

x(nx    ,:,:) = x(nx    ,:,:) - (0.5_kind_real/deltax)*v(1   ,:,:)
x(1:nx-1,:,:) = x(1:nx-1,:,:) - (0.5_kind_real/deltax)*v(2:nx,:,:)
x(1     ,:,:) = x(1   ,:,:) + (0.5_kind_real/deltax)*v(nx    ,:,:) 
x(2:nx  ,:,:) = x(2:nx,:,:) + (0.5_kind_real/deltax)*v(1:nx-1,:,:)
v(:,:,:) = 0.0_kind_real

end subroutine meridional_wind_ad
