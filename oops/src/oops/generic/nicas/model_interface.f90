!----------------------------------------------------------------------
! Module: model_interface.f90
!> Purpose: model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module model_interface

use model_aro, only: model_aro_coord,model_aro_read,model_aro_write
use model_arp, only: model_arp_coord,model_arp_read,model_arp_write
use model_gem, only: model_gem_coord,model_gem_read,model_gem_write
use model_geos, only: model_geos_coord,model_geos_read,model_geos_write
use model_gfs, only: model_gfs_coord,model_gfs_read,model_gfs_write
use model_ifs, only: model_ifs_coord,model_ifs_read,model_ifs_write
use model_mpas, only: model_mpas_coord,model_mpas_read,model_mpas_write
use model_nemo, only: model_nemo_coord,model_nemo_read,model_nemo_write
use model_oops, only: model_oops_write
use model_wrf, only: model_wrf_coord,model_wrf_read,model_wrf_write
use module_namelist, only: nam
use netcdf
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msvalr,msr
use tools_nc, only: ncfloat,ncerr
use type_mpl, only: mpl,mpl_bcast
use type_sdata, only: sdatatype
implicit none

private
public :: model_coord,model_read,model_write

contains

!----------------------------------------------------------------------
! Subroutine: model_coord
!> Purpose: get coordinates
!----------------------------------------------------------------------
subroutine model_coord(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! TODO: change that one day
sdata%nl0 = nam%nl

! Select model
if (trim(nam%model)=='aro') then
   call model_aro_coord(sdata)
elseif (trim(nam%model)=='arp') then
   call model_arp_coord(sdata)
elseif (trim(nam%model)=='gem') then
   call model_gem_coord(sdata)
elseif (trim(nam%model)=='geos') then
   call model_geos_coord(sdata)
elseif (trim(nam%model)=='gfs') then
   call model_gfs_coord(sdata)
elseif (trim(nam%model)=='ifs') then
   call model_ifs_coord(sdata)
elseif (trim(nam%model)=='mpas') then
   call model_mpas_coord(sdata)
elseif (trim(nam%model)=='nemo') then
   call model_nemo_coord(sdata)
elseif (trim(nam%model)=='oops') then
   ! Specific interface, see module_oops.f90
elseif (trim(nam%model)=='wrf') then
   call model_wrf_coord(sdata)
else
   call msgerror('wrong model')
end if

end subroutine model_coord

!----------------------------------------------------------------------
! Subroutine: model_read
!> Purpose: read model field
!----------------------------------------------------------------------
subroutine model_read(filename,varname,sdata,fld)

implicit none

! Passed variables
character(len=*),intent(in) :: filename                 !< File name
character(len=*),intent(in) :: varname                  !< Variable name
type(sdatatype),intent(in) :: sdata                     !< Sampling data
real(kind_real),intent(out) :: fld(sdata%nc0,sdata%nl0) !< Read field

! Local variables
integer :: ncid
character(len=1024) :: subr = 'model_read'

! Open file
call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid))

! Select model
if (trim(nam%model)=='aro') then
   call model_aro_read(ncid,varname,sdata,fld)
elseif (trim(nam%model)=='arp') then
   call model_arp_read(ncid,varname,sdata,fld)
elseif (trim(nam%model)=='gem') then
   call model_gem_read(ncid,varname,sdata,fld)
elseif (trim(nam%model)=='geos') then
   call model_geos_read(ncid,varname,sdata,fld)
elseif (trim(nam%model)=='gfs') then
   call model_gfs_read(ncid,varname,sdata,fld)
elseif (trim(nam%model)=='ifs') then
   call model_ifs_read(ncid,varname,sdata,fld)
elseif (trim(nam%model)=='mpas') then
   call model_mpas_read(ncid,varname,sdata,fld)
elseif (trim(nam%model)=='nemo') then
   call model_nemo_read(ncid,varname,sdata,fld)
elseif (trim(nam%model)=='oops') then
   ! Specific interface, see module_oops.f90 TODO
elseif (trim(nam%model)=='wrf') then
   call model_wrf_read(ncid,varname,sdata,fld)
else
   call msgerror('wrong model')
end if

! Close file
call ncerr(subr,nf90_close(ncid))

end subroutine model_read

!----------------------------------------------------------------------
! Subroutine: model_write
!> Purpose: write model field
!----------------------------------------------------------------------
subroutine model_write(filename,varname,sdata,fld)

implicit none

! Passed variables
character(len=*),intent(in) :: filename                !< File name
character(len=*),intent(in) :: varname                 !< Variable name
type(sdatatype),intent(in) :: sdata                    !< Sampling data
real(kind_real),intent(in) :: fld(sdata%nc0,sdata%nl0) !< Written field

! Local variables
integer :: ierr
integer :: ncid
character(len=1024) :: subr = 'model_write_grid'

if (mpl%main) then
   ! Check if the file exists
   ierr = nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_noclobber,nf90_64bit_offset),ncid)
   if (ierr/=nf90_noerr) then
      call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_write,ncid))
      call ncerr(subr,nf90_redef(ncid))
      call ncerr(subr,nf90_put_att(ncid,nf90_global,'_FillValue',msvalr))
   end if
   call ncerr(subr,nf90_enddef(ncid))

   ! Select model
   if (trim(nam%model)=='aro') then
      call model_aro_write(ncid,varname,sdata,fld)
   elseif (trim(nam%model)=='arp') then
      call model_arp_write(ncid,varname,sdata,fld)
   elseif (trim(nam%model)=='gem') then
      call model_gem_write(ncid,varname,sdata,fld)
   elseif (trim(nam%model)=='geos') then
      call model_geos_write(ncid,varname,sdata,fld)
   elseif (trim(nam%model)=='gfs') then
      call model_gfs_write(ncid,varname,sdata,fld)
   elseif (trim(nam%model)=='ifs') then
      call model_ifs_write(ncid,varname,sdata,fld)
   elseif (trim(nam%model)=='mpas') then
      call model_mpas_write(ncid,varname,sdata,fld)
   elseif (trim(nam%model)=='nemo') then
      call model_nemo_write(ncid,varname,sdata,fld)
   elseif (trim(nam%model)=='oops') then
      ! Specific interface, see module_oops.f90
      call model_oops_write(ncid,varname,sdata,fld)
   elseif (trim(nam%model)=='wrf') then
      call model_wrf_write(ncid,varname,sdata,fld)
   else
      call msgerror('wrong model')
   end if

   ! Close file
   call ncerr(subr,nf90_close(ncid))
end if

end subroutine model_write

end module model_interface