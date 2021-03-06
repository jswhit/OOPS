! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module handling observations for the QG model

module c_qg_obs_data

use qg_obs_data
use qg_obs_vectors
use qg_obsoper_mod
use iso_c_binding
use config_mod
use fckit_log_module, only : fckit_log
use string_f_c_mod

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine qg_obs_setup(c_key_self, c_conf) bind(c,name='qg_obsdb_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf !< configuration

type(obs_data), pointer :: self
character(len=max_string) :: fin, fout
character(len=max_string+30) :: record

if (config_element_exists(c_conf,"ObsData.ObsDataIn")) then
  fin  = config_get_string(c_conf,max_string,"ObsData.ObsDataIn.obsfile")
else
  fin  = ""
endif
write(record,*)'qg_obs_setup: file in =',trim(fin)
call fckit_log%info(record)

fout = config_get_string(c_conf,max_string,"ObsData.ObsDataOut.obsfile")
write(record,*)'qg_obs_setup: file out=',trim(fout)
call fckit_log%info(record)

call obs_data_registry%init()
call obs_data_registry%add(c_key_self)
call obs_data_registry%get(c_key_self, self)
call obs_setup(trim(fin), trim(fout), self)

end subroutine qg_obs_setup

! ------------------------------------------------------------------------------

subroutine qg_obs_delete(c_key_self) bind(c,name='qg_obsdb_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(obs_data), pointer :: self

call obs_data_registry%get(c_key_self, self)
call obs_delete(self)
call obs_data_registry%remove(c_key_self)

end subroutine qg_obs_delete

! ------------------------------------------------------------------------------

subroutine qg_obs_get(c_key_self, lreq, c_req, lcol, c_col, c_key_ovec) bind(c,name='qg_obsdb_get_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lreq, lcol
character(kind=c_char,len=1), intent(in) :: c_req(lreq+1), c_col(lcol+1)
integer(c_int), intent(in) :: c_key_ovec

type(obs_data), pointer :: self
type(obs_vect), pointer :: ovec
character(len=lreq) :: req
character(len=lcol) :: col

call obs_data_registry%get(c_key_self, self)
call qg_obs_vect_registry%get(c_key_ovec,ovec)
call c_f_string(c_req, req)
call c_f_string(c_col, col)

call obs_get(self, trim(req), trim(col), ovec)

end subroutine qg_obs_get

! ------------------------------------------------------------------------------

subroutine qg_obs_put(c_key_self, lreq, c_req, lcol, c_col, c_key_ovec) bind(c,name='qg_obsdb_put_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lreq, lcol
character(kind=c_char,len=1), intent(in) :: c_req(lreq+1), c_col(lcol+1)
integer(c_int), intent(in) :: c_key_ovec

type(obs_data), pointer :: self
type(obs_vect), pointer :: ovec
character(len=lreq) :: req
character(len=lcol) :: col

call obs_data_registry%get(c_key_self, self)
call qg_obs_vect_registry%get(c_key_ovec,ovec)
call c_f_string(c_req, req)
call c_f_string(c_col, col)

call obs_put(self, trim(req), trim(col), ovec)

end subroutine qg_obs_put

! ------------------------------------------------------------------------------

end module c_qg_obs_data
