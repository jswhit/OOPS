
!#include <fms_platform.h>

! ------------------------------------------------------------------------------

subroutine fv3gfs_setup(c_conf) bind(c,name='fv3gfs_setup_f')

use iso_c_binding
use config_mod
use mpp_mod,         only: mpp_init
use mpp_domains_mod, only: mpp_domains_init
use mpp_domains_mod, only: mpp_domains_set_stack_size
use fms_io_mod,      only: fms_io_init

implicit none

type(c_ptr), intent(in) :: c_conf
integer :: stackmax = 4000000

call mpp_init
call mpp_domains_init
call fms_io_init

if (config_element_exists(c_conf,"stackmax")) stackmax = config_get_int(c_conf,"stackmax")
call mpp_domains_set_stack_size(stackmax)

end subroutine fv3gfs_setup

! ------------------------------------------------------------------------------

subroutine fv3gfs_finalize() bind(c,name='fv3gfs_finalize_f')

use mpp_mod,         only: mpp_exit
use mpp_domains_mod, only: mpp_domains_exit
use fms_io_mod,      only: fms_io_exit

implicit none

call fms_io_exit
call mpp_domains_exit
call mpp_exit

end subroutine fv3gfs_finalize

! ------------------------------------------------------------------------------
