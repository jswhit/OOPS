program fv3gfs_fields

    #include <fms_platform.h>

    use mpp_mod,         only: mpp_pe, mpp_npes, mpp_init, mpp_exit
    use mpp_mod,         only: stdout, mpp_error, FATAL, NOTE
    use mpp_mod,         only: input_nml_file
    use mpp_domains_mod, only: domain2d, EAST, WEST, NORTH, CENTER, SOUTH,&
                               CORNER
    use mpp_domains_mod, only: mpp_domains_init, mpp_domains_exit
    use mpp_domains_mod, only: mpp_domains_set_stack_size
    use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain
    use mpp_io_mod,      only: mpp_open, mpp_close, MPP_ASCII, MPP_RDONLY
    use fms_io_mod,      only: fms_io_init, fms_io_exit, get_tile_string, &
                               restart_file_type, register_restart_field, &
                               save_restart, restore_state, file_exist, &
                               set_domain, nullify_domain, set_filename_appendix, &
                               get_mosaic_tile_file, get_instance_filename, & 
                               save_restart_border, restore_state_border, free_restart_type, &
                               field_exist

    use fv3gfs_mod, only: setup_geom

    implicit none

    type fv_atmos_type

       logical :: allocated = .false.
       logical :: hydrostatic = .false.  ! nonhydrostatic fields in restart?
       logical :: agrid_vel_rst = .false. ! agrid winds in restart?

!-----------------------------------------------------------------------
! Five prognostic state variables for the f-v dynamics
!-----------------------------------------------------------------------
! dyn_state:
! D-grid prognostatic variables: u, v, and delp (and other scalars)
!
!     o--------u(i,j+1)----------o
!     |           |              |
!     |           |              |
!  v(i,j)------scalar(i,j)----v(i+1,j)
!     |           |              |
!     |           |              |
!     o--------u(i,j)------------o
!
! The C grid component is "diagnostic" in that it is predicted every time step
! from the D grid variables.
      real, allocatable :: u(:,:,:)      ! D grid zonal wind (m/s)
      real, allocatable :: v(:,:,:)      ! D grid meridional wind (m/s)
      real, allocatable :: pt(:,:,:)     ! temperature (K)
      real, allocatable :: delp(:,:,:)   ! pressure thickness (pascal)
      real, allocatable :: q(:,:,:,:)    ! specific humidity and prognostic constituents

!----------------------
! non-hydrostatic state:
!----------------------------------------------------------------------
      real, allocatable ::     w(:,:,:)    ! cell center vertical wind (m/s)
      real, allocatable ::  delz(:,:,:)    ! layer thickness (meters)

!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
      real, allocatable :: phis(:,:)       ! Surface geopotential (g*Z_surf)
      real, allocatable :: ua(:,:,:)       ! (ua, va) are mostly used as the A grid winds
      real, allocatable :: va(:,:,:)     

      real, allocatable :: ak(:)  
      real, allocatable :: bk(:)  

    end type fv_atmos_type

    integer :: sizex_latlon_grid = 144
    integer :: sizey_latlon_grid = 90
    integer :: size_cubic_grid = 48
    integer :: nz = 10
    integer :: nq = 3
    logical :: hydrostatic = .false.  ! nonhydrostatic fields in restart?
    logical :: agrid_vel_rst = .false. ! agrid winds in restart?
    integer :: halo = 1
    integer :: stackmax = 4000000
    integer :: layout_cubic(2)  = (/4,2/)
    integer :: layout_latlon(2) = (/4,2/)
    integer :: io_layout(2) = (/1,1/) ! set ndivs_x and ndivs_y to divide each tile into io_layout(1)*io_layout(2)
                                      ! group and write out data from the root pe of each group.
    namelist /fv3gfs_geom_nml/ sizex_latlon_grid, sizey_latlon_grid, size_cubic_grid, &
                               nq, nz, hydrostatic, halo, stackmax, layout_cubic, layout_latlon, io_layout, &
                               agrid_vel_rst

    type(domain2D), save :: domain_cubic
    type(restart_file_type) :: Fv_restart
    type(fv_atmos_type) :: Atm
    character(len=64)    :: filename
    integer, parameter :: ntile_cubic = 6
    integer :: pe, npes, isc, isd, iec, ied, jsc, jsd, jec, jed, id_restart
    integer :: nmlunit, outunit, io_status

! initialization.
    call mpp_init
    call mpp_domains_init
    call fms_io_init

    pe = mpp_pe()
    npes = mpp_npes()

! read namelist.
    if (file_exist('input.nml') )then
       call mpp_open(nmlunit, 'input.nml', form=MPP_ASCII, action=MPP_RDONLY)
       read(nmlunit,fv3gfs_geom_nml,iostat=io_status)
       call mpp_close(nmlunit)
    endif
    if (io_status > 0) then
       call mpp_error(FATAL, '=>fv3gfs_geom: Error reading fv3gfs_geom_nml')
    endif

    outunit = stdout()
    write(outunit, fv3gfs_geom_nml)
    call mpp_domains_set_stack_size(stackmax)

! set up domain geometry.
    call setup_geom(domain_cubic, "cubic_grid", size_cubic_grid, size_cubic_grid, ntile_cubic, layout_cubic, io_layout, halo)
    if (mod(npes,ntile_cubic) == 0) &
        call mpp_error(NOTE, "fv3gfs_geom: setup_geom is done for cubic_grid")

! get compute and data domains (used to allocate data arrays).
    call mpp_get_compute_domain(domain_cubic, isc, iec, jsc, jec)
    call mpp_get_data_domain(domain_cubic, isd, ied, jsd, jed)
    print *,pe,'compute domain',isc, iec, jsc, jec
    print *,pe,'data domain',isd, ied, jsd, jed

! allocate fv_atmos structure.
    call allocate_fv_atmos_type(Atm, isd, ied, jsd, jed, isc, iec, jsc, jec, nz, nq, &
                                hydrostatic, agrid_vel_rst)

! read restart fields from fv_core.
    filename = 'fv_core.res.nc'
    id_restart = register_restart_field(Fv_restart, filename, 'ak', Atm%ak(:), no_domain=.true.)
    call restore_state(Fv_restart, id_restart, directory='INPUT')
    id_restart = register_restart_field(Fv_restart, filename, 'bk', Atm%bk(:), no_domain=.true.)
    call restore_state(Fv_restart, id_restart, directory='INPUT')
    call free_restart_type(Fv_restart)
    if (pe == 0) print *,'ak=',Atm%ak
    if (pe == 0) print *,'bk=',Atm%bk
    id_restart = register_restart_field(Fv_restart, filename, 'phis', Atm%phis, &
                 domain=domain_cubic)
    call restore_state(Fv_restart, id_restart, directory='INPUT')
    print *,'pe,shape,minval,maxval for phis',pe,shape(Atm%phis),minval(Atm%phis),maxval(Atm%phis)
    id_restart = register_restart_field(Fv_restart, filename, 'T', Atm%pt, &
                 domain=domain_cubic)
    call restore_state(Fv_restart, id_restart, directory='INPUT')
    print *,'pe,shape,minval,maxval for pt',pe,shape(Atm%pt),minval(Atm%pt),maxval(Atm%pt)
! FIXME:  u and v don't work, get this error
!FATAL from PE     2: fms_io(setup_one_field): data should be on either compute
!domain or data domain when domain is present for field u of file fv_core.res.nc
!   id_restart = register_restart_field(Fv_restart, filename, 'u', Atm%u, &
!                domain=domain_cubic,position=NORTH)
!   call restore_state(Fv_restart, id_restart, directory='INPUT')
!   id_restart = register_restart_field(Fv_restart, filename, 'v', Atm%v, &
!                domain=domain_cubic,position=EAST)
!   call restore_state(Fv_restart, id_restart, directory='INPUT')
    id_restart = register_restart_field(Fv_restart, filename, 'DELP', Atm%delp, &
                 domain=domain_cubic)
    call restore_state(Fv_restart, id_restart, directory='INPUT')
    if (Atm%agrid_vel_rst) then
        id_restart =  register_restart_field(Fv_restart, filename, 'ua', Atm%ua, &
                      domain=domain_cubic)
        call restore_state(Fv_restart, id_restart, directory='INPUT')
        id_restart =  register_restart_field(Fv_restart, filename, 'va', Atm%va, &
                      domain=domain_cubic)
        call restore_state(Fv_restart, id_restart, directory='INPUT')
    endif
    if (.not. Atm%hydrostatic) then
        id_restart =  register_restart_field(Fv_restart, filename, 'W', Atm%w, &
                      domain=domain_cubic)
        call restore_state(Fv_restart, id_restart, directory='INPUT')
        id_restart =  register_restart_field(Fv_restart, filename, 'DZ', Atm%delz, &
                      domain=domain_cubic)
        call restore_state(Fv_restart, id_restart, directory='INPUT')
    endif

! read tracers (should use field_table for this).
    call free_restart_type(Fv_restart)
    filename = 'fv_tracer.res.nc'
    id_restart = register_restart_field(Fv_restart, filename, 'sphum', Atm%q(:,:,:,1), &
                 domain=domain_cubic)
    call restore_state(Fv_restart, id_restart, directory='INPUT')
    id_restart = register_restart_field(Fv_restart, filename, 'o3mr', Atm%q(:,:,:,2), &
                 domain=domain_cubic)
    call restore_state(Fv_restart, id_restart, directory='INPUT')
    id_restart = register_restart_field(Fv_restart, filename, 'liq_wat', Atm%q(:,:,:,3), &
                 domain=domain_cubic)
    call restore_state(Fv_restart, id_restart, directory='INPUT')

! clean up and exit.
    call fms_io_exit
    call mpp_domains_exit
    call mpp_exit

contains

  subroutine allocate_fv_atmos_type(Atm, isd, ied, jsd, jed, isc, iec, jsc, jec, &
                                    nz, nq, hydrostatic, agrid_vel_rst)

    !WARNING: Before calling this routine, be sure to have set up the
    ! proper domain parameters.

    implicit none
    type(fv_atmos_type), intent(INOUT), target :: Atm
    logical, intent(IN) :: hydrostatic, agrid_vel_rst
    integer, intent(IN) :: isd, ied, jsd, jed, isc, iec, jsc, jec, nz, nq

    Atm%hydrostatic = hydrostatic
    Atm%agrid_vel_rst = agrid_vel_rst

    if (Atm%allocated) return

    allocate (    Atm%u(isd:ied  ,jsd:jed+1,nz) )
    allocate (    Atm%v(isd:ied+1,jsd:jed  ,nz) )
    allocate (   Atm%pt(isd:ied  ,jsd:jed  ,nz) )
    allocate ( Atm%delp(isd:ied  ,jsd:jed  ,nz) )
    allocate (    Atm%q(isd:ied  ,jsd:jed  ,nz, nq) )
    allocate ( Atm%phis(isd:ied  ,jsd:jed  ) )
    allocate ( Atm%ak(nz+1) )
    allocate ( Atm%bk(nz+1) )

    !--- include agrid winds in restarts for use in data assimilation 
    if (Atm%agrid_vel_rst) then
       allocate ( Atm%ua(isd:ied  ,jsd:jed  ,nz) )
       allocate ( Atm%va(isd:ied  ,jsd:jed  ,nz) )
    endif

    !--------------------------
    ! Non-hydrostatic dynamics:
    !--------------------------
    if (.not. Atm%hydrostatic ) then
       allocate (    Atm%w(isd:ied, jsd:jed  ,nz  ) )
       allocate ( Atm%delz(isd:ied, jsd:jed  ,nz) )
    endif

  end subroutine allocate_fv_atmos_type

end program fv3gfs_fields
