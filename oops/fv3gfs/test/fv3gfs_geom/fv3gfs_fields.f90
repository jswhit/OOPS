program fv3gfs_cubic_test

    #include <fms_platform.h>

    use mpp_mod,         only: mpp_pe, mpp_npes, mpp_init, mpp_exit
    use mpp_mod,         only: stdout, mpp_error, FATAL, NOTE
    use mpp_mod,         only: input_nml_file
    use mpp_domains_mod, only: domain2D
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

    integer :: sizex_latlon_grid = 144
    integer :: sizey_latlon_grid = 90
    integer :: size_cubic_grid = 48
    integer :: nz = 10
    integer :: halo = 1
    integer :: stackmax = 4000000
    integer :: layout_cubic(2)  = (/4,2/)
    integer :: layout_latlon(2) = (/4,2/)
    integer :: io_layout(2) = (/1,1/) ! set ndivs_x and ndivs_y to divide each tile into io_layout(1)*io_layout(2)
                                      ! group and write out data from the root pe of each group.
    namelist /fv3gfs_geom_nml/ sizex_latlon_grid, sizey_latlon_grid, size_cubic_grid, &
                               nz, halo, stackmax, layout_cubic, layout_latlon, io_layout

    type(domain2D), save :: domain_cubic
    type(restart_file_type) :: Fv_restart
    character(len=64)    :: filename
    integer, parameter :: ntile_cubic = 6
    integer :: pe, npes, isc, isd, iec, ied, jsc, jsd, jec, jed, id_restart
    integer :: nmlunit, outunit, io_status
    real, allocatable, dimension(:,:) :: phis
    real, allocatable, dimension(:,:,:) :: pt

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

! register 2d restart field.
    filename = 'fv_core.res.nc'
    allocate ( phis(isd:ied  ,jsd:jed ) )
    id_restart = register_restart_field(Fv_restart, filename, 'phis', phis, &
                 domain=domain_cubic)
! read 2d restart field.
    call restore_state(Fv_restart, id_restart, directory='INPUT')
    print *,'pe,shape,minval,maxval for phis',pe,shape(phis),minval(phis),maxval(phis)

! register 3d restart field.
    allocate ( pt(isd:ied  ,jsd:jed, nz ) )
    id_restart = register_restart_field(Fv_restart, filename, 'T', pt, &
                 domain=domain_cubic)
! read 2d restart field.
    call restore_state(Fv_restart, id_restart, directory='INPUT')
    print *,'pe,shape,minval,maxval for pt',pe,shape(pt),minval(pt),maxval(pt)

! clean up and exit.
    call fms_io_exit
    call mpp_domains_exit
    call mpp_exit

end program fv3gfs_cubic_test
