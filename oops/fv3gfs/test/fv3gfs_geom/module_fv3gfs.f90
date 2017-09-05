module fv3gfs_mod

    #include <fms_platform.h>

    use mpp_mod,         only: mpp_pe, mpp_npes, mpp_init, mpp_exit
    use mpp_mod,         only: stdout, mpp_error, FATAL, NOTE
    use mpp_mod,         only: input_nml_file
    use mpp_domains_mod, only: domain2D, mpp_define_layout, mpp_define_mosaic
    use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain
    use mpp_domains_mod, only: mpp_domains_init, mpp_domains_exit
    use mpp_domains_mod, only: mpp_domains_set_stack_size, mpp_define_io_domain
    use mpp_io_mod,      only: mpp_open, mpp_close, MPP_ASCII, MPP_RDONLY
    use fms_io_mod,      only: fms_io_init, fms_io_exit
    use fms_io_mod,      only: file_exist

    implicit none

    private

    public :: setup_geom

contains

subroutine setup_geom(domain, gtype, nx, ny, ntiles, layout_in, io_layout_in, halo)

    implicit none

    type(domain2D),   intent(inout) :: domain
    character(len=*), intent(in)    :: gtype
    integer,          intent(in)    :: nx, ny, ntiles
    integer,          intent(in)    :: layout_in(:), io_layout_in(:)
    integer,          intent(in)    :: halo

    integer                              :: pe, npes, npes_per_tile, tile
    integer                              :: num_contact
    integer                              :: n, layout(2)
    integer, allocatable, dimension(:,:) :: global_indices, layout2D
    integer, allocatable, dimension(:)   :: pe_start, pe_end
    integer, allocatable, dimension(:)   :: tile1, tile2
    integer, allocatable, dimension(:)   :: istart1, iend1, jstart1, jend1
    integer, allocatable, dimension(:)   :: istart2, iend2, jstart2, jend2
    integer, allocatable :: tile_id(:)
    integer                              :: isc, iec, jsc, jec
    integer                              :: isd, ied, jsd, jed
    logical is_symmetry

    pe = mpp_pe()
    npes = mpp_npes()

    if (mod(npes,ntiles) /= 0) then
       call mpp_error(NOTE, "setup_geom: npes can not be divided by ntiles, no test will be done for "//trim(gtype))
       return
    endif
    npes_per_tile = npes/ntiles
    tile = pe/npes_per_tile + 1

    if (layout_in(1)*layout_in(2) == npes_per_tile) then
       layout = layout_in
    else
       call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
    endif

    if (io_layout_in(1) <1 .or. io_layout_in(2) <1) call mpp_error(FATAL, &
            "setup_geom: both elements of variable io_layout_in must be positive integer")
    if (mod(layout(1), io_layout_in(1)) /= 0 ) call mpp_error(FATAL, &
         "setup_geom: layout(1) must be divided by io_layout_in(1)")
    if (mod(layout(2), io_layout_in(2)) /= 0 ) call mpp_error(FATAL, &
         "setup_geom: layout(2) must be divided by io_layout_in(2)")

    allocate(global_indices(4,ntiles), layout2D(2,ntiles), pe_start(ntiles), pe_end(ntiles) )
    do n = 1, ntiles
       global_indices(:,n) = (/1,nx,1,ny/)
       layout2D(:,n)       = layout
       pe_start(n)         = (n-1)*npes_per_tile
       pe_end(n)           = n*npes_per_tile-1
    enddo

    if (gtype == 'cubic_grid') then
        num_contact = 12
        allocate(tile1(num_contact), tile2(num_contact) )
        allocate(tile_id(ntiles))
        allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) )
        allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) )
        !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
        tile1(1) = 1; tile2(1) = 2
        istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
        istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
        !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
        tile1(2) = 1; tile2(2) = 3
        istart1(2) = 1;  iend1(2) = nx; jstart1(2) = ny; jend1(2) = ny
        istart2(2) = 1;  iend2(2) = 1;  jstart2(2) = ny; jend2(2) = 1
        !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
        tile1(3) = 1; tile2(3) = 5
        istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
        istart2(3) = nx; iend2(3) = 1;  jstart2(3) = ny; jend2(3) = ny
        !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
        tile1(4) = 1; tile2(4) = 6
        istart1(4) = 1;  iend1(4) = nx; jstart1(4) = 1;  jend1(4) = 1
        istart2(4) = 1;  iend2(4) = nx; jstart2(4) = ny; jend2(4) = ny
        !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
        tile1(5) = 2; tile2(5) = 3
        istart1(5) = 1;  iend1(5) = nx; jstart1(5) = ny; jend1(5) = ny
        istart2(5) = 1;  iend2(5) = nx; jstart2(5) = 1;  jend2(5) = 1
        !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
        tile1(6) = 2; tile2(6) = 4
        istart1(6) = nx; iend1(6) = nx; jstart1(6) = 1;  jend1(6) = ny
        istart2(6) = nx; iend2(6) = 1;  jstart2(6) = 1;  jend2(6) = 1
        !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
        tile1(7) = 2; tile2(7) = 6
        istart1(7) = 1;  iend1(7) = nx; jstart1(7) = 1;  jend1(7) = 1
        istart2(7) = nx; iend2(7) = nx; jstart2(7) = ny; jend2(7) = 1
        !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
        tile1(8) = 3; tile2(8) = 4
        istart1(8) = nx; iend1(8) = nx; jstart1(8) = 1;  jend1(8) = ny
        istart2(8) = 1;  iend2(8) = 1;  jstart2(8) = 1;  jend2(8) = ny
        !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
        tile1(9) = 3; tile2(9) = 5
        istart1(9) = 1;  iend1(9) = nx; jstart1(9) = ny; jend1(9) = ny
        istart2(9) = 1;  iend2(9) = 1;  jstart2(9) = ny; jend2(9) = 1
        !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
        tile1(10) = 4; tile2(10) = 5
        istart1(10) = 1;  iend1(10) = nx; jstart1(10) = ny; jend1(10) = ny
        istart2(10) = 1;  iend2(10) = nx; jstart2(10) = 1;  jend2(10) = 1
        !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
        tile1(11) = 4; tile2(11) = 6
        istart1(11) = nx; iend1(11) = nx; jstart1(11) = 1;  jend1(11) = ny
        istart2(11) = nx; iend2(11) = 1;  jstart2(11) = 1;  jend2(11) = 1
        !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
        tile1(12) = 5; tile2(12) = 6
        istart1(12) = nx; iend1(12) = nx; jstart1(12) = 1;  jend1(12) = ny
        istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;  jend2(12) = ny
        is_symmetry = .true.
    else if (gtype == 'latlon_grid') then
        num_contact = 8
        allocate(tile1(num_contact), tile2(num_contact), tile_id(1))
        allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) )
        allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) )
        !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
        tile1(1) = 1; tile2(1) = 2
        istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
        istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
        !--- Contact line 2, between tile 1 (SOUTH) and tile 3 (NORTH)  --- cyclic
        tile1(2) = 1; tile2(2) = 3
        istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
        istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
        !--- Contact line 3, between tile 1 (WEST) and tile 2 (EAST) --- cyclic
        tile1(3) = 1; tile2(3) = 2
        istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
        istart2(3) = nx; iend2(3) = nx; jstart2(3) = 1;  jend2(3) = ny
        !--- Contact line 4, between tile 1 (NORTH) and tile 3 (SOUTH)
        tile1(4) = 1; tile2(4) = 3
        istart1(4) = 1;  iend1(4) = nx; jstart1(4) = ny;  jend1(4) = ny
        istart2(4) = 1;  iend2(4) = nx; jstart2(4) = 1;   jend2(4) = 1
        !--- Contact line 5, between tile 2 (SOUTH) and tile 4 (NORTH) --- cyclic
        tile1(5) = 2; tile2(5) = 4
        istart1(5) = 1;  iend1(5) = nx; jstart1(5) = 1;  jend1(5) = 1
        istart2(5) = 1;  iend2(5) = nx; jstart2(5) = ny; jend2(5) = ny
        !--- Contact line 6, between tile 2 (NORTH) and tile 4 (SOUTH)
        tile1(6) = 2; tile2(6) = 4
        istart1(6) = 1;  iend1(6) = nx; jstart1(6) = ny;  jend1(6) = ny
        istart2(6) = 1;  iend2(6) = nx; jstart2(6) = 1;   jend2(6) = 1
        !--- Contact line 7, between tile 3 (EAST) and tile 4 (WEST)
        tile1(7) = 3; tile2(7) = 4
        istart1(7) = nx; iend1(7) = nx; jstart1(7) = 1;  jend1(7) = ny
        istart2(7) = 1;  iend2(7) = 1;  jstart2(7) = 1;  jend2(7) = ny
        !--- Contact line 8, between tile 3 (WEST) and tile 4 (EAST) --- cyclic
        tile1(8) = 3; tile2(8) = 4
        istart1(8) = 1;  iend1(8) = 1;  jstart1(8) = 1;  jend1(8) = ny
        istart2(8) = nx; iend2(8) = nx; jstart2(8) = 1;  jend2(8) = ny
        is_symmetry = .false.
    else
       call mpp_error(NOTE, "setup_geom: unknown gtype "//trim(gtype))
       return
    endif
    do n = 1, ntiles
       tile_id(n) = n
    enddo

    call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                           istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                           pe_start, pe_end, whalo=halo, ehalo=halo, shalo=halo, nhalo=halo,    &
                           symmetry=is_symmetry,tile_id=tile_id,name=gtype)

    if (io_layout_in(1) /= 1 .or. io_layout_in(2) /= 1) &
        call mpp_define_io_domain(domain, io_layout_in)

    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(domain, isd, ied, jsd, jed)

    deallocate(pe_start,pe_end)
    deallocate(layout2D, global_indices)
    deallocate(tile1, tile2, tile_id)
    deallocate(istart1, iend1, jstart1, jend1)
    deallocate(istart2, iend2, jstart2, jend2)

end subroutine setup_geom

end module fv3gfs_mod
