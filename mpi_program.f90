Program mpi_program
  ! ****************************************************************************
  !
  !  An MPI wrapper that calls a program multiple times in parallel.
  !
  ! ****************************************************************************

  use mpi
  
  implicit none
  
  ! MPI parameters
  integer                 , parameter :: ndim = 2
  logical, dimension(ndim), parameter :: periodic = (/.true.,.true./)
  integer, dimension(MPI_STATUS_SIZE) :: status
  
  ! MPI communicators
  integer :: ncpu
  integer :: cpuid
  integer :: CPU_MASTER
  integer :: ierr
  integer :: CART
  
  ! MPI local  variables
  integer, dimension(ndim) :: coord_l
  integer, dimension(ndim) :: dims
  integer :: x1, x2, dxl
  integer :: y1, y2, dyl

  ! MPI global variables
  integer, dimension(:,:), allocatable :: coord_g
  integer, dimension(:)  , allocatable :: xmin, xmax, dxg
  integer, dimension(:)  , allocatable :: ymin, ymax, dyg

  ! program parameters
  character(len=11), parameter :: sname    = 'MPI Program'
  
  ! program variables
  integer                    :: nx, ny, nz
  integer                    :: ix, iy
  character(len=250)         :: outstr
  
  ! ****************************************************************************

  ! initialize MPI
  call init_MPI

  ! This should be set by your input files
  nx = 10
  ny = 10
  nz = 250
  
  ! sets up the simple parallel scheme for the horizontal grid points input
  call setup_MPI

  write(outstr,'("computes X =",X,I3.3,X,"to",X,I3.3,", Y =",&
       &X,I3.3,X,"to",X,I3.3)') x1,x2,y1,y2
  call print(outstr,cpuid)

  do ix = x1, x2
     do iy = y1, y2
        ! This is where the call to the code code be made
        write(outstr,'("Running x =",X,I3.3,X,", y =",X,I3.3)') ix,iy
        call print(outstr,cpuid)
     end do
  end do

  ! End all communicators
  call MPI_Finalize(ierr)
  
contains
  ! ----------------------------------------------------------------------------
  !
  ! ----------------------------------------------------------------------------
  subroutine init_MPI
    ! **************************************************************************
    !
    !  This routine initialises MPI procedures.
    !
    ! **************************************************************************

    implicit none

    ! **************************************************************************

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,  ncpu, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, cpuid, ierr)
    CPU_MASTER = 0

  end subroutine init_MPI
  ! ----------------------------------------------------------------------------
  !
  ! ----------------------------------------------------------------------------
  subroutine setup_MPI
    ! **************************************************************************
    !
    !  This routine initialises all MPI procedures.
    !
    ! **************************************************************************

    implicit none

    integer :: i

    ! temporary subdomain ranges in master CPU
    integer :: ixr
    integer :: iyr

    ! **************************************************************************

    ! set up subdomain distribution
    dims = 0
    call MPI_DIMS_CREATE(ncpu, ndim, dims, ierr)

    ! get coordinate system
    call MPI_CART_CREATE(MPI_COMM_WORLD, ndim, dims, periodic, .true., &
                         CART, ierr)

    ! get the unique CPU ID number that distinguishes it from other processes
    ! This will be used for most of the initialization and all screen printing
    call MPI_COMM_RANK(CART, cpuid, ierr)

    call MPI_CART_COORDS(CART, cpuid, ndim, coord_l, ierr)
    
    call MPI_COMM_RANK(CART, cpuid, ierr)
    CPU_MASTER=0

    ! Gather data from every rank, combine and redistribute to all processors
    allocate(coord_g(2*ndim,0:ncpu-1))
    call MPI_ALLGATHER(coord_l          , ndim, MPI_INTEGER, &
                       coord_g(1:ndim,:), ndim, MPI_INTEGER, &
                       CART, ierr)

    ! set global and local grid point indices

    ! x-grid gridpoints
    allocate(xmin(0:ncpu-1)) ; xmin  = 0
    allocate(xmax(0:ncpu-1)) ; xmax  = 0
    allocate( dxg(0:ncpu-1)) ; dxg   = 0

    ! y-grid gridpoints
    allocate(ymin(0:ncpu-1)) ; ymin  = 0
    allocate(ymax(0:ncpu-1)) ; ymax  = 0
    allocate( dyg(0:ncpu-1)) ; dyg   = 0

    ixr = int(nx / (1.0*dims(1)) + 0.5)
    iyr = int(ny / (1.0*dims(2)) + 0.5)

    do i = 0, ncpu-1

       ! minimum indices for nx, ny, and nsnap
       xmin(i) = max(1, ixr * coord_g(1,i) + 1)
       ymin(i) = max(1, iyr * coord_g(2,i) + 1)

       ! maximum indices for nx, ny, and nsnap
       xmax(i) = min(nx, ixr * (coord_g(1,i) + 1))
       ymax(i) = min(ny, iyr * (coord_g(2,i) + 1))

       ! limiting maximum index by parameter
       if (coord_g(1,i) == (dims(1)-1)) xmax(i) = nx
       if (coord_g(2,i) == (dims(2)-1)) ymax(i) = ny

       ! limiting minimum index to maximum if min > max
       if (xmin(i) .gt. xmax(i)) xmin(i) = xmax(i)
       if (ymin(i) .gt. ymax(i)) ymin(i) = ymax(i)

       dxg(i) = xmax(i) - xmin(i) + 1
       dyg(i) = ymax(i) - ymin(i) + 1

    end do

    ! assign indices to each processor
    x1  = xmin(cpuid)
    x2  = xmax(cpuid)

    y1  = ymin(cpuid)
    y2  = ymax(cpuid)

    ! determine range of each subdomain
    dxl = x2 - x1  + 1
    dyl = y2 - y1  + 1

    ! syncronize indices on all subdomains
    call MPI_ALLGATHER(x1 , 1, MPI_INTEGER, xmin, 1, MPI_INTEGER, CART, ierr)
    call MPI_ALLGATHER(x2 , 1, MPI_INTEGER, xmax, 1, MPI_INTEGER, CART, ierr)
    call MPI_ALLGATHER(y1 , 1, MPI_INTEGER, ymin, 1, MPI_INTEGER, CART, ierr)
    call MPI_ALLGATHER(y2 , 1, MPI_INTEGER, ymax, 1, MPI_INTEGER, CART, ierr)
    call MPI_ALLGATHER(dxl, 1, MPI_INTEGER, dxg , 1, MPI_INTEGER, CART, ierr)
    call MPI_ALLGATHER(dyl, 1, MPI_INTEGER, dyg , 1, MPI_INTEGER, CART, ierr)

  end subroutine setup_MPI
  ! ----------------------------------------------------------------------------
  !
  ! ----------------------------------------------------------------------------
  subroutine print(text,pid)
    ! **************************************************************************
    !
    !  Prints text to master CPU. This avoids the message being printed N_CPU
    !  times
    !
    ! **************************************************************************

    implicit none

    ! dummy arguments
    character(len=*)          , intent(IN) :: text
    integer         , optional, intent(IN) :: pid

    ! local variables
    character(len=5) :: spid
    integer          :: strl

    ! **************************************************************************

    strl = len_trim(text)

    if (present(pid)) then
       write(spid,'(I4.4)') pid
       if (cpuid .eq. pid) then
          write(*,'(A6,X,A5,A,X,A)') 'CPU ID',spid,':',text(1:strl)
       endif
    elseif (cpuid .eq. CPU_MASTER) then
       write(*,'(A)') text(1:strl)
    end if

  end subroutine print
  ! ----------------------------------------------------------------------------
  !
  ! ----------------------------------------------------------------------------
  subroutine stop(text)
    ! **************************************************************************
    !
    !  Stops Linfor3D (with an optional message) and terminates the MPI
    !  execution environment
    !
    ! **************************************************************************

    implicit none
    
    ! dummy arguments
    character(len=*), optional, intent(IN) :: text

    ! local variables
    integer :: error

    ! **************************************************************************

    if (present(text)) call print(text, cpuid)

    call MPI_ABORT(MPI_COMM_WORLD, error, ierr)

    stop

  end subroutine stop
  ! ----------------------------------------------------------------------------
  !
  ! ----------------------------------------------------------------------------
end Program mpi_program
