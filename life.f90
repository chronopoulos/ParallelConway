! Parallel Game of Life
! Usage:
!  mpirun -np 4 ./life conway-900x900.pgm -blockrow -niter=10000 -count=1000
!  mpirun -np 6 ./life conway-900x900.pgm -blockcol -niter=5000 -count=1000
!  mpirun -np 9 ./life conway-900x900.pgm -checker -niter=5000 -count=1000
!  mpirun -np 9 ./life conway-900x900.pgm -checker -verbose
!  mpirun -np 9 ./life conway-900x900.pgm -checker -io=1 -niter=10 -verbose

program life

   use mpi
   use pgm
   implicit none

   integer :: np, sqrtnp, myrank, ierr
   character(len=64) :: arg, initfile, dbgmsg, catcommand, outfile
   character(len=5) :: rankstr
   character(len=6) :: tstring

   logical :: cnt=.false., verbose=.false., io=.false.
   integer :: t, niter=1000, cntfreq=1000, iofreq=1
   integer :: part=0   ! partition mode, default is blockrow

   double precision :: tstart, tstop   ! timing

   integer :: i, j, drow, dcol
   integer :: mybugs, allbugs, live_neighbors

   integer :: nrows, ncols, mycol, myrow, local_height, local_width
   integer :: gheight, gwidth ! height and width of the global array
   integer*1, allocatable :: mygrid(:,:), mychange(:,:)

   integer:: stat(MPI_STATUS_SIZE)

   logical :: is_self
   

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   ! Parse arguments

   do i=1,command_argument_count()
      call get_command_argument(i,arg)
      if (arg(1:9)=='-blockrow') then
         part=0 ! block row partition mode
      else if (arg(1:9)=='-blockcol') then
         part=1 ! block column partition mode
      else if (arg(1:8)=='-checker') then
         part=2 ! checkerboard partition mode
      else if (arg(1:7)=='-niter=') then
         read(arg(8:),*) niter
      else if (arg(1:7)=='-count=') then
         read(arg(8:),*) cntfreq
         cnt=.true.
      else if (arg(1:8)=='-verbose') then
         verbose=.true.
      else if (arg(1:4)=='-io=') then
         read(arg(5:),*) iofreq
         io=.true.
      else
         initfile=arg
      end if
   end do

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   ! Initialize MPI

   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Assign subdomains, read in PGM file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if (part .eq. 0) then ! block row
      nrows=np
      ncols=1
      myrow=myrank+1
      mycol=1
   else if (part .eq. 1) then ! block column
      nrows=1
      ncols=np
      myrow=1
      mycol=myrank+1
   else if (part .eq. 2) then ! checkerboard
      if (isPerfectSquare(np)) then
         sqrtnp=sqrt(real(np))
         nrows=sqrtnp
         ncols=sqrtnp
         ! this is a row-major decomposition
         myrow=myrank/sqrtnp+1
         mycol=mod(myrank,sqrtnp)+1
      else
         stop 'for checkerboard, np must be perfect square'
      end if
   end if

   call read_pgm(trim(initfile)//char(0), nrows, ncols, mycol, myrow, &
      gheight, gwidth, local_height, local_width, mygrid)

   allocate (mychange(0:local_height+1,0:local_width+1))
   mychange=0  ! possibly unnecessary

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Count Initial Bugs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if (cnt) then

      mybugs = 0
      do j=1,local_width
      do i=1,local_height
         if (mygrid(i,j) > 0) then
            mybugs = mybugs + 1
         endif
      enddo
      enddo

      call MPI_Allreduce(mybugs, allbugs, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! output
      if (myrank == 0) then
         print*, 'Total Bugs at Initialization: ', allbugs
      endif

   end if

   ! Write out initial state
   if (io) then
      write(outfile,'(A6,I6.6)') 'output', 0
      if (myrank==0) print*, 'Writing out to '//outfile
      call parallel_io(outfile)
   end if

   stop 'CHECKPOINT'

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Timestepping
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   tstart=MPI_Wtime()

   do t=1,niter

      ! Compute Change (Rules)
      do j=1,local_width
      do i=1,local_height

         ! count live neighbors
         live_neighbors=0
         do drow=-1,1
         do dcol=-1,1
            is_self = (drow==0) .and. (dcol==0)
            if (.not. is_self) live_neighbors = live_neighbors + mygrid(i+drow,j+dcol)
         end do
         end do

         ! assign changes accordingly
         if (mygrid(i,j)==1) then ! if alive
            if (live_neighbors .le. 1 .or. live_neighbors .ge. 4) then
               mychange(i,j)=-1 ! cell dies
            else
               mychange(i,j)=0  ! cell survives
            end if
         else ! if dead
            if (live_neighbors==3) then
               mychange(i,j)=1  ! cell is born
            else
               mychange(i,j)=0  ! cell remains dead
            end if
         end if
      
      end do
      end do

      ! Apply changes to internal field
      mygrid=mygrid+mychange

      ! Update ghost points
      if (np /= 1) then ! if np=1, skip this part

         if (part .eq. 0) then ! block row
            call ghostpoints_blockrow(mygrid)
         else if (part .eq. 1) then ! block column
            call ghostpoints_blockcol(mygrid)
         else if (part .eq. 2) then ! checkerboard
            call ghostpoints_checker(mygrid)
         end if

      end if

      ! Print Timestep
      if ((verbose).and.(myrank==0)) print*, 'Timestep: ', t

      ! Count Bugs
      if ((cnt).and.(mod(t,cntfreq)==0)) then

         mybugs = 0
         do j=1,local_width
         do i=1,local_height
            if (mygrid(i,j) > 0) then
               mybugs = mybugs + 1
            endif
         enddo
         enddo

         call MPI_Allreduce(mybugs, allbugs, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

         ! output
         if (myrank == 0) then
            print*, 'Total Bugs after ', t, ' timesteps: ', allbugs
         endif

      end if

      ! Parallel IO
      if ((io).and.(mod(t,iofreq)==0)) then
         write(outfile,'(A6,I6.6)') 'output', t
         if (myrank==0) print*, 'Writing out to '//outfile
         call parallel_io(outfile)
         ! now append header (tooooo slow)
         !call system('cat header.txt '//trim(outfile)//' > '//trim(outfile)//'.pgm')
      end if


   end do ! timestepping

   tstop=MPI_Wtime()
   if (myrank==0) print*, 'Execution Time in Seconds: ', tstop-tstart

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   call MPI_Finalize(ierr)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

   logical function isPerfectSquare(n)
      implicit none
      integer, intent(in) :: n
      real :: root

      root=sqrt(real(n))
      if (root==aint(root)) then
         isPerfectSquare = .true.
      else
         isPerfectSquare = .false.
      end if
   end function isPerfectSquare


   subroutine ghostpoints_blockrow(mygrid)

      implicit none
      integer*1, intent(inout) :: mygrid(0:local_height+1,0:local_width+1)
      integer*1 :: myghostrow(local_width)

      ! first send bottom rows down
      if (myrank==0) then ! rank 0 only sends
         myghostrow=mygrid(local_height,1:local_width)
         call MPI_Send(myghostrow, local_width, MPI_INTEGER1, myrank+1, 0, MPI_COMM_WORLD, ierr)
      else if (myrank==np-1) then ! rank np-1 only recv's
         call MPI_Recv(myghostrow, local_width, MPI_INTEGER1, myrank-1, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(0,1:local_width)=myghostrow
      else ! all other ranks send first, then recv
         myghostrow=mygrid(local_height,1:local_width)
         call MPI_Send(myghostrow, local_width, MPI_INTEGER1, myrank+1, 0, MPI_COMM_WORLD, ierr)
         call MPI_Recv(myghostrow, local_width, MPI_INTEGER1, myrank-1, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(0,1:local_width)=myghostrow
      end if

      ! next send top rows up
      if (myrank==0) then ! rank 0 only recvs
         call MPI_Recv(myghostrow, local_width, MPI_INTEGER1, myrank+1, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(local_height+1,1:local_width)=myghostrow
      else if (myrank==np-1) then ! rank np-1 only sends
         myghostrow=mygrid(1,1:local_width)
         call MPI_Send(myghostrow, local_width, MPI_INTEGER1, myrank-1, 0, MPI_COMM_WORLD, ierr)
      else ! all other ranks send first, then recv
         myghostrow=mygrid(1,1:local_width)
         call MPI_Send(myghostrow, local_width, MPI_INTEGER1, myrank-1, 0, MPI_COMM_WORLD, ierr)
         call MPI_Recv(myghostrow, local_width, MPI_INTEGER1, myrank+1, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(local_height+1,1:local_width)=myghostrow
      end if

   end subroutine ghostpoints_blockrow


   subroutine ghostpoints_blockcol(mygrid)

      implicit none
      integer*1, intent(inout) :: mygrid(0:local_height+1,0:local_width+1)
      integer*1 :: myghostcol(local_height)

      ! first send rightmost rows to the right
      if (myrank==0) then ! rank 0 only sends
         myghostcol=mygrid(1:local_height,local_width)
         call MPI_Send(myghostcol, local_height, MPI_INTEGER1, myrank+1, 0, MPI_COMM_WORLD, ierr)
      else if (myrank==np-1) then ! rank np-1 only recv's
         call MPI_Recv(myghostcol, local_height, MPI_INTEGER1, myrank-1, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(1:local_height,0)=myghostcol
      else ! all other ranks send first, then recv
         myghostcol=mygrid(1:local_height,local_width)
         call MPI_Send(myghostcol, local_height, MPI_INTEGER1, myrank+1, 0, MPI_COMM_WORLD, ierr)
         call MPI_Recv(myghostcol, local_height, MPI_INTEGER1, myrank-1, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(1:local_height,0)=myghostcol
      end if

      ! next send leftmost rows to the left
      if (myrank==0) then ! rank 0 only recvs
         call MPI_Recv(myghostcol, local_height, MPI_INTEGER1, myrank+1, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(1:local_height,local_width+1)=myghostcol
      else if (myrank==np-1) then ! rank np-1 only sends
         myghostcol=mygrid(1:local_height,1)
         call MPI_Send(myghostcol, local_height, MPI_INTEGER1, myrank-1, 0, MPI_COMM_WORLD, ierr)
      else ! all other ranks send first, then recv
         myghostcol=mygrid(1:local_height,1)
         call MPI_Send(myghostcol, local_height, MPI_INTEGER1, myrank-1, 0, MPI_COMM_WORLD, ierr)
         call MPI_Recv(myghostcol, local_height, MPI_INTEGER1, myrank+1, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(1:local_height,local_width+1)=myghostcol
      end if

   end subroutine ghostpoints_blockcol


   subroutine ghostpoints_checker(mygrid)

      implicit none
      integer*1, intent(inout) :: mygrid(0:local_height+1,0:local_width+1)
      integer*1 :: myghostcol(local_height), myxrow(0:local_width+1)

      ! step one: move internal columns across the rows of procs

      if (mod(myrank,sqrtnp)==0) then ! only send
         myghostcol=mygrid(1:local_height,local_width)
         call MPI_Send(myghostcol, local_height, MPI_INTEGER1, myrank+1, 0, MPI_COMM_WORLD, ierr)
      else if (mod(myrank,sqrtnp)==sqrtnp-1) then ! only recv
         call MPI_Recv(myghostcol, local_height, MPI_INTEGER1, myrank-1, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(1:local_height,0)=myghostcol
      else ! all other ranks send first, then recv
         myghostcol=mygrid(1:local_height,local_width)
         call MPI_Send(myghostcol, local_height, MPI_INTEGER1, myrank+1, 0, MPI_COMM_WORLD, ierr)
         call MPI_Recv(myghostcol, local_height, MPI_INTEGER1, myrank-1, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(1:local_height,0)=myghostcol
      end if

      if (mod(myrank,sqrtnp)==0) then ! only recv
         call MPI_Recv(myghostcol, local_height, MPI_INTEGER1, myrank+1, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(1:local_height,local_width+1)=myghostcol
      else if (mod(myrank,sqrtnp)==sqrtnp-1) then ! only send
         myghostcol=mygrid(1:local_height,1)
         call MPI_Send(myghostcol, local_height, MPI_INTEGER1, myrank-1, 0, MPI_COMM_WORLD, ierr)
      else ! all other ranks send first, then recv
         myghostcol=mygrid(1:local_height,1)
         call MPI_Send(myghostcol, local_height, MPI_INTEGER1, myrank-1, 0, MPI_COMM_WORLD, ierr)
         call MPI_Recv(myghostcol, local_height, MPI_INTEGER1, myrank+1, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(1:local_height,local_width+1)=myghostcol
      end if

      ! step two: move extended rows (internal plus ghost points) across cols of procs

      if (myrank/sqrtnp==0) then ! only send
         myxrow=mygrid(local_height,:)
         call MPI_Send(myxrow, local_width+2, MPI_INTEGER1, myrank+sqrtnp, 0, MPI_COMM_WORLD, ierr)
      else if (myrank/sqrtnp==sqrtnp-1) then ! only recv
         call MPI_Recv(myxrow, local_width+2, MPI_INTEGER1, myrank-sqrtnp, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(0,:)=myxrow
      else ! all other ranks send first, then recv
         myxrow=mygrid(local_height,:)
         call MPI_Send(myxrow, local_width+2, MPI_INTEGER1, myrank+sqrtnp, 0, MPI_COMM_WORLD, ierr)
         call MPI_Recv(myxrow, local_width+2, MPI_INTEGER1, myrank-sqrtnp, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(0,:)=myxrow
      end if

      if (myrank/sqrtnp==0) then ! only recv
         call MPI_Recv(myxrow, local_width+2, MPI_INTEGER1, myrank+sqrtnp, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(local_height+1,:)=myxrow
      else if (myrank/sqrtnp==sqrtnp-1) then ! only send
         myxrow=mygrid(1,:)
         call MPI_Send(myxrow, local_width+2, MPI_INTEGER1, myrank-sqrtnp, 0, MPI_COMM_WORLD, ierr)
      else ! all other ranks send first, then recv
         myxrow=mygrid(1,:)
         call MPI_Send(myxrow, local_width+2, MPI_INTEGER1, myrank-sqrtnp, 0, MPI_COMM_WORLD, ierr)
         call MPI_Recv(myxrow, local_width+2, MPI_INTEGER1, myrank+sqrtnp, 0, MPI_COMM_WORLD, stat, ierr)
         mygrid(local_height+1,:)=myxrow
      end if

   end subroutine ghostpoints_checker

   subroutine parallel_io(filename)

      implicit none
      character(len=64), intent(in) :: filename

      ! local variables
      integer :: lun
      integer :: gsizes(2), distribs(2), dargs(2), psizes(2)
      integer :: filetype
      integer*1, allocatable :: buff(:)
      integer :: buffsize
      integer(kind=MPI_OFFSET_KIND) :: disp


      ! open file
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, lun, ierr)

      ! create filetype darray
      gsizes=(/ gheight, gwidth /)
      distribs=(/ MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK /)
      dargs=(/ MPI_DISTRIBUTE_DFLT_DARG, MPI_DISTRIBUTE_DFLT_DARG /)
      psizes=(/ nrows, ncols /)
      call MPI_Type_Create_darray(np, myrank, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_FORTRAN, MPI_INTEGER1, filetype, ierr)
      call MPI_Type_Commit(filetype, ierr)

      ! set view to filetype darray
      disp=myrank*local_height*local_width*1 ! one-byte ints
      call MPI_FILE_SET_VIEW(lun, disp, MPI_INTEGER1, filetype, 'native', MPI_INFO_NULL, ierr) 

      ! write to file
      buffsize=local_height*local_width
      allocate (buff(buffsize))
      do i=1,local_height
      do j=1,local_width
         buff((i-1)*local_height+j)=mygrid(i,j)
      end do
      end do

      ! call MPI_FILE_WRITE(lun, mygrid(1:local_height,1:local_width), buffsize, MPI_INTEGER1, MPI_STATUS_IGNORE, ierr) 
      call MPI_FILE_WRITE(lun, buff, buffsize, MPI_INTEGER1, MPI_STATUS_IGNORE, ierr) 

      ! close file
      call MPI_FILE_CLOSE(lun, ierr) 

   end subroutine parallel_io


end program life






















