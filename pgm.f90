module pgm

contains

   subroutine read_pgm(filename, nrows, ncols, my_col, my_row, &
      height, width, local_height, local_width, mygrid)
      !
      ! nrows, ncols : how to divide up the board
      ! my_col, my_row : the col number and row num I start at
      ! local_height, local_width : dimensions of the game board (minus 
      !                             ghost cells)
      ! field_a, field_b : two copies of the game tile (plus ghost rows),
      !                    one for reading and one for updating

      implicit none

      ! arguments
      character (len=*), intent(in) :: filename
      integer, intent(in):: nrows, ncols, my_col, my_row
      integer, intent(out):: local_height, local_width
      integer*1, intent(out), allocatable :: mygrid(:,:)
      integer, intent(out) :: height, width

      ! local
      integer :: i, j
      character (len=64) :: the_header
      character (len=64) :: aheight,awidth, adepth, foo
      character :: bar(20)
      character :: last, mychar
      integer :: depth
      integer :: ierror, newval
      integer :: start_x, start_y
      integer :: lx, ly

      ! read the header info
      call cgetheader(filename, width, height, depth)

      if ( MOD(width, ncols) /= 0 ) then
         print*, 'width=', width
         print*, 'ncols=', ncols
         stop 'Width does not divide evenly.'
      end if

      if ( MOD(height, nrows) /= 0) then
         print*, 'height=', height
         print*, 'nrows=', nrows
         stop 'Height does not divide evenly.'
      end if

      ! Divide the total image among the local processors
      local_width = width / ncols
      local_height = height / nrows

      ! Find out where my starting range is
      ! my_col should range from 1 to np (not 0 to np-1) for col decomp
      start_x = local_width*(my_col-1)+1
      start_y = local_height*(my_row-1)+1

      ! allocate data with room for ghost points
      ! these arrays go from index 0 to n+1
      ! indices 1 through n are the local set
      ! indices 0 and n+1 are the ghost points
      allocate (mygrid(0:local_height+1,0:local_width+1))
      mygrid=0 ! maybe not necessary, want make sure true boundaries set to zero

      ! loop through all points
      do j=1,height
      do i=1,width

         ! get the next value from the file
         ! black cells (b=0) are bugs, all else background
         call cfgetc(newval)
         if (newval .eq. 0) then
            newval = 1
         else
            newval = 0
         end if

         ! If the character is local, then save it!
         if ((i .ge. start_x-1) .AND. (i .le. (start_x + local_width)) .AND. &
            (j .ge. start_y-1) .AND. (j .le. (start_y + local_height)) ) then

            ! Calculate the local pixels (+1 for ghost row,col)
            !  -> this will work so long as field_a,b(0:local_height+1,0:local_width+1)
            lx = i - start_x + 1
            ly = j - start_y + 1
            mygrid(ly,lx) = newval

         end if
       
      enddo
      enddo

      call cfclose()

   end subroutine read_pgm

   subroutine pval_read_pgm(filename, nrows, ncols, my_col, my_row, &
      height, width, local_height, local_width, mygrid)
      ! same as read_pgm but preserves the value of each cell, instead of projecting into 0 or 1

      implicit none

      ! arguments
      character (len=*), intent(in) :: filename
      integer, intent(in):: nrows, ncols, my_col, my_row
      integer, intent(out):: local_height, local_width
      integer*1, intent(out), allocatable :: mygrid(:,:)
      integer, intent(out) :: height, width

      ! local
      integer :: i, j
      character (len=64) :: the_header
      character (len=64) :: aheight,awidth, adepth, foo
      character :: bar(20)
      character :: last, mychar
      integer :: depth
      integer :: ierror, newval
      integer :: start_x, start_y
      integer :: lx, ly

      ! read the header info
      call cgetheader(filename, width, height, depth)

      if ( MOD(width, ncols) /= 0 ) then
         print*, 'width=', width
         print*, 'ncols=', ncols
         stop 'Width does not divide evenly.'
      end if

      if ( MOD(height, nrows) /= 0) then
         print*, 'height=', height
         print*, 'nrows=', nrows
         stop 'Height does not divide evenly.'
      end if

      ! Divide the total image among the local processors
      local_width = width / ncols
      local_height = height / nrows

      ! Find out where my starting range is
      ! my_col should range from 1 to np (not 0 to np-1) for col decomp
      start_x = local_width*(my_col-1)+1
      start_y = local_height*(my_row-1)+1

      ! allocate data with room for ghost points
      ! these arrays go from index 0 to n+1
      ! indices 1 through n are the local set
      ! indices 0 and n+1 are the ghost points
      allocate (mygrid(0:local_height+1,0:local_width+1))
      mygrid=0 ! maybe not necessary, want make sure true boundaries set to zero

      ! loop through all points
      do j=1,height
      do i=1,width

         ! get the next value from the file
         ! black cells (b=0) are bugs, all else background
         call cfgetc(newval)

         ! If the character is local, then save it!
         if ((i .ge. start_x-1) .AND. (i .le. (start_x + local_width)) .AND. &
            (j .ge. start_y-1) .AND. (j .le. (start_y + local_height)) ) then

            ! Calculate the local pixels (+1 for ghost row,col)
            !  -> this will work so long as field_a,b(0:local_height+1,0:local_width+1)
            lx = i - start_x + 1
            ly = j - start_y + 1
            mygrid(ly,lx) = newval

         end if
       
      enddo
      enddo

      call cfclose()

   end subroutine pval_read_pgm

end module pgm
