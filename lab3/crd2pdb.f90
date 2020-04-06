program crd2pdb

  !============================================================
  !  produces pdb file for visualisation from dpd coord file
  !                      J.Ilnytskyi
  !                        May 2007
  !============================================================

  implicit none
  integer, parameter :: si = 4
  integer, parameter :: dp = 8

  !filenames
  character(256) :: inp_coord_fn
  character(256) :: inp_pdb_fn
  character(256) :: buf, buf1, buf2

  !molecule info
  real(dp) :: scale

  !dimensions
  integer(si) :: st_n, nbea, dummy
  real(dp), dimension(3) :: box
  integer(si) :: ibea, type
  real(dp), dimension(3) :: r, v, e, u

  !constants and vars
  integer(si) :: u1=10,u2=11, i, imol, ios

  !=======================================================================
  !get cmd line arguments
  call getarg(1,inp_coord_fn)
  call getarg(2,inp_pdb_fn)
  if (trim(inp_coord_fn)=="".or.trim(inp_pdb_fn)=="")&
    call err ("Syntax: crd2pdb <file>.coord <file>.pdb")

  !read inp.coord file
  open (u1,file=inp_coord_fn,status="old",iostat=ios)
    if (ios/=0) call err ("Error: No inp.coord file")
  read (u1,fmt=*,iostat=ios) st_n, nbea, dummy
    if (ios/=0) call err ("Error in inp.coord: st_n, nbea, ngb")
  read (u1,fmt=*,iostat=ios) box(:)
    if (ios/=0) call err ("Error in inp.coord: box")

  scale = 5.0_dp
  !load original coords and write colour coded pdb file
  open (u2,file=inp_pdb_fn,status="unknown",iostat=ios)
  do ibea = 1, nbea
    read (u1,fmt=*,iostat=ios) dummy, type
    read (u1,fmt=*,iostat=ios) r(:)
    read (u1,fmt=*,iostat=ios) v(:)
    if (type==1) then
      write(u2,fmt='("ATOM  ",i5," ",a4," NONE   ",i2,"    ",3f8.3,f6.2,f6.2)')&
            ibea,'N   ',1,scale*r(:), 0.0,0.0
    elseif (type==2) then
      write(u2,fmt='("ATOM  ",i5," ",a4," NONE   ",i2,"    ",3f8.3,f6.2,f6.2)')&
            ibea,'S   ',2,scale*r(:), 0.0,0.0
    endif
  enddo

  !colour records
  write(u2,fmt='(a)') 'COLO  ################## 1####   0.000   0.500   0.500  3.00'
  write(u2,fmt='(a)') 'COLO  ################## 2####   1.000   0.000   0.000  3.00'
  write(u2,fmt='(a)') 'COLO  ################## 3####   0.000   1.000   0.000  3.00'
  write(u2,fmt='(a)') 'COLO  ################## 4####   0.500   0.500   0.500  3.00'

  close (u1)
  close (u2)

contains

  !-------------------------------------------------
  subroutine err (str)
  !-------------------------------------------------

    character(len=*) :: str

    write (*,*) str
    stop

  end subroutine err

end program crd2pdb
