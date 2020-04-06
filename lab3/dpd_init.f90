!=========================================
!                 DPD init
! Create initial topology and coord file
!           for DPD simulations
!               J.Ilnytskyi
!                 May 2007
!========================================

  !==============================================
  program dpd_init
  !==============================================

  implicit none

  !Fortran-90 data types (compiler dependent)
  integer, parameter :: si = 4
  integer, parameter :: dp = 8
  real(dp) :: pi = 3.141592653589793_dp

  !file units
  integer(si) :: u1=20, u2=21

  !system pars
  integer(si) :: ntop, nmol, nbea, nbon
  real(dp)    :: den, vol, kt 
  real(dp), dimension(3) :: box

  !arrays
  real(dp), allocatable, dimension(:,:)  :: r, v
  integer(si), allocatable, dimension(:) :: t
  real(dp), allocatable, dimension(:)    :: m
  integer(si), allocatable, dimension(:) :: b1, b2
  real(dp), allocatable, dimension(:)    :: b0, bk

  !chain parameters
  integer(si) :: len, n1, n2, type1, type2
  real(dp)    :: frac
  real(dp)    :: sep, mass
  real(dp)    :: bond0, bondk

  !running
  integer(si) :: idum, i, k, l
  integer(si) :: kbea, kbon
  integer(si) :: imol, ibea, ibon
  real(dp), dimension(3) :: r0, e0
  character(256) :: str

    !------------------------------------------------
    !          melt of diblock copolymers
    !------------------------------------------------

    !read command line argument(s)
    call getarg(1,str)
    if (trim(str)=="") call err ("Syntax: dpd_init <fA>")
    read (str,*) frac

    !selected params
    box(:) = 12.0_dp
    len = 10

    !fixed params
    den = 3.0_dp
    kt = 1.0_dp
    sep = 0.9_dp
    mass = 1.0_dp
    bond0 = 0.0_dp
    bondk = 4.0_dp
    type1 = 1
    type2 = 2

    !get nmol
    vol = product(box(:))
    nbea = den*vol
    nmol = nbea/len

    !get n1, n2
    n1 = nint(frac*real(len,dp))
    n2 = len - n1

    !allocate arrays
    kbea = len
    kbon = len-1
    nbea = nmol*kbea
    nbon = nmol*kbon
    allocate(r(3,nbea),v(3,nbea),t(nbea),m(nbea))
    allocate(b1(nbon),b2(nbon),b0(nbon),bk(nbon))

    !fill-in beads and bonds of all chains

    ibea = 0
    ibon = 0
    do imol = 1, nmol
      r0(1) = ran3(idum)*box(1)
      r0(2) = ran3(idum)*box(2)
      r0(3) = ran3(idum)*box(3)
      call rnd_unit_vec(e0)
      call mk_chain(r0, e0)
    enddo

    !write topol file
    write(str,fmt='(i8.8,"_diblock_fA=",f4.2,".top")') 0,frac 
    open (u1, file=str, status="unknown")
    write (u1, fmt='(3i8)') 1,nmol,nbea
    write (u1, fmt='(3i8)') nbon,0,0
    write (u1, fmt='(3i8)') 0,0,0
    write (u1,fmt='(a)') "-------------------------------------------"
    write (u1,fmt='(2i6)') 1, nmol
    write (u1,fmt='(4i6)') kbea, kbon, 0, 0
    do ibea = 1, kbea
      write(u1,fmt='(i6,i4,f12.3)') ibea,t(ibea),m(ibea)
    enddo
    do ibon = 1, kbon
      write(u1,fmt='(2i6,2f8.3)') b1(ibon),b2(ibon),b0(ibon),bk(ibon)
    enddo

    write (u1,fmt='(a)') "-------------------------------------------"
    close(u1)

    !write coord file
    write(str,fmt='(i8.8,"_diblock_fA=",f4.2,".coord")') 0,frac 
    open (u2, file=str, status="unknown")
    write (u2, fmt='(3i8)')      0, nbea, 0
    write (u2, fmt='(3e16.8)')   box(:)
    do ibea = 1, nbea
      write (u2, fmt='(2i6)')    ibea, t(ibea)
      write (u2, fmt='(3e16.8)') r(:,ibea)
      write (u2, fmt='(3e16.8)') v(:,ibea)
    enddo
    write (u2,fmt='(3(g14.6))')  box(:)
    write (u2,fmt='(3(g14.6))')  0.0_dp,0.0_dp,0.0_dp
    write (u2,fmt='(3(g14.6))')  0.0_dp,0.0_dp,0.0_dp
    close(u2)

    deallocate(r,v,t,m)
    deallocate(b1,b2,b0,bk)


contains

  !-------------------------------------------------------------
  subroutine mk_chain (r0, e0)
  !-------------------------------------------------------------

  implicit none
  real(dp),dimension(3) :: r0, e0  !1st bead coords, chain orient.
  integer(si) :: i

  do i = 1, len
    ibea = ibea + 1
    r(:,ibea) = r0(:) + (i-1)*sep*e0(:)
    r(:,ibea) = merge(r(:,ibea)+box(:),r(:,ibea),r(:,ibea)<0.0_dp)
    r(:,ibea) = merge(r(:,ibea)-box(:),r(:,ibea),r(:,ibea)>=box(:))
    call rnd_unit_vec(v(:,ibea))
    v(:,ibea) = sqrt(3.0_dp*kt/mass)*v(:,ibea)
    if (i<=n1) then
      t(ibea) = type1
    else
      t(ibea) = type2
    endif
    m(ibea)   = mass
    if (i>1) then
      ibon = ibon + 1
      b1(ibon) = ibea-1
      b2(ibon) = ibea
      b0(ibon) = bond0
      bk(ibon) = bondk
    endif      
  enddo

  end subroutine mk_chain


  !---------------------------------------------------------
  subroutine rnd_unit_vec (vec)
  !---------------------------------------------------------
  !generate vector uniformly distributed on a sphere by rejection method

    implicit none
    real(dp),dimension(3) :: vec
    real(dp) :: len

    do
      vec(1) = 2.0_dp*(ran3(idum)-0.5_dp)
      vec(2) = 2.0_dp*(ran3(idum)-0.5_dp)
      vec(3) = 2.0_dp*(ran3(idum)-0.5_dp)
      len = sum(vec(:)*vec(:))
      if (len < 1.0_dp) exit
    enddo
    len = sqrt(len)
    vec(:) = vec(:)/len

  end subroutine rnd_unit_vec


  !---------------------------------------------------------
  function ran3(idum)
  !---------------------------------------------------------

    integer(si) :: idum
    integer(si) :: mbig,mseed,mz
    real(dp)    :: ran3,fac
    parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.d0/mbig)
    integer(si) :: i,iff,ii,inext,inextp,k
    integer(si) :: mj,mk,ma(55)
    save iff,inext,inextp,ma
    data iff /0/

    if (idum.lt.0.or.iff.eq.0) then
      iff = 1
      mj = mseed-iabs(idum)
      mj = mod(mj,mbig)
      ma(55) = mj
      mk = 1
      do i = 1, 54
        ii = mod(21*i,55)
        ma(ii) = mk
        mk = mj-mk
        if (mk.lt.mz) mk = mk+mbig
        mj = ma(ii)
      enddo
     do k = 1, 4
        do i = 1, 55
         ma(i) = ma(i)-ma(1+mod(i+30,55))
          if (ma(i).lt.mz) ma(i) = ma(i)+mbig
        enddo
      enddo
      inext = 0
      inextp = 31
      idum = 1
    endif
    inext = inext+1
    if (inext.eq.56) inext = 1
    inextp = inextp+1
    if (inextp.eq.56) inextp=1
    mj = ma(inext)-ma(inextp)
    if (mj.lt.mz) mj = mj+mbig
    ma (inext) = mj
    ran3 = mj*fac

  end function ran3

  !-------------------------------------------------
  subroutine err (str)
  !-------------------------------------------------

    character(len=*) :: str

    write (*,*) str
    stop

  end subroutine err


  !==============================================
  end program dpd_init
  !==============================================
