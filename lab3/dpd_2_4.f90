!=====================================================================================
!                DPD_2_4 PROGRAM
!
! new in flags for:
! attraction, reactivity, grafts, flat and cylindrical surfaces
!
!                 J.Ilnytskyi
!                 July, 2011
!
!* NVT, NPT, NPxPyPzT
!* spherical DPD beads
!* repulsive and attractive interactions
!* surfaces along OZ
!* extended for the case of blank topology
!*       changes made into 2.2 version:
!* found bug in eval_bond and eval_grf due to not accounting for PBC in X,Y
!* account for walls reactive force and distribute it over the system
!* dynamic crosslinking:
!* (includes corrected issues with "CNT" due to clearing reactive flag, a(:) array
!  and with bond arrays allocation for more than one topology case from dpd_2_2c)
!  two sorts of reactiveness for reactive beads, s(:) array
!* cylindrical cavities made as inpenetrable curved walls
!* unified subroutines are introduced:
!  eval_cdr_uni, eval_attr_uni, where surface imaginary bead is jbea=0 of type=0
!  arrays r,v,t,for are extended to contain 0 element
!  cbonmax cangmax ctormax are employed to define dimensions for crosslinked arrays
!* grafting interactions change: save/read 0,1 (bottom ot top wall) for z coordinate
!* unified subroutine is introduced:
!  eval_bond_uni
!* unification of yes/no surface and yes/no cylinders cases
!* corrected linked cells yes/no PBC in Z by adding if (surface)...
!  ====================================================================================
!  bead type coding:
!  XXXX XXXX   XXXX XXXX
!  ||||    |      | ||||
!  ||||    |      | tttt - bead type for repulsion  [iand(000F,t)]
!  ||||    |      e ------ special (chain end, etc) [iand(0010,t)]
!  ||||    a-------------- attraction bead          [iand(0100,t)]  
!  |||c------------------- reactive bead            [iand(1000,t)]
!  ddd-------------------- reactiveness sort        [ishift(iand(E000,t),-13)]
!  ================================================
!  quick formation of bead type:
!  reactive type 1:  3000x  or
!  reactive type 2:  5000x
!       add to
!  repulsive type 1: 0001x  or
!  repulsive type 2: 0010x
!  ====================================================================================
!  reactive atoms coding:
!  upon reading topology file:
!    initially:            a(ibea) = 0 (non-reactive)
!    if iand(1000,t)/=0:   a(ibea) = 1 (reactive)
!    if in crosslink list: a(ibea) = 2 (reactive but already crosslinked)
!  upon succesfull crosslinking:
!    a(ibea) = 2 (crosslinked) 
!  upon writing to topology file:
!    initial type t(ibea) is written for each bead (reflecting initial ability to crosslink)
!    crosslinked ones can be retrieved from the crosslinks list
!  upon writing to coords file:
!    initial type t(ibea) is written for each bead (the same as in topology)
!!=======================================================================================


  !==============================================
  program dpd_2_4
  !==============================================

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXX SHARED VARIABLES XXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  implicit none

  !Fortran-90 data types (compiler dependent)
  integer, parameter :: si = 4
  integer, parameter :: dp = 8

  !inp and out filenames
  integer(si), parameter :: fnlen = 64
  real(dp),    parameter :: pi = 3.141592654_dp
  character(fnlen), parameter :: fn_input = "dpd.inp"
  character(fnlen), parameter :: fn_grf   = "grf.inp"
  character(fnlen), parameter :: fn_cyl   = "cyl.inp"
  character(fnlen), parameter :: fn_rea   = "react.inp"
  character(fnlen), parameter :: fn_list  = "dpd.lst"
  character(fnlen), parameter :: fn_mom   = "mom.lst"

  !bead coding masks
  integer(si) :: btype=15    !000FH
  integer(si) :: bspec=16    !0010H
  integer(si) :: battr=256   !0100H
  integer(si) :: breac=4096  !1000H
  integer(si) :: brsrt=57344 !E000H
  integer(si) :: brshf=-13

  !max functionality (# reactive beads per molecule)
  integer(si) :: maxfunct=3

  !read-in run parameters
  !options strings
  character(3) :: ropt         !{"NEW","CNT"}
  character(4) :: topt         !{"TDPD","TDSE"}
  character(4) :: popt         !{"NONE", "BARO","ABAR","BRSE","ABSE", "BERE","ABER","BESE","AESE"}
  integer(si)  :: st_nmx       !max MD step
  !temperature control data
  real(dp)     :: ktfix, gamma, vvext
  real(dp)     :: ktfixs, ktfixe
  integer(si)  :: ktfixn
  !pressure control data
  character(3) :: isoani
  real(dp)     :: pfix, ptau
  real(dp)     :: pfixs, pfixe
  integer(si)  :: pfixn
  real(dp), dimension(3) :: pxyz, pxyzs, pxyze
  !filenames
  character(fnlen) :: fn_tread, fn_twrit, fn_cread, fn_cwrit
  !import and b/a/t factors
  integer(si) :: cbonmax, cangmax, ctormax
  !frequencies for list and coords writing
  integer(si) :: st_lst, st_cwrit
  !timestep and cutoff
  real(dp) :: dt, rc
  !soft repulsion params parameters
  real(dp), dimension (3,0:3) :: prep
  !attractive force range and cutoff
  real(dp) :: xi, qc
  !additional attractive part 
  real(dp), dimension (3,0:3) :: peps
  !attractive and reactive flags
  character(3) :: aflg, rflg  !{"ATR","NON"}, {"REA","NON"}  
  !surface, grafts and cylinders flags
  character(3) :: sflg, gflg, cflg  !{"SRF","NON"}, {"GRF","NON"}, {"CYL","NON"}    

  !reactive data
  character(3) :: rtype        !{"NON","SQW",..}
  integer(si) :: rfreq
  real(dp) :: rprob, rb0, rbk, cc

  !flags
  logical :: barostat   !true if popt = "BARO","ABAR","BRSE","ABSE"
  logical :: berendsen  !true if popt = "BERE","ABER","BESE","AESE"
  logical :: attraction !true if aflg = "ATR" and there are beads with attractive bit set
  logical :: reactive   !true if rflg = "REA" and there are beads with reactive bit set, fn_rea file is required
  logical :: surface    !true if sflg = "SRF"
  logical :: grafts     !true if gflg = "GRF", fn_grf file is required
  logical :: cylinders  !true if cflg = "CYL", fn_cyl file is required

  !derivative values
  integer(si) :: dtop  !=2 for surface, =3 for no surface
  real(dp) :: sigma, rc_2
  real(dp) :: qc_2, twodxi
  real(dp) :: cc_2

  !# of topol, mols, beads, bonds, angles, torsions ----
  integer(si) :: ntop, nmol, nbea, nbon, nang, ntor  !current #
  integer(si) ::       cbon, cang, ctor              !extra (crosslinks)
  integer(si) :: nbontop, nangtop, ntortop           !allocated top #
  integer(si) :: ngrf                                !grafts #
  integer(si) :: ncyl                                !cylinders #

  !Arrays
  !molecule data  
  integer(si),allocatable,dimension (:) :: top       !mol host topol
  integer(si),allocatable,dimension (:) :: nfir,nlst !mol first/last atom
  integer(si),allocatable,dimension (:) :: bfir,blst !mol first/last bond
  integer(si),allocatable,dimension (:) :: afir,alst !mol first/last angle
  integer(si),allocatable,dimension (:) :: tfir,tlst !mol first/last torsion
  !bead data
  integer(si),allocatable,dimension (:) :: mol     !bead host molecule
  integer(si),allocatable,dimension (:) :: t,a,s   !bead type
  real(dp), allocatable,dimension (:,:) :: r,v,vsave !bead coord, veloc, saved veloc
  real(dp), allocatable, dimension  (:) :: m       !bead mass
  real(dp), allocatable,dimension (:,:) :: for     !bead force
  !bonds data
  integer(si),allocatable,dimension (:) :: b1, b2
  real(dp), allocatable, dimension  (:) :: b0, bk
  !angles data
  integer(si),allocatable,dimension (:) :: a1, a2, a3
  real(dp), allocatable, dimension  (:) :: a0, ak
  !torsions
  integer(si),allocatable,dimension (:) :: t1, t2, t3, t4
  real(dp), allocatable,dimension (:,:) :: tc
  !grafts
  integer(si),allocatable,dimension (:) :: grft, grfi
  real(dp), allocatable,dimension (:,:) :: grfr
  integer(si), allocatable,dimension(:) :: grfw
  real(dp),allocatable,dimension    (:) :: grfb0, grfbk
  !cylinders
  real(dp), allocatable,dimension (:,:) :: rcylc
  real(dp), allocatable,dimension (:)   :: cylr

  !instant thermodynamics
  integer(si) :: st_n0, st_n  !DPD start and current steps
  integer(si) :: idum         !ran3 integer
  real(dp) :: mass, vol, dens, kin, kt, pres
  real(dp), dimension (3) :: box, kinv, virv, presv
  real(dp), dimension (3) :: boxpv, xippv, xipv
  real(dp), dimension (3) :: smom, tmom   !surface and total momenta

  !list arrays
  integer(si) :: lst_dim, lst_ptr
  real(dp),allocatable,dimension(:) :: &
    lst_dens, lst_kt,   lst_pres, &
    lst_boxx, lst_boxy, lst_boxz, &
    lst_pxx,  lst_pyy,  lst_pzz,  &
    lst_momx, lst_momy, lst_momz

  !work
  integer(si) :: i
  character(256) :: str

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXX MAIN ROUTINE XXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !read data
  st_n = 0
  call read_run_pars
  call read_topol
  call read_coords
  call remove_com
  if (reactive)  call read_react
  if (grafts)    call read_grf
  if (cylinders) call read_cyl

  !initial evaluations 
  call eval_for_vir
  call eval_kin_kt_pres

  !perform MD loop
  do st_n = st_n0, st_nmx

    !set kT, P for slow control options 
    call set_ktfix_pfix

    !integrate eqs of motion
    if (barostat .or. berendsen) then
      call integr_vv_npt
    else
      call integr_vv_nvt
    endif

    !write(*,*) st_n, kt, pres

    !store/save data
    if (mod (st_n,st_lst)==0)   call add_to_lst
    if (mod (st_n,st_cwrit)==0) then
      call flush_lst
      call write_topol
      call write_coords
      call remove_com
    endif

  enddo

contains


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXX SUBROUTINES XXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !==============================================
  !================== READ-IN ===================
  !==============================================

  !----------------------------------------------
  subroutine read_run_pars
  !----------------------------------------------
  ! read dpd.inp file
  ! check for absence of outout files for new run
  ! allocate lists arrays
  !----------------------------------------------
  !format of dpd.inp file
  !"job title"
  !{"NEW","CNT"} {"TDPD","TDSE"} {"NONE","BARO","ABAR","BRSE","ABSE",
  !                                      "BERE","ABER","BESE","AESE"} st_nmx
  !{ 
  !  {ktfix gamma vvext},                        //for TDPD
  !  {ktfixs ktfixe ktfixn gamma vvext}          //for TDSE
  !{ {-no line present-}                         //for NONE
  !  {"ISO","ANI"} pfix ptau,                    //for BARO, BERE
  !  {"ISO","ANI"} pxyz(:) ptau,                 //for ABAR, ABER
  !  {"ISO","ANI"} pfixs pfixe pfixn ptau,       //for BRSE, BESE
  !  {"ISO","ANI"} pxyzs(:) pxyze(:) pfixn ptau} //for ABSE, AESE
  !  <read_topol_filename>
  !  <read_coord_filename>
  !  <write_topol_filename>
  !  <write_coord_filename>
  !  bonfct angfct torfct
  !  st_lst  st_cwrit
  !  dt   rc                          //timestep, repulsion cutoff
  !  prep(1,1) prep(1,2) prep(1,3)    //pairwise repulsion parameters matrix
  !  prep(2,1) prep(2,2) prep(2,3)
  !  prep(3,1) prep(3,2) prep(3,3)
  !  prep(1,0) prep(2,0) prep(3,0)     //surface repulsion parameters array
  !  xi qc                             //pairwise attraction range, cutoff
  !  peps(1,1) peps(1,2) peps(1,3)     //pairwise attraction epsilon matrix
  !  peps(2,1) peps(2,2) peps(2,3)
  !  peps(3,1) peps(3,2) peps(3,3)
  !  peps(1,0) peps(2,0) peps(3,0)     //surface attraction epsilon array          
  !  {"ATR", "NON"} {"REA","NON"}      //attracion and reactivness switches
  !  {"SRF","NON"} {"GRF","NON"} {"CYL","NON"} //surface, grafts, cylinders switches
  !----------------------------------------------
 
    implicit none
    integer(si) :: inp_u=10, lst_u=21
    integer(si) :: ios
    character(fnlen) :: str1, str2, str3, str4, str5
    character(256) :: title

    open (inp_u, file=fn_input, status="old", iostat=ios)
    if (ios /= 0) call error_exit ("Error: No input file")

    !read header
    read (inp_u, fmt=*, iostat=ios) title
    if (ios/=0) call error_exit ("Error: Invalid inp file, header")

    !read options line
    read (inp_u, fmt=*, iostat=ios) str1, str2, str3, st_nmx
    if (ios/=0) call error_exit ("Error: Invalid inp file, options line")
    ropt   = str1(1:3)
    topt   = str2(1:4)
    popt   = str3(1:4)
    !check for valid syntax
    if ((ropt/="NEW".and.ropt/="CNT").or.(topt/="TDPD".and.topt/="TDSE").or.&
        (popt/="NONE".and.&
         popt/="BARO".and.popt/="ABAR".and.popt/="BRSE".and.popt/="ABSE".and.&
         popt/="BERE".and.popt/="ABER".and.popt/="BESE".and.popt/="AESE").or.&
         st_nmx<1) then
      call error_exit ("Error: Invalid inp file, options line")
    endif

    !read temperature control line
    if (topt=="TDPD") then
      read (inp_u, fmt=*, iostat=ios) ktfix, gamma, vvext
      if (ios/=0.or.ktfix<0.0) call error_exit &
        ("Error: Invalid inp file, dpd temperature parameters control line")
      sigma = sqrt(2.0_dp*gamma*ktfix)
    elseif (topt=="TDSE") then
      read (inp_u, fmt=*, iostat=ios) ktfixs, ktfixe, ktfixn, gamma, vvext
      if (ios/=0.or.ktfix<0.0) call error_exit &
        ("Error: Invalid inp file, dpd temperature parameters control line")
      sigma = sqrt(2.0_dp*gamma*ktfix)
    endif

    !read pressure control line, if any
    if (popt=="BARO".or.popt=="BERE") then
      read (inp_u, fmt=*, iostat=ios) str1, pfix, ptau
      if (ios/=0) call error_exit ("Error: Invalid inp file, pressure control line")
      isoani = str1(1:3)
      pxyz(:) = pfix/3.0_dp
      if ((isoani/="ISO".and.isoani/="ANI").or.ptau<0.0) call error_exit &
        ("Error: Invalid inp file, pressure control line")
    elseif (popt=="ABAR".or.popt=="ABER") then
      read (inp_u, fmt=*, iostat=ios) str1, pxyz(1:3), ptau
      if (ios/=0) call error_exit ("Error: Invalid inp file, pressure control line")
      isoani = str1(1:3)
      if ((isoani/="ISO".and.isoani/="ANI").or.ptau<0.0) call error_exit &
        ("Error: Invalid inp file, pressure control line")
    elseif (popt=="BRSE".or.popt=="BESE") then
      read (inp_u, fmt=*, iostat=ios) str1, pfixs, pfixe, pfixn, ptau
      if (ios/=0) call error_exit ("Error: Invalid inp file, pressure control line")
      isoani = str1(1:3)
      pxyzs(:) = pfixs/3.0_dp
      pxyze(:) = pfixe/3.0_dp
      if ((isoani/="ISO".and.isoani/="ANI").or.ptau<0.0) call error_exit &
        ("Error: Invalid inp file, pressure control line")
    elseif (popt=="ABSE".or.popt=="AESE") then
      read (inp_u, fmt=*, iostat=ios) str1, pxyzs(1:3), pxyze(1:3), pfixn, ptau
      if (ios/=0) call error_exit ("Error: Invalid inp file, pressure control line")
      isoani = str1(1:3)
      if ((isoani/="ISO".and.isoani/="ANI").or.ptau<0.0) call error_exit &
        ("Error: Invalid inp file, pressure control line")
    endif
    barostat  = (popt=="BARO".or.popt=="ABAR".or.popt=="BRSE".or.popt=="ABSE")
    berendsen = (popt=="BERE".or.popt=="ABER".or.popt=="BESE".or.popt=="AESE")

    !read filenames
    read (inp_u, fmt=*, iostat=ios) fn_tread
    if (ios/=0) call error_exit &
      ("Error: missing/invalid input topology filename")
    read (inp_u, fmt=*, iostat=ios) fn_cread
    if (ios/=0) call error_exit &
      ("Error: missing/invalid input coord file filename")
    read (inp_u, fmt=*, iostat=ios) fn_twrit
    if (ios/=0) call error_exit &
      ("Error: missing/invalid output topology base filename")
    read (inp_u, fmt=*, iostat=ios) fn_cwrit
    if (ios/=0) call error_exit &
      ("Error: missing/invalid output coord base filename")

    !read max # of crosslinked bonds, angles, torsions
    read (inp_u, fmt=*, iostat=ios) cbonmax, cangmax, ctormax
    if (ios/=0) call error_exit &
       ("Error: Invalid inp file: cbonmx cangmax ctormax")

    !read step flags
    read (inp_u, fmt=*, iostat=ios) st_lst, st_cwrit
    if (st_lst==0)   st_lst   = 100000000
    if (st_cwrit==0) st_cwrit = 100000000
    if (mod(st_cwrit,st_lst)/=0 .or. ios/=0) call error_exit &
      ("Error: Invalid input file: st_lst should be divider of st_cwrit")

    !read timestep and cutoff
    read (inp_u, fmt=*, iostat=ios) dt, rc
    if (ios/=0) call error_exit ("Error: Invalid inp file: dt, rcut")
    if (dt < 0.0 .or. rc < 0.0) call error_exit ("Error: Invalid dt or rc")
    rc_2 = rc**2

    !read repulsion parameters matrix
    read (inp_u, fmt=*, iostat=ios) prep(1,1:3)
    read (inp_u, fmt=*, iostat=ios) prep(2,1:3)
    read (inp_u, fmt=*, iostat=ios) prep(3,1:3)
    if (prep(1,2)/=prep(2,1).or.prep(2,3)/=prep(3,2).or.prep(3,1)/=prep(1,3))&
      call error_exit ("Error: Invalid input file: prep(:,:) is non-symmetric")
    read (inp_u, fmt=*, iostat=ios) prep(1:3,0)

    !read attractive force range and cutoff
    read (inp_u, fmt=*, iostat=ios) xi, qc
    if (ios/=0) call error_exit ("Error: Invalid inp file: xi, qc")
    if (xi <= 0.0 .or. qc <= 0.0) call error_exit ("Error: Invalid xi or qc")
    qc_2 = qc**2
    twodxi = 2.0_dp/xi

    !read attractive epsilon parameters matrix
    read (inp_u, fmt=*, iostat=ios) peps(1,1:3)
    read (inp_u, fmt=*, iostat=ios) peps(2,1:3)
    read (inp_u, fmt=*, iostat=ios) peps(3,1:3)
    if (peps(1,2)/=peps(2,1).or.peps(2,3)/=peps(3,2).or.peps(3,1)/=peps(1,3))&
      call error_exit ("Error: Invalid input file: peps(:,:) is non-symmetric")
    read (inp_u, fmt=*, iostat=ios) peps(1:3,0)

    !read flags strings, set flags
    read (inp_u, fmt=*, iostat=ios) str1, str2
    if (ios/=0) call error_exit ("Error: Invalid inp file, flags: aflg, rflg")
    read (inp_u, fmt=*, iostat=ios) str3, str4, str5
    if (ios/=0) call error_exit ("Error: Invalid inp file, flags: sflg, gflg, cflg")
    aflg   = str1(1:3)
    rflg   = str2(1:3)
    sflg   = str3(1:3)
    gflg   = str4(1:3)
    cflg   = str5(1:3)
    if ((aflg/="ATR".and.aflg/="NON").or.(rflg/="REA".and.rflg/="NON").or.&
        (sflg/="SRF".and.sflg/="NON").or.(gflg/="GRF".and.gflg/="NON").or.(cflg/="CYL".and.cflg/="NON")) then
      call error_exit ("Error: Invalid inp file, flags: aflg, rflg, sflg, gflg or cflg")
    endif
    attraction = merge(.true.,.false.,aflg=="ATR")  !may be set to false if no beads of this kind
    reactive   = merge(.true.,.false.,rflg=="REA")  !may be set to false if no beads of this kind
    surface    = merge(.true.,.false.,sflg=="SRF")
    grafts     = merge(.true.,.false.,gflg=="GRF")
    cylinders  = merge(.true.,.false.,cflg=="CYL")

    close (inp_u)

    !set dtop for periodic boundary conditions
    dtop = merge(2,3,surface)

    !prevent accidental rewriting of list file
    if (ropt=="NEW") then
      open (lst_u, file=fn_list, status="old", iostat=ios)
      if (ios == 0) then
        close(lst_u)
        call error_exit ("Error: New run but dpd.lst already exists")
      endif
    endif

    !allocate list buffers
    lst_dim = st_cwrit/st_lst
    allocate(lst_dens(lst_dim),lst_kt(lst_dim),lst_pres(lst_dim),&
             lst_boxx(lst_dim),lst_boxy(lst_dim),lst_boxz(lst_dim),&
             lst_pxx(lst_dim),lst_pyy(lst_dim),lst_pzz(lst_dim),&
	     lst_momx(lst_dim),lst_momy(lst_dim),lst_momz(lst_dim),stat=ios)
    if (ios /= 0) call error_exit("Error: cannot allocate list arrays")
    !init list pointer
    lst_ptr = 0

  end subroutine read_run_pars


  !----------------------------------------------
  subroutine read_topol
  !----------------------------------------------
  !format of topol file
  !  ntop nmol nbea                     //total # in the melt
  !  nbon nang ntor                     //total regular bonding
  !  cbon cang ctor                     //total crosslinked bonding
  !{
  !  ----------------------------
  !  itop kmol                          //topol index
  !  kbea kbon kang ktor                //# of mols and one mol #
  !  [i t(i) m(i)]                      //beads types and masses
  !  [b1(i) b2(i) b0(i) bk(i)]          //regular bonds
  !  [a1(i) a2(i) a3(i) a0(i) ak(i)     //regular angles
  !  [t1(i) t2(i) t3(i) t4(i) tc(:,i)]  //regular torsions
  !}
  !  ----------------------------
  !  [b1(i) b2(i) b0(i) bk(i)]          //crosslinked bonds
  !  [a1(i) a2(i) a3(i) a0(i) ak(i)     //crosslinked angles
  !  [t1(i) t2(i) t3(i) t4(i) tc(:,i)]  //crosslinked torsions
  !----------------------------------------------
  

    implicit none
    integer(si) :: top_u = 11
    integer(si) :: itop, i, j, ios
    integer(si) :: i1, i2, i3, i4, ofs
    integer(si) :: imol,ibea,ibon,iang,itor
    integer(si) :: kmol,kbea,kbon,kang,ktor

    !open file
    open (top_u, file=fn_tread, status="old", iostat=ios)
    if (ios /= 0) call error_exit ("Error: no topology file")
    !read # of initial totals
    read (top_u, fmt=*, iostat=ios) ntop,nmol,nbea
    if (ios /= 0) call error_exit&
      ("Error: invalid topology: ntop,nmol,nbea")
    read (top_u, fmt=*, iostat=ios) nbon,nang,ntor
    if (ios /= 0) call error_exit&
      ("Error: invalid topology: nbon,nang,ntor")
    read (top_u, fmt=*, iostat=ios) cbon,cang,ctor
    if (ios /= 0) call error_exit&
      ("Error: invalid topology: cbon,cang,ctor")
    !allocate arrays
    !polymer bonds = nbon
    !crosslinks    = cbon
    !total possible crosslinks from maxfunct mols = (maxfunct*nmol)/2
    if (.not. reactive) then
      nbontop = nbon+cbon
      nangtop = nang+cang
      ntortop = ntor+ctor
    else
      nbontop = nbon+cbonmax
      nangtop = nang+cangmax
      ntortop = ntor+ctormax
    endif
    !allocate arrays
    allocate (top(nmol),nfir(nmol),nlst(nmol),bfir(nmol),blst(nmol),&
                        afir(nmol),alst(nmol),tfir(nmol),tlst(nmol),&
             r(3,0:nbea),v(3,0:nbea),t(0:nbea),for(3,0:nbea),&  !ibea=0 corresponds to surfaces
             a(nbea),s(nbea),m(nbea),mol(nbea),vsave(3,nbea),&
             b1(nbontop),b2(nbontop),b0(nbontop),bk(nbontop),&
             a1(nangtop),a2(nangtop),a3(nangtop),a0(nangtop),ak(nangtop),&
             t1(ntortop),t2(ntortop),t3(ntortop),t4(ntortop),tc(3,ntortop),stat=ios)
    if (ios /= 0) call error_exit ("Error: cannot allocate main set of arrays")

    !init pointers
    imol = 0
    ibea = 0
    ibon = 0
    iang = 0
    itor = 0
    !loop over topologies
    do itop = 1, ntop
      !read-in data for current topology
      read (top_u, fmt=*, iostat=ios) str
      read (top_u, fmt=*, iostat=ios) j, kmol
      if (j /= itop) call error_exit&
        ("Error: invalid topology: unsorted topologies")
      if (kmol==0) then
        read (top_u, fmt=*, iostat=ios) kbea,kbon,kang,ktor
        do i = 1, kbea
          read (top_u, fmt=*, iostat=ios) str
        enddo
        do i = 1, kbon
          read (top_u, fmt=*, iostat=ios) str
        enddo
        do i = 1, kang
          read (top_u, fmt=*, iostat=ios) str
        enddo
        do i = 1, ktor
          read (top_u, fmt=*, iostat=ios) str
        enddo
        cycle
      endif
      read (top_u, fmt=*, iostat=ios) kbea,kbon,kang,ktor
      if (ios /= 0) call error_exit&
        ("Error: invalid topology: kbea,kbon,kang,ktor")
      imol = imol + 1
      !fill-in data for the first molecule of current topology
      do i = 1, kbea
        ibea = ibea + 1
        read (top_u, fmt=*, iostat=ios) i1,t(ibea),m(ibea)
        if (ios /= 0) call error_exit ("Error: invalid topology: bead #,type,mass")
        if (i1 /= i)   call error_exit ("Error: invalid topology: bead #,type,mass")
        a(ibea) = 0  !default: non-reactive bead
        s(ibea) = 0  !no reactiveness sort
        if (iand(breac,t(ibea))/=0) then
          a(ibea) = 1                          !initially active reactive bead
          s(ibea) = ishft(iand(brsrt,t(ibea)),-13)  !activeness sort
        endif
        top(imol) = itop
        if (i==1)    nfir(imol)=ibea
        if (i==kbea) nlst(imol)=ibea
        mol(ibea) = imol
      enddo
      !read-in bonded interactions
      ofs = nfir(imol)-1
      if (kbon==0) then
        bfir(imol) = 0
        blst(imol) = 0
      endif
      do i = 1, kbon
        ibon = ibon + 1
        if (i==1)    bfir(imol) = ibon
        if (i==kbon) blst(imol) = ibon
        read (top_u, fmt=*, iostat=ios) i1,i2,b0(ibon),bk(ibon)
        if (ios /= 0) call error_exit ("Error: invalid topology: bond data")
        b1(ibon)=ofs+i1
        b2(ibon)=ofs+i2
      enddo
      if (kang==0) then
        afir(imol) = 0
        alst(imol) = 0
      endif
      do i = 1, kang
        iang = iang + 1
        if (i==1)    afir(imol) = iang
        if (i==kang) alst(imol) = iang
        read (top_u, fmt=*, iostat=ios) i1,i2,i3,a0(iang),ak(iang)
        if (ios /= 0) call error_exit ("Error: invalid topology: angle data")
        a1(iang)=ofs+i1
        a2(iang)=ofs+i2
        a3(iang)=ofs+i3
      enddo
      if (ktor==0) then
        tfir(imol) = 0
        tlst(imol) = 0
      endif
      do i = 1, ktor
        itor = itor + 1
        if (i==1)    tfir(imol) = itor
        if (i==ktor) tlst(imol) = itor
        read (top_u, fmt=*, iostat=ios) i1,i2,i3,i4,tc(:,itor)
        if (ios /= 0) call error_exit ("Error: invalid topology: torsion data")
        t1(itor)=ofs+i1
        t2(itor)=ofs+i2
        t3(itor)=ofs+i3
        t4(itor)=ofs+i4
      enddo
      !fill-in data for the rest mols of this topology
      !(chain-like filling, the data taken from the prev mol)
      do j = 2, kmol
        imol = imol + 1
        do i = 1, kbea
          ibea = ibea + 1
          top(imol) = top(imol-1)
          if (i==1)    nfir(imol) = ibea
          if (i==kbea) nlst(imol) = ibea
          t(ibea) = t(ibea-kbea)
          a(ibea) = a(ibea-kbea)
          s(ibea) = s(ibea-kbea)
          m(ibea) = m(ibea-kbea)
          mol(ibea) = imol
        enddo
        do i = 1, kbon
          ibon = ibon + 1
          if (i==1)    bfir(imol) = ibon
          if (i==kbon) blst(imol) = ibon
          b1(ibon) = b1(ibon-kbon)+kbea
          b2(ibon) = b2(ibon-kbon)+kbea
          b0(ibon) = b0(ibon-kbon)
          bk(ibon) = bk(ibon-kbon)
        enddo
        do i = 1, kang
          iang = iang + 1
          if (i==1)    afir(imol) = iang
          if (i==kang) alst(imol) = iang
          a1(iang) = a1(iang-kang)+kbea
          a2(iang) = a2(iang-kang)+kbea
          a3(iang) = a3(iang-kang)          
          a0(iang) = a0(iang-kang)
          ak(iang) = ak(iang-kang)
        enddo
        do i = 1, ktor
          itor = itor + 1
          if (i==1)    tfir(imol) = itor
          if (i==ktor) tlst(imol) = itor
          t1(itor) = t1(itor-ktor)+kbea
          t2(itor) = t2(itor-ktor)+kbea
          t3(itor) = t3(itor-ktor)          
          t4(itor) = t4(itor-ktor)          
          tc(:,itor) = tc(:,itor-ktor)          
        enddo
      enddo
    enddo
    t(0) = 0  !set type 0 for ibea=0 (surface beads)
    read (top_u, fmt=*, iostat=ios) str

    !check # consistency
    if (imol/=nmol) call error_exit&
       ("Error: invalid topology: total nmol doesn't match")
    if (ibea/=nbea) call error_exit&
       ("Error: invalid topology: total nbea doesn't match")
    if (ibon/=nbon) call error_exit&
       ("Error: invalid topology: total nbon doesn't match")
    if (iang/=nang) call error_exit&
       ("Error: invalid topology: total nang doesn't match")
    if (itor/=ntor) call error_exit&
       ("Error: invalid topology: total ntor doesn't match")

    !read-in extra crosslinking information
    !crosslink bonds, angles, torsions are appended to
    !the bonds, angles, torsions arrays 
    do i = 1, cbon
      ibon = ibon + 1
      read (top_u, fmt=*, iostat=ios) b1(ibon),b2(ibon),b0(ibon),bk(ibon)
      if (ios /= 0) call error_exit ("Error: invalid topology: crosslink bond data")
      if (iand(breac,t(b1(ibon)))==0.or.iand(breac,t(b2(ibon)))==0) &
        call error_exit ("Error: invalid topology: non-reactive atom listed as crosslinked")
      if (a(b1(ibon))/=1.or.a(b2(ibon))/=1) &
        call error_exit ("Error: invalid topology: crosslinked atom reactive but not marked active")
      a(b1(ibon)) = 2  !initially reactive atom but already linked (inactive)
      a(b2(ibon)) = 2  !initially reactive atom but already linked (inactive)
    enddo
    do i = 1, cang
      iang = iang + 1
      read (top_u, fmt=*, iostat=ios) a1(iang),a2(iang),a3(iang),a0(iang),ak(iang)
      if (ios /= 0) call error_exit ("Error: invalid topology: crosslink angle data")
    enddo
    do i = 1, ctor
      itor = itor + 1
      read (top_u, fmt=*, iostat=ios) t1(itor),t2(itor),t3(itor),t4(itor),tc(:,itor)
      if (ios /= 0) call error_exit ("Error: invalid topology: crosslink torsion data")
    enddo

    close(top_u)

    !get total mass
    mass = sum(m(1:nbea))

    !set attraction and reactive flags to false if no appropriate beads or all reacted (for reactive flag)
    if (attraction .and. count(iand(battr,t(:))/=0)==0) attraction = .false.
    if (reactive   .and. (count(iand(breac,t(:))/=0)==0.or.count(a(:)==1)==0)) reactive = .false.

  end subroutine read_topol


  !----------------------------------------------
  subroutine read_coords
  !----------------------------------------------
  ! format of coord file:
  !  st_n, nbea      //DPD step when this file was flushed
  !  box(:)          //box dimensions
  !{
  !  ibea            //bead index
  !  r(:,ibea)       //bead coords
  !  v(:,ibea)       //bead veloc
  !}
  ! boxpv(:)         //box dim on prev DPD step
  ! xippv(:)         //xip 2 DPD steps before
  ! xipv(:)          //xip on prev DPD steps
  !----------------------------------------------

    implicit none

    integer(si), parameter :: coord_u = 13
    integer(si) :: i, ibea, idx, nbea2, dummy, ios

    !open coord file and read-in header
    open (coord_u, file=fn_cread, status="old", iostat=ios)
    if (ios /= 0) call error_exit ("Error: No input coord file")
    read (coord_u, fmt=*, iostat=ios) st_n0, nbea2, dummy
    if (nbea2 /= nbea) call error_exit &
      ("Error: # of beads in topol and coord files does not match")
    st_n0 = st_n0 + 1
    if (ropt=="NEW") st_n0 = 1
    if (st_nmx < st_n0) call error_exit &
      ("Error: Nothing to do, job until st_nmx completed before")
    read (coord_u, fmt=*, iostat=ios) box(:)
    vol = box(1)*box(2)*box(3)
    dens = nbea/vol
    !read-in idx/coord/veloc
    !here we ignore bead type as it is read-in before from the topol file
    !this opens up the possibility to use old coords with redefined types
    !loaded from some new topology file
    do ibea = 1, nbea
      read (coord_u, fmt=*, iostat=ios) idx, dummy
      if (idx /= ibea) call error_exit &
        ("Error: coord file: wrong bead index")
      read (coord_u, fmt=*, iostat=ios) r(:,ibea)
      if (any(r(:,ibea)<0).or.any(r(:,ibea)>box(:))) call error_exit &
        ("Error: coord: bead out of box")
      read (coord_u, fmt=*, iostat=ios) v(:,ibea)
    enddo

    !read barostat parameters
    if (barostat) then
      read (coord_u, fmt=*, iostat=ios) boxpv(:)
      if (ios /= 0) call error_exit ("Error: coord: no required barostat record")
      read (coord_u, fmt=*, iostat=ios) xippv(:)
      if (ios /= 0) call error_exit ("Error: coord: no required barostat record")
      read (coord_u, fmt=*, iostat=ios) xipv(:)
      if (ios /= 0) call error_exit ("Error: coord: no required barostat record")
    endif

    close (coord_u) 

  end subroutine read_coords


  !----------------------------------------------
  subroutine remove_com
  !----------------------------------------------
  !remove centre of mass motion
  !----------------------------------------------

    implicit none
    integer(si) :: i

    !total momenta
    do i = 1, 3
      tmom(i) = sum(m(1:nbea)*v(i,1:nbea))
    enddo
    !remove com motion
    tmom(:) = tmom(:)/real(nbea,dp)
    do i = 1, 3
      v(i,1:nbea) = v(i,1:nbea)-tmom(i)/m(1:nbea)
    enddo

  end subroutine remove_com
  

  !----------------------------------------------
  subroutine read_react
  !----------------------------------------------
  !format of react file
  !  rtype                       //reactivity option, "SQW", ...
  !  for "SQW" option:
  !    rfreq, rprob, rb0, rbk, cc //frequency, probability, lbond, kbond, cutoff
  !----------------------------------------------

    implicit none
    integer(si) :: i, ios, react_u = 12
    character(256) :: str

    !open file
    open (react_u, file=fn_rea, status="old", iostat=ios)
    if (ios /= 0) call error_exit ("Error: no react file")
    !read rtype
    read (react_u, fmt=*, iostat=ios) str
    if (ios/=0) call error_exit ("Error: Invalid react file, reactivity option line")
    rtype   = str(1:3)
    !check for valid syntax
    if (rtype/="SQW") then
      call error_exit ("Error: Invalid react file, unknown reactivity option")
    endif
    if (rtype=="SQW") then
      read (react_u, fmt=*, iostat=ios) rfreq, rprob, rb0, rbk, cc
      if (rfreq<1 .or. rprob<0.0 .or. rb0<0.0 .or. cc<=0.0) &
        call error_exit ("Error: Invalid react file, rprob, rb0, rbk, cc")
    !reserved for other future options
    !elseif ()
    !elseif ()
    endif
    cc_2 = cc**2
    close (react_u)
 
  end subroutine read_react


  !----------------------------------------------
  subroutine read_grf
  !----------------------------------------------
  !format of grf file
  !  ngrf                     //total # of grafts
  !  [grft(i)={"fix","flt"} grfi(i) grfx(i) grfy(i) grfz(i) grfb0(i) grfbk(i)] 
  !  [grft(i)={"fix","flt"} grfi(i) grfx(i) grfy(i) grfz(i) grfb0(i) grfbk(i)] 
  !  [grft(i)={"fix","flt"} grfi(i) grfx(i) grfy(i) grfz(i) grfb0(i) grfbk(i)] 
  !  \____ graft type ____/  bead   \___ graft coords ____/ \_ bond params _/
  !----------------------------------------------

    implicit none
    integer(si) :: i, ios, grf_u = 12
    character(256) :: str

    !open file
    open (grf_u, file=fn_grf, status="old", iostat=ios)
    if (ios /= 0) call error_exit ("Error: no grf file")
    !read # of grafts
    read (grf_u, fmt=*, iostat=ios) ngrf
    if (ios /= 0) call error_exit&
      ("Error: invalid grf: ngrf")
    allocate (grft(ngrf),grfi(ngrf),grfr(2,ngrf),grfw(ngrf),grfb0(ngrf),grfbk(ngrf),stat=ios)
    if (ios /= 0) call error_exit("Error: cannot allocate grf arrays")
    do i = 1, ngrf
      read (grf_u, fmt=*, iostat=ios) str,grfi(i),grfr(1:2,i),grfw(i),grfb0(i),grfbk(i)
      if (ios /= 0) call error_exit&
        ("Error: invalid grf: graft params line")
      if (str(1:3)=="fix") then
        grft(i) = 1
      elseif (str(1:3)=="flt") then
        grft(i) = 2
      else
        call error_exit ("Error: invalid grf: valid graft types are 'fix' and 'flt'")
      endif
      if (.not.(grfw(i)==0) .and. .not.(grfw(i)==1))&
        call error_exit ("Error: invalid grfw: valid options are 0 or 1")
    enddo
    close(grf_u)
 
  end subroutine read_grf

  !----------------------------------------------
  subroutine read_cyl
  !----------------------------------------------
  !format of cyl file
  !  ncyl                              //total # of cylindrical pores
  !  rcylc(1,1)     rcylc(2,1)     cylr(1) 
  !  .......................... ............
  !  rcylc(1,ncyl)  rcylc(2,ncyl)  cylr(ncyl) 
  !  \___xpos___/   \___ypos___/   \_radius_/
  !----------------------------------------------

    implicit none
    integer(si) :: i, ios, cyl_u = 12
    character(256) :: str

    !open file
    open (cyl_u, file=fn_cyl, status="old", iostat=ios)
    if (ios /= 0) call error_exit ("Error: no cyl file")
    !read # of cyls
    read (cyl_u, fmt=*, iostat=ios) ncyl
    if (ios /= 0) call error_exit&
      ("Error: invalid cyl: ncyl")
    allocate (rcylc(3,ncyl),cylr(ncyl),stat=ios)
    if (ios /= 0) call error_exit("Error: cannot allocate cyl arrays")
    do i = 1, ncyl
      read (cyl_u, fmt=*, iostat=ios) rcylc(1,i), rcylc(2,i), cylr(i)
      if (ios /= 0) call error_exit&
        ("Error: invalid cyl: cyl params line")
      if (rcylc(1,i)<=0.0_dp .or. rcylc(1,i)>=box(1) .or.&
          rcylc(2,i)<=0.0_dp .or. rcylc(2,i)>=box(2) .or.&
          cylr(i)<0.0_dp)&
        call error_exit ("Error: invalid cyl: cylinder center outside the box or invalid radius")
    enddo
    close(cyl_u)
 
  end subroutine read_cyl


  !==============================================
  !================ WRITE FILES  ================
  !==============================================

  !----------------------------------------------
  subroutine add_to_lst
  !----------------------------------------------
  ! save instant data into list arrays
  !----------------------------------------------

    implicit none

    lst_ptr = lst_ptr + 1

    lst_dens(lst_ptr) = dens
    lst_kt(lst_ptr)   = kt
    lst_pres(lst_ptr) = pres

    lst_boxx(lst_ptr) = box(1)
    lst_boxy(lst_ptr) = box(2)
    lst_boxz(lst_ptr) = box(3)

    lst_pxx(lst_ptr)  = presv(1)
    lst_pyy(lst_ptr)  = presv(2)
    lst_pzz(lst_ptr)  = presv(3)

    lst_momx(lst_ptr) = tmom(1)
    lst_momy(lst_ptr) = tmom(2)
    lst_momz(lst_ptr) = tmom(3)


  end subroutine add_to_lst

  
  !----------------------------------------------
  subroutine flush_lst
  !----------------------------------------------
  ! flush list arrays
  ! format of list file:
  !  "header string with properties names listed"
  !{
  !  step dens kT pres Lx Ly Lz Px Py Pz
  !} 
  !----------------------------------------------
  
  implicit none
  integer(si) :: lst_u=40, i

  !open list file
  if (ropt=="NEW".and.st_n==st_cwrit) then
    open (lst_u,file=fn_list,position="rewind")
    write(lst_u,fmt='("   step     dens     kT     pres         ",&
                  "Lx      Ly      Lz          Px       Py       Pz         px        py        pz")')
  else
    open (lst_u, file=fn_list, position="append")
  endif
  
  !flush data buffer
  do i = 1, lst_dim
    !         step   dens  kt  pres     Lxyz       Pxyz         px        py        pz
     write (lst_u,fmt='(i9,f8.4,f8.4,f9.4,3X,3(f8.3),3X,3(f9.4),3X,3(g10.2))') &
      st_n-st_cwrit+i*st_lst,lst_dens(i),lst_kt(i),lst_pres(i),&
      lst_boxx(i),lst_boxy(i),lst_boxz(i),&
      lst_pxx(i),lst_pyy(i),lst_pzz(i),&
      lst_momx(i),lst_momy(i),lst_momz(i)
  enddo
  close(lst_u)

  !init list pointer
  lst_ptr = 0

  end subroutine flush_lst


  !----------------------------------------------
  subroutine write_topol
  !----------------------------------------------
  !format of topol file
  !  ntop nmol nbea                     //total # in the melt
  !  nbon nang ntor                     //total regular bonding
  !  cbon cang ctor                     //total crosslinked bonding
  !{
  !  ----------------------------
  !  itop kmol                          //topol index
  !  kbea kbon kang ktor                //# of mols and one mol #
  !  [i t(i) m(i)]                      //beads types and masses
  !  [b1(i) b2(i) b0(i) bk(i)]          //regular bonds
  !  [a1(i) a2(i) a3(i) a0(i) ak(i)     //regular angles
  !  [t1(i) t2(i) t3(i) t4(i) tc(:,i)]  //regular torsions
  !}
  !  ----------------------------
  !  [b1(i) b2(i) b0(i) bk(i)]          //crosslinked bonds
  !  [a1(i) a2(i) a3(i) a0(i) ak(i)     //crosslinked angles
  !  [t1(i) t2(i) t3(i) t4(i) tc(:,i)]  //crosslinked torsions
  !----------------------------------------------

    implicit none
    character(fnlen) :: filename
    integer(si) :: top_u = 11
    integer(si) :: i, itop, kmol, ofs, ios
    integer(si) :: imol,ibea,ibon,iang,itor
    integer(si) ::      kbea,kbon,kang,ktor

    !form filename and open file
    write (filename, fmt='(i8.8,"_",a)') st_n, trim(fn_twrit)
    open (top_u, file=filename, status="unknown")

    !write header
    write (top_u, fmt='(3i8)') ntop,nmol,nbea
    write (top_u, fmt='(3i8)') nbon,nang,ntor
    write (top_u, fmt='(3i8)') cbon,cang,ctor

    !loop over topologies
    do itop = 1, ntop
      !find kmol and first molecule imol for topology itop
      kmol = 0
      do i = 1, nmol
        if (top(i)==itop) then
          kmol = kmol + 1
          if (kmol==1) imol=i
        endif
      enddo
      if (kmol>0) then
        !save topology information using molecule imol
        write (top_u,fmt='(a)') "-------------------------------------------"
        !write itop, kmol
        write (top_u,fmt='(2i8)') itop, kmol
        !get and write kbea,kbon,kang,ktor
        kbea = nlst(imol)-nfir(imol)+1  !there is always at least one bead
        kbon = 0
        kang = 0
        ktor = 0
        if (bfir(imol)>0) kbon=blst(imol)-bfir(imol)+1
        if (afir(imol)>0) kang=alst(imol)-afir(imol)+1
        if (tfir(imol)>0) ktor=tlst(imol)-tfir(imol)+1
        write (top_u,fmt='(4i6)') kbea, kbon, kang, ktor
        !write beads data 
        ofs = nfir(imol)-1
        do ibea = nfir(imol), nlst(imol)
          write(top_u,fmt='(i6,i8,f12.3)') ibea-ofs,t(ibea),m(ibea)
        enddo
        !write bonds
        do i = 1, kbon
          ibon = bfir(imol)+i-1
          write(top_u,fmt='(2i6,2f8.3)') b1(ibon)-ofs,b2(ibon)-ofs,b0(ibon),bk(ibon)
        enddo
        !write angles
        do i = 1, kang
          iang = afir(imol)+i-1
          write(top_u,fmt='(3i6,2f8.3)') a1(iang)-ofs,a2(iang)-ofs,a3(iang)-ofs,&
                                         a0(iang),ak(iang)
        enddo
        !write torsions
        do i = 1, ktor
          itor = tfir(imol)+i-1
          write(top_u,fmt='(4i6,3f8.3)') t1(itor)-ofs,t2(itor)-ofs,t3(itor)-ofs,&
                                         t4(itor)-ofs,tc(:,itor)
        enddo
      else
        !save empty topology
        write (top_u,fmt='(a)') "-------------------------------------------"
        write (top_u,fmt='(2i8)') itop, 0
        write (top_u,fmt='(4i6)') 1, 0, 0, 0
        write (top_u,fmt='(i6,i8,f12.3)') 1,t(1),m(1)
      endif
    enddo
    write (top_u,fmt='(a)') "-------------------------------------------"

    !write crosslink data
    !write bonds
    do ibon = nbon+1, nbon+cbon
      write(top_u,fmt='(2i8,2f8.3)') b1(ibon),b2(ibon),b0(ibon),bk(ibon)
    enddo
    !write angles
    do iang = nang+1, nang+cang
      write(top_u,fmt='(3i8,2f8.3)') a1(iang),a2(iang),a3(iang),a0(iang),ak(iang)
    enddo
    !write torsionss
    do itor = ntor+1, ntor+ctor
      write(top_u,fmt='(4i8,3f8.3)') t1(itor),t2(itor),t3(itor),t4(itor),tc(:,itor)
    enddo
    
    close(top_u)

  end subroutine write_topol


  !----------------------------------------------
  subroutine write_coords
  !----------------------------------------------
  ! format of coord file:
  !  st_n, nbea      //DPD step when this file was flushed
  !  box(:)          //box dimensions
  !{
  !  ibea            //bead index
  !  r(:,ibea)       //bead coords
  !  v(:,ibea)       //bead veloc
  !}
  ! boxpv(:)         //box dim on prev DPD step
  ! xippv(:)         //xip 2 DPD steps before
  ! xipv(:)          //xip on prev DPD steps
  !----------------------------------------------

    implicit none
    character(fnlen) :: filename
    integer(si) :: coord_u=20, bar_u=21
    integer(si) :: ibea
    character(32) :: str1,str2,str3

    !form filename and open file
    write (filename, fmt='(i8.8,"_",a)') st_n, trim(fn_cwrit)
    open (coord_u, file=filename, status="unknown")

    !write header
    !write (coord_u, fmt='(3i8)') st_n, nbea, 0
    write (str1,*) st_n
    write (str2,*) nbea
    write (str3,*) 0
    write (coord_u,fmt='(a,X,a,X,a)') &
           trim(adjustl(str1)), trim(adjustl(str2)), trim(adjustl(str3))
    !write (coord_u, fmt='(3e16.8)') box(:)
    write (str1,fmt='(f16.6)') box(1)
    write (str2,fmt='(f16.6)') box(2)
    write (str3,fmt='(f16.6)') box(3)
    write (coord_u,fmt='(a,X,a,X,a)') &
           trim(adjustl(str1)), trim(adjustl(str2)), trim(adjustl(str3))

    !write beads information
    do ibea = 1, nbea
      !write (coord_u, fmt='(2i6)')    ibea, t(ibea)
      write (str1,*) ibea
      write (str2,*) t(ibea)
      write (coord_u,fmt='(a,X,a)') trim(adjustl(str1)), trim(adjustl(str2))
      !write (coord_u, fmt='(3e16.8)') r(:,ibea)
      write (str1,fmt='(f16.6)') r(1,ibea)
      write (str2,fmt='(f16.6)') r(2,ibea)
      write (str3,fmt='(f16.6)') r(3,ibea)
      write (coord_u,fmt='(a,X,a,X,a)') &
             trim(adjustl(str1)), trim(adjustl(str2)), trim(adjustl(str3))
      !write (coord_u, fmt='(3e16.8)') v(:,ibea)
      write (str1,fmt='(f16.6)') v(1,ibea)
      write (str2,fmt='(f16.6)') v(2,ibea)
      write (str3,fmt='(f16.6)') v(3,ibea)
      write (coord_u,fmt='(a,X,a,X,a)') &
             trim(adjustl(str1)), trim(adjustl(str2)), trim(adjustl(str3))
    enddo

    !append barostat data
    if (barostat) then
      !write (coord_u,fmt='(3(g14.6))') boxpv(:)
      write (str1,fmt='(f16.6)') boxpv(1)
      write (str2,fmt='(f16.6)') boxpv(2)
      write (str3,fmt='(f16.6)') boxpv(3)
      write (coord_u,fmt='(a,X,a,X,a)') &
             trim(adjustl(str1)), trim(adjustl(str2)), trim(adjustl(str3))
      !write (coord_u,fmt='(3(g14.6))') xippv(:)
      write (str1,fmt='(g16.6)') xippv(1)
      write (str2,fmt='(g16.6)') xippv(2)
      write (str3,fmt='(g16.6)') xippv(3)
      write (coord_u,fmt='(a,X,a,X,a)') &
             trim(adjustl(str1)), trim(adjustl(str2)), trim(adjustl(str3))
      !write (coord_u,fmt='(3(g14.6))') xipv(:)
      write (str1,fmt='(g16.6)') xipv(1)
      write (str2,fmt='(g16.6)') xipv(2)
      write (str3,fmt='(g16.6)') xipv(3)
      write (coord_u,fmt='(a,X,a,X,a)') &
             trim(adjustl(str1)), trim(adjustl(str2)), trim(adjustl(str3))
    endif

    close(coord_u)

  end subroutine write_coords


  !===============================================
  !============= FORCES AND VIRIAL ===============
  !===============================================


  !----------------------------------------------
  subroutine eval_for_vir
  !----------------------------------------------

    implicit none
    integer(si) :: ibea, jbea
    integer(si) :: ncell, icell, jcell
    integer(si), allocatable, dimension(:) :: head, list
    integer(si), dimension(3) :: ncxyz, ixyz, jxyz
    real(dp), dimension(3) :: cdim, rofs, rijv, vijv
    integer(si), dimension(3,14) :: jofs = reshape ( &
       (/ 0,0,0,  1, 0,0, 1, 1,0, -1, 1,0,  0,1,0, 0,0,1, -1,0,1,    &
          1,0,1, -1,-1,1, 0,-1,1,  1,-1,1, -1,1,1, 0,1,1,  1,1,1 /), &
       (/3,14/))
    integer(si) :: i1, i2, i3, nn
    integer(si) :: ibon, iang, itor, igrf, icyl
    real(dp) :: rij, rij_2, dth, cut, z, z_2, bz
    real(dp), dimension(2) :: dxy
    integer(si) :: ios

    !clear arrays
    for(:,:) = 0.0_dp
    virv(:)  = 0.0_dp

    !pairwise repulsive forces (all pairs within rc distance are evaluated)
    !linked cells params
    ncxyz(:) = int (box(:)/rc+0.0001,si)
    if (any(ncxyz(:)<3)) call error_exit &
      ("Error: repulsive forces, box too small, less than 3 cells in one of dimensions")
    cdim(:) = box(:)/ncxyz(:)
    ncell = ncxyz(1)*ncxyz(2)*ncxyz(3)
    allocate (head(ncell), list(nbea), stat=ios)
    if (ios /= 0) call error_exit("Error: cannot allocate head(), list() for repulsive")

    !build linked cells lists
    head(:) = 0
    do ibea = 1, nbea
      ixyz(:) = int(r(:,ibea)/cdim(:),si)
      icell = 1 + ixyz(1) + ncxyz(1)*(ixyz(2) + ncxyz(2)*ixyz(3))
      if (icell > ncell) then
        write (str,*) "Error: repulsive forces, cell idx",icell,", max =",ncell,", bead",ibea
        call error_exit (str)
      endif
      jbea = head(icell)
      head(icell) = ibea
      list(ibea) = jbea
    enddo

    !use linked cells to evaluate pairwise repulsive forces
    !loop over i-th cells using lattice sites
    do i1 = 0, ncxyz(1)-1
    do i2 = 0, ncxyz(2)-1
    do i3 = 0, ncxyz(3)-1
      ixyz(1) = i1
      ixyz(2) = i2
      ixyz(3) = i3
      icell = 1 + ixyz(1) + ncxyz(1)*(ixyz(2) + ncxyz(2)*ixyz(3))
      if (head(icell) == 0) cycle
      !loop over j-th cells as i-th cell neighbours
      do nn = 1, 14
        jxyz(:) = ixyz(:) + jofs(:,nn)
        if (surface) then
          if (jxyz(3)==-1 .or. jxyz(3)==ncxyz(3)) cycle  !skip PBC in Z for jcell
        endif
        !coords offset for j-th cell beads
        rofs(:) = merge(-box(:), 0.0_dp,  jxyz(:)==-1)
        rofs(:) = merge(+box(:), rofs(:), jxyz(:)==ncxyz(:))
        !cells PBC
        jxyz(:) = merge(ncxyz(:)-1, jxyz(:), jxyz(:)==-1)
        jxyz(:) = merge(0         , jxyz(:), jxyz(:)==ncxyz(:))
        jcell = 1 + jxyz(1) + ncxyz(1)*(jxyz(2) + ncxyz(2)*jxyz(3))
        if (head(jcell) == 0) cycle
        !primary cell head
        ibea = head (icell)
        !loop over beads in a primary cell
        do
          if (ibea == 0) exit
          !secondary cell head
          jbea = head (jcell)
          if (jcell == icell) jbea = list(ibea)
          !loop over beads in a secondary cell
          do
            if (jbea == 0) exit
            if (jbea == ibea) call error_exit ("Error: repulsive forces, i=j when processing linked lists")
            !check beads separation
            rijv(:) = r(:,ibea) - (r(:,jbea)+rofs(:))
            rij_2 = rijv(1)*rijv(1)+rijv(2)*rijv(2)+rijv(3)*rijv(3)
            if (rij_2 <= rc_2) then
              vijv(:) =  v(:,ibea) - v(:,jbea)
              call eval_cdr_uni (ibea, jbea, rijv, vijv, rij_2)
            endif
            jbea = list(jbea)
          enddo
          ibea = list(ibea)
        enddo
      enddo
    enddo
    enddo
    enddo

    deallocate (head, list)

    !pairwise attractive forces (interacting beads are marked by 'battr' bit in their type)
    if (attraction) then

      !linked cells params
      ncxyz(:) = int (box(:)/qc+0.0001,si)
      if (any(ncxyz(:)<3)) call error_exit &
        ("Error: attractive forces, box too small, less than 3 cells detected in one of dimensions")
      cdim(:) = box(:)/ncxyz(:)
      ncell = ncxyz(1)*ncxyz(2)*ncxyz(3)
      allocate (head(ncell), list(nbea), stat=ios)
        if (ios /= 0) call error_exit ("Error: cannot allocate head(), list() for attractive")

      !build linked cells lists
      head(:) = 0
      do ibea = 1, nbea
        if (iand(battr,t(ibea))==0) cycle
        ixyz(:) = int(r(:,ibea)/cdim(:),si)
        icell = 1 + ixyz(1) + ncxyz(1)*(ixyz(2) + ncxyz(2)*ixyz(3))
        if (icell > ncell) then
          write (str,*) "Error: attractive forces, cell idx",icell,", max =",ncell,", bead",ibea
          call error_exit (str)
        endif
        jbea = head(icell)
        head(icell) = ibea
        list(ibea) = jbea
      enddo

      !use linked cells to evaluate pairwise attractive forces
      !loop over i-th cells using lattice sites
      do i1 = 0, ncxyz(1)-1
      do i2 = 0, ncxyz(2)-1
      do i3 = 0, ncxyz(3)-1
        ixyz(1) = i1
        ixyz(2) = i2
        ixyz(3) = i3
        icell = 1 + ixyz(1) + ncxyz(1)*(ixyz(2) + ncxyz(2)*ixyz(3))
        if (head(icell) == 0) cycle
        !loop over j-th cells as i-th cell neighbours
        do nn = 1, 14
          jxyz(:) = ixyz(:) + jofs(:,nn)
          if (surface) then
            if (jxyz(3)==-1 .or. jxyz(3)==ncxyz(3)) cycle  !skip PBC in Z for jcell
          endif
          !coords offset for j-th cell beads
          rofs(:) = merge(-box(:), 0.0_dp,  jxyz(:)==-1)
          rofs(:) = merge(+box(:), rofs(:), jxyz(:)==ncxyz(:))
          !cells PBC
          jxyz(:) = merge(ncxyz(:)-1, jxyz(:), jxyz(:)==-1)
          jxyz(:) = merge(0         , jxyz(:), jxyz(:)==ncxyz(:))
          jcell = 1 + jxyz(1) + ncxyz(1)*(jxyz(2) + ncxyz(2)*jxyz(3))
          if (head(jcell) == 0) cycle
          !primary cell head
          ibea = head (icell)
          !loop over beads in a primary cell
          do
            if (ibea == 0) exit
            !secondary cell head
            jbea = head (jcell)
            if (jcell == icell) jbea = list(ibea)
            !loop over beads in a secondary cell
            do
              if (jbea == 0) exit
              if (jbea == ibea) call error_exit ("Error: attractive forces, i=j when processing linked lists")
              !check beads separation
              rijv(:) = r(:,ibea) - (r(:,jbea)+rofs(:))
              rij_2 = rijv(1)*rijv(1)+rijv(2)*rijv(2)+rijv(3)*rijv(3)
              if (rij_2 >= rc_2 .and. rij_2 <= qc_2) call eval_attr_uni (ibea, jbea, rijv, rij_2)
              jbea = list(jbea)
            enddo
            ibea = list(ibea)
          enddo
        enddo
      enddo
      enddo
      enddo

      deallocate (head, list)
    endif

    !reactive forces (only active beads with a(ibea)=1 within cc distance are checked)
    if (reactive .and. mod(st_n,rfreq)==0) then

      !linked cells params for reactivity
      ncxyz(:) = int (box(:)/cc+0.0001,si)
      if (any(ncxyz(:)<3)) call error_exit &
        ("Error: reactivity, box too small, less than 3 cells detected in one of dimensions")
      cdim(:) = box(:)/ncxyz(:)
      ncell = ncxyz(1)*ncxyz(2)*ncxyz(3)
      allocate (head(ncell), list(nbea))

      !build linked cells lists for currently available reactive sites
      head(:) = 0
      do ibea = 1, nbea
        if (a(ibea)==1) then
          ixyz(:) = int(r(:,ibea)/cdim(:),si)
          icell = 1 + ixyz(1) + ncxyz(1)*(ixyz(2) + ncxyz(2)*ixyz(3))
          if (icell > ncell) then
            write (str,*) "Error: reactivity, cell idx",icell,", max =",ncell,", bead",ibea
            call error_exit (str)
          endif
          jbea = head(icell)
          head(icell) = ibea
          list(ibea) = jbea
        endif
      enddo

      !use linked cells to evaluate reactivity
      !(some reactive sites will loose reactivity within the loop)
      !loop over i-th cells using lattice sites
      do i1 = 0, ncxyz(1)-1
      do i2 = 0, ncxyz(2)-1
      do i3 = 0, ncxyz(3)-1
        ixyz(1) = i1
        ixyz(2) = i2
        ixyz(3) = i3
        icell = 1 + ixyz(1) + ncxyz(1)*(ixyz(2) + ncxyz(2)*ixyz(3))
        if (head(icell) == 0) cycle
        !loop over j-th cells as i-th cell neighbours
        do nn = 1, 14
          jxyz(:) = ixyz(:) + jofs(:,nn)
          if (surface) then
            if (jxyz(3)==-1 .or. jxyz(3)==ncxyz(3)) cycle  !skip PBC in Z for jcell
          endif
          !coords offset for j-th cell beads
          rofs(:) = merge(-box(:), 0.0_dp,  jxyz(:)==-1)
          rofs(:) = merge(+box(:), rofs(:), jxyz(:)==ncxyz(:))
          !cells PBC
          jxyz(:) = merge(ncxyz(:)-1, jxyz(:), jxyz(:)==-1)
          jxyz(:) = merge(0         , jxyz(:), jxyz(:)==ncxyz(:))
          jcell = 1 + jxyz(1) + ncxyz(1)*(jxyz(2) + ncxyz(2)*jxyz(3))
          if (head(jcell) == 0) cycle
          !primary cell head
          ibea = head (icell)
          !loop over beads in a primary cell
          do
            if (ibea == 0) exit
            if (a(ibea)==1) then
              !secondary cell head
              jbea = head (jcell)
              if (jcell == icell) jbea = list(ibea)
              !loop over beads in a secondary cell
              do
                if (jbea == 0) exit
                if (jbea == ibea) call error_exit ("Error: reactivity, i=j when processing linked lists")
                if (a(jbea)==1) then
                  !avoid crosslinking of beads from the same molecule
                  !avoid crosslinking beads of the same reactivity sort
                  !if (mol(jbea)/=mol(ibea).and.s(jbea)/=s(ibea)) then
                  if (mol(jbea)/=mol(ibea)) then
                    !check beads separation
                    rijv(:) = r(:,ibea) - (r(:,jbea)+rofs(:))
                    rij_2 = rijv(1)*rijv(1)+rijv(2)*rijv(2)+rijv(3)*rijv(3)
                    if (rij_2 <= cc_2) then
                      !if (a(ibea)/=1.or.a(jbea)/=1) write (*,*) "   error:",ibea,jbea,t(ibea),t(jbea)
                      if (nbon+cbon==nbontop) call error_exit("No room to store the next crosslink bond")
                      call attempt_crosslink (ibea, jbea)
                      !write(51,*) ibea, jbea, mol(ibea), mol(jbea), s(ibea), s(jbea)
                    endif
                  endif
                endif
                if (a(ibea)/=1) exit  !ibea was crosslinked inside the jbea loop
                jbea = list(jbea)
              enddo
            endif
            ibea = list(ibea)
          enddo
        enddo
      enddo
      enddo
      enddo
      deallocate (head, list)
    endif

    !bonded forces (intramolecular + crosslinked)
    do ibon = 1, nbon+cbon
      call eval_bond_uni (b1(ibon), b2(ibon), b0(ibon), bk(ibon), 1)
    enddo
    do iang = 1, nang+cang
      call eval_ang (a1(iang), a2(iang), a3(iang), a0(iang), ak(iang))
    enddo
    do itor = 1, ntor+ctor
      call eval_tors (t1(itor), t2(itor), t3(itor), t4(itor), tc(:,itor))
    enddo

    !surface forces: each bead (of repulsive type 1,2,3) is repulsed from the surface with given prep(1:3,0) if z<rc
    !and if attraction then each bead (of repulsive type 1,2,3) is attracted to it with given peps(1:3,0) if rc<z<qc
    if (surface) then
      bz = box(3)
      do ibea = 1, nbea
        z = r(3,ibea)
        !check for repulsive interaction
        !imaginary bead has the same (x,y) and z=0 or z=box(3)
        if (z<rc) then
          rijv(1:2) = 0.0_dp  !bottom surface
          rijv(3) = z
          vijv(:) = v(:,ibea)
          rij_2 = z**2
          call eval_cdr_uni (ibea, 0, rijv, vijv, rij_2)
        elseif ((bz-z)<rc) then
          rijv(1:2) = 0.0_dp  !upper surface
          rijv(3) = z-bz
          vijv(:) = v(:,ibea)
          rij_2 = (bz-z)**2
          call eval_cdr_uni (ibea, 0, rijv, vijv, rij_2)
        endif
        !check for attractive interaction
        if (attraction .and. peps(iand(15,t(ibea)),0)>0.0_dp) then
          if (z>rc.and.z<qc) then
            rijv(1:2) = 0.0_dp  !bottom surface
            rijv(3) = z
            rij_2 = z**2
            call eval_attr_uni (ibea, 0, rijv, rij_2)
          elseif ((bz-z)>rc.and.(bz-z)<qc) then
            rijv(1:2) = 0.0_dp  !upper surface
            rijv(3) = bz-z
            rij_2 = (bz-z)**2
          call eval_attr_uni (ibea, 0, rijv, rij_2)
          endif
        endif
      enddo
    endif

    !graft forces
    if (grafts) then
      do igrf = 1, ngrf
        call eval_bond_uni (grfi(igrf), igrf, grfb0(igrf), grfbk(igrf), 2)
      enddo
    endif

    !cylinders: bead of type 1,2,3 is repulsed from the cylinder surface with given prep(1:3,0) if z<rc
    !and attracted to it with given peps(1:3,0) if rc<z<qc
    if (cylinders) then
      do icyl = 1, ncyl
        do ibea = 1, nbea
          !imaginary bead is on the surface of the cylinder
          rijv(1:2) = r(1:2,ibea)-rcylc(1:2,icyl)  !vector from cylinder's center
          rijv(1:2) = rijv(1:2) - box(1:2)*nint(rijv(1:2)/box(1:2))  !PBC
          rij = sqrt(sum(rijv(1:2)**2))  !distance to cyl center
          !assign rijv and rij_2
          rijv(1:2) = rijv(1:2)*(1.0_dp - cylr(icyl)/rij) !vector from imaginary bead to ibea
          rijv(3) = 0.0_dp  !we work in XY plane
          rij_2 = sum(rijv(1:2)**2)
          if (rij_2<rc_2) then
            vijv(:) = v(:,ibea)
            call eval_cdr_uni (ibea, 0, rijv, vijv, rij_2)
          endif
          !check for attractive interaction
          if (attraction .and. peps(iand(15,t(ibea)),0)>0.0_dp) then
            if (rij_2>rc_2.and.rij_2<qc_2) then
              call eval_attr_uni (ibea, 0, rijv, rij_2)
            endif
          endif
        enddo
      enddo
    endif

    !account for the cumulative force imposed on walls and cylinders' surfaces
    for(:,0) = for(:,0)/real(nbea,dp)  !contribution per each bead
    do ibea = 1, nbea
      for(:,ibea) = for(:,ibea) + for(:,0)
    enddo

    !replace forces  by f/m*dt/2,
    dth = dt/2.0_dp
    do ibea = 1, nbea
      for(:,ibea) = for(:,ibea)/m(ibea)*dth
    enddo

    !we use Allen-Tildesley definition for the virial W=1/3*sum(r*f)
    virv(:) = virv(:)/3.0_dp

  end subroutine eval_for_vir


  !----------------------------------------------
  subroutine eval_kin_kt_pres
  !----------------------------------------------
  ! evaluate components of kinetic energy, kinetic energy, kT
  !------------------------------------------------- 

    implicit none
    integer(si) :: ibea

    !kinetic energy
    kinv(:) = 0.0_dp
    do ibea = 1, nbea
      kinv(:) = kinv(:) + m(ibea)*v(:,ibea)**2
    enddo
    kinv(:) = kinv(:)/2.0_dp
    kin = kinv(1)+kinv(2)+kinv(3)

    !kT
    kt = 2.0_dp*kin/(3.0_dp*(nbea-1))

    !pressure
    presv(:) = (2.0_dp*kinv(:)+3.0_dp*virv(:))/(3.0_dp*vol)
    pres = presv(1)+presv(2)+presv(3)

  end subroutine eval_kin_kt_pres


  !===============================================
  !================= INTEGRATOR ==================
  !===============================================

  !----------------------------------------------
  subroutine set_ktfix_pfix
  !----------------------------------------------

    implicit none

    if (topt=="TDSE") then
      ktfix=merge(ktfixs+(ktfixe-ktfixs)/(ktfixn-1)*(st_n-1),ktfixe,st_n<=ktfixn)
    endif

    if ((barostat .and.(popt=="BRSE".or.popt=="ABSE")).or.&
        (berendsen.and.(popt=="BESE".or.popt=="AESE"))) then
      pxyz(:)=merge(pxyzs(:)+(pxyze(:)-pxyzs(:))/(pfixn-1)*(st_n-1),pxyze(:),st_n<=pfixn)
    endif

  end subroutine set_ktfix_pfix


  !-------------------------------------------------
  subroutine integr_vv_nvt
  !-------------------------------------------------
  ! velocity-Verlet algorithm, NVT
  ! forces at r(t+dt), v(t+dt/2+vvext*dt/2)
  ! after Groot and Warren, J.Chem.Phys.107, 4423 (1997)
  !-------------------------------------------------

    implicit none
    real(dp), dimension(3) :: rijv, drhv, tailv, cenintv, ceninthv, flipv, dvxy
    real(dp), dimension(2) :: root

    integer(si) :: ibea, icyl
    real(dp) :: dist, rij, dvz, drlen, b, c, discr, tail
    logical :: reflected

    !force array: f(t)/m*dth (calculated on the previous step)
    !preliminary advances (forces at time t)
    smom(:) = 0.0_dp

    do ibea = 1, nbea
      !velocity and coordinates primary advances
      v(:,ibea) = v(:,ibea) + for(:,ibea)  !veloc: v(t+dt/2)
      r(:,ibea) = r(:,ibea) + v(:,ibea)*dt !coord: r(t+dt) based on v(t+dt/2)
      !PBC, dtop=2 for surface, =3 otherwise
      if (any (r(1:dtop,ibea)<0.0_dp .or. r(1:dtop,ibea)>box(1:dtop))) then
        r(1:dtop,ibea) = merge(r(1:dtop,ibea)+box(1:dtop),r(1:dtop,ibea),r(1:dtop,ibea)<0.0_dp)
        r(1:dtop,ibea) = merge(r(1:dtop,ibea)-box(1:dtop),r(1:dtop,ibea),r(1:dtop,ibea)>=box(1:dtop))
      endif

      !check for reflections out of cylinders and walls
      do
        !check for crossing cylinders surfaces
        reflected = .false.
        if (cylinders) then
          do icyl = 1, ncyl
            !check for crossing cylinder's surface
            rijv(1:2) = r(1:2,ibea)-rcylc(1:2,icyl)  !2D vector from the circle center to r(t+dt)
            rijv(1:2) = rijv(1:2) - box(1:2)*nint(rijv(1:2)/box(1:2))  !PBC
            rij = sqrt(sum(rijv(1:2)**2))  !2D distance to the circle center
            if (rij < cylr(icyl)) then
              !apply reflection out of a tangential plane at the intersection point
              drlen = sqrt(sum(v(1:2,ibea)**2))
              drhv(1:2) = v(1:2,ibea)/drlen  !2D unit vector along delta_r=r(t+dt)-r(t) [=v(t+dt/2)*dt]
              b = -2.0_dp*sum(drhv(1:2)*rijv(1:2))   !coefficients of quadratic equation
              c = sum(rijv(1:2)**2) - cylr(icyl)**2  !for delta_r tail inside a cylinder
              discr = b**2 - 4.0_dp*c                !discriminant
              if (discr<0.0_dp) call error_exit("Error in reflecting off cylinder: D<0")
              root(1) = (-b-sqrt(discr))/2.0_dp  !roots
              root(2) = (-b+sqrt(discr))/2.0_dp
              if (all(root(1:2)>0.0_dp)) call error_exit("Error in reflecting off cylinder: both roots>0")
              tail = merge(root(1),root(2),root(1)>0.0_dp)
              tailv(1:2) = tail*drhv(1:2)  !2D vector from 2D intersection point to r(t+dt)
              cenintv(1:2) = rijv(1:2) - tailv(1:2)  !2D vector from the circle center to intersection point
              drlen = sqrt(sum(cenintv(1:2)**2))
              ceninthv(1:2) = cenintv(1:2)/drlen  !2D unit vector along cenintv(1:2)
              flipv(1:2) = -2.0_dp*(sum(ceninthv(1:2)*tailv(1:2)))*ceninthv(1:2)  !2D flipping vector
              !apply reflection off the cylinder
              r(1:2,ibea) = r(1:2,ibea) + flipv(1:2)  !flip tail, XY coord
              dvxy(1:2) = flipv(1:2)/dt  !contribution to the XY velocity
              v(1:2,ibea) = v(1:2,ibea) + dvxy(1:2)  !velocity update
              smom(1:2) = smom(1:2)-m(ibea)*dvxy(1:2)  !momentum asquired by the cylinder
              reflected = .true.
            endif
          enddo
        endif
        !check for the crossing the walls
        if (surface) then
          if (r(3,ibea)<0.0_dp) then
            dist = -r(3,ibea)          !distance from the bottom wall
            dvz = 2.0_dp*dist/dt       !vel shift on z when flipping
            r(3,ibea) = dist           !reflection of z coord
            v(3,ibea) = v(3,ibea)+dvz  !v(t+dt/2) after reflection
            smom(3) = smom(3)-m(ibea)*dvz  !momentum asquired by the wall
            reflected = .true.
          elseif (r(3,ibea)>box(3)) then
            dist = r(3,ibea)-box(3)    !distance from the upper wall
            dvz = -2.0_dp*dist/dt      !vel shift on z when flipping
            r(3,ibea) = box(3) - dist  !reflection of z coord
            v(3,ibea) = v(3,ibea)+dvz  !v(t+dt/2) after reflection
            smom(3) = smom(3)-m(ibea)*dvz  !momentum asquired by the wall
            reflected = .true.
          endif
        endif
        if (.not.reflected) exit  !no reflection occurred from both cylinders and surfaces 
      enddo

      vsave(:,ibea) = v(:,ibea)    !save advanced (and possibly reflected) v(t+dt/2)
      !approx velocity v(t+dt) for forces evaluation
      !(advanced and possibly reflected v(t+dt/2) plus correction assuming old f(t)/m*dth)
      v(:,ibea) = v(:,ibea)+vvext*for(:,ibea)
    enddo

    !calculate forces and gorques at r(t+dt), approx.v(t+dt)
    call eval_for_vir

    !complete velocities advance (now forces at time t+dt)
    smom(:) = smom(:)/real(nbea,dp)  !surface: contribution to each bead from the wall momentum, =0 for no surface
    do ibea = 1, nbea
      v(:,ibea) = vsave(:,ibea) + for(:,ibea) + smom(:)/m(ibea) !veloc: v(t+dt)
    enddo

    call eval_kin_kt_pres

    !evaluate total momentum
    tmom(:) = 0.0_dp
    do ibea = 1, nbea
      tmom(:) = tmom(:) + m(ibea)*v(:,ibea)
    enddo

  end subroutine integr_vv_nvt


  !----------------------------------------------
  subroutine integr_vv_npt
  !----------------------------------------------
  ! velocity-Verlet algorithm, NPxxPyyPzzT
  ! forces at r(t+dt), v(t+dt/2+vvext*dt/2)
  ! extension of NVT algorithm of 
  ! Groot and Warren, J.Chem.Phys.107, 4423 (1997)
  !----------------------------------------------

    implicit none
    integer(si) :: ibea, icyl
    real(dp) :: dth, qpmass, volp, voln
    real(dp) :: dist, rij, dvz, drlen, b, c, discr, tail
    real(dp), dimension(3) :: xipdthv
    real(dp), dimension(3) :: shft
    real(dp), dimension(3) :: xipnv, boxnv, boxsclv, velsclv
    real(dp), dimension(3) :: rijv, drhv, tailv, cenintv, ceninthv, flipv, dvxy
    real(dp), dimension(2) :: root
    logical :: reflected

    if (barostat) then
      !Qp mass is 1/3 of the one for isotropic volume NPT formalism
      qpmass = real((nbea+1),dp)*ktfix*(ptau)**2
    else
      xipv(:) = 0.0_dp
    endif

    !shorthands
    dth = dt/2.0_dp
    xipdthv(:) = xipv(:)*dth

    !v(t) -> v(t+dt/2)
    !for: f(t)/m*dth calculated on the previous step
    do ibea = 1, nbea
      v(:,ibea) = v(:,ibea)*(1.0_dp-xipdthv(:)) + for(:,ibea)
    enddo

    !box scale parameters due to barostat
    if (barostat) then
      xipnv(:) = xippv(:) + 6.0_dp*vol/qpmass*(presv(:)-pxyz(:))*dt
      if (isoani=="ANI") then
        boxnv(:) = boxpv(:) + 2.0_dp*box(:)*xipv(:)*dt
        boxsclv(:) = boxnv(:)/box(:)
        velsclv(:) = 2.0_dp*boxnv(:)/(box(:)+boxnv(:))
      else
        volp = boxpv(1)*boxpv(2)*boxpv(3)
        voln = volp + 6.0_dp*vol*xipv(1)*dt
        boxsclv(:) = (voln/vol)**(1.0_dp/3.0_dp)
        boxnv(:) = box(:)*boxsclv(:)
        velsclv(:) = (2.0_dp*voln/(vol+voln))**(1.0_dp/3.0_dp)
      endif
    elseif (berendsen) then
      if (isoani=="ANI") then
        boxsclv(:) = 1.0_dp - 1.0_dp/ptau*(pxyz(:)-presv(:))*dt
      else
        boxsclv(:) = 1.0_dp - 1.0_dp/ptau*(pfix/3.0_dp-pres/3.0_dp)*dt
      endif
      velsclv(:) = 1.0_dp
    endif

    !transform to universal shifted coords -box/2..+box/2
    !these must be scaled, not 0..box coords
    shft(:) = -box(:)/2.0_dp
    do ibea = 1, nbea
      r(:,ibea) = r(:,ibea) + shft(:)
    enddo

    !change box dimensions
    if (barostat) then
      boxpv(:) = box(:)
      box(:) = boxnv(:)
    elseif (berendsen) then
      box(:) = box(:)*boxsclv(:)
    endif
    !new shift to be applied after full step for coords
    shft(:) = -box(:)/2.0_dp

    !preliminary advances (forces at time t): r(t) -> r(t+dt)
    smom(:) = 0.0_dp  !clear surface momenta
    do ibea = 1, nbea
      !apply eqs for box coords and transform to new 0..box coords
      r(:,ibea) = r(:,ibea)*boxsclv(:) + v(:,ibea)*velsclv(:)*dt - shft(:)
      !PBC, dtop=2 for surface, =3 otherwise
      if (any (r(1:dtop,ibea)<0.0_dp .or. r(1:dtop,ibea)>box(1:dtop))) then
        r(1:dtop,ibea) = merge(r(1:dtop,ibea)+box(1:dtop),r(1:dtop,ibea),r(1:dtop,ibea)<0.0_dp)
        r(1:dtop,ibea) = merge(r(1:dtop,ibea)-box(1:dtop),r(1:dtop,ibea),r(1:dtop,ibea)>=box(1:dtop))
      endif

      !check for reflections out of cylinders and walls
      do
        !check for crossing cylinders surfaces
        reflected = .false.
        if (cylinders) then
          do icyl = 1, ncyl
            !check for crossing cylinder's surface
            rijv(1:2) = r(1:2,ibea)-rcylc(1:2,icyl)  !2D vector from the circle center to r(t+dt)
            rijv(1:2) = rijv(1:2) - box(1:2)*nint(rijv(1:2)/box(1:2))  !PBC
            rij = sqrt(sum(rijv(1:2)**2))  !2D distance to the circle center
            if (rij < cylr(icyl)) then
              !apply reflection out of a tangential plane at the intersection point
              drlen = sqrt(sum(v(1:2,ibea)**2))
              drhv(1:2) = v(1:2,ibea)/drlen  !2D unit vector along delta_r=r(t+dt)-r(t) [=v(t+dt/2)*dt]
              b = -2.0_dp*sum(drhv(1:2)*rijv(1:2))   !coefficients of quadratic equation
              c = sum(rijv(1:2)**2) - cylr(icyl)**2  !for delta_r tail inside a cylinder
              discr = b**2 - 4.0_dp*c                !discriminant
              if (discr<0.0_dp) call error_exit("Error in reflecting off cylinder: D<0")
              root(1) = (-b-sqrt(discr))/2.0_dp  !roots
              root(2) = (-b+sqrt(discr))/2.0_dp
              if (all(root(1:2)>0.0_dp)) call error_exit("Error in reflecting off cylinder: both roots>0")
              tail = merge(root(1),root(2),root(1)>0.0_dp)
              tailv(1:2) = tail*drhv(1:2)  !2D vector from 2D intersection point to r(t+dt)
              cenintv(1:2) = rijv(1:2) - tailv(1:2)  !2D vector from the circle center to intersection point
              drlen = sqrt(sum(cenintv(1:2)**2))
              ceninthv(1:2) = cenintv(1:2)/drlen  !2D unit vector along cenintv(1:2)
              flipv(1:2) = -2.0_dp*(sum(ceninthv(1:2)*tailv(1:2)))*ceninthv(1:2)  !2D flipping vector
              !apply reflection off the cylinder
              r(1:2,ibea) = r(1:2,ibea) + flipv(1:2)  !flip tail, XY coord
              dvxy(1:2) = flipv(1:2)/dt  !contribution to the XY velocity
              v(1:2,ibea) = v(1:2,ibea) + dvxy(1:2)  !velocity update
              smom(1:2) = smom(1:2)-m(ibea)*dvxy(1:2)  !momentum asquired by the cylinder
              reflected = .true.
            endif
          enddo
        endif
        !check for the crossing the walls
        if (surface) then
          if (r(3,ibea)<0.0_dp) then
            dist = -r(3,ibea)          !distance from the bottom wall
            dvz = 2.0_dp*dist/dt       !vel shift on z when flipping
            r(3,ibea) = dist           !reflection of z coord
            v(3,ibea) = v(3,ibea)+dvz  !v(t+dt/2) after reflection
            smom(3) = smom(3)-m(ibea)*dvz  !momentum asquired by the wall
            reflected = .true.
          elseif (r(3,ibea)>box(3)) then
            dist = r(3,ibea)-box(3)    !distance from the upper wall
            dvz = -2.0_dp*dist/dt      !vel shift on z when flipping
            r(3,ibea) = box(3) - dist  !reflection of z coord
            v(3,ibea) = v(3,ibea)+dvz  !v(t+dt/2) after reflection
            smom(3) = smom(3)-m(ibea)*dvz  !momentum asquired by the wall
            reflected = .true.
          endif
        endif
        if (.not.reflected) exit  !no reflection occurred from both cylinders and surfaces 
      enddo  !end for reflections check
    enddo  !end for beads loop

    !rename Verlet-Stoermer integrated vars for the next step
    if (barostat) then
      xippv(:) = xipv(:)
      xipv(:)  = xipnv(:)
    endif

    !update volume and density
    if (barostat.or.berendsen) then
      vol  = box(1)*box(2)*box(3)
      dens = nbea/vol
    endif

    !v(t+dt/2) -> approx. v(t+dt)
    do ibea = 1, nbea
      vsave(:,ibea) = v(:,ibea)                 !save   v(t+dt/2)
      v(:,ibea) = v(:,ibea) + vvext*for(:,ibea) !approx.v(t+dt)
    enddo

    !calculate forces and gorques on running step st_n
    call eval_for_vir

    !update xipdth
    xipdthv(:) = xipv(:)*dth

    !v(t+dt/2) -> v(t+dt)
    !for: f(t+dt)/m*dth just calculated
    smom(:) = smom(:)/real(nbea,dp)  !surface: contribution to each bead from the wall momentum
    !TRY TO CHECK HERE DO WE SCALE OR NOT WALLS MOMENTUM CONTRUBUTION
    do ibea = 1, nbea
      v(:,ibea) = (vsave(:,ibea) + for(:,ibea))/(1.0_dp+xipdthv(:)) + smom(:)/m(ibea)
    enddo

    call eval_kin_kt_pres

    !evaluate total momentum
    tmom(:) = 0.0_dp
    do ibea = 1, nbea
      tmom(:) = tmom(:) + m(ibea)*v(:,ibea)
    enddo

  end subroutine integr_vv_npt


  !===============================================
  !=============== INTERACTIONS ==================
  !===============================================


  !----------------------------------------------
  subroutine eval_cdr_uni (ibea, jbea, rijv, vijv, rij_2)
  !----------------------------------------------
  !unifies eval_beabea, eval_beasurf, eval_beacyl
  !ibea - current bead
  !jbea - >0: regular bead
  !       =0: surface imaginary bead
  !           (on wall or cylinder surface)
  !---------------------------------------------- 

    implicit none
    integer(si) :: ibea, jbea
    real(dp), dimension(3) :: rijv, vijv
    real(dp) :: rij_2

    integer(si) :: i, j
    real(dp) :: rij, wr, wd, theta
    real(dp), dimension(3) :: rhij, fcij, fdij, frij, fij
    logical :: same

    !evaluate required data
    rij = sqrt(rij_2)
    rhij(:) = rijv(:)/rij
    wr  = 1.0_dp - rij/rc
    wd  = wr*wr
    call randg (theta)

    !conservative force
    !rightmost 4 bits of type are used for repulsive type coding, surface bead got type 0
    i = iand(15,t(ibea))
    j = iand(15,t(jbea))
    fcij(:) = prep(i,j)*wr*rhij(:)

    !dissipative force
    fdij(:) = -gamma*wd*sum(rhij(:)*vijv(:))*rhij(:)

    !random force
    frij(:) = sigma*wr*theta/sqrt(dt)*rhij(:)

    !total pair force
    fij(:) = fcij(:) + fdij(:) + frij(:)

    !add contributions to the forces on i,j and virial
    for(:,ibea) = for(:,ibea) + fij(:)
    for(:,jbea) = for(:,jbea) - fij(:)  !for(:,0) accumulates force on surfaces
    virv(:) = virv(:) + fcij(:)*rijv(:)

  end subroutine eval_cdr_uni


  !----------------------------------------------------
  subroutine eval_attr_uni (ibea, jbea, rijv, rij_2)
  !----------------------------------------------------

    implicit none
    integer(si) :: ibea, jbea
    real(dp), dimension(3) :: rijv
    real(dp) :: rij_2

    integer(si) :: i, j
    real(dp) :: rij, rm1
    real(dp), dimension(3) :: rhij, faij

    !evaluate required data
    rij = sqrt(rij_2)
    rm1 = rij - 1.0_dp
    rhij(:) = rijv(:)/rij

    !rightmost 4 bits of type are used for repulsive type coding, surface bead got type 0
    !attractive force
    i = iand(15,t(ibea))
    j = iand(15,t(jbea))
    faij(:) = -peps(i,j)*twodxi*rm1*exp(-rm1*rm1/xi)*rhij(:)

    !add contributions to the forces on i,j and virial
    for(:,ibea) = for(:,ibea) + faij(:)
    for(:,jbea) = for(:,jbea) - faij(:)  !for(:,0) accumulates force on surfaces
    virv(:) = virv(:) + faij(:)*rijv(:)

  end subroutine eval_attr_uni


  !----------------------------------------------------
  subroutine attempt_crosslink (ibea, jbea)
  !----------------------------------------------------

    implicit none
    integer(si) :: ibea, jbea
    integer(si) :: bi

    !add crosslinking bond if criteria is met
    if (rtype=="SQW") then
      !square well, link is formed with even probability "rprob" within [0;cc]
      if (rprob > ran3(idum)) then
        cbon = cbon + 1
        bi = nbon + cbon
        b1(bi) = ibea
        b2(bi) = jbea
        b0(bi) = rb0
        bk(bi) = rbk
        !change flag for ibea and jbea for being linked
        a(ibea) = 2        
        a(jbea) = 2
      endif
    endif

  end subroutine attempt_crosslink


  !--------------------------------------------------
  subroutine eval_bond_uni (ibea, jbea, b0, bk, mode)
  !--------------------------------------------------
  ! ibea: first bead
  ! jbea: second bead (for mode==1), igrf (for mode==2)
  ! b0, bk  - harmonic spring parameters
  !--------------------------------------------------

    implicit none
    integer(si) :: ibea, jbea
    real(dp) :: b0, bk
    integer(si) :: mode

    real(dp), dimension(3) :: rijv, rfix, fij
    real(dp) :: rij

    !set rijv(:)
    if (mode==1) then
      rijv(:) = r(:,ibea) - r(:,jbea)  !normal bond
    else
      rfix(1:2) = grfr(1:2,jbea)  !grafting bond
      rfix(3) = merge(0.0_dp, box(3), grfw(jbea)==0)
      rijv(:) = r(:,ibea) - rfix(:)
      if (grft(jbea)==2) rijv(1:2) = 0.0_dp
    endif

    !PBC and rij
    rijv(1:dtop) = rijv(1:dtop) - box(1:dtop)*nint(rijv(1:dtop)/box(1:dtop))
    rij = sqrt(sum(rijv(:)**2))

    !force on i
    fij(:) = -bk*(rij-b0)*rijv(:)/rij

    !add contribution to the forces and virial
    for(:,ibea) = for(:,ibea) + fij(:)
    if (mode==1) then
      for(:,jbea) = for(:,jbea) - fij(:)
    else
      for(:,0) = for(:,0) - fij(:)
    endif
    virv(:) = virv(:) + fij(:)*rijv(:)
   
  end subroutine eval_bond_uni


  !----------------------------------------------
  subroutine eval_ang (ibea, jbea, kbea, a0, ak)
  !----------------------------------------------
  ! mode=1: angles, mode=2: z-angles
  !----------------------------------------------

    implicit none
    integer(si) :: ibea, jbea, kbea
    real(dp) :: a0, ak

    real(dp), dimension(3) :: rjiv, rkjv
    real(dp) :: rji_2,  rkj_2
    real(dp) :: angle,  angfac, en
    real(dp) :: cjiji,  ckjkj,  ckjji
    real(dp) :: cfact1, cfact2, cfact3
    real(dp), dimension(3) :: fi, fj, fk

    !get vectors
    rjiv(:) = r(:,jbea) - r(:,ibea)
    rkjv(:) = r(:,kbea) - r(:,jbea)
    !rjiv(:) = rjiv(:) - box(:)*nint (rjiv(:)/box(:))
    !rkjv(:) = rkjv(:) - box(:)*nint (rkjv(:)/box(:))
    !surface
    rjiv(1:dtop) = rjiv(1:dtop) - box(1:dtop)*nint (rjiv(1:dtop)/box(1:dtop))
    rkjv(1:dtop) = rkjv(1:dtop) - box(1:dtop)*nint (rkjv(1:dtop)/box(1:dtop))
    !vector's lengths
    rji_2 = sum (rjiv(:)**2)
    rkj_2 = sum (rkjv(:)**2)
    !cos of angle
    angle = (sum (rjiv(:)*rkjv(:)))/sqrt(rji_2*rkj_2)
    !correct possible rounding effects
    if (abs(angle) > 1.0_dp) angle = angle - sign (0.0000001_dp, angle)
    angle = acos (angle)
    !get bending energy
    en = ak * (angle - a0)**2 / 2.0_dp

    !set angfac by resolving possible 0/0 situations
    if (abs (angle)<0.000001 .and. a0<0.0001) then
      angfac = - ak
    elseif (abs (angle - pi)<0.000001 .and. abs(a0)>3.14_dp) then
      angfac = ak
    else
      angfac = (-ak * (angle - a0)) / sin(angle)
    endif     

    !shorthands
    cjiji = rji_2
    ckjkj = rkj_2
    ckjji = sum (rkjv(:)*rjiv(:))
    cfact1 = angfac / sqrt (ckjkj * cjiji)
    cfact2 = ckjji / ckjkj
    cfact3 = ckjji / cjiji

    !forces acting on each bead
    fk(:) = cfact1 * (cfact2 * rkjv(:) - rjiv(:))
    fi(:) = cfact1 * (rkjv(:) - cfact3 * rjiv(:))
    fj(:) = -fi(:) - fk(:)

    !add contributions to the forces
    for(:,ibea) = for(:,ibea) + fi(:)
    for(:,jbea) = for(:,jbea) + fj(:)
    for(:,kbea) = for(:,kbea) + fk(:)

    !add contributions to the virial
    ! virv(:) = virv(:) + r(:,ibea)*fi(:) + &
    !          (r(:,ibea)+rjiv(:))*fj(:) + (r(:,ibea)+rjiv(:)+rkjv(:))*fk(:)

  end subroutine eval_ang


  !----------------------------------------------
  subroutine eval_tors (ibea, jbea, kbea, lbea, tc)
  !----------------------------------------------
  ! option 3 of gbmoldd is left only 
  !----------------------------------------------

    implicit none
    integer(si) :: ibea, jbea, kbea, lbea
    real(dp), dimension(:) :: tc

    real(dp), dimension(3) :: rjiv, rkjv, rlkv
    real(dp) :: rji_2, rkj_2, rlk_2
    real(dp) :: rji,   rkj,   rlk
    real(dp), dimension(3) :: plkkjv, pkjjiv
    real(dp) :: plkkj, pkjji
    real(dp) :: cosphi, cosphi_2, cosphi_3
    real(dp) :: de
    real(dp) :: clklk, clkkj, clkji, ckjkj, ckjji, cjiji
    real(dp) :: dlkkj, dkjji
    real(dp) :: cfact1, dfact1, dfact2, dfact3
    real(dp), dimension(3) :: fi, fj, fk, fl, ft

    !rji, rkj, rlk vectors
    rjiv(:) = r(:,jbea) - r(:,ibea)
    rkjv(:) = r(:,kbea) - r(:,jbea)
    rlkv(:) = r(:,lbea) - r(:,kbea)
    !rjiv(:) = rjiv(:) - box(:) * nint (rjiv(:)/box(:))
    !rkjv(:) = rkjv(:) - box(:) * nint (rkjv(:)/box(:))
    !rlkv(:) = rlkv(:) - box(:) * nint (rlkv(:)/box(:))
    !surface
    rjiv(1:dtop) = rjiv(1:dtop) - box(1:dtop) * nint (rjiv(1:dtop)/box(1:dtop))
    rkjv(1:dtop) = rkjv(1:dtop) - box(1:dtop) * nint (rkjv(1:dtop)/box(1:dtop))
    rlkv(1:dtop) = rlkv(1:dtop) - box(1:dtop) * nint (rlkv(1:dtop)/box(1:dtop))
    !vector's lengths
    rji_2 = sum (rjiv(:)**2)
    rkj_2 = sum (rkjv(:)**2)
    rlk_2 = sum (rlkv(:)**2)
    rji = sqrt (rji_2)
    rkj = sqrt (rkj_2)
    rlk = sqrt (rlk_2)

    !vector products
    plkkjv(1) = rlkv(2)*rkjv(3) - rlkv(3)*rkjv(2)
    plkkjv(2) = rlkv(3)*rkjv(1) - rlkv(1)*rkjv(3)
    plkkjv(3) = rlkv(1)*rkjv(2) - rlkv(2)*rkjv(1)
    pkjjiv(1) = rkjv(2)*rjiv(3) - rkjv(3)*rjiv(2)
    pkjjiv(2) = rkjv(3)*rjiv(1) - rkjv(1)*rjiv(3)
    pkjjiv(3) = rkjv(1)*rjiv(2) - rkjv(2)*rjiv(1)
    !their lenghts
    plkkj = sqrt (sum (plkkjv(:)**2))
    pkjji = sqrt (sum (pkjjiv(:)**2))

    !cosphi as in Allen-Tildesley (0 for all-trans)
    cosphi = - (sum (plkkjv(:)*pkjjiv(:))) / (plkkj*pkjji)
    !correct possible rounding effects
    if (abs(cosphi) > 1.0_dp) cosphi = cosphi - sign (0.0000001_dp, cosphi)
    cosphi_2 = cosphi**2
    cosphi_3 = cosphi_2*cosphi

    !calculate the torsional energy
    !en = tor + tc(1)*cosphi + tc(2)*cosphi_2 + tc(3)*cosphi_3

    !calculate force parameters
    clklk = rlk_2
    clkkj = sum (rlkv(:)*rkjv(:))
    clkji = sum (rlkv(:)*rjiv(:))
    ckjkj = rkj_2
    ckjji = sum (rkjv(:)*rjiv(:))
    cjiji = rji_2
    dlkkj = clklk*ckjkj - clkkj*clkkj
    dkjji = ckjkj*cjiji - ckjji*ckjji

    de = - (tc(1) + 2.0_dp*tc(2)*cosphi + 3.0_dp*tc(3)*cosphi_2)

    cfact1 = (clkkj*ckjji - clkji*ckjkj)
    dfact1 = -de / sqrt (dlkkj*dkjji)
    dfact2 = -cfact1/dkjji
    dfact3 = -cfact1/dlkkj

    !get contributions to the torsion forces acting on i,j,k,l
    fl(:) = dfact1 * ( ckjji*rkjv(:) - ckjkj*rjiv(:) + &
                       dfact3*(ckjkj*rlkv(:) - clkkj*rkjv(:)) )
    fi(:) = dfact1 * ( -clkkj*rkjv(:) + ckjkj*rlkv(:) + &
                       dfact2*(-ckjkj*rjiv(:) + ckjji*rkjv(:)) )
    ft(:) = dfact1 * ( ckjji*rlkv(:) + clkkj*rjiv(:) - 2.0_dp*clkji*rkjv(:) + &
                       dfact2*(cjiji*rkjv(:) - ckjji*rjiv(:)) + &
                       dfact3*(clklk*rkjv(:) - clkkj*rlkv(:)))
    fj(:) = -ft(:) - fi(:)
    fk(:) =  ft(:) - fl(:)

    !add contributions to the torsion forces acting on i,j,k,l
    for(:,ibea) = for(:,ibea) + fi(:)
    for(:,jbea) = for(:,jbea) + fj(:)
    for(:,kbea) = for(:,kbea) + fk(:)
    for(:,lbea) = for(:,lbea) + fl(:)

    !add contribution to the virial
    !virv(:) = virv(:) + r(:,ibea)*fi(:) + &
    !    (r(:,ibea)+rjiv(:))*fj(:) + (r(:,ibea)+rjiv(:)+rkjv(:))*fk(:) + &
    !    (r(:,ibea)+rjiv(:)+rkjv(:)+rlkv(:))*fl(:)

  end subroutine eval_tors


  !===============================================
  !===================== MISC ====================
  !===============================================

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

  !---------------------------------------------------------
  subroutine randg (x)
  !routine to generate random number 'x' distributed according to
  !Gaussian distribution with mean 0 and standard deviation 1
  !---------------------------------------------------------

    implicit none
    real(dp) :: x
    real(dp) :: v1, v2, r

    do
      v1 = 2.0*ran3(idum)-1.0
      v2 = 2.0*ran3(idum)-1.0
      r = v1*v1+v2*v2
      if (r < 1.0) exit
    enddo
    x = v1*sqrt(-2.0*log(r)/r)

  end subroutine randg

  !---------------------------------------------------------
  subroutine randg_ext(x,mean,sig)
  !routine to generate random number 'x' distributed according to
  !Gaussian distribution with mean 'mean' and standard deviation 'sig'
  !---------------------------------------------------------

    implicit none
    real(dp) :: x, mean, sig
    real(dp) :: v1, v2, r

    do
      v1 = 2.0*ran3(idum)-1.0
      v2 = 2.0*ran3(idum)-1.0
      r = v1*v1+v2*v2
      if (r < 1.0) exit
    enddo
    x = v1*sqrt(-2.0*log(r)/r)
    x = mean + x*sig 

  end subroutine randg_ext

  !----------------------------------------------
  subroutine error_exit (str)
  !----------------------------------------------

    character(len=*) :: str

    write (*,fmt='(a,i8,a)') "Step",st_n,": ",trim(str)
    call exit()

  end subroutine error_exit

  end program dpd_2_4
