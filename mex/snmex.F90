!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! mxsnWork module for SNOPT Fortran mex
!
! 13 Sep 2013: First version.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "fintrf.h"

module mxsnWork
  implicit none
  public

  ! SNOPT workspace
  integer            :: leniw = 5000, lenrw = 5000, lencw = 500
  integer,          allocatable :: iw(:), iw0(:)
  double precision, allocatable :: rw(:), rw0(:)
  character*8,      allocatable :: cw(:), cw0(:)


  ! SNOPT mex variables
  logical              :: firstCall = .true.,  &
                          memCall   = .false., &
                          printOpen = .false., &
                          summOpen  = .false., &
                          screenON  = .false.

  integer              :: callType = 0

  mwPointer            :: fgHandle, HxHandle
  mwPointer            :: objHandle, conHandle ! for fmincon-style
  mwPointer            :: stopHandle

  integer              :: neG


  integer, parameter   :: iPrint     = 9, iSpecs   = 4, iSumm    = 55, &
                          systemCall = 0, userCall = 1
  integer, parameter   :: snSolve    = 1,  &
                          snSetXX    = 2,  &
                          snSetIX    = 3,  &
                          snSetRX    = 4,  &
                          snGetXX    = 5,  &
                          snGetCX    = 6,  &
                          snGetIX    = 7,  &
                          snGetRX    = 8,  &
                          snSpecs    = 9,  &
                          snOpenP    = 10, &
                          snOpenS    = 11, &
                          snClosP    = 12, &
                          snClosS    = 13, &
                          snSetWork  = 14, &
                          snscrnON   = 15, &
                          snscrnOff  = 16, &
                          snFindJac  = 17, &
                          snSolveN   = 18, &
                          snEnd      = 999

  ! SNOPT arrays
  integer,          allocatable :: xstate(:), Fstate(:)
  double precision, allocatable :: x(:), xmul(:), xlow(:), xupp(:), &
                                   F(:), Fmul(:), Flow(:), Fupp(:)

  double precision, allocatable :: rtmpa(:), G1(:)

  integer,          allocatable :: iAfun(:), jAvar(:)
  double precision, allocatable :: riAfun(:), rjAvar(:), A(:)

  integer,          allocatable :: iGfun(:), jGvar(:)
  double precision, allocatable :: riGfun(:), rjGvar(:)

  integer,          allocatable :: tFstate(:)
  double precision, allocatable :: tF(:), tFmul(:), tFlow(:), tFupp(:)


  ! SQOPT arrays
  integer,          allocatable :: hEtype(:), hs(:), indA(:), locA(:)
  double precision, allocatable :: cObj(:), pi(:), rc(:), bl(:), bu(:), &
                                   valA(:), rlocA(:), rindA(:)

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine resetSNOPT
    !---------------------------------------------------------------------------
    ! resetSNOPT for new problem.
    ! Reset variables, deallocate all arrays.
    ! (Also registered with mexAtExit to deallocate workspace and close files)
    !---------------------------------------------------------------------------

    !    if (printOpen) close(iPrint)
    close(iPrint)
    printOpen = .false.

    !    if (summOpen) close(iSumm)
    close(iSumm)
    summOpen  = .false.

    close(iSpecs)

    firstCall = .true.
    memCall   = .false.

    leniw     = 5000
    lenrw     = 5000
    lencw     = 500

    call deallocSQOPT

    call deallocSNOPT
    call deallocA
    call deallocG
    call deallocF

    if (allocated(cw)) deallocate(cw)
    if (allocated(iw)) deallocate(iw)
    if (allocated(rw)) deallocate(rw)

    if (allocated(cw0)) deallocate(cw0)
    if (allocated(iw0)) deallocate(iw0)
    if (allocated(rw0)) deallocate(rw0)

  end subroutine resetSNOPT


  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine allocSNOPT( n, nF )
    integer, intent(in) :: n, nF
    !---------------------------------------------------------------------------
    ! Allocate space for SNOPT solve.
    !---------------------------------------------------------------------------

    call deallocSNOPT

    allocate(rtmpa(max(n,nF)))
    allocate(x(n),  xlow(n),  xupp(n),  xmul(n),  xstate(n))
    allocate(F(nF), Flow(nF), Fupp(nF), Fmul(nF), Fstate(nF))

  end subroutine allocSNOPT

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine deallocSNOPT
    !---------------------------------------------------------------------------
    ! Deallocate x,F arrays involved in solve routine.
    !---------------------------------------------------------------------------

    call deallocR(rtmpa)

    call deallocR(x)
    call deallocR(xmul)
    call deallocR(xlow)
    call deallocR(xupp)
    call deallocI(xstate)

    call deallocR(F)
    call deallocR(Fmul)
    call deallocR(Flow)
    call deallocR(Fupp)
    call deallocI(Fstate)

  end subroutine deallocSNOPT

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine allocA( lenA )
    integer, intent(in) :: lenA
    !---------------------------------------------------------------------------
    ! Allocate space for A
    !---------------------------------------------------------------------------

    call deallocA

    allocate(iAfun(lenA), jAvar(lenA), A(lenA))
    allocate(riAfun(lenA), rjAvar(lenA))

  end subroutine allocA

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine deallocA
    !---------------------------------------------------------------------------
    ! Deallocate space for A
    !---------------------------------------------------------------------------

    call deallocI(iAfun)
    call deallocI(jAvar)
    call deallocR(A)

    call deallocR(riAfun)
    call deallocR(rjAvar)

  end subroutine deallocA

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine allocG( lenG )
    integer, intent(in) :: lenG
    !---------------------------------------------------------------------------
    ! Allocate space for G
    !---------------------------------------------------------------------------

    call deallocG

    allocate( G1(lenG))
    allocate( iGfun(lenG),  jGvar(lenG))
    allocate(riGfun(lenG), rjGvar(lenG))

  end subroutine allocG

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine deallocG
    !---------------------------------------------------------------------------
    ! Deallocate space for G
    !---------------------------------------------------------------------------

    call deallocR(G1)

    call deallocI(iGfun)
    call deallocI(jGvar)

    call deallocR(riGfun)
    call deallocR(rjGvar)

  end subroutine deallocG

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine allocJac( n, lenA, lenG )
    integer, intent(in) :: n, lenA, lenG
    !---------------------------------------------------------------------------
    ! Allocate space for snJac.
    !---------------------------------------------------------------------------

    call deallocR(x)
    call deallocR(xlow)
    call deallocR(xupp)
    allocate( x(n), xlow(n), xupp(n) )

    call allocA(lenA)
    call allocG(lenG)

  end subroutine allocJac

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine deallocJac
    !---------------------------------------------------------------------------
    ! Deallocate space for snJac.
    !---------------------------------------------------------------------------

    call deallocR(x)
    call deallocR(xlow)
    call deallocR(xupp)

    call deallocI(iAfun)
    call deallocI(jAvar)

    call deallocI(iGfun)
    call deallocI(jGvar)

  end subroutine deallocJac

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine allocF( nF )
    integer, intent(in) :: nF
    !---------------------------------------------------------------------------
    ! Allocate F space for snSTOP.
    !---------------------------------------------------------------------------

    call deallocF
    allocate( tF(nF), tFstate(nF), tFmul(nF), tFlow(nF), tFupp(nF) )

  end subroutine allocF

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine deallocF
    !---------------------------------------------------------------------------
    ! Deallocate F space for snSTOP.
    !---------------------------------------------------------------------------

    if (allocated(tF))      deallocate(tF)
    if (allocated(tFmul))   deallocate(tFmul)
    if (allocated(tFlow))   deallocate(tFlow)
    if (allocated(tFupp))   deallocate(tFupp)
    if (allocated(tFstate)) deallocate(tFstate)

  end subroutine deallocF

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine deallocI( array )
    integer, allocatable :: array(:)

    if (allocated(array)) deallocate(array)

  end subroutine deallocI

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine deallocR( array )
    double precision, allocatable :: array(:)

    if (allocated(array)) deallocate(array)

  end subroutine deallocR

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine deallocSQOPT

    if (allocated(cObj))   deallocate(cObj)

    if (allocated(x))      deallocate(x)
    if (allocated(pi))     deallocate(pi)
    if (allocated(rc))     deallocate(rc)
    if (allocated(bl))     deallocate(bl)
    if (allocated(bu))     deallocate(bu)
    if (allocated(hs))     deallocate(hs)
    if (allocated(hEtype)) deallocate(hEtype)

    if (allocated(indA))   deallocate(indA)
    if (allocated(locA))   deallocate(locA)
    if (allocated(valA))   deallocate(valA)

    if (allocated(rindA))  deallocate(rindA)
    if (allocated(rlocA))  deallocate(rlocA)

  end subroutine deallocSQOPT

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine checkCol ( pm, n, name )
    mwPointer     :: pm
    integer       :: n
    character*(*) :: name
    !---------------------------------------------------------------------------
    ! Check column dimension of pm is equal to n.
    !---------------------------------------------------------------------------
    character*80 :: str
    mwSize       :: m, mxGetN

    m = mxGetN(pm)
    if ( m /= n ) then
       write(str,100) name, m, n
       call mexErrMsgTxt ( str )
    end if

    return

100 format ( a, ' has incorrect column dimension ', i5, &
                '.  Should be length ', i5 )

  end subroutine checkCol

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine checkRow ( pm, n, name )
    character*(*) :: name
    mwPointer     :: pm
    integer       :: n
    !---------------------------------------------------------------------------
    ! Check row dimension of pm is equal to n.
    !---------------------------------------------------------------------------
    character*80 :: str
    mwSize       :: m, mxGetM

    m = mxGetM(pm)
    if ( m /= n ) then
       write(str,100) name, m, n
       call mexErrMsgTxt ( str )
    end if

    return

100 format ( a, ' has incorrect row dimension ', i5, &
                '.  Should be length ', i5 )

  end subroutine checkRow

end module mxsnWork
