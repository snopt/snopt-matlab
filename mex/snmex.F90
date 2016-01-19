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
  double precision, allocatable :: riGfun(:), rjGvar(:)

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

contains

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
