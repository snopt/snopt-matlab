!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File snoptmex.F90
! Mex function for SNOPT7.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "fintrf.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine mexFunction ( nlhs, plhs, nrhs, prhs )
  use mxsnWork

  implicit none

  integer*4  :: nlhs, nrhs
  mwPointer  :: prhs(*), plhs(*)
  !=====================================================================
  ! Mex function for SNOPT7
  !
  ! Option      Action
  !    1        Solve the problem (fmincon-style)
  !    2        Set option
  !    3        Set integer option
  !    4        Set real option
  !    5        Get option
  !    6        Get character option
  !    7        Get integer option
  !    8        Get real option
  !    9        Read specs file
  !    10       Openprint file
  !    11       Opensummary file
  !    12       Closeprint file
  !    13       Closesummary file
  !    14
  !    15       Screen on
  !    16       Screen off
  !    17       Compute the Jacobian structure via snJac
  !
  ! 18 Sep 2013: Current version.
  ! 01 May 2015: Added ability to modify initial amount of workspace
  !=====================================================================
  ! Matlab
  mwPointer        :: mxGetN, mxGetPr
  integer          :: mxIsChar
  double precision :: mxGetScalar

  ! SNOPT
  character        :: filename*80
  integer          :: info, iOpt, strlen
  double precision :: rOpt, rleniw, rlenrw
  external         :: snInit

  ! Get option.
  if ( nrhs < 1 ) call mexErrMsgTxt( 'Need an option input argument' )
  rOpt = mxGetScalar(prhs(1))
  iOpt = rOpt


  ! Deal with on/off screen, print, and summary files first.
  if ( iOpt == snOpenP ) then

     if ( nrhs /= 2 ) call mexErrMsgTxt( 'Wrong number of input arguments' )

     info = mxIsChar(prhs(2))
     if ( info /= 1 ) call mexErrMsgTxt( 'Need a filename string' )

     strlen = mxGetN(prhs(2))
     if ( strlen > 80 ) call mexErrMsgTxt( 'Print filename is too long' )

     if ( strlen > 0 ) then
        call mxGetString( prhs(2), filename, strlen )
     else
        call mexErrMsgTxt( 'Empty print filename' )
     end if

     if ( printOpen) close( iPrint )

     open( iPrint, file=filename, status='unknown' )
     printOpen= .true.
     return

  else if ( iOpt == snOpenS ) then

     if ( nrhs /= 2 ) call mexErrMsgTxt( 'Wrong number of input arguments' )

     info = mxIsChar(prhs(2))
     if ( info /= 1 ) call mexErrMsgTxt( 'Need a filename string' )

     strlen = mxGetN(prhs(2))
     if ( strlen > 80 ) call mexErrMsgTxt( 'Summary filename is too long' )

     if ( strlen > 0 ) then
        call mxGetString( prhs(2), filename, strlen )
     else
        call mexErrMsgTxt( 'Empty summary filename' )
     end if

     if ( summOpen) close( iSumm )

     open( iSumm, file=filename, status='unknown' )
     summOpen= .true.
     return

  else if ( iOpt == snClosP ) then
     if ( printOpen) close( iPrint )
     printOpen= .false.
     return

  else if ( iOpt == snClosS ) then
     if ( summOpen) close( iSumm )
     summOpen= .false.
     return

  else if ( iOpt == snscrnOn ) then
     screenOn = .true.
     return

  else if ( iOpt == snscrnOff ) then
     screenOn = .false.
     return

  else if ( iOpt == snsetwork ) then
     rleniw = mxGetScalar(prhs(2))
     rlenrw = mxGetScalar(prhs(3))
     leniw  = rleniw
     lenrw  = rlenrw

     if ( leniw < 500 .or. lenrw < 500 ) &
          call mexErrMsgTxt('Workspace size must be at least 500')
     return
  end if

  if ( firstCall ) then
     allocate( cw(lencw), iw(leniw), rw(lenrw) )

     callType = userCall
     call snInit ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )
     callType = systemCall

     memCall   = .false.
     firstCall = .false.
  end if

  ! Do whatever we need to do.
  if      ( iOpt == snSolve ) then

     callType = userCall
     call snmxSolve ( nlhs, plhs, nrhs, prhs )
     callType = systemCall

  else if ( iOpt == snFindJac ) then

     callType = userCall
     call snmxFindJac ( nlhs, plhs, nrhs, prhs )
     callType = systemCall

  else if ( iOpt == snSetXX .or. &
            iOpt == snSetIX .or. &
            iOpt == snSetRX .or. &
            iOpt == snGetXX .or. &
            iOpt == snGetCX .or. &
            iOpt == snGetIX .or. &
            iOpt == snGetRX ) then

     callType = userCall
     call snmxOptions ( iOpt, nlhs, plhs, nrhs, prhs )
     callType = systemCall

  else if ( iOpt == snSpecs ) then

     callType = userCall
     call snmxSpecs ( nlhs, plhs, nrhs, prhs )
     callType = systemCall

  else if ( iOpt == snEnd ) then
     if ( printOpen) close( iPrint )
     printOpen= .false.

     if ( summOpen) close( iSumm )
     summOpen = .false.
     firstCall = .true.
     memCall   = .false.

     if ( allocated(cw) ) deallocate( cw )
     if ( allocated(iw) ) deallocate( iw )
     if ( allocated(rw) ) deallocate( rw )

  end if

  return

end subroutine mexFunction

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snmxSolve ( nlhs, plhs, nrhs, prhs )
  use mxsnWork

  implicit none
  integer*4  :: nlhs, nrhs
  mwPointer  :: prhs(*), plhs(*)
  !---------------------------------------------------------------------
  ! Solve the problem
  ! The matlab call is
  !   [x, F, xmul, Fmul, info ] =
  !       snoptmex ( solveopt, x, xlow, xupp, xmul, xstate,
  !                            F, Flow, Fupp, Fmul, Fstate, ObjAdd,
  !                  ObjRow, A, iAfun, jAvar, iGfun, jGvar, userfun )
  !---------------------------------------------------------------------
  ! Matlab
  mwPointer        :: mxDuplicateArray, mxGetM, mxGetPr, &
                      mxCreateDoubleMatrix, mxCreateDoubleScalar
  integer          :: mxIsClass
  double precision :: mxGetScalar

  ! SNOPT
  character*8      :: probName
  integer          :: info, Errors
  integer          :: Start, ObjRow, n, nF, lenA, lenG, neA, tmp, &
                      nxname, nFname, mincw, miniw, minrw, nInf, nS
  double precision :: rinfo, ObjAdd, sInf
  external         :: snMemA, snoptA, matlabFG

  double precision, parameter   :: infBnd = 1.0d+20

  character*8,      allocatable :: xname(:), Fname(:)
  integer,          allocatable :: xstate(:), Fstate(:), &
                                   iAfun(:), jAvar(:), iGfun(:), jGvar(:)
  double precision, allocatable :: x(:), xmul(:), xlow(:), xupp(:), &
                                   F(:), Fmul(:), Flow(:), Fupp(:), A(:), &
                                   rtmp(:), riA(:), rjA(:)


  ! Check number of input and output arguments.
  if ( nrhs /= 18 .and. nrhs /= 19 ) call mexErrMsgTxt( 'Wrong number of input arguments' )


  !---------------------------------------------------------------------
  ! Compute number of variables, constraints, etc
  ! Check dimension of some of the input
  !---------------------------------------------------------------------
  ! Get number of variables and constraints
  n  = mxGetM(prhs(2))
  nF = mxGetM(prhs(7))

  allocate( x(n),  xlow(n),  xupp(n),  xmul(n),  xstate(n) )
  allocate( F(nF), Flow(nF), Fupp(nF), Fmul(nF), Fstate(nF) )
  allocate( rtmp(max(n,nF)) )

  ! Get initial x
  tmp = mxGetM(prhs(2))
  if ( tmp > 0 ) then
     call checkCol( prhs(2), 1, 'x' )
     call mxCopyPtrToReal8( mxGetPr(prhs(2)), x, n )
  else
     x(1:n) = 0.0
  end if

  ! Get lower bounds
  tmp = mxGetM(prhs(3))
  if ( tmp > 0 ) then
     call checkRow( prhs(3), n, 'xlow' )
     call checkCol( prhs(3), 1, 'xlow' )
     call mxCopyPtrToReal8( mxGetPr(prhs(3)), xlow, n )
  else
     xlow = -infBnd
  end if

  ! Get upper bounds
  tmp = mxGetM(prhs(4))
  if ( tmp > 0 ) then
     call checkRow( prhs(4), n, 'xupp' )
     call checkCol( prhs(4), 1, 'xupp' )
     call mxCopyPtrToReal8( mxGetPr(prhs(4)), xupp, n )
  else
     xupp = infBnd
  end if

  ! Get initial multipliers
  tmp = mxGetM(prhs(5))
  if ( tmp > 0 ) then
     call checkRow( prhs(5), n, 'xmul' )
     call checkCol( prhs(5), 1, 'xmul' )
     call mxCopyPtrToReal8( mxGetPr(prhs(5)), xmul, n )
  else
     xmul = 0
  end if

  ! Get initial states
  tmp = mxGetM(prhs(6))
  if ( tmp > 0 ) then
     call checkRow( prhs(6), n, 'xstate' )
     call checkCol( prhs(6), 1, 'xstate' )
     call mxCopyPtrToReal8( mxGetPr(prhs(6)), rtmp, n )
     xstate = rtmp(1:n)
  else
     xstate = 0
  end if

  tmp = mxGetM(prhs(7))
  if ( tmp > 0 ) then
     call checkCol( prhs(7),   1, 'Flow' )
     call mxCopyPtrToReal8( mxGetPr(prhs(7)),  Flow, nF )
  else
     Flow = -infBnd
  end if

  tmp = mxGetM(prhs(8))
  if ( tmp > 0 ) then
     call checkRow( prhs(8),  nF, 'Fupp' )
     call checkCol( prhs(8),   1, 'Fupp' )
     call mxCopyPtrToReal8( mxGetPr(prhs(8)),  Fupp, nF )
  else
     Fupp = infBnd
  end if

  tmp = mxGetM(prhs(9))
  if ( tmp > 0 ) then
     call checkRow( prhs(9),  nF, 'Fmul' )
     call checkCol( prhs(9),   1, 'Fmul' )
     call mxCopyPtrToReal8( mxGetPr(prhs(9)),  Fmul, nF )
  else
     Fmul = 0.0
  end if

  tmp = mxGetM(prhs(10))
  if ( tmp > 0 ) then
     call checkRow( prhs(10), nF, 'Fstate' )
     call checkCol( prhs(10),  1, 'Fstate' )
     call mxCopyPtrToReal8( mxGetPr(prhs(10)), rtmp, nF )
     Fstate = rtmp(1:nF)
  else
     Fstate = 0
  end if

  ObjAdd  = mxGetScalar(prhs(11))
  rtmp(1) = mxGetScalar(prhs(12))
  ObjRow  = rtmp(1)

  !---------------------------------------------------------------------
  ! Get the Jacobian structure (linear and nonlinear)
  !---------------------------------------------------------------------
  neA  = mxGetM ( prhs(13) )
  lenA = neA

  if ( lenA > 0 ) then
     call checkCol( prhs(13),   1, 'A' )

     call checkRow( prhs(14), neA, 'iAfun' )
     call checkCol( prhs(14),   1, 'iAfun' )
     call checkRow( prhs(15), neA, 'jAvar' )
     call checkCol( prhs(15),   1, 'jAvar' )

     allocate( iAfun(lenA), jAvar(lenA), A(lenA) )
     allocate( riA(lenA), rjA(lenA) )

     call mxCopyPtrToReal8( mxGetPr(prhs(13)),   A, neA )
     call mxCopyPtrToReal8( mxGetPr(prhs(14)), riA, neA )
     call mxCopyPtrToReal8( mxGetPr(prhs(15)), rjA, neA )

     iAfun = riA(1:neA)
     jAvar = rjA(1:neA)
     deallocate( riA, rjA )

  else
     neA  = 0
     lenA = 0
     allocate( iAfun(1), jAvar(1), A(1) )
  end if


  neG  = mxGetM ( prhs(16) )
  lenG = neG
  if ( neG > 0 ) then
     call checkCol( prhs(16), 1, 'iGfun' )

     call checkRow( prhs(17), neG, 'jGvar' )
     call checkCol( prhs(17),   1, 'jGvar' )

     allocate(  iGfun(lenG),  jGvar(lenG) )
     allocate( riGfun(lenG), rjGvar(lenG) )

     call mxCopyPtrToReal8( mxGetPr(prhs(16)), riGfun, neG )
     call mxCopyPtrToReal8( mxGetPr(prhs(17)), rjGvar, neG )

     iGfun = riGfun(1:neG)
     jGvar = rjGvar(1:neG)

  else
     neG  = 0
     lenG = 0
     allocate( iGfun(1), jGvar(1) )
  end if


  !---------------------------------------------------------------------
  ! Set userfun function
  !---------------------------------------------------------------------
  if ( nrhs == 18 ) then
     info = mxIsClass( prhs(18), 'function_handle' )
     if ( info /= 1 ) call mexErrMsgTxt( 'Wrong input type for userfg' )
     fgHandle  = mxDuplicateArray( prhs(18) )
     objHandle = 0
     conHandle = 0

  else
     fgHandle  = 0

     info = mxIsClass( prhs(18), 'function_handle' )
     if ( info /= 1 ) call mexErrMsgTxt( 'Wrong input type for funobj' )
     objHandle = mxDuplicateArray( prhs(18) )

     info = mxIsClass( prhs(19), 'function_handle' )
     if ( info /= 1 ) call mexErrMsgTxt( 'Wrong input type for nonlcon' )
     conHandle = mxDuplicateArray( prhs(19) )
  end if

  !---------------------------------------------------------------------
  ! Allocateother space for SNOPT
  !---------------------------------------------------------------------
  nfname   = 1
  nxname   = 1

  allocate( xname(nxname), Fname(nFname) )


  !---------------------------------------------------------------------
  ! Set workspace
  !---------------------------------------------------------------------
100 if ( .not. memCall ) then
     call snMemA &
          ( INFO, nF, n, nxname, nfname, neA, neG, &
            mincw, miniw, minrw, &
            cw, lencw, iw, leniw, rw, lenrw )
     memCall = .true.

     if ( leniw .le. miniw ) then
        ! Not enough integer space
        leniw = miniw
        allocate( iw0(leniw) )
        iw0(1:500) = iw(1:500)

        call move_alloc( from=iw0, to=iw )

     end if

     if ( lenrw .le. minrw ) then
        ! Not enough real space
        lenrw = minrw
        allocate( rw0(lenrw) )
        rw0(1:500) = rw(1:500)

        call move_alloc( from=rw0, to=rw )
     end if

     if ( lencw .le. mincw ) then
        ! Not enough character space
        lencw = mincw
        allocate( cw0(lencw) )
        cw0(1:500) = cw(1:500)

        call move_alloc( from=cw0, to=cw )
     end if

     call snSeti &
          ( 'Total character workspace', lencw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw )
     call snSeti &
          ( 'Total integer   workspace', leniw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw )
     call snSeti &
          ( 'Total real      workspace', lenrw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw )
  end if

  !---------------------------------------------------------------------
  ! Solve the problem
  !---------------------------------------------------------------------
  Start    = 0  ! cold start
  probName = 'matlabMx'

  call snoptA &
       ( Start, nF, n, nxname, nFname, ObjAdd, ObjRow, &
         probName, matlabFG, &
         iAfun, jAvar, lenA, neA, A, &
         iGfun, jGvar, lenG, neG, &
         xlow, xupp, xname, Flow, Fupp, Fname, &
         x, xstate, xmul, F, Fstate, Fmul, &
         INFO, mincw, miniw, minrw, nS, nInf, sInf, &
         cw, lencw, iw, leniw, rw, lenrw, &
         cw, lencw, iw, leniw, rw, lenrw )

  if ( INFO == 82 .or. INFO == 83 .or. INFO == 84 ) then
     memCall = .false.
     go to 100
  end if


  !---------------------------------------------------------------------
  ! Set output
  !---------------------------------------------------------------------
  if ( nlhs > 0 ) then
     plhs(1) = mxCreateDoubleMatrix ( n, 1, 0 )
     call mxCopyReal8ToPtr( x, mxGetPr(plhs(1)), n )
  end if


  ! Constraints
  if ( nlhs > 1 ) then
     plhs(2) = mxCreateDoubleMatrix ( nF, 1, 0 )
     call mxCopyReal8ToPtr( F, mxGetPr(plhs(2)), nF )
  end if

  ! Exit flag
  rinfo = info
  if ( nlhs > 2 ) plhs(3) = mxCreateDoubleScalar ( rinfo )

  ! Multipliers for bounds
  if ( nlhs > 3 ) then
     plhs(4) = mxCreateDoubleMatrix ( n, 1, 0 )
     call mxCopyReal8ToPtr( xmul, mxGetPr(plhs(4)), n )
  end if


  ! Multipliers for linear inequalities
  if ( nlhs > 4 ) then
     plhs(5) = mxCreateDoubleMatrix ( nF, 1, 0 )
     call mxCopyReal8ToPtr( Fmul, mxGetPr(plhs(5)), nF )
  end if


  ! Deallocate memory
  if ( fgHandle  /= 0 ) call mxDestroyArray( fgHandle )
  if ( objHandle /= 0 ) call mxDestroyArray( objHandle )
  if ( conHandle /= 0 ) call mxDestroyArray( conHandle )

  fgHandle  = 0
  objHandle = 0
  conHandle = 0

  if ( allocated(rtmp) )   deallocate( rtmp )
  if ( allocated(x) )      deallocate( x )
  if ( allocated(xmul) )   deallocate( xmul )
  if ( allocated(xlow) )   deallocate( xlow )
  if ( allocated(xupp) )   deallocate( xupp )
  if ( allocated(xname) )  deallocate( xname )
  if ( allocated(xstate) ) deallocate( xstate )

  if ( allocated(F) )      deallocate( F )
  if ( allocated(Fmul) )   deallocate( Fmul )
  if ( allocated(Flow) )   deallocate( Flow )
  if ( allocated(Fupp) )   deallocate( Fupp )
  if ( allocated(Fname) )  deallocate( Fname )
  if ( allocated(Fstate) ) deallocate( Fstate )

  if ( allocated(iAfun) )  deallocate( iAfun )
  if ( allocated(jAvar) )  deallocate( jAvar )
  if ( allocated(A) )      deallocate( A )
  if ( allocated(iGfun) )  deallocate( iGfun )
  if ( allocated(jGvar) )  deallocate( jGvar )

  if ( allocated(riGfun) ) deallocate( riGfun )
  if ( allocated(rjGvar) ) deallocate( rjGvar )

end subroutine snmxSolve

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snmxFindJac ( nlhs, plhs, nrhs, prhs )
  use mxsnWork

  implicit none
  integer*4  :: nlhs, nrhs
  mwPointer  :: prhs(*), plhs(*)
  !---------------------------------------------------------------------
  ! Find the structure of the Jacobian matrix.
  ! The matlab call is
  !  [A,iAfun,jAvar,iGfun,jGvar] = snJac( usrfun, x0, xlow, xupp, nF )
  !
  ! 25 Sep 2013: First Fortran version.
  !---------------------------------------------------------------------
  mwPointer        :: mxGetM, mxGetPr, mxCreateDoubleMatrix, &
                      mxDuplicateArray
  integer          :: mxIsClass
  double precision :: mxGetScalar

  integer          :: iExit, Errors, info, lenA, lenG, n, neA, nF, &
                      nfname, nxname, mincw, miniw, minrw
  double precision :: rtmp

  integer,          allocatable :: iAfun(:), jAvar(:), iGfun(:), jGvar(:)
  double precision, allocatable :: riAfun(:), rjAvar(:)
  double precision, allocatable :: A(:), x(:), xlow(:), xupp(:)

  external         :: snMemA, snJac, matlabFG

  if ( nlhs /= 5 ) call mexErrMsgTxt( 'Wrong number of output variables' )
  if ( nrhs /= 6 ) call mexErrMsgTxt( 'Wrong number of input arguments' )


  ! Get userfg function
  info = mxIsClass( prhs(2), 'function_handle')
  if ( info /= 1 ) call mexErrMsgTxt( 'Wrong input type for userfg' )
  fgHandle = mxDuplicateArray( prhs(2) )


  ! Set the number of variables
  n = mxGetM ( prhs(3) )
  allocate( x(n), xlow(n), xupp(n) )


  ! Get initial point x0
  call mxCopyPtrToReal8( mxGetPr(prhs(3)), x, n )


  ! Get lower and upper bounds
  call mxCopyPtrToReal8( mxGetPr(prhs(4)), xlow, n )
  call mxCopyPtrToReal8( mxGetPr(prhs(5)), xupp, n )


  ! Get number of constraints
  rtmp = mxGetScalar(prhs(6))
  nF   = rtmp


  ! Allocatespace
  lenA = n*nF
  neA  = lenA
  lenG = lenA
  neG  = lenA

  allocate( iAfun(lenA), jAvar(lenA), A(lenA) )
  allocate( iGfun(lenG), jGvar(lenG) )

  nfname = 1
  nxname = 1

  !---------------------------------------------------------------------
  ! Set workspace
  !---------------------------------------------------------------------
100 if ( .not. memCall ) then
     call snMemA &
          ( INFO, nF, n, nxname, nfname, neA, neG, &
            mincw, miniw, minrw, &
            cw, lencw, iw, leniw, rw, lenrw )

     memCall = .true.

     if ( leniw .le. miniw ) then
        ! Not enough integer space
        leniw = miniw
        allocate( iw0(leniw) )
        iw0(1:500) = iw(1:500)

        call move_alloc( from=iw0, to=iw )
     end if

     if ( lenrw .le. minrw ) then
        ! Not enough real space
        lenrw = minrw
        allocate( rw0(lenrw) )
        rw0(1:500) = rw(1:500)

        call move_alloc( from=rw0, to=rw )
     end if

     if ( lencw .le. mincw ) then
        ! Not enough character space
        lencw = mincw
        allocate( cw0(lencw) )
        cw0(1:500) = cw(1:500)

        call move_alloc( from=cw0, to=cw )
     end if

     call snSeti &
          ( 'Total character workspace', lencw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw )
     call snSeti &
          ( 'Total integer   workspace', leniw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw )
     call snSeti &
          ( 'Total real      workspace', lenrw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw )
  end if


  !---------------------------------------------------------------------
  ! Compute the structure of the Jacobian
  !---------------------------------------------------------------------
  call snJac( iExit, nF, n, matlabFG, &
              iAfun, jAvar, lenA, neA, A, &
              iGfun, jGvar, lenG, neG, &
              x, xlow, xupp, mincw, miniw, minrw, &
              cw, lencw, iw, leniw, rw, lenrw, &
              cw, lencw, iw, leniw, rw, lenrw )

  if ( INFO == 82 .or. INFO == 83 .or. INFO == 84 ) then
     memCall = .false.
     go to 100
  end if


  allocate( riAfun(neA), rjAvar(neA), riGfun(neG), rjGvar(neG) )
  riAfun = iAfun(1:neA)
  rjAvar = jAvar(1:neA)
  riGfun = iGfun(1:neG)
  rjGvar = jGvar(1:neG)


  ! Set output [A, iAfun, jAvar, iGfun, jGvar ]
  plhs(1) = mxCreateDoubleMatrix( neA, 1, 0 )
  plhs(2) = mxCreateDoubleMatrix( neA, 1, 0 )
  plhs(3) = mxCreateDoubleMatrix( neA, 1, 0 )

  plhs(4) = mxCreateDoubleMatrix( neG, 1, 0 )
  plhs(5) = mxCreateDoubleMatrix( neG, 1, 0 )

  call mxCopyReal8ToPtr( A(1:neA), mxGetPr(plhs(1)), neA )
  call mxCopyReal8ToPtr( riAfun,   mxGetPr(plhs(2)), neA )
  call mxCopyReal8ToPtr( rjAvar,   mxGetPr(plhs(3)), neA )

  call mxCopyReal8ToPtr( riGfun,   mxGetPr(plhs(4)), neG )
  call mxCopyReal8ToPtr( rjGvar,   mxGetPr(plhs(5)), neG )


  ! Destroy arrays
  deallocate( x, xlow, xupp )
  deallocate(  iAfun,  jAvar,  iGfun,  jGvar, A )
  deallocate( riAfun, rjAvar, riGfun, rjGvar )

  call mxDestroyArray( fgHandle )
  fgHandle = 0

end subroutine snmxFindJac

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snmxOptions( iOpt, nlhs, plhs, nrhs, prhs )
  use mxsnWork

  implicit none
  integer*4  :: iOpt, nlhs, nrhs
  mwPointer  :: prhs(*), plhs(*)
  !---------------------------------------------------------------------
  ! Set/get options.
  !---------------------------------------------------------------------
  ! Matlab
  mwPointer        :: mxGetN, mxGetPr, mxCreateDoubleScalar, &
                      mxCreateString
  double precision :: mxGetScalar

  character        :: buffer*50, cvalue*8
  integer          :: Errors, ivalue, strlen
  double precision :: rvalue

  integer          :: snGet
  external         :: snSet, snSetI, snSetR, &
                      snGet, snGetC, snGetI, snGetR


  if ( iOpt == snSetIX .or. iOpt == snSetRX ) then
     if ( nrhs /= 3 ) call mexErrMsgTxt( 'Wrong number of input arguments' )
  else
     if ( nrhs /= 2 ) call mexErrMsgTxt( 'Wrong number of input arguments' )
  end if


  ! Get string
  strlen = mxGetN(prhs(2))
  if ( strlen > 50 ) call mexErrMsgTxt( 'Option string is too long' )

  if ( strlen > 0 ) then
     call mxGetString( prhs(2), buffer, strlen )
  else
     call mexErrMsgTxt( 'Empty option string' )
  end if


  if      ( iOpt == snSetXX ) then
     call snSet( buffer, iPrint, iSumm, Errors, &
                 cw, lencw, iw, leniw, rw, lenrw )

  else if ( iOpt == snSetIX ) then

     rvalue = mxGetScalar(prhs(3))
     ivalue = rvalue

     call snSetI( buffer, ivalue, iPrint, iSumm, Errors, &
                  cw, lencw, iw, leniw, rw, lenrw )

  else if ( iOpt == snSetRX ) then

     rvalue = mxGetScalar(prhs(3))

     call snSetR( buffer, rvalue, iPrint, iSumm, Errors, &
                  cw, lencw, iw, leniw, rw, lenrw )

  else if ( iOpt == snGetXX ) then

     ivalue  = snGet( buffer, Errors, cw, lencw, iw, leniw, rw, lenrw )

     rvalue  = ivalue
     plhs(1) = mxCreateDoubleScalar( rvalue )

  else if ( iOpt == snGetCX ) then

     call snGetC( buffer, cvalue, Errors, &
                  cw, lencw, iw, leniw, rw, lenrw )

     plhs(1) = mxCreateString(cvalue)

  else if ( iOpt == snGetIX ) then

     call snGetI( buffer, ivalue, Errors, &
                  cw, lencw, iw, leniw, rw, lenrw )

     rvalue = ivalue
     plhs(1) = mxCreateDoubleScalar( rvalue )

  else if ( iOpt == snGetRX ) then

     call snGetR( buffer, rvalue, Errors, &
                  cw, lencw, iw, leniw, rw, lenrw )

     plhs(1) = mxCreateDoubleScalar( rvalue )

  end if

end subroutine snmxOptions

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snmxSpecs( nlhs, plhs, nrhs, prhs )
  use mxsnWork

  implicit none
  integer*4  :: nlhs, nrhs
  mwPointer  :: prhs(*), plhs(*)
  !---------------------------------------------------------------------
  ! Read specs file.
  !---------------------------------------------------------------------
  ! Matlab
  mwPointer        :: mxCreateDoubleScalar, mxGetN

  character        :: filename*120
  integer          :: info, strlen
  double precision :: rvalue

  external         :: snSpec


  if ( nrhs /= 2 ) call mexErrMsgTxt( 'Wrong number of input arguments' )
  if ( nlhs /= 1 ) call mexErrMsgTxt( 'Wrong number of output arguments' )


  strlen = mxGetN(prhs(2))
  if ( strlen > 120 ) call mexErrMsgTxt( 'Specs filename is too long' )

  if ( strlen > 0 ) then
     call mxGetString( prhs(2), filename, strlen )
  else
     call mexErrMsgTxt( 'Empty spc filename' )
  end if

  open( iSpecs, file=filename, status='unknown' )
  call snSpec( iSpecs, info, cw, lencw, iw, leniw, rw, lenrw )
  rewind ( iSpecs )
  close( iSpecs )

  ! snSpec will return info == 101 or 107 if successful
  ! The matlab version returns 0 if successful
  if ( info == 101 .or. info == 107 ) then
     rvalue = 0
  else
     rvalue = 1
  end if

  plhs(1) = mxCreateDoubleScalar( rvalue )

end subroutine snmxSpecs

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine matlabFG( Status, n, x, needF, nF, F, needG, lenG, G, &
                     cu, lencu, iu, leniu, ru, lenru )
  use mxsnWork
  implicit none

  integer          :: Status, n, nF, needF, needG, lenG, &
                      lencu, leniu, lenru, iu(leniu)
  double precision :: F(nF), G(lenG), x(n), ru(lenru)
  character*8      :: cu(lencu)

  !---------------------------------------------------------------------
  ! Matlab callback to evaluate objective function and gradient at x.
  !---------------------------------------------------------------------
  integer*4        :: nlhs, nrhs, nlhs1, nrhs1
  double precision :: tmp
  double precision, allocatable :: G1(:)

  mwIndex          :: j

  mwPointer        :: prhs(6), plhs(2), prhs1(3), plhs1(1)
  mwPointer        :: mxCreateDoubleMatrix, mxCreateDoubleScalar, &
                      mxDuplicateArray, mxGetN, mxGetPr
  integer          :: mxIsNaN

  if ( needF == 0 .and. needG == 0 ) return

  nlhs    = 2
  nrhs    = 4

  prhs(1) = mxCreateDoubleMatrix( n, 1, 0 )
  call mxCopyReal8ToPtr( x, mxGetPr(prhs(1)), n )

  tmp = nF
  prhs(3) = mxCreateDoubleScalar( tmp )

  tmp = n
  prhs(4) = mxCreateDoubleScalar( tmp )


  if ( fgHandle /= 0 ) then
     prhs(2) = mxDuplicateArray( fgHandle )

  else
     prhs(2) = mxCreateDoubleScalar( tmp )

     if ( objHandle /= 0 ) then
        nrhs = 5
        prhs(5) = mxDuplicateArray( objHandle )

        if ( conHandle /= 0 ) then
           nrhs = 6
           prhs(6) = mxDuplicateArray( conHandle )
        end if
     end if
  end if

  ! Call Matlab: [F,G] = snwrapper(x,userfun,nF,n)
  call mexCallMatlab( nlhs, plhs, nrhs, prhs(1:nrhs), 'snwrapper' )

  ! Objective and constraints
  if ( needF > 0 ) then
     call checkRow( plhs(1), nF, 'F' )
     call checkCol( plhs(1),  1, 'F' )

     call mxCopyPtrToReal8( mxGetPr(plhs(1)), F, nF )
  end if


  ! Gradient and Jacobian
  if ( needG > 0 ) then

     if ( mxGetN(plhs(2)) == 0 ) then
        ! No derivatives are defined

     else
        if ( mxGetN(plhs(2)) > 1 ) then
           ! Dense G given
           call checkRow( plhs(2), nF, 'G' )
           call checkCol( plhs(2),  n, 'G' )

           ! Convert to vector format for SNOPTA
           nlhs1    = 1
           nrhs1    = 3

           prhs1(1) = mxCreateDoubleMatrix( neG, 1, 0 )
           prhs1(2) = mxCreateDoubleMatrix( neG, 1, 0 )
           prhs1(3) = mxDuplicateArray(plhs(2))

           call mxCopyReal8ToPtr( riGfun(1:neG), mxGetPr(prhs1(1)), neG )
           call mxCopyReal8ToPtr( rjGvar(1:neG), mxGetPr(prhs1(2)), neG )

           ! Call Matlab: [G] = snfindG(iGfun,jGvar,Gdense)
           call mexCallMatlab( nlhs1, plhs1, nrhs1, prhs1, 'snfindG' )

           ! Check dimensions of G
           call checkRow( plhs1(1), neG, 'G' )
           call checkCol( plhs1(1),   1, 'G' )


           ! Copy non-NaN entries of G.
           allocate( G1(neG) )
           call mxCopyPtrToReal8( mxGetPr(plhs1(1)), G1, neG )

           do j = 1, neG
              if ( mxIsNaN(G1(j)) == 0 ) then
                 G(j) = G1(j)
              end if
           end do

           deallocate( G1 )


           ! Destroy array
           call mxDestroyArray( prhs1(1) )
           call mxDestroyArray( prhs1(2) )
           call mxDestroyArray( prhs1(3) )

        else
           ! Sparse G given
           ! Copy non-NaN entries of G.
           allocate( G1(neG) )
           call mxCopyPtrToReal8( mxGetPr(plhs(2)), G1, neG )

           do j = 1, neG
              if ( mxIsNaN(G1(j)) == 0 ) then
                 G(j) = G1(j)
              end if
           end do

           deallocate( G1 )
        end if
     end if
  end if


  ! Destroy arrays
  call mxDestroyArray( prhs(1) )
  call mxDestroyArray( prhs(2) )
  call mxDestroyArray( prhs(3) )
  call mxDestroyArray( prhs(4) )
  if ( nrhs >= 5 ) call mxDestroyArray( prhs(5) )
  if ( nrhs >= 6 ) call mxDestroyArray( prhs(6) )

end subroutine matlabFG

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
