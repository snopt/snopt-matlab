!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File snoptmex.F90
! Mex function for SNOPT7.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "fintrf.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
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
  ! 09 Nov 2015: Added states, iteration counts as output
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
  if (nrhs < 1) call mexErrMsgIdAndTxt('SNOPT:InputArgs','Need an option input argument')
  rOpt = mxGetScalar(prhs(1))
  iOpt = rOpt

  ! Register exit function.
  call mexAtExit( resetSNOPT )


  ! Deal with on/off screen, print, and summary files first.
  if (iOpt == snOpenP) then

     if (nrhs /= 2) call mexErrMsgIdAndTxt('SNOPT:InputArgs','Wrong number of input arguments')

     info = mxIsChar(prhs(2))
     if (info /= 1) call mexErrMsgIdAndTxt('SNOPT:FileArgs','Need a filename string')

     strlen = mxGetN(prhs(2))
     if (strlen > 80) call mexErrMsgIdAndTxt('SNOPT:FileArg','Print filename is too long')

     if (strlen > 0) then
        call mxGetString(prhs(2), filename, strlen)
     else
        call mexErrMsgIdAndTxt('SNOPT:FileArg','Empty print filename')
     end if

     if (printOpen) close(iPrint)

     open(iPrint, file=filename, status='unknown')
     printOpen= .true.
     return

  else if (iOpt == snOpenS) then

     if (nrhs /= 2) call mexErrMsgIdAndTxt('SNOPT:InputArgs','Wrong number of input arguments')

     info = mxIsChar(prhs(2))
     if (info /= 1) call mexErrMsgIdAndTxt('SNOPT:FileArgs','Need a filename string')

     strlen = mxGetN(prhs(2))
     if (strlen > 80) call mexErrMsgIdAndTxt('SNOPT:FileArgs','Summary filename is too long')

     if (strlen > 0) then
        call mxGetString(prhs(2), filename, strlen)
     else
        call mexErrMsgIdAndTxt('SNOPT:FileArg','Empty print filename')
     end if

     if (summOpen) close(iSumm)

     open(iSumm, file=filename, status='unknown')
     summOpen= .true.
     return

  else if (iOpt == snClosP) then
     if (printOpen) close(iPrint)
     printOpen= .false.
     return

  else if (iOpt == snClosS) then
     if (summOpen) close(iSumm)
     summOpen= .false.
     return

  else if (iOpt == snscrnOn) then
     screenOn = .true.
     return

  else if (iOpt == snscrnOff) then
     screenOn = .false.
     return

  else if (iOpt == snsetwork) then
     rleniw = mxGetScalar(prhs(2))
     rlenrw = mxGetScalar(prhs(3))
     leniw  = rleniw
     lenrw  = rlenrw

     if (leniw < 500 .or. lenrw < 500) &
          call mexErrMsgIdAndTxt('SNOPT:Workspace','Workspace size must be at least 500')
     return
  end if

  ! What calls get to this point:
  !  snSolve, snGet, snSet, snSpecs, snJac
  if (firstCall) then
     allocate(cw(lencw), iw(leniw), rw(lenrw))

     callType = userCall
     call snInit (iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw)
     callType = systemCall

     memCall   = .false.
     firstCall = .false.
  end if

  ! Do whatever we need to do.
  if      (iOpt == snSolve) then

     callType = userCall
     call snmxSolve (nlhs, plhs, nrhs, prhs)
     callType = systemCall

  else if (iOpt == snFindJac) then

     callType = userCall
     call snmxFindJac (nlhs, plhs, nrhs, prhs)
     callType = systemCall

  else if (iOpt == snSetXX .or. &
           iOpt == snSetIX .or. &
           iOpt == snSetRX .or. &
           iOpt == snGetXX .or. &
           iOpt == snGetCX .or. &
           iOpt == snGetIX .or. &
           iOpt == snGetRX) then

     callType = userCall
     call snmxOptions (iOpt, nlhs, plhs, nrhs, prhs)
     callType = systemCall

  else if (iOpt == snSpecs) then

     callType = userCall
     call snmxSpecs (nlhs, plhs, nrhs, prhs)
     callType = systemCall

  else if (iOpt == snEnd) then

     call resetSNOPT

  end if

  return

end subroutine mexFunction

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snmxSolve (nlhs, plhs, nrhs, prhs)
  use mxsnWork

  implicit none
  integer*4  :: nlhs, nrhs
  mwPointer  :: prhs(*), plhs(*)
  !---------------------------------------------------------------------
  ! Solve the problem
  ! The matlab call is
  !   [x, F, xmul, Fmul, info ] =
  !       snoptmex (solveopt, x, xlow, xupp, xmul, xstate,
  !                            F, Flow, Fupp, Fmul, Fstate, ObjAdd,
  !                  ObjRow, A, iAfun, jAvar, iGfun, jGvar, userfun)
  !---------------------------------------------------------------------
  ! Matlab
  mwPointer        :: mxDuplicateArray, mxGetM, mxGetN, mxGetPr, &
                      mxCreateDoubleMatrix, mxCreateDoubleScalar
  integer          :: mxIsChar, mxIsClass, mxIsEmpty, mxIsNumeric
  double precision :: mxGetScalar

  ! SNOPT
  character*8      :: probName
  integer          :: info, Errors, tmp
  integer          :: Start, ObjRow, n, nF, lenA, lenG, neA, &
                      mincw, miniw, minrw, nInf, nS
  double precision :: rinfo, ObjAdd, sInf
  external         :: snMemA, snKerA, matlabFG, matlabSTOP
  external         :: snLog, snLog2, sqLog

  double precision, parameter   :: infBnd = 1.0d+20

  integer,          parameter   :: nxname = 1, nFname = 1
  character*8                   :: xname(1), Fname(1)

  ! Check number of input and output arguments.
  if (nrhs /= 20 .and. nrhs /= 21) &
       call mexErrMsgIdAndTxt('SNOPT:InputArgs','Wrong number of input arguments')


  !---------------------------------------------------------------------
  ! Compute number of variables, constraints, etc
  ! Check dimension of some of the input
  !---------------------------------------------------------------------
  ! Get number of variables and constraints
  n  = mxGetM(prhs(2))
  nF = mxGetM(prhs(7))

  call allocSNOPT(n, nF)

  ! Get initial x
  if (mxIsEmpty(prhs(2)) > 0 ) then
     x(1:n) = 0.0
  else
     call checkRow(prhs(2), n, 'x0')
     call checkCol(prhs(2), 1, 'x0')
     call mxCopyPtrToReal8(mxGetPr(prhs(2)), x, n)
  end if

  ! Get lower bounds
  if (mxIsEmpty(prhs(3)) > 0) then
     xlow(1:n) = -infBnd
  else
     call checkRow(prhs(3), n, 'xlow')
     call checkCol(prhs(3), 1, 'xlow')
     call mxCopyPtrToReal8(mxGetPr(prhs(3)), xlow, n)
  end if

  ! Get upper bounds
  if (mxIsEmpty(prhs(4)) > 0) then
     xupp(1:n) = infBnd
  else
     call checkRow(prhs(4), n, 'xupp')
     call checkCol(prhs(4), 1, 'xupp')
     call mxCopyPtrToReal8(mxGetPr(prhs(4)), xupp, n)
  end if

  ! Get initial multipliers
  if (mxIsEmpty(prhs(5)) > 0) then
     xmul(1:n) = 0
  else
     call checkRow(prhs(5), n, 'xmul')
     call checkCol(prhs(5), 1, 'xmul')
     call mxCopyPtrToReal8(mxGetPr(prhs(5)), xmul, n)
  end if

  ! Get initial states
  if (mxIsEmpty(prhs(6)) > 0) then
     xstate(1:n) = 0
  else
     call checkRow(prhs(6), n, 'xstate')
     call checkCol(prhs(6), 1, 'xstate')
     call mxCopyPtrToReal8(mxGetPr(prhs(6)), rtmpa, n)
     xstate(1:n) = rtmpa(1:n)
  end if

  ! Get lower bounds
  if (mxIsEmpty(prhs(7)) > 0) then
     Flow(1:nF) = -infBnd
  else
     call checkRow(prhs(7),  nF, 'Fupp')
     call checkCol(prhs(7),   1, 'Flow')
     call mxCopyPtrToReal8(mxGetPr(prhs(7)),  Flow, nF)
  end if

  ! Get upper bounds
  if (mxIsEmpty(prhs(8)) > 0) then
     Fupp(1:nF) = infBnd
  else
     call checkRow(prhs(8),  nF, 'Fupp')
     call checkCol(prhs(8),   1, 'Fupp')
     call mxCopyPtrToReal8(mxGetPr(prhs(8)),  Fupp, nF)
  end if

  ! Get multipliers
  if (mxIsEmpty(prhs(9)) > 0) then
     Fmul(1:nF) = 0.0
  else
     call checkRow(prhs(9),  nF, 'Fmul')
     call checkCol(prhs(9),   1, 'Fmul')
     call mxCopyPtrToReal8(mxGetPr(prhs(9)),  Fmul, nF)
  end if

  ! Get states
  if (mxIsEmpty(prhs(10)) > 0) then
     Fstate(1:nF) = 0
  else
     call checkRow(prhs(10), nF, 'Fstate')
     call checkCol(prhs(10),  1, 'Fstate')
     call mxCopyPtrToReal8(mxGetPr(prhs(10)), rtmpa, nF)
     Fstate(1:nF) = rtmpa(1:nF)
  end if

  ObjAdd   = mxGetScalar(prhs(11))
  rtmpa(1) = mxGetScalar(prhs(12))
  ObjRow   = rtmpa(1)

  !---------------------------------------------------------------------
  ! Get the Jacobian structure (linear and nonlinear)
  !---------------------------------------------------------------------
  if (mxIsEmpty(prhs(13)) > 0) then
     neA  = 0
     lenA = 0
  else
     neA  = mxGetM(prhs(13))
     lenA = neA
  end if

  if (neA > 0) then
     call checkCol(prhs(13),   1, 'A')

     call checkRow(prhs(14), neA, 'iAfun')
     call checkCol(prhs(14),   1, 'iAfun')
     call checkRow(prhs(15), neA, 'jAvar')
     call checkCol(prhs(15),   1, 'jAvar')

     call allocA( lenA )

     call mxCopyPtrToReal8(mxGetPr(prhs(13)),      A, neA)
     call mxCopyPtrToReal8(mxGetPr(prhs(14)), riAfun, neA)
     call mxCopyPtrToReal8(mxGetPr(prhs(15)), rjAvar, neA)

     iAfun(1:neA) = riAfun(1:neA)
     jAvar(1:neA) = rjAvar(1:neA)

     deallocate(riAfun, rjAvar)
  else
     neA  = 0
     lenA = 0

     call allocA( 1 )
  end if


  if (mxIsEmpty(prhs(16)) > 0) then
     neG  = 0
     lenG = 0
  else
     neG  = mxGetM(prhs(16))
     lenG = neG
  end if

  if (neG > 0) then
     call checkCol(prhs(16),   1, 'iGfun')

     call checkRow(prhs(17), neG, 'jGvar')
     call checkCol(prhs(17),   1, 'jGvar')

     call allocG( lenG )

     call mxCopyPtrToReal8(mxGetPr(prhs(16)), riGfun, neG)
     call mxCopyPtrToReal8(mxGetPr(prhs(17)), rjGvar, neG)

     iGfun(1:neG) = riGfun(1:neG)
     jGvar(1:neG) = rjGvar(1:neG)

  else
     neG  = 0
     lenG = 0
     call allocG( 1 )
  end if


  !---------------------------------------------------------------------
  ! Set userfun, problem name, and snSTOP functions
  !---------------------------------------------------------------------
  probName = ''
  if (nrhs == 20) then
     ! snopt.m
     info = mxIsClass(prhs(18), 'function_handle')
     if (info /= 1) call mexErrMsgIdAndTxt('SNOPT:FunArg','Wrong input type for userfg')
     fgHandle  = mxDuplicateArray(prhs(18))
     objHandle = 0
     conHandle = 0

     ! Problem name
     info = mxIsChar(prhs(19))
     if (info /= 1) call mexErrMsgIdAndTxt('SNOPT:InputArg','Wrong input type for problem name')
     tmp = mxGetN(prhs(19))
     call mxGetString(prhs(19), probName, min(8,tmp))

     ! Check for STOP function
     stopHandle = 0
     info = mxIsNumeric(prhs(20))
     if ( info /= 1 ) then
        ! Check if STOP function is actually a function
        info = mxIsClass(prhs(20), 'function_handle')
        if (info /= 1) call mexErrMsgIdAndTxt('SNOPT:FunArg','Wrong input type for snSTOP')
        stopHandle = mxDuplicateARray(prhs(20))
     end if

  else
     ! snsolve.m
     fgHandle  = 0

     info = mxIsClass(prhs(18), 'function_handle')
     if (info /= 1) call mexErrMsgIdAndTxt('SNOPT:FunArg','Wrong input type for funobj')
     objHandle = mxDuplicateArray(prhs(18))

     info = mxIsClass(prhs(19), 'function_handle')
     if (info /= 1) call mexErrMsgIdAndTxt('SNOPT:FunArg','Wrong input type for nonlcon')
     conHandle = mxDuplicateArray(prhs(19))

     ! Problem name
     info = mxIsChar(prhs(20))
     if (info /= 1) call mexErrMsgIdAndTxt('SNOPT:InputArg','Wrong input type for problem name')
     tmp = mxGetN(prhs(20))
     call mxGetString(prhs(20), probName, min(8,tmp))

     ! Check for STOP function
     stopHandle = 0
     info = mxIsNumeric(prhs(21))
     if ( info /= 1 ) then
        ! Check if STOP function is actually a function
        info = mxIsClass(prhs(21), 'function_handle')
        if (info /= 1) call mexErrMsgIdAndTxt('SNOPT:FunArg','Wrong input type for snSTOP')
        stopHandle = mxDuplicateARray(prhs(21))
     end if

  end if


  !---------------------------------------------------------------------
  ! Set workspace
  !---------------------------------------------------------------------
100 if (.not. memCall) then
     call snMemA &
          (INFO, nF, n, nxname, nfname, neA, neG, &
            mincw, miniw, minrw, &
            cw, lencw, iw, leniw, rw, lenrw)
     memCall = .true.

     if (leniw .le. miniw) then
        ! Not enough integer space
        leniw = miniw
        allocate(iw0(leniw))
        iw0(1:500) = iw(1:500)

        call move_alloc(from=iw0, to=iw)

     end if

     if (lenrw .le. minrw) then
        ! Not enough real space
        lenrw = minrw
        allocate(rw0(lenrw))
        rw0(1:500) = rw(1:500)

        call move_alloc(from=rw0, to=rw)
     end if

     if (lencw .le. mincw) then
        ! Not enough character space
        lencw = mincw
        allocate(cw0(lencw))
        cw0(1:500) = cw(1:500)

        call move_alloc(from=cw0, to=cw)
     end if

     call snSeti &
          ('Total character workspace', lencw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw)
     call snSeti &
          ('Total integer   workspace', leniw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw)
     call snSeti &
          ('Total real      workspace', lenrw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw)
  end if

  !---------------------------------------------------------------------
  ! Solve the problem
  !---------------------------------------------------------------------
  Start    = 0  ! cold start

  call snKerA                                        &
       ( Start, nF, n, nxname, nFname,               &
         objAdd, objRow, probName,                   &
         matlabFG, snLog, snLog2, sqLog, matlabSTOP, &
         iAfun, jAvar, lenA, neA, A,                 &
         iGfun, jGvar, lenG, neG,                    &
         xlow, xupp, xname, Flow, Fupp, Fname,       &
         x, xstate, xmul, F, Fstate, Fmul,           &
         INFO, mincw, miniw, minrw,                  &
         nS, nInf, sInf,                             &
         cw, lencw, iw, leniw, rw, lenrw,            &
         cw, lencw, iw, leniw, rw, lenrw )

  if (INFO == 82 .or. INFO == 83 .or. INFO == 84) then
     memCall = .false.
     go to 100
  end if


  !---------------------------------------------------------------------
  ! Set output
  !---------------------------------------------------------------------
  if (nlhs > 0) then
     plhs(1) = mxCreateDoubleMatrix (n, 1, 0)
     call mxCopyReal8ToPtr(x, mxGetPr(plhs(1)), n)
  end if


  ! Constraints
  if (nlhs > 1) then
     plhs(2) = mxCreateDoubleMatrix (nF, 1, 0)
     call mxCopyReal8ToPtr(F, mxGetPr(plhs(2)), nF)
  end if

  ! Exit flag
  rinfo = info
  if (nlhs > 2) plhs(3) = mxCreateDoubleScalar (rinfo)

  ! Multipliers for bounds
  if (nlhs > 3) then
     plhs(4) = mxCreateDoubleMatrix (n, 1, 0)
     call mxCopyReal8ToPtr(xmul, mxGetPr(plhs(4)), n)
  end if


  ! Multipliers for linear inequalities
  if (nlhs > 4) then
     plhs(5) = mxCreateDoubleMatrix (nF, 1, 0)
     call mxCopyReal8ToPtr(Fmul, mxGetPr(plhs(5)), nF)
  end if

  ! State variables
  if ( nlhs > 5 ) then
     plhs(6) = mxCreateDoubleMatrix(n,1,0)
     call mxCopyReal8ToPtr(real(xstate,8), mxGetPr(plhs(6)),n)
  end if

  ! State variables
  if ( nlhs > 6 ) then
     plhs(7) = mxCreateDoubleMatrix(nF,1,0)
     call mxCopyReal8ToPtr(real(Fstate,8), mxGetPr(plhs(7)),nF)
  end if

  ! Number of total iterations, major itns
  rinfo = iw(421)
  if (nlhs > 7) plhs(8) = mxCreateDoubleScalar(rinfo)

  rinfo = iw(422)
  if (nlhs > 8) plhs(9) = mxCreateDoubleScalar(rinfo)


  ! Deallocate memory
  if (fgHandle    /= 0) call mxDestroyArray(fgHandle)
  if (objHandle   /= 0) call mxDestroyArray(objHandle)
  if (conHandle   /= 0) call mxDestroyArray(conHandle)
  if (stopHandle  /= 0) call mxDestroyArray(stopHandle)

  fgHandle   = 0
  objHandle  = 0
  conHandle  = 0
  stopHandle = 0

  call deallocSNOPT
  call deallocA
  call deallocG

end subroutine snmxSolve

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snmxFindJac (nlhs, plhs, nrhs, prhs)
  use mxsnWork

  implicit none
  integer*4  :: nlhs, nrhs
  mwPointer  :: prhs(*), plhs(*)
  !---------------------------------------------------------------------
  ! Find the structure of the Jacobian matrix.
  ! The matlab call is
  !  [A,iAfun,jAvar,iGfun,jGvar] = snJac(usrfun, x0, xlow, xupp, nF)
  !
  ! 25 Sep 2013: First Fortran version.
  !---------------------------------------------------------------------
  mwPointer        :: mxGetM, mxGetPr, mxCreateDoubleMatrix, &
                      mxDuplicateArray
  integer          :: mxIsClass
  double precision :: mxGetScalar

  integer          :: iExit, Errors, info, lenA, lenG, n, neA, nF, &
                       mincw, miniw, minrw
  double precision :: rtmp

  external         :: snMemA, snJac, matlabFG

  integer, parameter :: nxname = 1, nFname = 1


  if (nlhs /= 5) call mexErrMsgIdAndTxt('SNOPT:InputArg','Wrong number of output variables')
  if (nrhs /= 6) call mexErrMsgIdAndTxt('SNOPT:InputArg','Wrong number of input arguments')


  ! Get userfg function
  info = mxIsClass(prhs(2), 'function_handle')
  if (info /= 1) call mexErrMsgIdAndTxt('SNOPT:FunArg','Wrong input type for userfg')
  fgHandle = mxDuplicateArray(prhs(2))


  ! Set the number of variables
  n = mxGetM (prhs(3))

  ! Get number of constraints
  rtmp = mxGetScalar(prhs(6))
  nF   = rtmp

  if (nF == 0 .or. n == 0) call mexErrMsgIdAndTxt('SNOPT:Input','Empty inputs to snJac')

  ! Allocate space
  lenA = n*nF
  neA  = lenA
  lenG = lenA
  neG  = lenA

  call allocJac( n, lenA, lenG )


  ! Get initial point x0
  call mxCopyPtrToReal8(mxGetPr(prhs(3)), x, n)

  ! Get lower and upper bounds
  call mxCopyPtrToReal8(mxGetPr(prhs(4)), xlow, n)
  call mxCopyPtrToReal8(mxGetPr(prhs(5)), xupp, n)


  !---------------------------------------------------------------------
  ! Set workspace
  !---------------------------------------------------------------------
100 if (.not. memCall) then
     call snMemA &
          (INFO, nF, n, nxname, nfname, neA, neG, &
            mincw, miniw, minrw, &
            cw, lencw, iw, leniw, rw, lenrw)

     memCall = .true.

     if (leniw .le. miniw) then
        ! Not enough integer space
        leniw = miniw
        allocate(iw0(leniw))
        iw0(1:500) = iw(1:500)

        call move_alloc(from=iw0, to=iw)
     end if

     if (lenrw .le. minrw) then
        ! Not enough real space
        lenrw = minrw
        allocate(rw0(lenrw))
        rw0(1:500) = rw(1:500)

        call move_alloc(from=rw0, to=rw)
     end if

     if (lencw .le. mincw) then
        ! Not enough character space
        lencw = mincw
        allocate(cw0(lencw))
        cw0(1:500) = cw(1:500)

        call move_alloc(from=cw0, to=cw)
     end if

     call snSeti &
          ('Total character workspace', lencw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw)
     call snSeti &
          ('Total integer   workspace', leniw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw)
     call snSeti &
          ('Total real      workspace', lenrw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw)
  end if


  !---------------------------------------------------------------------
  ! Compute the structure of the Jacobian
  !---------------------------------------------------------------------
  call snJac(iExit, nF, n, matlabFG, &
              iAfun, jAvar, lenA, neA, A, &
              iGfun, jGvar, lenG, neG, &
              x, xlow, xupp, mincw, miniw, minrw, &
              cw, lencw, iw, leniw, rw, lenrw, &
              cw, lencw, iw, leniw, rw, lenrw)

  if (INFO == 82 .or. INFO == 83 .or. INFO == 84) then
     memCall = .false.
     go to 100
  end if

  riAfun(1:neA) = iAfun(1:neA)
  rjAvar(1:neA) = jAvar(1:neA)
  riGfun(1:neG) = iGfun(1:neG)
  rjGvar(1:neG) = jGvar(1:neG)

  ! Set output [A, iAfun, jAvar, iGfun, jGvar ]
  plhs(1) = mxCreateDoubleMatrix(neA, 1, 0)
  plhs(2) = mxCreateDoubleMatrix(neA, 1, 0)
  plhs(3) = mxCreateDoubleMatrix(neA, 1, 0)

  plhs(4) = mxCreateDoubleMatrix(neG, 1, 0)
  plhs(5) = mxCreateDoubleMatrix(neG, 1, 0)

  call mxCopyReal8ToPtr(A(1:neA), mxGetPr(plhs(1)), neA)
  call mxCopyReal8ToPtr(riAfun,   mxGetPr(plhs(2)), neA)
  call mxCopyReal8ToPtr(rjAvar,   mxGetPr(plhs(3)), neA)

  call mxCopyReal8ToPtr(riGfun,   mxGetPr(plhs(4)), neG)
  call mxCopyReal8ToPtr(rjGvar,   mxGetPr(plhs(5)), neG)

  ! Destroy arrays
  call deallocJac

  call mxDestroyArray(fgHandle)
  fgHandle = 0

end subroutine snmxFindJac

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snmxOptions(iOpt, nlhs, plhs, nrhs, prhs)
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


  if (iOpt == snSetIX .or. iOpt == snSetRX) then
     if (nrhs /= 3) call mexErrMsgIdAndTxt('SNOPT:InputArg','Wrong number of input arguments')
  else
     if (nrhs /= 2) call mexErrMsgIdAndTxt('SNOPT:InputArg','Wrong number of input arguments')
  end if


  ! Get string
  strlen = mxGetN(prhs(2))
  if (strlen > 50) call mexErrMsgIdAndTxt('SNOPT:InputArg','Option string is too long')

  if (strlen > 0) then
     call mxGetString(prhs(2), buffer, strlen)
  else
     call mexErrMsgIdAndTxt('SNOPT:InputArg','Empty option string')
  end if


  if      (iOpt == snSetXX) then
     call snSet(buffer, iPrint, iSumm, Errors, &
                 cw, lencw, iw, leniw, rw, lenrw)

  else if (iOpt == snSetIX) then

     rvalue = mxGetScalar(prhs(3))
     ivalue = rvalue

     call snSetI(buffer, ivalue, iPrint, iSumm, Errors, &
                  cw, lencw, iw, leniw, rw, lenrw)

  else if (iOpt == snSetRX) then

     rvalue = mxGetScalar(prhs(3))

     call snSetR(buffer, rvalue, iPrint, iSumm, Errors, &
                  cw, lencw, iw, leniw, rw, lenrw)

  else if (iOpt == snGetXX) then

     ivalue  = snGet(buffer, Errors, cw, lencw, iw, leniw, rw, lenrw)

     rvalue  = ivalue
     plhs(1) = mxCreateDoubleScalar(rvalue)

  else if (iOpt == snGetCX) then

     call snGetC(buffer, cvalue, Errors, &
                  cw, lencw, iw, leniw, rw, lenrw)

     plhs(1) = mxCreateString(cvalue)

  else if (iOpt == snGetIX) then

     call snGetI(buffer, ivalue, Errors, &
                  cw, lencw, iw, leniw, rw, lenrw)

     rvalue = ivalue
     plhs(1) = mxCreateDoubleScalar(rvalue)

  else if (iOpt == snGetRX) then

     call snGetR(buffer, rvalue, Errors, &
                  cw, lencw, iw, leniw, rw, lenrw)

     plhs(1) = mxCreateDoubleScalar(rvalue)

  end if

end subroutine snmxOptions

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snmxSpecs(nlhs, plhs, nrhs, prhs)
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


  if (nrhs /= 2) call mexErrMsgIdAndTxt('SNOPT:InputArg','Wrong number of input arguments')
  if (nlhs /= 1) call mexErrMsgIdAndTxt('SNOPT:InputArg','Wrong number of output arguments')


  strlen = mxGetN(prhs(2))
  if (strlen > 120) call mexErrMsgIdAndTxt('SNOPT:InputArg','Specs filename is too long')

  if (strlen > 0) then
     call mxGetString(prhs(2), filename, strlen)
  else
     call mexErrMsgIdAndTxt('SNOPT:InputArg','Empty spc filename')
  end if

  open(iSpecs, file=filename, status='unknown')
  call snSpec(iSpecs, info, cw, lencw, iw, leniw, rw, lenrw)
  rewind (iSpecs)
  close(iSpecs)

  ! snSpec will return info == 101 or 107 if successful
  ! The matlab version returns 0 if successful
  if (info == 101 .or. info == 107) then
     rvalue = 0
  else
     rvalue = 1
  end if

  plhs(1) = mxCreateDoubleScalar(rvalue)

end subroutine snmxSpecs

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine matlabFG(Status, n, x, needF, nF, F, needG, lenG, G, &
                     cu, lencu, iu, leniu, ru, lenru)

  use mxsnWork, only : checkRow, checkCol, &
                       fgHandle, objHandle, conHandle, &
                       neG, riGfun, rjGvar, G1

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
  character(20) :: str
  mwIndex          :: j

  mwPointer        :: prhs(8), plhs(2), prhs1(3), plhs1(1)
  mwPointer        :: mxCreateDoubleMatrix, mxCreateDoubleScalar, &
                      mxDuplicateArray, mxGetN, mxGetPr
  integer          :: mxIsNaN

  if (needF == 0 .and. needG == 0) return

  nlhs    = 2
  nrhs    = 6

  prhs(1) = mxCreateDoubleMatrix(n, 1, 0)
  call mxCopyReal8ToPtr(x, mxGetPr(prhs(1)), n)

  tmp = nF
  prhs(3) = mxCreateDoubleScalar(tmp)

  tmp = n
  prhs(4) = mxCreateDoubleScalar(tmp)

  tmp = needF
  prhs(5) = mxCreateDoubleScalar(tmp)

  tmp = needG
  prhs(6) = mxCreateDoubleScalar(tmp)


  if (fgHandle /= 0) then
     prhs(2) = mxDuplicateArray(fgHandle)
  else
     prhs(2) = mxCreateDoubleScalar(tmp)

     if (objHandle /= 0) then
        nrhs = 7
        prhs(7) = mxDuplicateArray(objHandle)

        if (conHandle /= 0) then
           nrhs = 8
           prhs(8) = mxDuplicateArray(conHandle)
        end if
     end if
  end if

  ! Call Matlab: [F,G] = snwrapper(x,userfun,nF,n)
  call mexCallMatlab(nlhs, plhs, nrhs, prhs(1:nrhs), 'snwrapper')

  ! Objective and constraints
  if (needF > 0) then
     call checkRow(plhs(1), nF, 'F')
     call checkCol(plhs(1),  1, 'F')

     call mxCopyPtrToReal8(mxGetPr(plhs(1)), F, nF)
  end if


  ! Gradient and Jacobian
  if (needG > 0) then

     if (mxGetN(plhs(2)) == 0) then
        ! No derivatives are defined

     else
        if (mxGetN(plhs(2)) > 1) then
           ! Dense G given
           call checkRow(plhs(2), nF, 'G')
           call checkCol(plhs(2),  n, 'G')

           ! Convert to vector format for SNOPTA
           nlhs1    = 1
           nrhs1    = 3

           prhs1(1) = mxCreateDoubleMatrix(neG, 1, 0)
           prhs1(2) = mxCreateDoubleMatrix(neG, 1, 0)
           prhs1(3) = mxDuplicateArray(plhs(2))

           call mxCopyReal8ToPtr(riGfun(1:neG), mxGetPr(prhs1(1)), neG)
           call mxCopyReal8ToPtr(rjGvar(1:neG), mxGetPr(prhs1(2)), neG)

           ! Call Matlab: [G] = snfindG(iGfun,jGvar,Gdense)
           call mexCallMatlab(nlhs1, plhs1, nrhs1, prhs1, 'snfindG')

           ! Check dimensions of G
           call checkRow(plhs1(1), neG, 'G')
           call checkCol(plhs1(1),   1, 'G')


           ! Copy non-NaN entries of G.
           call mxCopyPtrToReal8(mxGetPr(plhs1(1)), G1, neG)

           do j = 1, neG
              if (mxIsNaN(G1(j)) == 0) then
                 G(j) = G1(j)
              end if
           end do


           ! Destroy array
           call mxDestroyArray(plhs1(1))
           call mxDestroyArray(prhs1(1))
           call mxDestroyArray(prhs1(2))
           call mxDestroyArray(prhs1(3))

        else
           ! Sparse G given
           ! Copy non-NaN entries of G.
           call mxCopyPtrToReal8(mxGetPr(plhs(2)), G1, neG)

           do j = 1, neG
              if (mxIsNaN(G1(j)) == 0) then
                 G(j) = G1(j)
              end if
           end do

        end if
     end if
  end if


  ! Destroy arrays
  call mxDestroyArray(plhs(1))
  call mxDestroyArray(plhs(2))

  call mxDestroyArray(prhs(1))
  call mxDestroyArray(prhs(2))
  call mxDestroyArray(prhs(3))
  call mxDestroyArray(prhs(4))
  call mxDestroyArray(prhs(5))
  call mxDestroyArray(prhs(6))
  if (nrhs >= 7) call mxDestroyArray(prhs(7))
  if (nrhs >= 8) call mxDestroyArray(prhs(8))

end subroutine matlabFG

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine matlabSTOP &
     ( iAbort,                                           &
       KTcond, mjrPrtlvl, minimize,                      &
       m, maxS, n, nb, nnCon0, nnCon, nnObj0, nnObj, nS, &
       itn, nMajor, nMinor, nSwap,                       &
       condZHZ, iObj, scaleObj, objAdd,                  &
       fObj, fMerit, penParm, step,                      &
       primalInf, dualInf, maxVi, maxViRel, hs,          &
       neJ, nlocJ, locJ, indJ, Jcol, negCon,             &
       scales, bl, bu, Fx, fCon, gCon, gObj,             &
       yCon, pi, rc, rg, x,                              &
       cu, lencu, iu, leniu, ru, lenru,                  &
       cw, lencw, iw, leniw, rw, lenrw )

  use mxsnWork, only : stopHandle, allocF, deallocF, &
                       tF, tFmul, tFlow, tFupp, tFstate
  implicit none

  logical, intent(in) :: KTcond(2)
  integer, intent(in) :: iObj, itn,                                   &
                         lencu, lencw, leniu, leniw, lenru, lenrw,    &
                         mjrPrtlvl, minimize, m, maxS, n, nb, neJ,    &
                         negCon, nlocJ, nnCon0, nnCon, nnObj0, nnObj, &
                         nMajor, nMinor, nS, nSwap,                   &
                         hs(nb), locJ(nlocJ), indJ(neJ),              &
                         iw(leniw)
  double precision, intent(in) ::                                     &
       condZHZ, scaleObj, objAdd, fObj, fMerit, penParm(4),           &
       maxViRel, maxVi, step, primalInf, dualInf,                     &
       scales(nb), bl(nb), bu(nb), Fx(nnCon0),                        &
       fCon(nnCon0), gCon(negCon), gObj(nnObj0), Jcol(neJ), pi(m),    &
       rc(nb), rg(maxS), yCon(nnCon0), x(nb), rw(lenrw)

  character(8), intent(in) :: cw(lencw)*8

  integer,          intent(inout) :: iu(leniu)
  double precision, intent(inout) :: ru(lenru)
  character(8),     intent(inout) :: cu(lencu)

  integer,          intent(out) :: iAbort

  !===================================================================
  ! snSTOP is called every major iteration.
  ! If iAbort > 0 on exit, the run is terminated.
  ! By specifying a custom version of snSTOP, the user can arrange for
  ! snopt to be terminated at any given major iteration.
  !
  ! 14 Oct 2004: First version of   snSTOP.
  ! 29 Aug 2007: Parameter list extended.
  ! 18 Dec 2015: New argument fObj.
  ! 22 Dec 2015: Matlab version
  !===================================================================
  integer*4, parameter :: nlhs = 1, nrhs = 22

  integer          :: i, iN, j, lkxN, nF, ObjRow, nkx
  double precision :: rtmp
  mwPointer        :: prhs(nrhs), plhs(nlhs)
  mwPointer        :: mxGetPr, mxDuplicateArray, mxCreateDoubleMatrix, &
                      mxCreateDoubleScalar
  double precision :: mxGetScalar

  iAbort = 0

  if ( stopHandle /= 0 ) then
     ! iAbort = snSTOP(...)
     prhs(1)  = mxDuplicateArray(stopHandle)

     ! Set up scalars
     rtmp     = itn
     prhs(2)  = mxCreateDoubleScalar(rtmp)

     rtmp     = nMajor
     prhs(3)  = mxCreateDoubleScalar(rtmp)

     rtmp     = nMinor
     prhs(4)  = mxCreateDoubleScalar(rtmp)

     prhs(5)  = mxCreateDoubleScalar(condZHZ)
     prhs(6)  = mxCreateDoubleScalar(fObj)
     prhs(7)  = mxCreateDoubleScalar(fMerit)
     prhs(8)  = mxCreateDoubleScalar(step)
     prhs(9)  = mxCreateDoubleScalar(primalInf)
     prhs(10) = mxCreateDoubleScalar(dualInf)
     prhs(11) = mxCreateDoubleScalar(maxVi)
     prhs(12) = mxCreateDoubleScalar(maxViRel)

     ! x, xlow, xupp, xmul, xstate
     prhs(13) = mxCreateDoubleMatrix(n,1,0)
     call mxCopyReal8ToPtr(x, mxGetPr(prhs(13)),n)

     prhs(14) = mxCreateDoubleMatrix(n,1,0)
     call mxCopyReal8ToPtr(bl, mxGetPr(prhs(14)),n)

     prhs(15) = mxCreateDoubleMatrix(n,1,0)
     call mxCopyReal8ToPtr(bu, mxGetPr(prhs(15)),n)

     prhs(16) = mxCreateDoubleMatrix(n,1,0)
     call mxCopyReal8ToPtr(rc, mxGetPr(prhs(16)),n)

     prhs(17) = mxCreateDoubleMatrix(n,1,0)
     call mxCopyReal8ToPtr(real(hs,8), mxGetPr(prhs(17)),n)

     ! F, Flow, Fupp, Fmul, Fstate,...
     nF     = iw(248)
     ObjRow = iw(103)
     nkx    = iw(247)
     lkxN   = iw(252) ! jN = kxN(j ) => col j of Jcol is variable jN

     call allocF( nF )

     do j  = n+1, nkx
        i  = j - n
        iN = iw(lkxN-1+j)
         if (iN == ObjRow) then
            if (iObj > 0) then
               tF(ObjRow) = fObj + x(n+iObj)
            else
               tF(ObjRow) = fObj
            end if
         else
            if (i <= nnCon) then
               tF(iN) = Fx(i)
            else
               tF(iN) =  x(j)
            end if
            tFstate(iN) = hs(j)
            tFmul(iN)   = rc(j)
            tFlow(iN)   = bl(j)
            tFupp(iN)   = bu(j)
         end if
      end do

     prhs(18) = mxCreateDoubleMatrix(nF,1,0)
     call mxCopyReal8ToPtr(tF, mxGetPr(prhs(18)), nF)

     prhs(19) = mxCreateDoubleMatrix(nF,1,0)
     call mxCopyReal8ToPtr(tFlow, mxGetPr(prhs(19)), nF)

     prhs(20) = mxCreateDoubleMatrix(nF,1,0)
     call mxCopyReal8ToPtr(tFupp, mxGetPr(prhs(20)), nF)

     prhs(21) = mxCreateDoubleMatrix(nF,1,0)
     call mxCopyReal8ToPtr(tFmul, mxGetPr(prhs(21)), nF)

     prhs(22) = mxCreateDoubleMatrix(nF,1,0)
     call mxCopyReal8ToPtr(real(tFstate,8), mxGetPr(prhs(22)), nF)

     call deallocF

     call mexCallMatlab(nlhs, plhs, nrhs, prhs, 'feval')

     iAbort = mxGetScalar(plhs(1))

     call mxDestroyArray(plhs(1))

     call mxDestroyArray(prhs(1))
     call mxDestroyArray(prhs(2))
     call mxDestroyArray(prhs(3))
     call mxDestroyArray(prhs(4))
     call mxDestroyArray(prhs(5))
     call mxDestroyArray(prhs(6))
     call mxDestroyArray(prhs(7))
     call mxDestroyArray(prhs(8))
     call mxDestroyArray(prhs(9))
     call mxDestroyArray(prhs(10))
     call mxDestroyArray(prhs(11))
     call mxDestroyArray(prhs(12))
     call mxDestroyArray(prhs(13))
     call mxDestroyArray(prhs(14))
     call mxDestroyArray(prhs(15))
     call mxDestroyArray(prhs(16))
     call mxDestroyArray(prhs(17))
     call mxDestroyArray(prhs(18))
     call mxDestroyArray(prhs(19))
     call mxDestroyArray(prhs(20))
     call mxDestroyArray(prhs(21))
     call mxDestroyArray(prhs(22))
  end if

end subroutine matlabSTOP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
