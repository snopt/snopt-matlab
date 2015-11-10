!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File sqoptmex.F90
! Mex function for SQOPT7.
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
  ! Mex function for SQOPT7
  !
  ! Option      Action
  !    1        Solve the problem (quadprog-style)
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
  !
  ! 13 Dec 2013: Current version.
  ! 01 May 2015: Added ability to modify initial amount of workspace
  ! 09 Nov 2015: Added states, iteration counts as output
  !=====================================================================
  ! Matlab
  mwPointer        :: mxGetN, mxGetPr
  integer          :: mxIsChar
  double precision :: mxGetScalar

  ! SQOPT
  character        :: filename*80
  integer          :: info, iOpt, strlen
  double precision :: rOpt, rleniw, rlenrw
  external         :: sqInit

  ! Get option.
  if (nrhs < 1) call mexErrMsgTxt('Need an option input argument')
  rOpt = mxGetScalar(prhs(1))
  iOpt = rOpt


  ! Deal with on/off screen, file, summary files first.
  if (iOpt == snOpenP) then

     if (nrhs /= 2) call mexErrMsgTxt('Wrong number of input arguments')

     info = mxIsChar(prhs(2))
     if (info /= 1) call mexErrMsgTxt('Need a filename string')

     strlen = mxGetN(prhs(2))
     if (strlen > 80) call mexErrMsgTxt('Print filename is too long')

     if (strlen > 0) then
        call mxGetString(prhs(2), filename, strlen)
     else
        call mexErrMsgTxt('Empty print filename')
     end if

     if (printOpen) close(iPrint)

     open(iPrint, file=filename, status='unknown')
     printOpen= .true.
     return

  else if (iOpt == snOpenS) then

     if (nrhs /= 2) call mexErrMsgTxt('Wrong number of input arguments')

     info = mxIsChar(prhs(2))
     if (info /= 1) call mexErrMsgTxt('Need a filename string')

     strlen = mxGetN(prhs(2))
     if (strlen > 80) call mexErrMsgTxt('Summary filename is too long')

     if (strlen > 0) then
        call mxGetString(prhs(2), filename, strlen)
     else
        call mexErrMsgTxt('Empty summary filename')
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
          call mexErrMsgTxt('Workspace size must be at least 500')
     return

  end if


  if (firstCall) then
     allocate(cw(lencw), iw(leniw), rw(lenrw))

     callType = userCall
     call sqInit(iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw)
     callType = systemCall

     memCall   = .false.
     firstCall = .false.
  end if

  ! Do whatever we need to do.
  if      (iOpt == snSolve) then

     callType = userCall
     call sqmxSolve(nlhs, plhs, nrhs, prhs)
     callType = systemCall

  else if (iOpt == snSetXX .or. &
            iOpt == snSetIX .or. &
            iOpt == snSetRX .or. &
            iOpt == snGetXX .or. &
            iOpt == snGetCX .or. &
            iOpt == snGetIX .or. &
            iOpt == snGetRX) then

     callType = userCall
     call snmxOptions(iOpt, nlhs, plhs, nrhs, prhs)
     callType = systemCall

  else if (iOpt == snSpecs) then

     callType = userCall
     call snmxSpecs(nlhs, plhs, nrhs, prhs)
     callType = systemCall

  else if (iOpt == snEnd) then
     if (printOpen) close(iPrint)
     printOpen= .false.

     if (summOpen) close(iSumm)
     summOpen  = .false.
     memCall   = .false.
     firstCall = .true.

     if (allocated(cw)) deallocate(cw)
     if (allocated(iw)) deallocate(iw)
     if (allocated(rw)) deallocate(rw)

  end if

  return

end subroutine mexFunction

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine sqmxSolve (nlhs, plhs, nrhs, prhs)
  use mxsnWork

  implicit none
  integer*4  :: nlhs, nrhs
  mwPointer  :: prhs(*), plhs(*)
  !---------------------------------------------------------------------
  ! Solve the problem
  ! The matlab call iss
  !   [x, fval, exit, itn, y ] =
  !       qpsolve (Hx, c, A, b, Aeq, beq, lb, ub, x0)
  !
  ! where
  !   Hx        is a user-defined subroutine to compute H*x
  !   c         is linear terms of the objective
  !   x0        is the initial point
  !   A, b      are the linear inequality constraints A*x <= b
  !   Aeq, beq  are the linear equality constraints Aeq*x = beq
  !   lb, ub    are the lower and upper bounds on x
  !---------------------------------------------------------------------
  ! Matlab
  mwPointer        :: mxDuplicateArray, mxGetM, mxGetN, mxGetPr, &
                      mxCreateDoubleMatrix, mxCreateDoubleScalar
  integer          :: mxIsChar, mxIsClass, mxIsEmpty
  double precision :: mxGetScalar

  ! SQOPT
  character*8      :: probName, Start
  integer          :: Errors, info, i1, i2, strlen
  integer          :: iObj, m, n, nnH, ncObj, neA, &
                      nNames, mincw, miniw, minrw, nInf, nS
  double precision :: rinfo, Obj, ObjAdd, sInf
  external         :: sqopt, matlabHx

  double precision, parameter   :: infBnd = 1.0d+20

  character*8,      allocatable :: Names(:)
  integer,          allocatable :: hEtype(:), hs(:), indA(:), locA(:)
  double precision, allocatable :: cObj(:), x(:), pi(:), rc(:), bl(:), bu(:), &
                                   valA(:), rlocA(:), rindA(:)


  ! Check number of input and output arguments.
  if (nrhs /= 15) call mexErrMsgTxt('Wrong number of input arguments')


  !---------------------------------------------------------------------
  ! Problem name
  !---------------------------------------------------------------------
  info = mxIsChar (prhs(2))
  if (info /= 1) call mexErrMsgTxt('Wrong input type for problem name')

  strlen = mxGetN(prhs(2))
  if (strlen > 8) strlen = 8

  probName = ''
  call mxGetString(prhs(2), probName, strlen)


  !---------------------------------------------------------------------
  ! Number of constraints and variables
  !---------------------------------------------------------------------
  m   = mxGetScalar(prhs(3))
  n   = mxGetScalar(prhs(4))
  nnH = n

  !---------------------------------------------------------------------
  ! Hessian matrix
  !---------------------------------------------------------------------
  info = mxIsClass(prhs(5), 'function_handle')
  if (info /= 1) call mexErrMsgTxt('Wrong input type for Hx')
  HxHandle = mxDuplicateArray(prhs(5))


  !---------------------------------------------------------------------
  ! Linear term of objective
  !---------------------------------------------------------------------
  if (mxIsEmpty(prhs(6)) > 0) then
     ncObj = 0
     allocate(cObj(1))
  else
     ncObj = mxGetM(prhs(6))
     allocate(cObj(ncObj))

     call checkRow(prhs(6), ncObj, 'cObj')
     call checkCol(prhs(6),     1, 'cObj')

     call mxCopyPtrToReal8(mxGetPr(prhs(6)), cObj(1:ncObj), ncObj)
  end if


  !---------------------------------------------------------------------
  ! Allocate SQOPT space
  !---------------------------------------------------------------------
  allocate(x(n+m), hs(n+m), hEtype(n+m), pi(m), rc(n+m))
  allocate(bl(n+m), bu(n+m))

  x      = 0.0
  hs     = 0
  hEtype = 0

  if (.not. mxIsEmpty(prhs(7)) > 0) then
     call checkRow(prhs(7), n, 'x0')
     call checkCol(prhs(7), 1, 'x0')
     call mxCopyPtrToReal8(mxGetPr(prhs(7)), x(1:n), n)
  end if


  !---------------------------------------------------------------------
  ! Set bounds on variables
  !---------------------------------------------------------------------
  ! Lower and upper bounds on x
  ! Check dimensions of lower and upper bounds on x
  if (mxIsEmpty(prhs(8)) > 0) then
     bl(1:n) = -infBnd
  else
     call checkRow(prhs(8), n, 'lb')
     call checkCol(prhs(8), 1, 'lb')
     call mxCopyPtrToReal8(mxGetPr(prhs(8)), bl(1:n), n)
  end if

  if (mxIsEmpty(prhs(9)) > 0) then
     bu(1:n) = infBnd
  else
     call checkRow(prhs(9), n, 'ub')
     call checkCol(prhs(9), 1, 'ub')
     call mxCopyPtrToReal8(mxGetPr(prhs(9)), bu(1:n), n)
  end if


  !---------------------------------------------------------------------
  ! Get the linear constraint matrix
  !---------------------------------------------------------------------
  neA = mxGetScalar(prhs(10))
  if (neA > 0) then
     call checkRow(prhs(11), neA, 'indA')
     call checkCol(prhs(11),   1, 'indA')

     call checkRow(prhs(12), n+1, 'locA')
     call checkCol(prhs(12),   1, 'locA')

     call checkRow(prhs(13), neA, 'valA')
     call checkCol(prhs(13),   1, 'valA')

     allocate(indA(neA), valA(neA), locA(n+1))
     allocate(rindA(neA), rlocA(n+1))

     call mxCopyPtrToReal8(mxGetPr(prhs(11)), rindA, neA)
     call mxCopyPtrToReal8(mxGetPr(prhs(12)), rlocA, n+1)
     call mxCopyPtrToReal8(mxGetPr(prhs(13)), valA, neA)

     indA(1:neA) = rindA(1:neA)
     locA(1:n+1) = rlocA(1:n+1)

     deallocate(rindA, rlocA)
  else
     neA  = 0
     allocate(indA(1), locA(1), valA(1))
  end if


  !---------------------------------------------------------------------
  ! Set the constraint bounds
  !---------------------------------------------------------------------
  if (m > 0) then
     i1 = 1+n
     i2 = m+n

     ! Check dimension of lA and uA
     if (mxIsEmpty(prhs(14)) > 0) then
        bl(i1:i2) = -infBnd
     else
        call checkRow(prhs(14), m, 'al')
        call checkCol(prhs(14), 1, 'al')
        call mxCopyPtrToReal8(mxGetPr(prhs(14)), bl(i1:i2), m)
     end if

     if (mxIsEmpty(prhs(15)) > 0) then
        bu(i1:i2) = infBnd
     else
        call checkRow(prhs(15), m, 'au')
        call checkCol(prhs(15), 1, 'au')
        call mxCopyPtrToReal8(mxGetPr(prhs(15)), bu(i1:i2), m)
     end if
  end if


  !---------------------------------------------------------------------
  ! Allocate other space for SNOPT
  !---------------------------------------------------------------------
  nNames = 1
  allocate(Names(nNames))

  iObj   = 0
  ObjAdd = 0.0


  !---------------------------------------------------------------------
  ! Set workspace
  !---------------------------------------------------------------------
100 if (.not. memCall) then
     call sqMem &
          (INFO, m, n, neA, ncObj, nnH, &
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

     call sqSeti &
          ('Total character workspace', lencw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw)
     call sqSeti &
          ('Total integer   workspace', leniw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw)
     call sqSeti &
          ('Total real      workspace', lenrw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw)
  end if


  !---------------------------------------------------------------------
  ! Solve the problem
  !---------------------------------------------------------------------
  Start  = 'Cold'  ! cold start

  call sqopt &
         (Start, matlabHx, m, n, neA, nNames, &
           ncObj, nnH, iObj, ObjAdd, probName, &
           valA, indA, locA, bl, bu, cObj, Names, &
           hEtype, hs, x, pi, rc, &
           INFO, mincw, miniw, minrw, nS, nInf, sInf, Obj, &
           cw, lencw, iw, leniw, rw, lenrw, &
           cw, lencw, iw, leniw, rw, lenrw)

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
  if (nlhs > 1) plhs(2) = mxCreateDoubleScalar (Obj)


  ! Exit flag
  rinfo = info
  if (nlhs > 2) plhs(3) = mxCreateDoubleScalar (rinfo)

  ! Multipliers for bounds
  if (nlhs > 3) then
     plhs(4) = mxCreateDoubleMatrix (m, 1, 0)
     call mxCopyReal8ToPtr(pi, mxGetPr(plhs(4)), m)
  end if


  ! Multipliers for linear inequalities
  if (nlhs > 4) then
     plhs(5) = mxCreateDoubleMatrix (n, 1, 0)
     call mxCopyReal8ToPtr(rc(1:n), mxGetPr(plhs(5)), n)
  end if

  ! States
  if ( nlhs > 5 ) then
     plhs(6) = mxCreateDoubleMatrix (n+m, 1, 0)
     call mxCopyReal8ToPtr(real(hs,8), mxGetPr(plhs(6)), n+m)
  end if

  ! Number of iterations
  rinfo = iw(421)
  if (nlhs > 6) plhs(7) = mxCreateDoubleScalar(rinfo)


  ! Deallocate memory
  if (HxHandle /= 0) call mxDestroyArray(HxHandle)
  HxHandle = 0

  if (allocated(x))      deallocate(x)
  if (allocated(pi))     deallocate(pi)
  if (allocated(rc))     deallocate(rc)
  if (allocated(bl))     deallocate(bl)
  if (allocated(bu))     deallocate(bu)
  if (allocated(Names))  deallocate(Names)
  if (allocated(hs))     deallocate(hs)
  if (allocated(hEtype)) deallocate(hEtype)

  if (allocated(indA))   deallocate(indA)
  if (allocated(locA))   deallocate(locA)
  if (allocated(valA))   deallocate(valA)

end subroutine sqmxSolve

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snmxOptions (iOpt, nlhs, plhs, nrhs, prhs)
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

  integer          :: sqGet
  external         :: sqSet, sqSetI, sqSetR, &
                      sqGet, sqGetC, sqGetI, sqGetR


  if (iOpt == snSetIX .or. iOpt == snSetRX) then
     if (nrhs /= 3) call mexErrMsgTxt('Wrong number of input arguments')
  else
     if (nrhs /= 2) call mexErrMsgTxt('Wrong number of input arguments')
  end if


  ! Get string
  strlen = mxGetN(prhs(2))
  if (strlen > 50) call mexErrMsgTxt('Option string is too long')

  if (strlen > 0) then
     call mxGetString(prhs(2), buffer, strlen)
  else
     call mexErrMsgTxt('Empty option string')
  end if


  if      (iOpt == snSetXX) then
     call sqSet (buffer, iPrint, iSumm, Errors, &
                  cw, lencw, iw, leniw, rw, lenrw)

  else if (iOpt == snSetIX) then

     rvalue = mxGetScalar(prhs(3))
     ivalue = rvalue

     call sqSetI (buffer, ivalue, iPrint, iSumm, Errors, &
                   cw, lencw, iw, leniw, rw, lenrw)

  else if (iOpt == snSetRX) then

     rvalue = mxGetScalar(prhs(3))

     call sqSetR (buffer, rvalue, iPrint, iSumm, Errors, &
                   cw, lencw, iw, leniw, rw, lenrw)

  else if (iOpt == snGetXX) then

     ivalue  = sqGet (buffer, Errors, cw, lencw, iw, leniw, rw, lenrw)

     rvalue  = ivalue
     plhs(1) = mxCreateDoubleScalar (rvalue)

  else if (iOpt == snGetCX) then

     call sqGetC (buffer, cvalue, Errors, &
                   cw, lencw, iw, leniw, rw, lenrw)

     plhs(1) = mxCreateString(cvalue)

  else if (iOpt == snGetIX) then

     call sqGetI (buffer, ivalue, Errors, &
                   cw, lencw, iw, leniw, rw, lenrw)

     rvalue = ivalue
     plhs(1) = mxCreateDoubleScalar (rvalue)

  else if (iOpt == snGetRX) then

     call sqGetR (buffer, rvalue, Errors, &
                   cw, lencw, iw, leniw, rw, lenrw)

     plhs(1) = mxCreateDoubleScalar (rvalue)

  end if

end subroutine snmxOptions

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snmxSpecs (nlhs, plhs, nrhs, prhs)
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

  external         :: sqSpec


  if (nrhs /= 2) call mexErrMsgTxt('Wrong number of input arguments')
  if (nlhs /= 1) call mexErrMsgTxt('Wrong number of output arguments')


  strlen = mxGetN(prhs(2))
  if (strlen > 120) call mexErrMsgTxt('Specs filename is too long')

  if (strlen > 0) then
     call mxGetString(prhs(2), filename, strlen)
  else
     call mexErrMsgTxt('Empty spc filename')
  end if

  open(iSpecs, file=filename, status='unknown')
  call sqSpec(iSpecs, info, cw, lencw, iw, leniw, rw, lenrw)
  rewind (iSpecs)
  close(iSpecs)

  ! sqSpec will return info == 101 or 107 if successful
  ! The matlab version returns 0 if successful
  if (info == 101 .or. info == 107) then
     rvalue = 0
  else
     rvalue = 1
  end if

  plhs(1) = mxCreateDoubleScalar (rvalue)

end subroutine snmxSpecs

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine matlabHx (nnH, x, Hx, Status, &
                      cu, lencu, iu, leniu, ru, lenru)
  use mxsnWork
  implicit none

  integer          :: Status, nnH, lencu, leniu, lenru, iu(leniu)
  double precision :: x(nnH), Hx(nnH), ru(lenru)
  character*8      :: cu(lencu)

  !---------------------------------------------------------------------
  ! Matlab callback to evaluate H*x.
  !---------------------------------------------------------------------
  integer*4 :: nlhs, nrhs, nlhs1, nrhs1
  mwPointer :: prhs(2), plhs(1)
  mwPointer :: mxCreateDoubleMatrix, mxGetPr, mxDuplicateArray

  nlhs = 1
  nrhs = 2

  prhs(1) = mxDuplicateArray(HxHandle)

  prhs(2) = mxCreateDoubleMatrix (nnH, 1, 0)
  call mxCopyReal8ToPtr(x, mxGetPr(prhs(2)), nnH)


  ! Call Matlab: [Hx] = userHx(x)
  call mexCallMatlab(nlhs, plhs, nrhs, prhs, 'feval')

  call checkRow(plhs(1), nnH, 'Hx')
  call checkCol(plhs(1),   1, 'Hx')

  call mxCopyPtrToReal8(mxGetPr(plhs(1)), Hx, nnH)

  ! Destroy arrays
  call mxDestroyArray(plhs(1))

  call mxDestroyArray(prhs(1))
  call mxDestroyArray(prhs(2))

end subroutine matlabHx

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
