function  [A,iAfun,jAvar,iGfun,jGvar] = snJac(userfun,x0,xlow,xupp,nF)
%function [A,iAfun,jAvar,iGfun,jGvar] = snJac(usrfun,x0,xlow,xupp,nF)
%         Finds the coordinate structure for the Jacobian.

findJacOption = 17;

if (ischar(userfun))
  userFG = str2func(userfun);
else
  userFG = userfun;
end
[A,iAfun,jAvar,iGfun,jGvar] = snoptmex(findJacOption,userFG,x0,xlow,xupp,nF);
