function  [A,iAfun,jAvar,iGfun,jGvar] = snJac(userfun,x0,xlow,xupp,nF)
%function [A,iAfun,jAvar,iGfun,jGvar] = snJac(usrfun,x0,xlow,xupp,nF)
%         Finds the coordinate structure for the Jacobian.


userFG = checkFun(userfun,'SNOPT','userfun');
[A,iAfun,jAvar,iGfun,jGvar] = snoptmex(17,userFG,x0,xlow,xupp,nF);
