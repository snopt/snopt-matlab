function [x,fval,exitflag,output,lambda,states] = solve_snfmincon(start, stopFun, ...
						  name, userfun, x0, xlow, ...
						  xupp, xmul, xstate, ...
						  Flow, Fupp, Fmul, Fstate, ...
						  ObjAdd, ObjRow, ...
						  Aij, iAfun, jAvar, iGfun,jGvar, ...
						  n, nonlin_ineq, nonlin_eq, ...
						  linear_ineq, linear_eq)
% function [x,fval,exitflag,output,lambda,states] = solve_snfmincon(start, stopFun, ...
% 						  name, userfun, x0, xlow, ...
% 						  xupp, xmul, xstate, ...
% 						  Flow, Fupp, Fmul, Fstate, ...
% 						  ObjAdd, ObjRow, ...
% 						  Aij, iAfun, jAvar, iGfun, ...
% 						  jGvar)


mexopt = 1;
[x,F,exitflag,xmul,Fmul,xstate,Fstate,itn,mjritn] = ...
    snoptmex(mexopt, ...
	     start, stopFun, name, ...
	     userfun, x0, ...
	     xlow, xupp, xmul, xstate, ...
	     Flow, Fupp, Fmul, Fstate, ...
	     ObjAdd, ObjRow, ...
	     Aij, iAfun, jAvar, iGfun, jGvar);


fval     = F(1);
zero     = zeros(n,1);
states.x = xstate(1:n);
lambda.x = xmul(1:n);

if nonlin_ineq > 0,
  i1 = 1+1; i2 = i1-1 + nonlin_ineq;
  lambda.ineqnonlin = Fmul(i1:i2);
  states.ineqnonlin = Fstate(i1:i2);
else
  lambda.ineqnonlin = [];
  states.ineqnonlin = [];
end

if nonlin_eq > 0,
  i1 = 1+nonlin_ineq; i2 = i1-1 + nonlin_eq;
  lambda.eqnonlin = Fmul(i1:i2);
  states.eqnonlin = Fmul(i1:i2);
else
  lambda.eqnonlin = [];
  states.eqnonlin = [];
end

if linear_ineq > 0,
  i1 = 1+nonlin_ineq+nonlin_eq; i2 = i1-1 + linear_ineq;
  lambda.ineqlin = Fmul(i1:i2);
  states.ineqlin = Fmul(i1:i2);
else
  lambda.ineqlin = [];
  states.ineqlin = [];
end

if linear_eq > 0,
  i1 = 1+nonlin_ineq+nonlin_eq+linear_ineq; i2 = i1-1 + linear_eq;
  lambda.eqlin = Fmul(i1:i2);
  states.eqlin = Fmul(i1:i2);
else
  lambda.eqlin = [];
  states.eqlin = [];
end

output.info       = exitflag;
output.iterations = itn;
output.majors     = mjritn;
