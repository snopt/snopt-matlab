function [x,Obj,INFO,lambda] = sqopt( name, Hx, c, x0, xl, xu, A, al, au )
% [x,Obj,INFO,pi,rc] = sqopt ( name, Hx, c, x0, xl, xu, A, al, au );
%
% This function solves the quadratic optimization problem:
%   minimize:
%              c'x + x'*H*x
%   subject to:
%            xl <=  x <= xu
%            al <= Ax <= au
% where:
%  x        is the column vector of initial values of the unknowns
%  xl, xu   are the lower and upper bounds of the variables
%  c        is the linear term of the objective
%  H        is the Hessian matrix of the objective
%  A        is the linear constraint matrix
%  al, au   are the lower and upper bounds of the linear constraints
%
% Calling sequences:
%  x = sqopt ( name, Hx, c, x0, xl, xu, A, al, au )
%
%  [x]                = sqopt ( name, Hx, c, x0, xl, xu, A, al, au )
%  [x,Obj]            = sqopt ( name, Hx, c, x0, xl, xu, A, al, au )
%  [x,Obj,INFO]       = sqopt ( name, Hx, c, x0, xl, xu, A, al, au )
%  [x,Obj,INFO,pi,rc] = sqopt ( name, Hx, c, x0, xl, xu, A, al, au )
%
%
% INPUT:
%  name              is the name of the problem
%
%  x0                is the initial guess for x
%
%  Hx                is a Matlab function that computes H*x for a given x.
%                    Hx can be a Matlab function handle or a string
%
%  c                 is the linear term of the quadratic objective
%
%  xl, xu            are the upper and lower bounds on x
%
%  A                 is the linear constraint matrix; A cannot be empty, but
%                    if there are no constraints, set A = zeros(0,n), where n
%                    is the number of variables
%
%  al, au            are the upper and lower bounds on the linear constraints A*x
%
% OUTPUT:
%  x                 is the final point
%
%  Obj               is the objective at the final point x
%
%  INFO              is the exit flag returned by SQOPT
%
%  lambda            is a structure containing the multipliers at the point x
%                    lambda.lb      are the lower bound multipliers
%                    lambda.ub      are the upper bound multipliers
%                    lambda.linear  are the linear constraint multipliers
%

solveopt = 1;

if ( nargin ~= 9 ),
  error('Wrong number of input arguments for sqopt');
end

if (ischar(Hx))
  userHx = str2func(Hx);
else
  userHx = Hx;
end

[m,n]                = size(A);
[neA,indA,locA,valA] = crd2spr(A);

[x,Obj,INFO,piA,y] = sqoptmex( solveopt, name, m, n, userHx, c, x0, xl, xu, ...
			       neA, indA, locA, valA, al, au );

if ( nargout >= 4 )
  n    = size(x);
  zero = zeros(n);

  lambda.lb         = max(y,zero);
  lambda.ub         = min(y,zero);
  lambda.linear     = piA;
end
