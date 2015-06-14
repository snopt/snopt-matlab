function [x,F,inform,xmul,Fmul] = snopt( x, xlow, xupp, xmul, xstate,...
					 Flow, Fupp, Fmul, Fstate,...
					 userfun, varargin );
% [x,F,inform,xmul,Fmul] = snopt( x, xlow, xupp, xmul, xstate, ...
%                                   Flow, Fupp, Fmul, Fstate, ...
%                                   userfun,varargin );
%
% This function solves the nonlinear optimization problem:
% minimize:
%            F(ObjRow) + ObjAdd
% subject to:
%            xlow <= x <= xupp
%            Flow <= F <= Fupp
% where:
%  x        is the column vector of initial values of the unknowns.
%  F        is a vector of objective and constraint functions specified
%           in the m-file userfun.
%  ObjRow   is the objective row of F (default ObjRow = 1).
%  userfun  is the handle for a Matlab funcion that defines the elements of F
%           and optionally, their derivatives.
%
% Additional arguments allow the specification of more detailed
% problem information.
%
% Calling sequence 1:
%  [x,F,inform,xmul,Fmul] = snopt ( x, xlow, xupp, xmul, xstate,...
%                                     Flow, Fupp, Fmul, Fstate, userfun )
%
% Calling sequence 2:
%  [x,F,inform,xmul,Fmul] = snopt ( x, xlow, xupp, xmul, xstate,...
%                                     Flow, Fupp, Fmul, Fstate, userfun,...
%                                     ObjAdd, ObjRow )
% Calling sequence 3:
%  [x,F,inform,xmul,Fmul] = snopt ( x, xlow, xupp, xmul, xstate,...
%                                     Flow, Fupp, Fmul, Fstate, userfun,...
%                                     A, iAfun, jAvar, iGfun, jGvar )
% Calling sequence 4:
%  [x,F,inform,xmul,Fmul] = snopt ( x, xlow, xupp, xmul, xstate,...
%                                     Flow, Fupp, Fmul, Fstate, userfun,...
%                                     ObjAdd, ObjRow,
%                                     A, iAfun, jAvar, iGfun, jGvar )
%
% INPUT:
%  x             is the initial guess for x.
%
%  xlow, xupp    are the upper and lower bounds on x.
%
%  xmul          contains the initial multipliers for x.
%
%  xstate        are the states of the variables x.
%
%  Flow, Fupp    are the upper and lower bounds on F.
%
%  Fmul          contains the initial multipliers for F.
%
%  Fstate        are the states of the constraints F.
%
%  userfun       is the user-defined Matlab function that computes the
%                objective and constraint functions and their corresponding
%                derivatives.
%                userfun can be either a function handle or a string.
%
%                WARNING: If arguments iAfun, jAvar, and A, are provided,
%                then it is crucial that the associated linear terms
%                are not included in the calculation of F in userfun.
%                **Example:
%                      >> [f,G] = feval(userfun,x);
%                      >> F  = sparse(iAfun,jAvar,A)*x + f
%                      >> DF = sparse(iAfun,jAvar,A) + sparse(iGfun,jGvar,G)
%                (where DF denotes F').
%                G must be either a dense matrix, or a vector of derivatives
%                in the same order as the indices iGfun and jGvar, i.e.,
%                if G is a vector, the Jacobian of f is given
%                by sparse(iGfun,jGvar,G).
%
%                **More details on G:
%                The user maybe define, all, some, or none of the
%                entries of G.  If the user does NOT intend to
%                supply ALL nonzero entries of G, it is imperative
%                that the proper derivative level is set prior to
%                a call to snopt():
%                     >> snseti("Derivative option",k);
%                where k = 0 or 1.  Meaning:
%                    1 -- Default.  All derivatives are provided.
%                    0 -- Some derivatives are not provided.
%                For the case k = 0 (ONLY), G may be returned as
%                an empty array [].
%                For the case k = 0, the user must denote
%                unknown NONZERO elements of G by NaN.
%                **Example: (vector case)
%                   >> G = [1, NaN, 3, NaN, -5]';
%                  or (full matrix case)
%                   >> G = [1, 0, NaN; 0 NaN 2; 0 0 3];
%
% ObjAdd         is the constant term in the objective function. The default
%                is 0.  If ObjRow /= 0, then ObjAdd is added to F(ObjRow).
%
% ObjRow         indicates which row of F acts as the objective function.
%                If not specified, the default is 1.
%
% A,iAfun,jAvar  hold the constant elements in the Jacobian of F in
%                coordinate form.
%                If i = iAfun(k) and j = jAvar(k), then A(k) is the (i,j)-th
%                element in A.
%
% iGfun, jGvar   hold the indices of the nonlinear elements in the Jacobian
%                of F.
%
% More IMPORTANT details:
%   1) The indices (iAfun,jAvar) must be DISJOINT from (iGfun,jGvar).
%      A nonzero element in F' must be either an element of G or an
%      element of A, but not the sum of the two.
%
%   2) If the user does not wish to provide iAfun, jAvar, iGfun,
%      jGvar, then snopt() will determine them by calling snJac().
%
%      WARNING: In this case, the derivative level will be set to zero
%      if constant elements exist.  This is because the linear
%      elements have not yet been deleted from the definition of
%      userfun.  Furthermore, if G is given in vector form, the
%      ordering of G may not necessarily correspond to (iGfun,jGvar)
%      computed by snJac().

n      = length(x);
ObjAdd = 0;
ObjRow = 1;

if (ischar(userfun))
  userFG = str2func(userfun);
else
  userFG = userfun;
end

F0 = userFG(x);
nF = length(F0);


if nargin == 10,

  % Calling sequence 1
  % Derivatives are estimated by differences.
  % Call snJac to estimate the pattern of nonzeros for the Jacobian.

  [A,iAfun,jAvar,iGfun,jGvar] = snJac(userFG,x,xlow,xupp,nF);
elseif nargin == 12

  % Calling sequence 2
  % Derivatives are estimated by differences.
  % Call snJac to estimate the pattern of nonzeros for the Jacobian.

  [A,iAfun,jAvar,iGfun,jGvar] = snJac(userFG,x,xlow,xupp,nF);
  ObjAdd = varargin{1};
  ObjRow = varargin{2};
elseif nargin == 15

  % Calling sequence 3
  % The user is providing derivatives.

  A      = varargin{1};
  iAfun  = varargin{2};
  jAvar  = varargin{3};
  iGfun  = varargin{4};
  jGvar  = varargin{5};
elseif ( nargin == 17 )

  % Calling sequence 4
  % The user is providing derivatives.

  ObjAdd = varargin{1};
  ObjRow = varargin{2};
  A      = varargin{3};
  iAfun  = varargin{4};
  jAvar  = varargin{5};
  iGfun  = varargin{6};
  jGvar  = varargin{7};
end

solveopt = 1;
[x,F,inform,xmul,Fmul] = snoptmex( solveopt, x, ...
				   xlow, xupp, xmul, xstate, ...
				   Flow, Fupp, Fmul, Fstate, ...
				   ObjAdd, ObjRow, A, iAfun, jAvar, ...
				   iGfun, jGvar, userFG );
