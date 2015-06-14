function [x,fval,exitflag,lambda] = snsolve(userobj,x0,A,b,varargin)

% A wrapper for snopt to make it look like fmincon.
%   [x,fval,exitflag,lambda] = snsolve(myobj,x0,A,b)
%   [x,fval,exitflag,lambda] = snsolve(myobj,x0,A,b,Aeq,beq)
%   [x,fval,exitflag,lambda] = snsolve(myobj,x0,A,b,Aeq,beq,xlow,xupp)
%   [x,fval,exitflag,lambda] = snsolve(myobj,x0,A,b,Aeq,beq,xlow,xupp,nonlcon)
%
% snsolve and fmincon assume problems are of the form:
%    minimize    f(x)
%   such that    c(x)   <= 0,
%                c_eq(x) = 0,
%                Ax     <= b,
%                A_eq x  = b_eq,
%                xlow   <= x <= xupp.
%
% 26 November 2012.


if isrow(x0),
  [xi,n] = size(x0);
elseif iscolumn(x0),
  [n,xi] = size(x0);
elseif isempty(x0);
  fprintf('Error: You must provide a non-empty starting point.\n');
  return
else
  fprintf('Error: You must provide a row or column vector for the starting point.\n');
  return
end

[mi,n0]  =  size(A);
nonlcon =  @dummyCon;

if (ischar(userobj))
  myobj = str2func(userobj);
else
  myobj = userobj;
end

if     nargin == 4,
  me       = 0;
  nli      = 0;
  nle      = 0;
  Aeq      = [];
  beq      = [];
  xlow     = [];
  xupp     = [];

elseif nargin == 6,
  Aeq      = varargin{1};
  beq      = varargin{2};

  [me,n0]  = size(Aeq);
  nli      = 0;
  nle      = 0;
  xlow     = [];
  xupp     = [];

elseif nargin == 8,
  Aeq      = varargin{1};
  beq      = varargin{2};
  xlow     = varargin{3};
  xupp     = varargin{4};

  [me,n0]  = size(Aeq);
  nli      = 0;
  nle      = 0;

elseif nargin == 9,
  Aeq      = varargin{1};
  beq      = varargin{2};
  xlow     = varargin{3};
  xupp     = varargin{4};
  nonlconU = varargin{5};

  if (ischar(nonlconU))
    nonlcon = str2func(nonlconU);
  else
    nonlcon = nonlconU;
  end

  [me,n0]  = size(Aeq);
  [c,ceq]  = feval(nonlcon,x0);
  nli      = size(c,1);
  nle      = size(ceq,1);

else
  error('Wrong number of input arguments')
end

nCon = 1 + nli + nle + mi + me;   % number of constraints (size of F(x))


% snoptA problem format:
%    minimize    F_obj (x)
%   such that  l_f <= F(x) <= u_f
%              l   <=   x  <= u
%
%           [ F_obj   ]            [ F_0'     ]   [  0    ]     1
%           [ c(x)    ]            [ c'(x)    ]   [  0    ]     nli
%           [ c_eq(x) ]            [ c_eq'(x) ]   [  0    ]     nle
%   F(x) =  [ Ax      ]    F'(x) = [ 0        ] + [  A    ]     mi
%           [ A_eq x  ]            [ 0        ]   [  A_eq ]     me
%                                    "G(x)"         "A"

[iAfun,jAvar,Aij] = find( [ zeros(1+nli+nle,n); A; Aeq ]);
[iGfun,jGvar,Gij] = find( [ ones(1+nli+nle,n); zeros(mi+me,n) ]);

xmul    =  zeros(n,1);
xstate  =  zeros(n,1);
Fmul    =  zeros(nCon,1);
Fstate  =  zeros(nCon,1);
ObjAdd  =  0;
ObjRow  =  1;
Flow    = [ -inf; -inf*ones(nli,1); zeros(nle,1); -inf*ones(mi,1); zeros(me,1) ];
Fupp    = [  inf; zeros(nli,1); zeros(nle,1); b; beq ];


% Solve the problem!
solveopt = 1;
[x,F,exitflag,xmul,Fmul] = snoptmex( solveopt, x0, ...
				     xlow, xupp, xmul, xstate, ...
				     Flow, Fupp, Fmul, Fstate, ...
				     ObjAdd, ObjRow, ...
				     Aij, iAfun, jAvar, ...
				     iGfun, jGvar, ...
				     myobj, nonlcon );

fval              = feval(myobj,x);
lambda.lower      = xmul;
lambda.upper      = xmul;
lambda.ineqnonlin = Fmul(2:nli+1);
lambda.eqnonlin   = Fmul(nli+2:nli+nle+1);
lambda.ineqlin    = Fmul(nli+nle+2:nli+nle+mi+1);
lambda.eqlin      = Fmul(nli+nle+mi+2:nCon);


function [c,ceq] = dummyCon(x)

c = []; ceq = [];
