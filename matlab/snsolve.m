function [x,fval,exitflag,lambda,states,output] = snsolve(userobj,x0,A,b,varargin)

% A wrapper for snopt to make it look like fmincon.
%   [...] = snsolve(myobj,x0,A,b)
%   [...] = snsolve(myobj,x0,A,b,Aeq,beq)
%   [...] = snsolve(myobj,x0,A,b,Aeq,beq,xlow,xupp)
%   [...] = snsolve(myobj,x0,A,b,Aeq,beq,xlow,xupp,nonlcon)
%   [...] = snsolve(myobj,x0,A,b,Aeq,beq,xlow,xupp,nonlcon,options)
%
% Output from snsolve:
%   [x,fval,exitflag,lambda,states,output] = snsolve(...)
%
%
% Output:
%   x                  solution
%
%   fval               final objective value at x
%
%   exitflag           exit condition from SNOPT
%
%   lambda.lower       final multipliers for lower bounds
%   lambda.upper       final multipliers for upper bounds
%   lambda.ineqnonlin  final multipliers for nonlinear inequalities
%   lambda.eqnonlin    final multipliers for nonlinear equalities
%   lambda.ineqlin     final multipliers for linear inequalities
%   lambda.eqlin       final multipliers for linear equalities
%
%   states.x           final state of variables
%   states.F           final state of slack (constraint) variables
%
%   output.info        exit code from SNOPT (same as exitflag)
%   output.iterations  the total number of minor iterations
%   output.majors      the total number of major iterations
%
%
% snsolve and fmincon assume problems are of the form:
%    minimize    f(x)
%   such that    c(x)   <= 0,
%                c_eq(x) = 0,
%                Ax     <= b,
%                A_eq x  = b_eq,
%                xlow   <= x <= xupp.
%

% Check for starting point x0.
x0 = colvec(x0,'x0',0,0);
n  = length(x0);

% Get user-defined functions.
nonlcon =  @dummyCon;
if (ischar(userobj))
  myobj = str2func(userobj);
else
  myobj = userobj;
end
F = myobj(x0);
if length(F) == 0,
  error('Error: userobj must return the objective function.');
end

if     nargin == 4,
  Aeq  = [];  beq  = [];
  xlow = [];  xupp = [];
  c    = [];  ceq  = [];

elseif nargin == 6,
  Aeq  = varargin{1};
  beq  = varargin{2};
  xlow = [];  xupp = [];
  c    = [];  ceq  = [];

elseif nargin == 8,
  Aeq      = varargin{1};
  beq      = varargin{2};
  xlow     = varargin{3};
  xupp     = varargin{4};
  c    = [];  ceq  = [];

elseif nargin == 9 || nargin == 10,
  Aeq      = varargin{1};
  beq      = varargin{2};
  xlow     = varargin{3};
  xupp     = varargin{4};
  nonlconU = varargin{5};

  if (ischar(nonlconU)),
    nonlcon = str2func(nonlconU);
  else
    nonlcon = nonlconU;
  end

  [c,ceq]  = feval(nonlcon,x0);

else
  error('Wrong number of input arguments')
end

% Options?
probName = '';
mlSTOP   = 0;

if nargin == 10,
  if isstruct(varargin{6}),
    options = varargin{6};
    if isfield(options,'name'),
      probName = options.name;
    end

    if isfield(options,'stop'),
      if (ischar(options.stop))
	mlSTOP = str2func(options.stop);
      else
	mlSTOP = options.stop;
      end
    end

  else
    error('Options struct error');
  end
end

% Check inputs
[mi,n0] =  size(A);
if ( ~isempty(A) && n0 ~= n ),
  error('Error: A has incorrect column dimension.');
end

b = colvec(b,'b',1,mi);
if isempty(b) && ~isempty(A),
  error('Error: b is empty, but A is not.');
end
if ~isempty(b) && isempty(A),
  error('Error: b is not empty, but A is.');
end

[me,n0] =  size(Aeq);
if ( ~isempty(Aeq) && n0 ~= n ),
  error('Error: Aeq has incorrect column dimension.');
end

beq = colvec(beq,'beq',1,me);
if isempty(beq) && ~isempty(Aeq),
  error('Error: beq is empty, but Aeq is not.');
end
if ~isempty(beq) && isempty(Aeq),
  error('Error: beq is not empty, but Aeq is.');
end

nli = length(c);
nle = length(ceq);


% Set total number of constraints (size of F(x)).
nCon = 1 + nli + nle + mi + me;

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

iAfun   = colvec(iAfun,'iAfun',1,0);
jAvar   = colvec(jAvar,'jAvar',1,0);
Aij     = colvec(Aij,'Aij',1,0);
if length(Aij) ~= length(iAfun) || ...
      length(Aij) ~= length(jAvar) || ...
      length(iAfun) ~= length(jAvar),
  error('Error: A, iAfun, jAvar must have the same length.');
end

iGfun   = colvec(iGfun,'iGfun',1,0);
jGvar   = colvec(jGvar,'jGvar',1,0);
if length(iGfun) ~= length(jGvar),
  error('Error: iGfun and jGvar must have the same length.');
end

xmul    =  zeros(n,1);
xstate  =  zeros(n,1);
Fmul    =  zeros(nCon,1);
Fstate  =  zeros(nCon,1);
ObjAdd  =  0;
ObjRow  =  1;
Flow    = [ -inf; -inf*ones(nli,1); zeros(nle,1); -inf*ones(mi,1); beq ];
Fupp    = [  inf; zeros(nli,1); zeros(nle,1); b; beq ];


% Solve the problem!
solveopt = 1;

[x,F,exitflag,xmul,Fmul,xstate,Fstate,itn,mjritn] = snoptmex( solveopt, x0, ...
						  xlow, xupp, xmul, xstate, ...
						  Flow, Fupp, Fmul, Fstate, ...
						  ObjAdd, ObjRow, ...
						  Aij, iAfun, jAvar, ...
						  iGfun, jGvar, ...
						  myobj, nonlcon, probName, mlSTOP );

fval              = feval(myobj,x);
lambda.lower      = xmul;
lambda.upper      = xmul;
lambda.ineqnonlin = Fmul(2:nli+1);
lambda.eqnonlin   = Fmul(nli+2:nli+nle+1);
lambda.ineqlin    = Fmul(nli+nle+2:nli+nle+mi+1);
lambda.eqlin      = Fmul(nli+nle+mi+2:nCon);

states.x          = xstate;
states.F          = Fstate;

output.info       = exitflag;
output.iterations = itn;
output.majors     = mjritn;


function [c,ceq] = dummyCon(x)

c = []; ceq = [];
