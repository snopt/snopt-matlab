function [x,fval,info,output,lambda,states] = snsolve(obj,x0,A,b,varargin)
% function [x,fval,info,output,lambda,states] = snsolve(obj,x0,A,b,varargin)
%
%   [...] = snsolve(obj,x0,A,b)
%   [...] = snsolve(obj,x0,A,b,options)
%
%   [...] = snsolve(obj,x0,A,b,Aeq,beq)
%   [...] = snsolve(obj,x0,A,b,Aeq,beq,options)
%
%   [...] = snsolve(obj,x0,A,b,Aeq,beq,xlow,xupp)
%   [...] = snsolve(obj,x0,A,b,Aeq,beq,xlow,xupp,options)
%
%   [...] = snsolve(obj,x0,A,b,Aeq,beq,xlow,xupp,nonlcon)
%   [...] = snsolve(obj,x0,A,b,Aeq,beq,xlow,xupp,nonlcon,options)
%
%   [...] = snsolve(obj,x0,A,b,Aeq,beq,xlow,xupp,nonlcon,lambda,states)
%   [...] = snsolve(obj,x0,A,b,Aeq,beq,xlow,xupp,nonlcon,lambda,states,options)
%
%
% Output from snsolve:
%   [x,fval,info,output,lambda,states] = snsolve(...)
%
%
% snsolve and fmincon assume problems are of the form:
%    minimize        f(x)
%   such that    c(x)   <= 0,
%                c_eq(x) = 0,
%                Ax     <= b,
%                A_eq x  = b_eq,
%                xlow   <= x <= xupp.
%
% Input:
% options       is an (optional) input argument of type struct.  SNOPT
%               options can be set using this structure by creating an entry with a
%               field name equal to the SNOPT keyword with spaces replaced by
%               underscores '_'.  For example,
%                      options.iterations_limit = 250;
%
%               Additional keywords include:
%
%               options.name        is the problem name
%
%               options.start       'Cold', 'Warm'
%
%               options.screen      is a string set to 'on' or 'off'.
%                                   Summary to the screen is controlled
%                                   by this option. (default 'on')
%
%               options.printfile   is a string denoting the print file.
%                                   By default, no print file is created.
%                                   Not setting this option or setting it to
%                                   '' turns off print output.
%
%               options.specsfile   is a string denoting the options
%                                   filename.
%
%               options.stop        is the "snSTOP" function called at every
%                                   major iteration.
%
%               options.iwork       is an integer defining the integer
%                                   SNOPT workspace length.
%
%               options.rwork       is an integer defining the real
%                                   SNOPT workspace length.
%
% Output:
%   x            solution
%
%   fval         final objective value at x
%
%   info         exit condition from SNOPT
%
%   output       is a struct with the following fields
%                output.info        exit code from SNOPT (same as info)
%                output.iterations  the total number of minor iterations
%                output.majors      the total number of major iterations
%
%   lambda       are the final multipliers
%                lambda.x            variables
%                lambda.ineqnonlin   nonlinear inequalities
%                lambda.eqnonlin     nonlinear equalities
%                lambda.ineqlin      linear inequalities
%                lambda.eqlin        linear equalities
%
%   states       are the final states
%                states.x            variables
%                states.ineqnonlin   nonlinear inequalities
%                states.eqnonlin     nonlinear equalities
%                states.ineqlin      linear inequalities
%                states.eqlin        linear equalities
%
%

name       = '';
istart     = 0;

printfile  = '';
screen     = 'on';
specsfile  = '';

iwork      = 0;
rwork      = 0;

stopFun    = 0;
optionsLoc = 0;


% Deal with options
if nargin == 5 || nargin == 7 || nargin == 9 || ...
	   nargin == 10 || nargin == 12,
  optionsLoc = nargin - 4;
  if isstruct(varargin{optionsLoc}),
    options = varargin{optionsLoc};

    % Name
    if isfield(options,'name'),
      name = options.name;
    end

    % Start
    if isfield(options,'start'),
      if strcmp(lower(options.start),'warm'),
	istart = 1;
      elseif strcmp(lower(options.start),'hot'),
	istart = 2;
      end
    end

    % Print output
    if isfield(options,'printfile'),
      if ischar(options.printfile),
	printfile = options.printfile;
      end
    end

    % Specs file
    if isfield(options,'specsfile'),
      if ischar(options.specsfile),
	specsfile = options.specsfile;
      end
    end

    % Screen
    if isfield(options,'screen'),
      if ischar(options.screen),
	screen = options.screen;
      end
    end

    % Stop function
    if isfield(options,'stop'),
      if ischar(options.stop),
	stopFun = str2func(options.stop);
      elseif isa(options.stop,'function_handle'),
	stopFun = options.stop;
      else
	error('SNOPT:InputArgs','options.stop should be a string or function handle');
      end
    end

    % iwork
    if isfield(options,'iwork'),
      if ischar(options.iwork),
	iwork = options.iwork;
      end
    end

    % rwork
    if isfield(options,'rwork'),
      if ischar(options.rwork),
	rwork = options.rwork;
      end
    end

  else
    optionsLoc = 0;
  end
end

% Set print, screen, workspace FIRST.
snprint(printfile);
snscreen(screen);
snsetwork(iwork,rwork);


% Read specsfile
if ~strcmp(specsfile,''),
  mexopt = 9;
  info = snspec(specsfile);

  if info ~= 101 && info ~= 107,
    x = []; xmul = []; xstate = [];
    F = []; Fmul = []; Fstate = [];
    output = [];

    end_snopt();
    return;
  end
end


% Handle other options
if (optionsLoc ~= 0),
  fields = fieldnames(options);
  for i = 1:numel(fields),
    if (ischar(fields{i})),
      keyword = strrep(fields{i}, '_', ' ');

      if ~strcmp(keyword,'screen') && ...
	    ~strcmp(keyword,'printfile') && ...
	    ~strcmp(keyword,'specsfile') && ...
	    ~strcmp(keyword,'name') && ...
	    ~strcmp(keyword,'iwork') && ...
	    ~strcmp(keyword,'rwork') && ...
	    ~strcmp(keyword,'stop') && ...
	    ~strcmp(keyword,'start'),

	option = options.(fields{i});

	if (isnumeric(option)),
	  option = num2str(option);
	end
	string = strjoin({keyword, option});

	snset(string);
      end
    end
  end
end


linear_ineq = size(A,1);
linear_eq   = 0;
nonlin_ineq = 0;
nonlin_eq   = 0;
nonlcon     = 0;
x0          = colvec(x0,'x0',0,0);
n           = length(x0);


% Check user-defined functions
myobj = checkFun(obj,'SNOPT','obj');

gotGrad = 0;
try
  [fobj,gobj] = myobj(x0);
  gotGrad = 1;
catch
  try
    fobj = myobj(x0);
    gotGrad = 0;
  catch
    error('SNOPT:InputArgs', ...
	  'Wrong number of output arguments for obj');
  end
end


if     nargin == 4 || nargin == 5,
  % snsolve(obj,x0,A,b)
  % snsolve(obj,x0,A,b,options)

  Aeq  = [];  beq  = [];
  xlow = [];  xupp = [];
  c    = [];  ceq  = [];
  xmul = [];  xstate = [];
  Fmul = [];  Fstate = [];


elseif nargin == 6 || nargin == 7,
  % snsolve(obj,x0,A,b,Aeq,beq)
  % snsolve(obj,x0,A,b,Aeq,beq,options)

  Aeq       = varargin{1};
  beq       = varargin{2};
  linear_eq = size(Aeq,1);

  xlow = [];  xupp = [];
  c    = [];  ceq  = [];
  xmul = [];  xstate = [];
  Fmul = [];  Fstate = [];

elseif nargin == 8 || (nargin == 9 && optionsLoc ~=0),
  % snsolve(obj,x0,A,b,Aeq,beq,xlow,xupp)
  % snsolve(obj,x0,A,b,Aeq,beq,xlow,xupp,options)

  Aeq       = varargin{1};
  beq       = varargin{2};
  xlow      = varargin{3};
  xupp      = varargin{4};
  linear_eq = size(Aeq,1);
  c    = [];  ceq  = [];
  xmul = [];  xstate = [];
  Fmul = [];  Fstate = [];

elseif nargin >= 9 && nargin <= 12,
  % snsolve(obj,x0,A,b,Aeq,beq,xlow,xupp,nonlcon)
  % snsolve(obj,x0,A,b,Aeq,beq,xlow,xupp,nonlcon,options)
  % snsolve(obj,x0,A,b,Aeq,beq,xlow,xupp,nonlcon,lambda,states)
  % snsolve(obj,x0,A,b,Aeq,beq,xlow,xupp,nonlcon,lambda,states,options)

  Aeq       = varargin{1};
  beq       = varargin{2};
  xlow      = varargin{3};
  xupp      = varargin{4};
  nonlc     = varargin{5};
  linear_eq = size(Aeq,1);
  xmul = [];  xstate = [];
  Fmul = [];  Fstate = [];

  if nargin == 11 || nargin == 12,
    lambda = varargin{6};
    states = varargin{7};

    xmul   = lambda.x;  xstate = states.x;
    Fmul   = [ 0; lambda.ineqnonlin; lambda.eqnonlin;
	       lambda.ineqlin; lambda.eqlin];
    Fstate = [ 0; states.ineqnonlin; states.eqnonlin;
	       states.ineqlin; states.eqlin];

  end

  nonlcon  = checkFun(nonlc,'SNOPT','nonlcon');

  gotDeriv = 0;
  try
    [c,ceq,J,Jeq] = nonlcon(x0);
    J = J';  Jeq = Jeq';
    gotDeriv = 1;
  catch
    try
      [c,ceq]  = nonlcon(x0);
      gotDeriv = 0;
    catch
      error('SNOPT:InputArgs', ...
	    'Wrong number of output arguments for nonlcon');
    end
  end
  nonlin_ineq = size(c,1);
  nonlin_eq   = size(ceq,1);

else
  error('SNOPT:InputArgs','Wrong number of input arguments for snsolve')
end


% Set total number of constraints (size of F(x)).
nCon = 1 + nonlin_ineq + nonlin_eq + linear_ineq + linear_eq;

% snoptA problem format:
%    minimize    F_obj (x)
%   such that  l_f <= F(x) <= u_f
%              l   <=   x  <= u
%
%           [ F_obj   ]            [ F_0'     ]   [  0    ]     1
%           [ c(x)    ]            [ c'(x)    ]   [  0    ]     nonlin_ineq
%           [ c_eq(x) ]            [ c_eq'(x) ]   [  0    ]     nonlin_eq
%   F(x) =  [ Ax      ]    F'(x) = [ 0        ] + [  A    ]     linear_ineq
%           [ A_eq x  ]            [ 0        ]   [  A_eq ]     linear_eq
%                                    "G(x)"         "A"

[iAfun,jAvar,Aij] = find([zeros(1+nonlin_ineq+nonlin_eq,n); A; Aeq]);
[iGfun,jGvar,~]   = find([ones(1+nonlin_ineq+nonlin_eq,n); zeros(linear_ineq+linear_eq,n)]);

iAfun   = colvec(iAfun,'iAfun',1,0);
jAvar   = colvec(jAvar,'jAvar',1,0);
Aij     = colvec(Aij,'Aij',1,0);
if length(Aij) ~= length(iAfun) || ...
      length(Aij) ~= length(jAvar) || ...
      length(iAfun) ~= length(jAvar),
  error('SNOPT:InputArgs','A, iAfun, jAvar must have the same length.');
end

iGfun   = colvec(iGfun,'iGfun',1,0);
jGvar   = colvec(jGvar,'jGvar',1,0);
if length(iGfun) ~= length(jGvar),
  error('SNOPT:InputArgs','iGfun and jGvar must have the same length.');
end

ObjAdd  =  0;
ObjRow  =  1;
Flow    = [ -inf;
	    -inf*ones(nonlin_ineq,1);
	    zeros(nonlin_eq,1);
	    -inf*ones(linear_ineq,1);
	    beq ];
Fupp    = [  inf;
	     zeros(nonlin_ineq,1);
	     zeros(nonlin_eq,1);
	     b;
	     beq ];

if isa(nonlcon,'function_handle'),
  [x,fval,info,output,lambda,states] = ...
      solve_snfmincon(istart, stopFun, name, ...
		      @(x,needF,needG)snfun(x,needF,needG,myobj,gotGrad, ...
					    nonlcon,gotDeriv,iGfun,jGvar), ...
		      x0, ...
		      xlow, xupp, xmul, xstate, ...
		      Flow, Fupp, Fmul, Fstate, ...
		      ObjAdd, ObjRow, ...
		      Aij, iAfun, jAvar, iGfun, jGvar, ...
		      n, nonlin_ineq, nonlin_eq, linear_ineq, linear_eq);
else
  [x,fval,info,output,lambda,states] = ...
      solve_snfmincon(istart, stopFun, name, ...
		      @(x,needF,needG)snfun(x,needF,needG,myobj,gotGrad), ...
		      x0, ...
		      xlow, xupp, xmul, xstate, ...
		      Flow, Fupp, Fmul, Fstate, ...
		      ObjAdd, ObjRow, ...
		      Aij, iAfun, jAvar, iGfun, jGvar, ...
  		      n, nonlin_ineq, nonlin_eq, linear_ineq, linear_eq);
end

% End
end_snopt();



function [F,G] = snfun(x,needF,needG,obj,gotGrad,varargin)
% Wrapper for obj and nonlcon in snsolve call.

% Compute objective function and gradients
fobj = []; gobj = [];

if needG > 0,
  if gotGrad,
    [fobj,gobj] = obj(x);
  else
    if needF > 0,
      fobj = obj(x);
    end
  end
else
  fobj = obj(x);
end

% Compute constraint functions and gradients
c    = []; ceq  = [];
J    = []; Jeq  = [];

if nargin == 9,
  nonlcon  = varargin{1};
  gotDeriv = varargin{2};
  iGfun    = varargin{3};
  jGvar    = varargin{4};

  if needG > 0,
    if gotDeriv,
      [c,ceq,J,Jeq] = nonlcon(x);
    else
      if needF > 0,
	[c,ceq] = nonlcon(x);
      end
    end
  else
    [c,ceq] = nonlcon(x);
  end
end

F = [  fobj; c; ceq ];
G = [ gobj'; J; Jeq ];

% Convert G to vector format to match SNOPTA and (iGfun,jGvar)
[~,n] = size(G);
if n > 1,
  G = snfindG(iGfun,jGvar,G);
end
