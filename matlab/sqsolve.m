function [x,fval,exitFlag,output,lambda] = sqsolve(H, f, varargin)
% function [x,fval,exitFlag,output,lambda] = sqsolve(H, f, varargin)
%
% This function interface is similar to the MATLAB function quadprog.
%
% Solve the given quadratic problem:
%       minimize        q(x) = half*x'*H*x + f'*x
%     subject to   lb <=  x  <= ub
%                       A*x  <= b
%                     Aeq*x   = beq
%
% Calling sequences:
%  x = sqsolve(H, f)
%  x = sqsolve(H, f, A, b)
%  x = sqsolve(H, f, A, b, Aeq, beq)
%  x = sqsolve(H, f, A, b, Aeq, beq, lb, ub)
%  x = sqsolve(H, f, A, b, Aeq, beq, lb, ub, x0)
%  x = sqsolve(H, f, A, b, Aeq, beq, lb, ub, x0, options)
%  x = sqsolve(H, f, A, b, Aeq, beq, lb, ub, x0, lambda, states, options)
%
%  [x,fval]                        = sqsolve(H, f, ...)
%  [x,fval,exitflag]               = sqsolve(H, f, ...)
%  [x,fval,exitflag,output]        = sqsolve(H, f, ...)
%  [x,fval,exitflag,output,lambda] = sqsolve(H, f, ...)
%
%
%   INPUT:
%    H         is a Matlab function (either a function handle or string)
%              that computes H*x for a given x or a matrix (dense or sparse).
%              If the problem is an LP (H = 0), then set H = 0 or H = []
%              (or call lpopt).
%
%     f        is the linear term of the objective
%
%     A, b     contain the linear inequality constraints A*x <= b
%
%     Aeq, beq(optional) contain the lineaer equality constraints Aeq*x <= beq
%
%     lb, ub  (optional) are the lower and upper bounds of x
%
%     x0       is the initial point x
%
%  options     is an (optional) input argument of type struct.  SNOPT
%              options can be set using this structure by creating an entry with a
%              field name equal to the SNOPT keyword with spaces replaced by
%              underscores '_'.  For example,
%                 options.iterations_limit = 250;
%
%              Additional keywords include:
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
%               options.iwork       is an integer defining the integer
%                                   SNOPT workspace length.
%
%               options.rwork       is an integer defining the real
%                                   SNOPT workspace length.
%
%   OUTPUT:
%     x        is the final point
%
%     fval     is the final objective value
%
%     exitFlag is the exit flag returned by SQOPT
%
%     output   is a structure containing run information --
%              output.iterations is the total number of iterations
%              output.funcCount   is the total number of function evaluations
%
%     lambda   is a structure containing the multipliers
%              lambda.x          are for the variables
%              lambda.linear     are for the linear constraints
%
%     states   is a structure containing the states
%              states.x          are for the variables
%              states.linear     are for the linear constraints
%

name     = '';
start    = 'Cold';

printfile  = '';
screen     = 'on';
specsfile  = '';

iwork      = 0;
rwork      = 0;

% Deal with options.
optionsLoc = 0;
if nargin == 10 || nargin == 12,
  optionsLoc = nargin - 9;
  if isstruct(varargin{optionsLoc}),
    options = varargin{optionsLoc};
    % Name
    if isfield(options,'name'),
      probName = options.name;
    end

    % Start
    if isfield(options,'start'),
      start = options.start;
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
sqprint(printfile);
sqscreen(screen);
sqsetwork(iwork,rwork);


% Read specsfile
if ~strcmp(specsfile,''),
  mexopt = 9;
  info = sqspecs(specsfile);

  if info ~= 101 && info ~= 107,
    x = []; obj = 0; output = []; lambda = []; states = [];

    end_sqopt();
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
	    ~strcmp(keyword,'start'),

	option = options.(fields{i});

	if (isnumeric(option)),
	  option = num2str(option);
	end
	string = strjoin({keyword, option});

	sqset(string);
      end
    end
  end
end

if isempty(H),
    warning('No Hessian detected: the problem is an LP');
    userHx = 0;
else
  if isnumeric(H),
    if H == 0,
      warning('No Hessian detected: the problem is an LP');
      userHx = 0;
    else
      userHx = @(x)myHx(H,x);
    end
  else
    userHx = checkFun(H,'SQOPT','H');
  end
end


if nargin == 2,
  % sqsolve(H, f)
  A   = [];  b   = [];
  Aeq = [];  beq = [];
  lb  = [];  ub  = [];
  x0  = [];

  xstate = []; xmul = [];
  astate = []; amul = [];

elseif nargin == 4,
  % sqsolve(H, f, A, b)
  A = varargin{1};
  b = varargin{2};
  Aeq = [];  beq = [];
  lb  = [];  ub  = [];
  x0  = [];

  xstate = []; xmul = [];
  astate = []; amul = [];


elseif nargin == 6,
  % sqsolve(H, f, A, b, Aeq, beq)
  A   = varargin{1};
  b   = varargin{2};
  Aeq = varargin{3};
  beq = varargin{4};
  lb  = [];  ub  = [];
  x0  = [];

  xstate = []; xmul = [];
  astate = []; amul = [];

elseif nargin == 8,
  % sqsolve(H, f, A, b, Aeq, beq, lb, ub)
  A   = varargin{1};
  b   = varargin{2};
  Aeq = varargin{3};
  beq = varargin{4};
  lb  = varargin{5};
  ub  = varargin{6};
  x0  = [];

  xstate = []; xmul = [];
  astate = []; amul = [];

elseif nargin == 9,
  % sqsolve(H, f, A, b, Aeq, beq, lb, ub, x0)
  A   = varargin{1};
  b   = varargin{2};
  Aeq = varargin{3};
  beq = varargin{4};
  lb  = varargin{5};
  ub  = varargin{6};
  x0  = varargin{7};

  xstate = []; xmul = [];
  astate = []; amul = [];

elseif nargin == 10 || nargin == 12,
  % sqsolve(H, f, A, b, Aeq, beq, lb, ub, x0, options)
  % sqsolve(H, f, A, b, Aeq, beq, lb, ub, x0, lambda, states, options)
  A   = varargin{1};
  b   = varargin{2};
  Aeq = varargin{3};
  beq = varargin{4};
  lb  = varargin{5};
  ub  = varargin{6};
  x0  = varargin{7};

  xstate = []; xmul = [];
  astate = []; amul = [];

  if nargin == 12,
    lambda = varargin{8};
    states = varargin{9};

    xstate = states.x;
    xmul   = lambda.x;
    astate = states.linear;
    amul   = lambda.linear;
  end

else
  error('SQOPT:InputArgs','Wrong number of input arguments for sqsolve');
end

ineq = size(A,1);
AA   = [                 A; Aeq ];
al   = [ -inf*ones(ineq,1); beq ];
au   = [                 b; beq ];

m    = size(AA,1);
n    = size(x0,1);

[x,fval,exitFlag,output,lambda] = solve_sqopt(start,name, m, n, ...
					      userHx, f, x0, lb, ub, xstate, xmul, ...
					      AA,  al, au, astate, amul);

end_sqopt();


function [Hx] = myHx(H,x)
  Hx = H*x;