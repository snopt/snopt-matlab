function [x,obj,info,output,lambda,states] = sqopt(H, c, x0, xl, xu, A, al, au, varargin)
% function [x,obj,info,output,lambda,states] = sqopt(H, c, x0, xl, xu, A, al, au, varargin)
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
%  [] = sqopt(H, c, x0, xl, xu, A, al, au)
%  [] = sqopt(H, c, x0, xl, xu, A, al, au, options)
%
%  [] = sqopt(H, c, x0, xl, xu, A, al, au, states, lambda)
%  [] = sqopt(H, c, x0, xl, xu, A, al, au, states, lambda, options)
%
%  [x,obj,info,output,lambda,states] = sqopt(...)
%
%
% INPUT:
%  x0       is the initial guess for x
%
%  H        is a Matlab function (either a function handle or string)
%           that computes H*x for a given x or a matrix (dense or sparse).
%           If the problem is an LP (H = 0), then set H = 0 or H = []
%           (or call lpopt).
%
%  c        is the linear term of the quadratic objective
%
%  xl, xu   are the upper and lower bounds on x
%
%  A        is the linear constraint matrix. A can be a structure, or a
%           dense or sparse matrix.
%           If A is a structure, then A is represented as a
%           sparse-by-column matrix and should  have fields:
%               A.loc -- column pointers
%               A.ind -- row indices
%               A.val -- matrix values
%
%  al, au   are the upper and lower bounds on the linear constraints A*x
%
%  options  is an (optional) input argument of type struct.  SNOPT
%           options can be set using this structure by creating an entry with a
%           field name equal to the SNOPT keyword with spaces replaced by
%           underscores '_'.  For example,
%              options.iterations_limit = 250;
%
%           Additional keywords include:
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
%
% OUTPUT:
%  x        is the final point
%
%  obj      is the final objective value
%
%  info     is the exit flag returned by the solver
%
%  output   is a structure containing run information --
%           output.iterations is the total number of iterations
%           output.funcCount   is the total number of function evaluations
%
%  lambda   is a structure containing the multipliers
%           lambda.x          are for the variables
%           lambda.linear     are for the linear constraints
%
%  states   is a structure
%           states.x          are for the variables
%           states.linear     are for the linear constraints
%

name       = '';
start      = 'Cold';

printfile  = '';
screen     = 'on';
specsfile  = '';

iwork      = 0;
rwork      = 0;

% Deal with options.
optionsLoc = 0;
if nargin == 9 || nargin == 11,
  optionsLoc = nargin - 8;
  if isstruct(varargin{optionsLoc}),
    options = varargin{optionsLoc};
    % Name
    if isfield(options,'name'),
      name = options.name;
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
  info = sqspec(specsfile);

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
    userHx = checkFun(H,'SQOPT','Hx');
  end
end


if nargin == 8 || nargin == 9,
  % sqopt(H, c, x0, xl, xu, A, al, au)
  % sqopt(H, c, x0, xl, xu, A, al, au, options)

  xstate = []; xmul = [];
  astate = []; amul = [];

elseif nargin == 10 || nargin == 11,
  % sqopt(H, c, x0, xl, xu, A, al, au, states, lambda)
  % sqopt(H, c, x0, xl, xu, A, al, au, states, lambda, options)

  states = varargin{1};
  lambda = varargin{2};

  xstate = []; xmul   = [];
  astate = []; amul   = [];

  if isfield(states,'x'),
    xstate = states.x;
  end

  if isfield(lambda,'x'),
    xmul = lambda.x;
  end

  if isfield(states,'linear'),
    astate = states.linear;
  end

  if isfield(lambda,'linear'),
    amul = lambda.linear;
  end

else
  error('SQOPT:InputArgs','Wrong number of arguments in SQOPT');
end

if isempty(A),
  % Setup fake constraint matrix and bounds
  warning('SQOPT:InputArgs','No linear constraints detected; dummy constraint created');

  m = 1;
  n = numel(x0);

  neA     = 1;
  indA(1) = 1;
  valA(1) = 1.0;

  locA    = zeros(n+1,1);
  locA(1) = 1;
  locA(2:n+1) = 2;
  al = [-inf]; au = [inf];

else
  if isstruct(A),
    if isfield(A,'ind') && isfield(A,'loc') && isfield(A,'val'),
      % In sparse-by-col form
      n    = numel(x0);
      locA = colvec(A.loc,'A.loc',1,n+1);
      indA = colvec(A.ind,'A.ind',1,0);
      valA = colvec(A.val,'A.val',1,0);
      m    = max(indA);
      neA  = numel(valA);
    else
      error('SQOPT:InputArgs','Matrix must have ind, loc, and val fields')
    end

  else
    [m,n]                = size(A);
    [neA,indA,locA,valA] = crd2spr(A);
  end
end

x0  = colvec(x0,'x0',1,n);
xl  = colvec(xl,'xl',1,n);
xu  = colvec(xu,'xu',1,n);
al  = colvec(al,'al',1,m);
au  = colvec(au,'au',1,m);
c   = colvec(c,'c',1,0);

[x,obj,info,output,lambda,states] = solve_sqopt(start, name, ...
						m, n, userHx, c, ...
						x0, xl, xu, xstate, xmul, ...
						neA, indA, locA, valA, ...
						al, au, astate, amul);

end_sqopt();


function [Hx] = myHx(H,x)
  Hx = H*x;