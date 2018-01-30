function [x,F,info,xmul,Fmul,xstate,Fstate,output] = snopt(x, xlow, xupp, xmul, xstate,...
						  Flow, Fupp, Fmul, Fstate,...
						  userfun, varargin);
% This function solves the nonlinear optimization problem:
% minimize:
%            F(ObjRow) + ObjAdd
% subject to:
%            xlow <= x <= xupp
%            Flow <= F <= Fupp
% where:
%  x        is the column vector of initial values of the unknowns.
%
%  F        is a vector of objective and constraint functions specified
%           in the m-file userfun.
%
%  ObjRow   is the objective row of F (default ObjRow = 1).
%
%  userfun  is the handle for a Matlab funcion that defines the elements of F
%           and optionally, their derivatives.
%
% Additional arguments allow the specification of more detailed problem information.
%
% Calling sequence 1:
%  [...] = snopt(x, xlow, xupp, xmul, xstate,
%                     Flow, Fupp, Fmul, Fstate, userfun,
%                     [options])
%
% Calling sequence 2:
%  [...] = snopt(x, xlow, xupp, xmul, xstate,...
%                     Flow, Fupp, Fmul, Fstate, userfun,...
%                     ObjAdd, ObjRow, [options])
%
% Calling sequence 3:
%  [...] = snopt(x, xlow, xupp, xmul, xstate,...
%                     Flow, Fupp, Fmul, Fstate, userfun,...
%                     ObjAdd, ObjRow,
%                     A, G, [options])
%
% Output from snopt:
%  [x,F,info,xmul,Fmul,xstate,Fstate,output] = snopt(...)
%
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
%                      >> [f,G] = userfun(x);
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
%                that the proper derivative level is set
%                     >> option.derivative_option = k;
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
% A              is either a struct type or a matrix (dense or sparse)
%                that defines the constant elements in the Jacobian of F.
%
%                If A is a struct, the structure is provided in
%                coordinate form with fields:
%                   A.row
%                   A.col
%                   A.val
%                If i = A.row(k) and j = A.col(k), then A.val(k) is the (i,j)-th
%                element in A.
%
% G              is a struct type or a matrix (dense or sparse) defining
%                the nonlinear elements in the Jacobian of F.
%
%                If G is a struct, the Jacobian structure is provided in
%                coordinate form with fields:
%                   G.row
%                   G.col
%                If G is a matrix, then nonzero elements in the matrix
%                denote the nonzero elements of the nonlinear elements of
%                the Jacobian.
%
%            More IMPORTANT details:
%            1) The indices (A.row,A.col) must be DISJOINT from (G.row,G.col).
%               A nonzero element in F' must be either an element of G or an
%               element of A, but not the sum of the two.
%
%            2) If the user does not wish to provide A.row, A.col, G.row,
%               G.col, then snopt() will determine them by calling snJac().
%
%               WARNING: In this case, the derivative level will be set to zero
%               if constant elements exist.  This is because the linear
%               elements have not yet been deleted from the definition of
%               userfun.  Furthermore, if G is given in vector form, the
%               ordering of G may not necessarily correspond to (iGfun,jGvar)
%               computed by snJac().
%
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
name       = '';
istart     = 0;

printfile  = '';
screen     = 'on';
specsfile  = '';

iwork      = 0;
rwork      = 0;

stopFun    = 0;
optionsLoc = 0;

% Deal with options first.
if nargin == 11 || nargin == 13 || nargin == 15,
  optionsLoc = nargin - 10;
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


% Check userfun
userFG = checkFun(userfun,'SNOPT','userfun');

gotDeriv = 0;
try
  [F,G] = userFG(x);
  gotDeriv = 1;
catch
  try
    F = userFG(x);
    gotDeriv = 0;
  catch
    error('SNOPT:InputArgs','userfun should return 1 or 2 arguments');
  end
end


ObjAdd  = 0;
ObjRow  = 1;

callJac = 0;

if nargin == 10 || nargin == 11,
  %  snopt(x, xlow, xupp, xmul, xstate,
  %           Flow, Fupp, Fmul, Fstate, userfun,
  %           [options])
  % User is not providing derivative structures.

  warning('SNOPT:Input','User is not providing SNOPT derivatives structures');

  F0 = userFG(x);
  nF = length(F0);
  n  = length(x);

  callJac = 1;
  if ~gotDeriv,
    % User is also not providing derivatives.
    % Call snJac to estimate the pattern of nonzeros for the Jacobian.
    warning('SNOPT:Input','Derivative structures estimated via snJac');
  else
    % User IS providing derivatives via userfun.
    warning('SNOPT:Input',['Derivatives provided but not structures: estimating' ...
		    ' structure via snJac.']);
  end

elseif nargin == 12 || nargin == 13,
  % snopt(x, xlow, xupp, xmul, xstate,...
  %          Flow, Fupp, Fmul, Fstate, userfun,...
  %       ObjAdd, ObjRow, [options])
  % User is not providing derivative structures.

  warning('SNOPT:Input','User is not providing SNOPT derivatives structures');

  F0 = userFG(x);
  nF = length(F0);
  n  = length(x);

  callJac = 1;
  if ~gotDeriv,
    % Call snJac to estimate the pattern of nonzeros for the Jacobian.
    warning('SNOPT:Input','Derivative structures estimated via snJac');
  else
    % User IS providing derivatives via userfun.
    warning('SNOPT:Input',['Derivatives provided but not structures: estimating' ...
		    ' structure via snJac.']);
  end

  ObjAdd = varargin{1};
  ObjRow = varargin{2};

elseif nargin == 14 || nargin == 15,
  % snopt(x, xlow, xupp, xmul, xstate,...
  %          Flow, Fupp, Fmul, Fstate, userfun,...
  %       ObjAdd, ObjRow,
  %       A, G, [options])
  % The user is providing derivatives.

  ObjAdd = varargin{1};
  ObjRow = varargin{2};
  A      = varargin{3};
  G      = varargin{4};

  if isstruct(A),
    if isfield(A,'row') && isfield(A,'col') && isfield(A,'val'),
      % In coordinate form
      iAfun = colvec(A.row,'A.row',1,0);
      jAvar = colvec(A.col,'A.col',1,0);
      valA  = colvec(A.val,'A.val',1,0);
    else
      error('SNOPT:InputArgs','Matrix must have row, col and val fields')
    end

  elseif isnumeric(A),
    % Dense or sparse
    [iAfun,jAvar,valA] = find(A);
  else
    error('SNOPT:InputArgs','Wrong input type for A')
  end

  if isstruct(G),
    % In coordinate form
    if isfield(G,'row') && isfield(G,'col'),
      % In coordinate form
      iGfun = colvec(G.row,'G.row',1,0);
      jGvar = colvec(G.col,'G.col',1,0);
    else
      error('SNOPT:InputArgs','Matrix must have row and col fields')
    end

  elseif isnumeric(G),
    % Dense or sparse
    [iGfun,jGvar] = find(G);
  else
    error('SNOPT:InputArgs','Wrong input type for G')
  end

else
  error('SNOPT:InputArgs','Wrong number of input arguments for SNOPT')
end

% Call snJac
if callJac > 0,
  [valA,iAfun,jAvar,iGfun,jGvar] = snjac(userFG,x,xlow,xupp,nF);
end

[x,F,info,xmul,Fmul, ...
 xstate,Fstate,output] = solve_snopt(istart, stopFun, name, ...
				     @(x,needF,needG)snfun(x,needF,needG,...
						  userFG,gotDeriv,iGfun,jGvar), ...
				     x, ...
				     xlow, xupp, xmul, xstate, ...
				     Flow, Fupp, Fmul, Fstate, ...
				     ObjAdd, ObjRow, ...
				     valA, iAfun, jAvar, ...
				     iGfun, jGvar);
end_snopt();


function [F,G] = snfun(x,needF,needG,userfun,gotDeriv,iGfun,jGvar)
% Wrapper for userfun
% Compute functions and gradients (if necessary)

F = []; G = [];

if needG > 0
  if gotDeriv,
    [F,G] = userfun(x);

    % Convert G to vector format to match SNOPTA and (iGfun,jGvar)
    [~,n] = size(G);
    if n > 1,
      G = snfindG(iGfun,jGvar,G);
    end
  else
    if needF > 0,
      F = userfun(x);
    end
  end

else
  if needF > 0,
    F = userfun(x);
  end
end
