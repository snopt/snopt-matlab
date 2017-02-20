function [myfun] = checkFun(fun,solver,varargin)
% Check function is a string or function handle.
% Optionally, check function has correct number of output arguments.
% Return the function handle.
%
% myfun = checkFun(fun,solver)
% myfun = checkFun(fun,solver,nouts)
%
errID = [solver ':InputArgs'];

if ischar(fun),
  myfun = str2func(fun);
elseif isa(fun,'function_handle'),
  myfun = fun;
else
  error(errID,'%s should be a function handle or string',inputname(fun));
end

if nargin == 3,
  nouts = varargin{1};
  if ~any(nargout(fun),nouts),
    error(errID,'%s must return the objective function.',inputname(fun));
  end
end
