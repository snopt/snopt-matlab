function info = sqspec(specsfile)
% function info = sqspec(specsfile)
%     Causes sqopt to read its optional parameters from the named file.
%     The format of this file is described in the snopt documentation.
%
%     Returns 101 or 107 if successful.

if ~strcmp(specsfile,''),
  mexopt = 9;
  info = sqoptmex(mexopt,specsfile);
end