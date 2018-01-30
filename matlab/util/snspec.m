function info = snspec(specsfile)
% function info = snspec(specsfile)
%     Causes snopt to read its optional parameters from the named file.
%     The format of this file is described in the snopt documentation.
%
%     Returns 101 or 107 if successful.

if ~strcmp(specsfile,''),
  mexopt = 9;
  info = snoptmex(mexopt,specsfile);
end