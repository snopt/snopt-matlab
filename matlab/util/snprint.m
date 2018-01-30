function snprint(printfile)
% function snprint(printfile)

if ~strcmp(printfile,''),
  mexopt = 10;
  snoptmex(mexopt,printfile);
end
