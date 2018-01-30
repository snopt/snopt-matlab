function sqprint(printfile)
% function sqprint(printfile)

if ~strcmp(printfile,''),
  mexopt = 10;
  sqoptmex(mexopt,printfile);
end
