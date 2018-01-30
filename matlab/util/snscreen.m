function snscreen(screen)
%function snscreen(screen)

if strcmp(screen,'on'),
  mexopt = 15;
  snoptmex(mexopt);
elseif strcmp(screen,'off')
  mexopt = 16;
  snoptmex(mexopt);
end
