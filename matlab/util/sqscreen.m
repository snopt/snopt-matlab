function sqscreen(screen)
%function sqscreen(screen)

if strcmp(screen,'on'),
  mexopt = 15;
  sqoptmex(mexopt);
elseif strcmp(screen,'off')
  mexopt = 16;
  sqoptmex(mexopt);
end
