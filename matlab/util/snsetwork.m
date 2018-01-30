function snsetwork(iwork,rwork)
% function snsetwork(leniw, lenrw)
%     Modify the initial amount of workspace for SNOPT.
%     Values must be at least 500.
%

if iwork > 0 && rwork > 0,
  mexopt = 14;
  snoptmex(mexopt,iwork,rwork);
end
