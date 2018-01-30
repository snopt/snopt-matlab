function sqsetwork(iwork,rwork)
% function sqsetwork(leniw, lenrw)
%     Modify the initial amount of workspace for SQOPT.
%     Values must be at least 500.
%

if iwork > 0 && rwork > 0,
  mexopt = 14;
  sqoptmex(mexopt,iwork,rwork);
end
