function [x,fval,exitFlag,output,lambda] = solve_sqprog(start, name, userHx, ...
						  f, x0, lb, ub, xstate, ...
						  xmul, A, al, au, astate, ...
						  amul)
%function [x,fval,exitFlag,output,lambda] = solve_sqprog(start, name, userHx, ...
%						  f, x0, lb, ub, xstate, ...
%						  xmul, A, al, au, astate, ...
%						  amul)

mexopt = 1;
[x,fval,exitFlag,itn,y,state] = sqoptmex(mexopt, start, name, ...
					 userHx, f, x0, lb, ub, xstate, xmul, ...
					 A,  al, au, astate, amul);

% Set output
output.iterations = itn;

m    = size( A,1);
n    = size(x0,1);
zero = zeros(n,1);

states.x      = state(1:n);
lambda.x      = y(1:n);
if m > 0,
  states.linear = state(n+1:n+m);
  lambda.linear = y(n+1:n+m);
end
