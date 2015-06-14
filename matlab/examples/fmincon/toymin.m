function [x,F,INFO]=toymin
% Mimics sntoyA.f in $SNOPT/examples
% Example of the 'fmincon'-style call to SNOPT.
%
%     Minimize      3*x(1) + (x(1) + x(2) + x(3))^2 + 5*x(4)
%
%     subject to             4*x(2)   + 2*x(3)               >= 0
%                     x(1) +   x(2)^2 +   x(3)^2              = 2
%                              x(2)^4 +   x(3)^4   +   x(4)   = 4
%
%                     x(1) >= 0,                       x(4) >= 0.
%
%

snscreen on;
snprint('toymin.out');  % By default, screen output is off;

sntoy.spc = which('sntoy.spc');
snspec (sntoy.spc);

snseti ('Major Iteration limit', 250);

x0     = ones(4,1);
A      = [ 0 -4 -2 0];
b      = [ 0 ];

Aeq    = [];
beq    = [];

lb     = [  0,-Inf,-Inf,   0]';
ub     = Inf*ones(4,1);


[x,F,INFO,lambda] = snsolve( @toyObj, x0, A, b, Aeq, beq, lb, ub, @toyCon);

snprint off;
snend;

