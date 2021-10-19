 function [x,F,INFO] = snoptmain2()
%function [x,F,INFO] = snoptmain2()
% Defines the NLP problem and calls the mex interface for snopt.
% First derivatives are provided.

options.name = 'snmain2';
options.screen = 'on';

options.specsfile = which('snoptmain2.spc');
options.printfile = 'snoptmain2.out';

%Get condensed data for the Hexagon problem.
[x,xlow,xupp,xmul,xstate, ...
   Flow,Fupp,Fmul,Fstate, ...
 A.val, A.row, A.col, G.row, G.col] = hexagon;

options.maximize = '';

ObjAdd = 0; ObjRow = 1;
[x,F,INFO] = snopt(x,xlow,xupp,xmul,xstate, ...
		   Flow,Fupp,Fmul,Fstate, ...
		   @snoptuserfun2, ...
		   ObjAdd, ObjRow, ...
		   A,G,options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [x,xlow,xupp,xmul,xstate, ...
	   Flow,Fupp,Fmul,Fstate, ...
	   A,iAfun,jAvar,iGfun,jGvar] = hexagon()
%function [x,xlow,xupp,Flow,Fupp,A,iAfun,jAvar,iGfun,jGvar] = hexagon()
%
% Defines the problem hexagon:
%   maximize F(1)  (the objective row)
%   subject to
%            xlow <=   x  <= xupp
%            Flow <= F(x) <= Fupp
%   where
%     F( 1)  =  x_2 x_6 - x_1 x_7 + x_3 x_7 + x_5 x_8 - x_4 x_9 - x_3 x_8
%     F( 2)  =    x_1^2 + x_6^2
%     F( 3)  =   (x_2   - x_1)^2  +  (x_7 - x_6)^2
%     F( 4)  =   (x_3   - x_1)^2  +   x_6^2
%     F( 5)  =   (x_1   - x_4)^2  +  (x_6 - x_8)^2
%     F( 6)  =   (x_1   - x_5)^2  +  (x_6 - x_9)^2
%     F( 7)  =    x_2^2 + x_7^2
%     F( 8)  =   (x_3   - x_2)^2  +   x_7^2
%     F( 9)  =   (x_4   - x_2)^2  +  (x_8 - x_7)^2
%     F(10)  =   (x_2   - x_5)^2  +  (x_7 - x_9)^2
%     F(11)  =   (x_4   - x_3)^2  +   x_8^2
%     F(12)  =   (x_5   - x_3)^2  +   x_9^2
%     F(13)  =    x_4^2 +  x_8^2
%     F(14)  =   (x_4   - x_5)^2 + (x_9 - x_8)^2
%     F(15)  =    x_5^2 + x_9^2
%     F(16)  =  -x_1 + x_2
%     F(17)  =        -x_2 + x_3
%     F(18)  =               x_3 - x_4
%     F(19)  =                     x_4 - x_5

neF    = 19;
n      =  9;
Obj    =  1; % The default objective row

% A dense Jacobian of F is provided in snoptuserfun2

[iGfun,jGvar] = find(ones(neF,n));

% Constant Jacobian elements are not treated specially.

iAfun = [];  jAvar = [];  A     = [];

% Multipliers and states

Fstate = zeros(neF,1);
Fmul   = zeros(neF,1);

% Ranges for F.

Flow = zeros(neF,1);    Fupp = ones (neF,1);

Flow( 1:15)  = -Inf;    Fupp(16:neF) =  Inf;

% The Objective row (row 1 by default) is free.

Flow(Obj)    = -Inf;    Fupp(Obj)    =  Inf;

% Ranges for x.

xlow = -Inf*ones(n,1);  xupp =  Inf*ones(n,1);

xlow(1) =  0;
xlow(3) = -1;
xlow(5) =  0;
xlow(6) =  0;
xlow(7) =  0;

xupp(3) =  1;
xupp(8) =  0;
xupp(9) =  0;


x   =  [ .1;
         .125;
         .666666;
         .142857;
         .111111;
         .2;
         .25;
        -.2;
        -.25 ];

% Multipliers and states

xstate = zeros(n,1);
xmul   = zeros(n,1);
