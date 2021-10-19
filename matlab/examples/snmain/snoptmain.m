 function [x,F,INFO] = snoptmain()
%function [x,F,INFO] = snoptmain()
% Defines the NLP problem and calls the simplest
% version of the mex interface for snopt.

options.name = 'snmain';
options.screen = 'on';

options.printfile = 'snoptmain.out';

%Get condensed data for the Hexagon problem.
[x,xlow,xupp,xmul,xstate,Flow,Fupp,Fmul,Fstate] = hexagon;

options.maximize = '';

[x,F,INFO,~,~,~,~,output] = snopt(x,xlow,xupp,xmul,xstate,...
				  Flow,Fupp,Fmul,Fstate,@snoptuserfun, ...
				  options);

itns   =  output.iterations;
majors =  output.majors;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [x,xlow,xupp,xmul,xstate,Flow,Fupp,Fmul,Fstate] = hexagon()
%function [x,xlow,xupp,Flow,Fupp] = hexagon()
%
% Defines the problem hexagon:
%   maximize F(1) (the objective row)
%   subject to
%            xlow <=   x  <= xupp
%            Flow <= F(x) <= Fupp
%   where
%     F( 1) =   x_2 x_6 - x_1 x_7 + x_3 x_7 + x_5 x_8 - x_4 x_9 - x_3 x_8
%     F( 2) =    x_1^2 + x_6^2
%     F( 3) =  (x_2   - x_1)^2  + (x_7 - x_6)^2
%     F( 4) =  (x_3   - x_1)^2  +   x_6^2
%     F( 5) =  (x_1   - x_4)^2  + (x_6 - x_8)^2
%     F( 6) =  (x_1   - x_5)^2  + (x_6 - x_9)^2
%     F( 7) =    x_2^2 + x_7^2
%     F( 8) =  (x_3   - x_2)^2  +   x_7^2
%     F( 9) =  (x_4   - x_2)^2  + (x_8 - x_7)^2
%     F(10) =  (x_2   - x_5)^2  + (x_7 - x_9)^2
%     F(11) =  (x_4   - x_3)^2  +   x_8^2
%     F(12) =  (x_5   - x_3)^2  +   x_9^2
%     F(13) =    x_4^2 +  x_8^2
%     F(14) =  (x_4   - x_5)^2 +(x_9 - x_8)^2
%     F(15) =    x_5^2 + x_9^2
%     F(16) =   -x_1 + x_2
%     F(17) =         -x_2 + x_3
%     F(18) =                x_3 - x_4
%     F(19) =                      x_4 - x_5

neF    = 19;
n      =  9;
Obj    =  1; % The default objective row

% Ranges for F.

Flow = zeros(neF,1);
Fupp = ones(neF,1);

Flow( 1:15)  = -Inf;
Fupp(16:neF) =  Inf;

% The Objective row is always free.

Flow(Obj)  = -Inf;
Fupp(Obj)  =  Inf;

% Multipliers and states
Fmul   = zeros(neF,1);
Fstate = zeros(neF,1);

% Ranges for x.

xlow = -Inf*ones(n,1);
xupp =  Inf*ones(n,1);

xlow(1) =  0;
xlow(3) = -1;
xlow(5) =  0;
xlow(6) =  0;
xlow(7) =  0;

xupp(3) =  1;
xupp(8) =  0;
xupp(9) =  0;

% Multipliers and states
xmul   = zeros(n,1);
xstate = zeros(n,1);

x   =  [ .1;
         .125;
         .666666;
         .142857;
         .111111;
         .2;
         .25;
        -.2;
        -.25 ];

