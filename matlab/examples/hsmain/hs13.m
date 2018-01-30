function [x,F,xmul,Fmul,INFO] = hs13()
% HS 13 (modified so CQ holds)
%
%     Minimize        x(1) + x(2)
%
%     subject to      x(1)^3 - x(2)  >= 0
%
%                              x(2) >= 1.
%
%

options.name = 'hs13';
options.screen = 'on';

options.printfile = 'hs13.out';  % By default, screen output is off;
options.specsfile = which('hs13.spc');

[x,xlow,xupp,xmul,xstate, ...
   Flow,Fupp,Fmul,Fstate, ...
 ObjAdd,ObjRow, ...
 A.val,A.row,A.col,G.row,G.col] = hs13data;

[x,F,INFO]= snopt(x, xlow, xupp, xmul, xstate, ...
		  Flow, Fupp, Fmul, Fstate, ...
		  'hs13userfun', ObjAdd, ObjRow, ...
		  A, G, options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,xlow,xupp,xmul,xstate, ...
	    Flow,Fupp,Fmul,Fstate, ...
	  ObjAdd,ObjRow,A,iAfun,jAvar,iGfun,jGvar] = hs13data()

ObjRow = 1;
ObjAdd = 0;

x      = [ 3 1]';
xlow   = [-Inf,  1 ]';
xupp   = [ Inf, Inf]';
xmul   = zeros(2,1);
xstate = zeros(2,1);

Flow   = [ -Inf,   0]';
Fupp   = [  Inf, Inf]';
Fmul   = zeros(2,1);
Fstate = zeros(2,1);

%     ------------------------------------------------------------------
%     The nonzero pattern of the Jacobian is as follows:
%
%             Column
%            | 1   2
%            +--------
%     row 1  | G   G     Objective row
%         2  | G   G
%
%

A     = [];
iAfun = [];
jAvar = [];

G = [ 1,  1;
      1,  2;
      2,  1;
      2,  2 ];

iGfun = G(:,1); jGvar = G(:,2);
