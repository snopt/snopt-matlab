function [x,F,info]=hs116_1()
% HS Problem 116 with explicit linear constraints.

options.printfile = '';
options.specsfile = which('hs116.spc');
options.system_information ='yes';

[x, xlow, xupp, xmul, xstate, ...
 Flow, Fupp, Fmul, Fstate, ...
 ObjAdd, ObjRow, ...
 A, G] = hs116data;
% A.val, A.row, A.col, ...
% G.row, G.col] = hs116data;
options.name = 'hs116';
options.printfile = 'hs116_1.out';

[x,F,info]= snopt(x, xlow, xupp, xmul, xstate,  ...
		  Flow, Fupp, Fmul, Fstate,     ...
		  @(x)hs116userfun(x), ObjAdd, ObjRow, ...
		  A, G, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,xlow,xupp,xmul,xstate, ...
          Flow,Fupp,Fmul,Fstate,   ...
          ObjAdd,ObjRow,A,G] = hs116data()
%     hs116data defines problem HS116.
%
%     Minimize      x(11) + x(12) + x(13)
%     subject to    linear and nonlinear constraints
    n      = 13;
    neF    = 15;
    Obj    =  1;
    ObjRow = Obj;

    %% parameters
    a = 0.002;
    b = 1.262626;
    c = 1.231059;
    d = 0.03475;
    e = 0.975;
    f = 0.00975;

    % HS Solution
    x = [   0.80377
            0.89999
            0.97095
            0.10000
            0.19081
            0.46057
          574.07758
           74.07758
          500.01615
            0.10000
           20.23309
           77.34768
            0.00673];

    %% initial x
    x = [   0.5;
            0.8;
            0.9;
            0.1;
            0.14;
            0.5;
            489;
             80;
            650;
            450;
            150;
            150;
            150   ];

    %% lower bounds on x
    xlow = [    0.1;
                0.1;
                0.1;
                1e-4;
                0.1;
                0.1;
                0.1;
                0.1;
                500;
                0.1;
                1;
                1e-4;
                1e-4    ];

    %% upper bounds on x
    xupp = [    1  ;
                1  ;
                1  ;
                0.1;
                0.9;
                0.9;
             1000;
             1000;
             1000;
              500;
              150;
              150;
              150   ];

    xstate = zeros(n,1);
    xmul   = zeros(n,1);

    ObjAdd = 0;

    %% Bounds on F

    Flow      = zeros(neF,1); Fupp      = Inf*ones(neF,1);
    Flow(Obj) = -Inf;         Fupp(Obj) = Inf;
    Flow(3)   = -Inf;         Fupp(3)   = 1;
    Flow(4)   =   50;         Fupp(4)   = 250;
    Flow(13)  = -Inf;         Fupp(13)  = 1;
    Flow(15)  =  0.9;         Fupp(15)  = Inf;

    Fmul      = zeros(neF,1);
    Fstate    = zeros(neF,1);

    A = zeros(neF, n);
    A(Obj,11) = 1;
    A(Obj,12) = 1;
    A(Obj,13) = 1;
    A(2,1)   = -1;
    A(2,2)   = 1;
    A(3,7)   = a;
    A(3,8)   = -a;
    A(4,11)  = 1;
    A(4,12)  = 1;
    A(4,13)  = 1;
    A(5,2)   = -1;
    A(5,3)   = 1;
    A(6,13)  = 1;
    A(10,12) = 1;
    A(11,11) = 1;

    G = zeros(neF,n);
    G(6,3) = c*x(10);
    G(6,10) = -b + c*x(3);
    G(7,2) = -d - e*x(5) + 2*f*x(2);
    G(7,5) =  1 - e*x(2);
    G(8,3) = -d - e*x(6) + 2*f*x(3);
    G(8,6) =  1 - e*x(3);
    G(9,1) = -d - e*x(4) + 2*f*x(1);
    G(9,4) =  1 - e*x(1);
    G(10,2) =   c*x(9);
    G(10,9) = -b + c*x(2);
    G(11,1) = c*x(8);
    G(11,8) =-b + c*x(1);
    G(12,1) =-x(8);
    G(12,4) =-x(7) + x(8);
    G(12,5) = x(7);
    G(12,7) = x(5) - x(4);
    G(12,8) =-x(1) + x(4);
    G(13,1) =-a*x(8);
    G(13,2) = a*x(9);
    G(13,5) = 1 + a*x(8);
    G(13,6) = 1 - a*x(9);
    G(13,8) = a*(x(5) - x(1));
    G(13,9) = a*(x(2) - x(6));
    G(14,2) =-500 + x(9) + x(10);
    G(14,3) =-x(10);
    G(14,6) = 500 - x(9);
    G(14,9) = x(2) - x(6);
    G(14,10) =-x(3) + x(2);
    G(15,2) = 1 - a*x(10);
    G(15,3) = a*x(10);
    G(15,10) = -a*(x(2) - x(3));
