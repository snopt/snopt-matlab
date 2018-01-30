function t1diet_lp()

options.screen = 'on';
options.printfile = 't1dietLP.out';
options.specsfile = which('t1diet.spc');


% Linear objective term
c = [3   24   13    9   20   19];

% Initial x and bounds on x
x0   = ones (6,1);
xlow = zeros(6,1);
xupp = [ 4
	 3
	 2
	 8
	 2
	 2 ];


% Linear constraint matrix and bounds
A = [ 110  205  160  160  420  260;
      4   32   13    8    4   14;
      2   12   54  285   22   80];

alow = [ 2000;
	 55;
	 800];
aupp =  Inf*ones(3,1);

options.name = 't1dietlp';

[x,obj,INFO,output,lambda,states] = lpopt(c, x0, xlow, xupp, ...
					  A, alow, aupp, options);
states
lambda
