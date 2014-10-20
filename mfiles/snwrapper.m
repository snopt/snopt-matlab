function [F,G] = snwrapper(x,userfun,nF,varargin)
%function [F,G] = snwrapper(x,userfun);
%   Wrapper to allow variable length arguments in
%   user-supplied functions.


if ( nargin == 3 ),
  try
    [F,G] = userfun(x);
  catch
    F     = userfun(x);
    G     = [];
  end

else
  if ( nargin == 4 ),
    myobj = varargin{1};

    try
      [obj,grad] = myobj(x);

      F          = [ obj; zeros(nF-1,1) ];
      G          = [ grad; zeros(nF-1,n) ];
    catch
      [obj] = myobj(x);

      F     = [ obj; zeros(nF-1,1) ];
      G     = [];
    end

  elseif ( nargin == 5 ),
    myobj   = varargin{1};
    nonlcon = varargin{2};

    try
      [obj,grad]      = myobj(x);
      [c,ceq,dc,dceq] = nonlcon(x);

      z = nF - 1 - size(c,1) - size(ceq,1);
      F = [ obj; c; ceq; zeros(z,1) ];
      G = [ grad; cd; ced; zeros(z,n) ];

    catch
      [obj]   = myobj(x);
      [c,ceq] = nonlcon(x);

      z = nF - 1 - size(c,1) - size(ceq,1);
      F = [ obj; c; ceq; zeros(z,1) ];
      G = [];
    end
  end
end
