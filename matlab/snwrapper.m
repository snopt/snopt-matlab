function [F,G] = snwrapper(x,userfun,nF,n,needF,needG,varargin)
%function [F,G] = snwrapper(x,userfun,n);
%   Wrapper to allow variable length arguments in
%   user-supplied functions.


if ( nargin == 6 ),
  % Normal SNOPTA userfun
  if ( needG > 0 ),
    try
      [F,G] = userfun(x);
    catch
      F     = userfun(x);
      G     = [];
    end
  else
    F     = userfun(x);
    G     = [];
  end

elseif ( nargin == 7 ),
  % 'fmincon'-style call with user-defined
  %  objective function
  myobj = varargin{1};

  if ( needG > 0 ),
    try
      [obj,grad] = myobj(x);
      F          = [ obj; zeros(nF-1,1) ];
      G          = [ grad; zeros(nF-1,n) ];
    catch
      [obj] = myobj(x);
      F     = [ obj; zeros(nF-1,1) ];
      G     = [];
    end
  else
    [obj] = myobj(x);
    F     = [ obj; zeros(nF-1,1) ];
    G     = [];
  end

elseif ( nargin == 8 ),
  % 'fmincon'-style call with user-defined
  %  objective and constraint functions

  myobj   = varargin{1};
  nonlcon = varargin{2};

  if ( needG > 0 ),
    try
      [obj,grad] = myobj(x);
      F = [ obj  ];
      G = [ grad ];
    catch
      [obj] = myobj(x);
      F = [ obj ];
      G = [];
    end

    try
      [c,ceq,dc,dceq] = nonlcon(x);
      z = nF - 1 - size(c,1) - size(ceq,1);
      F = [ F; c; ceq; zeros(z,1) ];
      G = [ G; dc; dceq; zeros(z,n) ];
    catch
      [c,ceq] = nonlcon(x);
      z = nF - 1 - size(c,1) - size(ceq,1);
      F = [ F; c; ceq; zeros(z,1) ];
      G = [];
    end

  else
    [obj] = myobj(x);
    F     = [ obj ];
    G     = [];

    [c,ceq] = nonlcon(x);
    z = nF - 1 - size(c,1) - size(ceq,1);
    F = [ F; c; ceq; zeros(z,1) ];
    G = [];
  end

end
