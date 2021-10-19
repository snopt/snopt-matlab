function sqsetr( option, rvalue )
% function sqsetr( option, rvalue )
%     Sets a optional real parameter of sqopt. The string "option" will be read
%     by sqopt. If the string contains a setting that sqopt understands,
%     sqopt will set internal parameters accordingly. For a description of
%     available parameters, please see the SQOPT documentation.
%
sqoptmex(4, option, rvalue);