function snset( option )
% function snset( option )
%     Sets a optional parameter of sqopt. The string "option" will be read
%     by sqopt. If the string contains a setting that sqopt understands,
%     sqopt will set internal parameters accordingly. For a description of
%     available parameters, please see the SNOPT documentation.
%
snoptmex(2, option);