function snseti( option, ivalue )
% function snset( option )
%     Sets a optional integer parameter of snopt. The string "option" will be read
%     by snopt. If the string contains a setting that snopt understands,
%     snopt will set internal parameters accordingly. For a description of
%     available parameters, please see the SNOPT documentation.
%
snoptmex(3, option, ivalue);