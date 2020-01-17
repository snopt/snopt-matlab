function snsetr( option, rvalue )
% function snset( option )
%     Sets a optional real parameter of snopt. The string "option" will be read
%     by snopt. If the string contains a setting that snopt understands,
%     snopt will set internal parameters accordingly. For a description of
%     available parameters, please see the SNOPT documentation.
%
snoptmex(4, option, rvalue);