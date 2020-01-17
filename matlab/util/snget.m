function [value] = snget( option )
% function snget( option )
%     Gets value of an option in snopt. The string "option" will be read
%     by snopt. If the string contains a setting that snopt understands,
%     snopt will retrive the parameter accordingly. For a description of
%     available parameters, please see the SNOPT documentation.
%
[value] = snoptmex(5, option);