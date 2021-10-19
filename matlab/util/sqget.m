function [value] = sqget( option )
% function sqget( option )
%     Gets value of an option in sqopt. The string "option" will be read
%     by sqopt. If the string contains a setting that sqopt understands,
%     sqopt will retrive the parameter accordingly. For a description of
%     available parameters, please see the SQOPT documentation.
%
[value] = sqoptmex(5, option);