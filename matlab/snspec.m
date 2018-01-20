% function inform = snspec( filename )
%     Causes snopt to read its optional parameters from the named file.
%     The format of this file is described in the snopt documentation.
%
%     Returns 101 or 107 if successful.
function inform = snspec( filename )

snoption = 9;
inform   = snoptmex( snoption, filename );
