% function inform = sqspec( filename )
%     Causes sqopt to read its optional parameters from the named file.
%     The format of this file is described in the sqopt documentation.
%
%     Returns 101 or 107 if successful.
function inform = sqspec( filename )

snoption = 9;
inform   = sqoptmex( snoption, filename );
