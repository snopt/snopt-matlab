%Test Script.

format compact;
setpath;  % defines the path

fprintf('\n============================================================= ');
fprintf('\nt1diet: Solving diet LP problem using SNOPT ... ');
t1diet;

fprintf('\n============================================================= ');
fprintf('\nsntoy: Solving toy problem using SNOPT ... ');
sntoy;

fprintf('\n============================================================= ');
fprintf('\ntoymin: Solving toy problem using fmincon-style SNOPT ... ');
toymin;

fprintf('\n============================================================= ');
fprintf('\nhsmain: snopt solves hs47 ... ');
hsmain;

fprintf('\n============================================================= ');
fprintf('\nsntoy2: snopt solves toy problem ... ');
sntoy2;

fprintf('\n============================================================= ');
fprintf('\nsnoptmain: snopt solves hexagon with no derivatives ... ');
snoptmain;

fprintf('\n============================================================= ');
fprintf('\nsnoptmain2: snopt solves hexagon with dense Jacobian ... ');
snoptmain2;

fprintf('\n============================================================= ');
fprintf('\nsnoptmain3: snopt solves hexagon with some derivatives ... ');
snoptmain3;

fprintf('\n============================================================= ');
fprintf('\nhs116: snopt solves hs116 ... ');
hs116;

fprintf('\n============================================================= ');
fprintf('\nspring: snopt solves spring ... ');
springa;
