function yi = interp1s(x,y,xi) 
% INTERP1S, linear interpolation, with constant extrapolation outside
%     slight modification of linear version of INTERP1, which
%     allows values xi values outside the range of x, and 
%     returns the closest values at such points.
% Inputs:
%     x  -  col. vector of x values of function (assumed increasing)
%     y  -  col. vector of y values of function
%     xi -  new col. vector of values to plug into function
% Output:
%     yi -  linearly interpolated approximate values of f(xi)
%
%     For more details, try "help interp1"

%    Copyright (c) J. S. Marron 1997


%  Initialize Output
%
nx = length(x) ;
ni = length(xi) ;
yi = zeros(ni,1) ;


%  Fill output with entries below left end
%
flagl = (xi < x(1)) ;
nflagl = sum(flagl) ;
if nflagl > 0 ;     %  If there are new arguments below the left end,
                        %  then need to use smallest value
  yi(flagl) = y(1) * ones(nflagl,1) ;
end ;


%  Fill output with entries above right end
%
flagr = (xi > x(nx)) ;
nflagr = sum(flagr) ;
if nflagr > 0 ;     %  If there are new arguments above the right end,
                        %  then need to use largest value
  yi(flagr) = y(nx) * ones(nflagr,1) ;
end ;


%  Fill rest of output vector
%
flag = (~flagl) & (~flagr) ;
nflag = sum(flag) ;
if nflag > 0 ;     %  If there some interior arguments,
                       %  then use interp1
  yi(flag) = interp1(x,y,xi(flag)) ;
end ;

