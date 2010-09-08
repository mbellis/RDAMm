function [msmooth,xgrid,vh] = gpfam(data,vgridp,eptflag,nh) 
% GPFAM, General Purpose FAMily approach to 1-d smoothing (kde or reg)
%     Does a family of kernel smooths (Gaussian kernel), for
%     either density estimation, or regression, using a binned 
%     implementation.  Family is 15 (usually) logarithmically
%     spaced smooths, centered at a good data driven bandwidth
%     (SJPI for kde, RSW for reg).  Spread is determined by peak
%     height.  Curves for bandwidths less than the binwidth, and
%     greater than twice the range are not shown (or computed)
%   Can use first 1, 2, 3 or 4 arguments.
% Inputs:
%     data   - either n x 1 column vector of density estimation data
%                  or n x 2 matrix of regression data:
%                            X's in first column,  Y's in second
%     vgridp - vector of parameters for grid to evaluate at:
%                  0 (or not specified)  -  use endpts of data and 401 bins
%                  [le; lr]  -  le is left end, re is right, 401 bins
%                         (get error message and no return if le > lr)
%                  [le; lr; nb] - le left, re right, and nb bins
%    eptflag - endpoint truncation flag (only has effect when imptyp = 0):
%                  0 (or not specified)  -  move data outside range to
%                                   nearest endpoint
%                  1  -  truncate data outside range
%    nh      - number of different bandwidths (family members).
%                  15 - default when no value specified
% Output:
%     (none)  -  Draws a graph of the result (in the current axes)
%     msmooth -  matrix of smooths
%     xgrid   -  col vector grid of points at which estimate(s) are 
%                    evaluted (useful for plotting), unless grid is input,
%                    can also get this from linspace(le,re,nb)'  
%          vh -  col vector of bandwidths, used in constructing the
%                    family
%
% Assumes path can find personal functions:
%    vec2mat.m
%    gplbinr.m
%    gpkde.m
%    bwsjpib.m
%    bwos.m
%    bwrswb.m

%    Copyright (c) J. S. Marron 1997



%  Set parameters and defaults according to number of input arguments
%
if nargin == 1 ;    %  only 1 argument input, use default vgridp
  ivgridp = 0 ;
else ;              %  xgrid was specified, use that
  ivgridp = vgridp ;
  if length(ivgridp) <= 2 ;  %  then need to add default number of bins
    ivgridp = [ivgridp(1); ivgridp(2); 401] ;
  end ;
end ;

if nargin <= 2 ;    %  at most 2 arguments input, use default endpt trunc
  ieptflag = 0 ;    %  Default
else ;
  ieptflag = eptflag ;    %  Endpt trunc was specified, so use it
end ;

if nargin <= 3 ;    %  at most 3 arguments input, use default # in family
  nhp = 15 ;
          %  Number of smooths in family grid
else ;
  nhp = nh ;
end ;


%  detect whether density or regression data
%
if size(data,2) == 1 ;    %  Then is density estimation
  xdat = data(:,1) ;
  idatyp = 1 ;
else ;                    %  Then assume regression ;
  xdat = data(:,1) ;
  ydat = data(:,2) ;
  idatyp = 2 ;
end ;
n = length(xdat) ;
if ivgridp == 0 ;   %  then use standard default grid
  ivgridp = [min(xdat),max(xdat),401] ;
end ;



%  Bin the data
%
if idatyp == 1 ;        %  Treating as density estimation
  bincts = gplbinr(xdat,ivgridp,ieptflag) ;
elseif idatyp == 2 ;    %  Treating as regression
  bincts = gplbinr([xdat ydat],ivgridp,ieptflag) ;
end ;



%  Get the Central Bandwidth
%
if idatyp == 1 ;        %  Treating as density estimation
  hcent = bwsjpib(bincts,ivgridp,0,-1,ieptflag) ;
          %   0 for default h grid
          %  -1 for using already binned data
elseif idatyp == 2 ;    %  Treating as regression
  hcent = bwrswb(data) ;
end ;



%  Get the upper bandwidth
%
if idatyp == 1 ;        %  Treating as density estimation

  fh = gpkde(bincts,hcent,ivgridp,-1,ieptflag) ;
          %  -1 for binned implementation, with given bincts
  [fmax,imax] = max(fh) ;
  hnew = hcent ;
  c = .6 ;   
          %  upper end is smoothed until height at peak is c * original
  ratio = 1 ;
  while ratio > c ;
    hnew = hnew * 1.05 ;
          %  Could tune 1.05 for: more accuracy - closer to 1
          %                       more speed    - larger
    fhnew = gpkde(bincts,hnew,ivgridp,-1,ieptflag) ;
    ratio = fhnew(imax) / fmax ;
  end ;

elseif idatyp == 2 ;    %  Treating as regression

  [fh, xgrid] = gpnpr(bincts,hcent,ivgridp,-1,1,ieptflag) ;
          %  -1 for binned implementation, with given bincts
          %   1 for local linear
  coeflin = polyfit(data(:,1),data(:,2),1) ;
          %  coefficients of linear fit
  flin = polyval(coeflin,xgrid) ;
          %  height of linear fit at each data point
  [fmaxdif,imaxdif] = max(abs(fh - flin)) ;
          %  largest distance from smooth to linear fit
  hnew = hcent ;
  c = .3 ;   
          %  upper end is smoothed until biggest difference 
          %  from linear is c * original
  ratio = 1 ;
  while ratio > c ;
    hnew = hnew * 1.05 ;
          %  Could tune 1.05 for: more accuracy - closer to 1
          %                       more speed    - larger
    fhnew = gpnpr(bincts,hnew,ivgridp,-1,1,ieptflag) ;
    ratio = abs(fhnew(imaxdif) - flin(imaxdif)) / fmaxdif ;
  end ;

end ;
  hmax = max(hnew, 3 * hcent) ;
          %  make sure range is big enough to be interesting



%  Get the bandwidth range
%
vh = 10.^linspace(log10(hcent^2 / hmax), log10(hmax), nhp) ;

hupper = ivgridp(2) - ivgridp(1) ;
          %  range of the data
hlower = hupper / ivgridp(3) ;
          %  lowest "reasonable h", the binwidth
hupper = 2 * hupper ;
          %  largest "reasonable h", twice the range

flag = vh < hlower ;    
if sum(flag) > 0 ;
  disp(['!!!   Warning from gpfam: truncated ' num2str(sum(flag)) ...
                   ' small h''s   !!!'] ) ;
  vh = vh(~flag) ;  
end ;

flag = vh > hupper ;    
if sum(flag) > 0 ;
  disp(['!!!   Warning from gpfam: truncated ' num2str(sum(flag)) ...
                   ' large h''s   !!!'] ) ;
  vh = vh(~flag) ;  
end ;




%  Get the matrix of smooths
%
if idatyp == 1 ;        %  Treating as density estimation
  [msmooth, xgrid] = gpkde(bincts,vh,ivgridp,-1,ieptflag) ;
elseif idatyp == 2 ;    %  Treating as regression
  [msmooth, xgrid] = gpnpr(bincts,vh,ivgridp,-1,1,ieptflag) ;
end ;



%  Make plots if no numerical output requested
%
if nargout == 0 ;  %  Then no numerical output, but make a plot
                   %  on the current axes
  plot(xgrid,msmooth,'g') ;
          %  Plot most curves in green
    vcurvh = get(gca,'Children') ;
          %  Vector of handles for curves
    set(vcurvh((nhp+1)/2),'LineWidth',2) ;
          %  use fatter line width for one in middle
    set(vcurvh((nhp+1)/2),'Color',[1 1 0]) ;
          %  use yellow color for one in middle
end ;

