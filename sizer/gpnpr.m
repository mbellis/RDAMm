%==================
% FUNCTION GPNR
%==================
%
% GPNPR, General Purpose NonParametric Regression (1-d, Local poly))
%     Does 1-d kernel local polynomial (usually linear) regression,
%     using binned (default), direct (either matrix, or loops for 
%     bigger data sets), or moving window (for higher degree)
%     implementations, with the bandwidth either user specified 
%     (can be vector), or data driven (Ruppert Sheather and Wand ROT
%      or DPI).  Kernel shape is Gaussian.
%   Can use first 1, 2, 3, 4, 5 or 6 arguments.
% Inputs:
%     data   - either n x 2 matrix of Xs (1st col) and Ys (2nd col)
%                or g x 2 matrix of bincts, when imptyp = -1
%     vh     - vector of bandwidths, or specifies data driven:
%                  0 (or not specified)  -  Ruppert, Sheather, Wand DPI
%                  -1  -  Ruppert, Sheather, Wand ROT
%                      Note: <= 0 only works for imptype = 0
%                  >0  -  Use input number (numbers if vector)
%                      Note: this (these) MUST be >0 for imptyp >= 1
%                                 (the direct implementations)
%     vxgrid - vector of parameters for, or values of, grid to evaluate at:
%                  0 (or not specified)  -  use endpts of data and 401 bins
%                  [le; lr]  -  le is left end, re is right, 401 bins
%                         (get error message and no return if le > lr)
%                  [le; lr; nb] - le left, re right, and nb bins
%                  xgrid  -  Use input values 
%                            Note:  need to have more than 3 entries,
%                                 and only works when imptyp = 1, 2 or 3
%    imptyp  - flag indicating implementation type:
%                 -1  -  binned version, and "data" is assumed to be
%                                   bincounts of prebinned data
%                  0 (or not specified)  -  linear binned version
%                                   and bin data here
%                  1  -  Direct matrix implementation
%                  2  -  Slow looped implementation (only useful when
%                            1 creates matrices that are too large)
%                  3  -  Moving window implementation (slow, but allows
%                            higher polynomial degrees than 1 (linear)
%    polydeg - scalar degree of local polynomial.  For imptyp < 3, can
%                  only use:    0 - local constant (ie. Nadaraya-Watson)
%                               1 - local linear (the default)
%                  For imptyp = 3, can be 0,1,2,...
%    eptflag - endpoint truncation flag (only has effect when imptyp = 0):
%                  0 (or not specified)  -  move data outside range to
%                                   nearest endpoint
%                  1  -  truncate data outside range
%  For equally spaced X's, and a return at those X's only,
%      use all 1's in 1st column of bincts, and Y's in the 2nd
% Output:
%     (none)  -  Draws a graph of the result (in the current axes)
%     npr     -  col vector of heights of kernel kernel density estimate,
%                    unless vh is a vector (then have matrix, with
%                    corresponding cols as density estimates)
%                    Note: when a grid is used where the data are
%                    too sparse for a given bandwidth, direct 
%                    implementations return an interpolation of
%                    the data, and binned implementations return
%                    an interpolation of the estimate
%     xgrid   -  col vector grid of points at which estimate(s) are 
%                    evaluted (useful for plotting), unless grid is input,
%                    can also get this from linspace(le,re,nb)'  
%     mker    -  matrix (vector) of kernel functions, evaluated at xgrid,
%                    which can be plotted to show "effective window 
%                    sizes" (currently scaled to have mass 0.05, may
%                    need some vertical rescaling)
%
% Used by:   gpfam.m
%
% Assumes path can find personal functions:
%    vec2mat.m
%    gplbinr.m
%    bwrswb.m
%
%    Copyright (c) J. S. Marron 1997

function [npr,xgrid,mker] = gpnpr(data,vh,vxgrid,imptyp,polydeg,eptflag) 

%  Set parameters and defaults according to number of input arguments
if nargin == 1 ;    %  only 1 argument input, use default bandwidth
  ivh = 0 ;      %  use default Ruppert Sheather Wand DPI
else ;              %  bandwidth was specified, use that
  ivh = vh ;
end ;

if nargin <= 2 ;    %  at most 2 arguments input, use default xgrid
  ivxgrid = 0 ;
else ;              %  xgrid was specified, use that
  ivxgrid = vxgrid ;
end ;

if nargin <= 3 ;    %  at most 3 arguments input, use default implementation
  iimptyp = 0 ;
else ;              %  implementation was specified, use that
  iimptyp = imptyp ;
end ;

if nargin <= 4 ;    %  at most 4 arguments input, use default poly degree
  ipolydeg = 1 ;    %  Default
else ;
  ipolydeg = polydeg ;    %  Poly degree was specified, so use it
end ;

if nargin <= 5 ;    %  at most 5 arguments input, use default endpt trunc
  ieptflag = 0 ;    %  Default
else ;
  ieptflag = eptflag ;    %  Endpt trunc was specified, so use it
end ;


%  Calculate local polynomial smooth
%
if iimptyp > 0 ;    %  Then do direct implementation

  if min(ivh) > 0 ;    %  Then have valid bandwidths, so proceed

    n = length(data) ;

    if length(ivxgrid) > 3 ;  %  Then use input grid
      xgrid = ivxgrid ;
      nbin = length(xgrid) ;
    else ;                    %  Need to generate a grid
      nbin = 401 ;         %  Default
      lend = min(data(:,1)) ;   %  Default
      rend = max(data(:,1)) ;   %  Default
      if length(ivxgrid) >= 2 ;      %  use input endpoints
        lend = ivxgrid(1) ;
        rend = ivxgrid(2) ;
      end ;
      if length(ivxgrid) == 3 ;      %  use number of grid points
        nbin = ivxgrid(3) ;
      end ;

      if lend > rend ;    %  Then bad range has been input
        disp('!!!   Error in gpnpr: invalid range chosen  !!!') ;
        xgrid = [] ;
      else ;
        xgrid = linspace(lend,rend,nbin)' ;
      end ;
    end ;


    %  Loop through bandwidths
    npr = [] ;
    for ih = 1:length(ivh) ;
      h = ivh(ih) ;

      if iimptyp ~= 2  &  iimptyp ~= 3 ;
                    %  Then do direct matrix implementation

        arg = vec2mat((data(:,1) ./ h),nbin) - vec2mat((xgrid' ./ h),n) ;
          %  efficient way to divide all dif's by h
        mwt = exp(-(arg .^2) / 2) ;
          %  exponential part of Gaussian density (constant not
          %         needed since these divide each other later)
        arg = arg * h ;
          %  put back on scale of data.

        s0 = sum(mwt)' ;
          %  sum part of s0, and make result a column vector
        sy0 = sum(mwt .* vec2mat(data(:,2),nbin))' ;
        if ipolydeg ~= 0 ;    %  Then need extra stuff for local linear
          s1 = sum(mwt .* arg)' ;
          s2 = sum(mwt .* arg.^2)' ;
          sy1 = sum(mwt .* arg .* vec2mat(data(:,2),nbin))' ;
        end ;

        arg = 0 ;
        mwt = 0 ;
          %  get rid of these huge matrices as soon as possible

        if ipolydeg == 0 ;    %  then do local constant (NW)
          denom = s0 ; 
        else ;                %  then do local linear
          denom = s2 .* s0 - s1 .* s1 ; 
        end ;

        if sum( (denom / max([denom; eps])) <= eps ) ;
                   %  If denominator has any entry that is effectively 0
          disp('!!!   Warning from gpnpr:  h is too small  !!!') ;
          disp('!!!      returning interpolant of data     !!!') ;

          %  sort data and interpolate
          [sxdat, vsind] = sort(data(:,1)) ;
          sydat = data(:,2) ;
          sydat = sydat(vsind) ;
          nprh = interp1s(sxdat,sydat,xgrid) ;
             %  specially modified version, that allows extrapolation

        else ;    %  then do usual calculations

          if ipolydeg == 0 ;    %  then do local constant (NW)
            nprh = sy0 ./ denom ;
          else ;                %  then do local linear
            nprh = (s2 .* sy0 - s1 .* sy1) ./ denom ;
          end ;

        end ;

        npr = [npr nprh] ;

      else ;   %  Do slower looped implementations
        nprh = [] ;
        for ixg = 1:nbin ;    %  Loop through grid points
          arg = (data(:,1) - xgrid(ixg)) / h ;
          vwt = exp(-(arg .^2) / 2) ;
          %  exponential part of Gaussian density (constant not
          %         needed since these divide each other later)

          if iimptyp == 2 ;   %  Then do formula type of linear fit
                              %  (direct modification of iimptyp = 1)
            arg = arg * h ;
          %  put back on scale of data.
            s0 = sum(vwt) ;
            sy0 = sum(vwt .* data(:,2)) ;
            if ipolydeg ~= 0 ;    %  Then need extra stuff for local linear
              s1 = sum(vwt .* arg)' ;
              s2 = sum(vwt .* arg.^2)' ;
              sy1 = sum(vwt .* arg .* data(:,2))' ;
            end ;

            if ipolydeg == 0 ;    %  then do local constant (NW)
              denom = s0 ; 
            else ;                %  then do local linear
              denom = s2 .* s0 - s1 .* s1 ; 
            end ;

            if denom <= eps ;
                   %  If denominator is effectively 0
              disp('!!!   Warning from gpnpr:  h is too small  !!!') ;
              disp('!!!      returning interpolant of data     !!!') ;

              %  sort data and interpolate
              [sxdat, vsind] = sort(data(:,1)) ;
              sydat = data(:,2) ;
              sydat = sydat(vsind) ;
              nprhx = interp1s(sxdat,sydat,xgrid(ixg)) ;

            else ;    %  then do usual calculations

              if ipolydeg == 0 ;    %  then do local constant (NW)
                nprhx = sy0 ./ denom ;
              else ;                %  then do local linear
                nprhx = (s2 .* sy0 - s1 .* sy1) ./ denom ;
              end ;

            end ;

          else ;   %  (iimptyp = 3)  Then do direct poly fit, in window
            mx = ones(size(data,1),1) ;
                  %  1st column of "polynomial fit design matrix"
            if ipolydeg >= 1 ;
              for ideg = 1:ipolydeg ;
                mx = [mx, (arg .* mx(:,size(mx,2)))] ;
                  %  next column of "polynomial fit design matrix"
              end ;
            end ;
            mwt = diag(vwt) ;

            xpwx = mx' * mwt * mx ;

            if (1 / cond(xpwx)) <= eps ;
                   %  If matrix is effectively singular
              disp('!!!   Warning from gpnpr:  h is too small  !!!') ;
              disp('!!!      returning interpolant of data     !!!') ;

              %  sort data and interpolate
              [sxdat, vsind] = sort(data(:,1)) ;
              sydat = data(:,2) ;
              sydat = sydat(vsind) ;
              nprhx = interp1s(sxdat,sydat,xgrid(ixg)) ;

            else ;    %  then do usual calculations

              polycoef = inv(xpwx) * mx' * mwt * data(:,2) ;
                  %  kernel weighted least squares fit of poly to Ys
              nprhx = polycoef(1) ;
                  %  1st entry is order 0 part
            end ;
          end ;

          nprh = [nprh; nprhx] ;
        end ;
        npr = [npr nprh] ;
      end ;
    end ;

  else ;    %  Have invalid bandwidths

    disp('!!!   Error in gpnpr: A bandwidth is invalid   !!!') ;
    disp('    (Note: cannot use data driven, with direct impl''s)') ;

  end ;

else ;     %  Then do binned implementation

  if iimptyp == -1 ;   %  Then data have already been binned

    if (length(ivxgrid) == 1) | (length(ivxgrid) > 3) ;
                         %  Then can't proceed because don't have bin ends
      disp('!!!   Error: gpnpr needs to know the endpoints   !!!') ;
      disp('!!!            to use this implementation        !!!') ;
      bincts = [] ;
    else ;
      n = sum(data(:,1)) ;
      bincts = data ;

      nbin = 401 ;    %  default value
      lend = ivxgrid(1) ;
      rend = ivxgrid(2) ;
      if length(ivxgrid) == 3 ;          %  then use number of grid points
        nbin = ivxgrid(3) ;
      end ;

      if nbin ~= length(bincts) ;    %  Then something is wrong
        disp('!!!   Warning: gpnpr was told the wrong number of bins   !!!') ;
        disp('!!!            will just use the number of counts.       !!!') ;
        nbin = rows(bincts) ;
      end ;
    end ;

  else ;               %  Then need to bin data

    if length(ivxgrid) > 3 ;  %  Then need to warn of change to default
      disp('!!!   Warning: gpnpr was given an xgrid, and also   !!!') ;
      disp('!!!       asked to bin; will bin and ignore xgrid   !!!') ;
    end ;

    %  Specify grid parameters
    nbin = 401 ;         %  Default
    lend = min(data(:,1)) ;   %  Default
    rend = max(data(:,1)) ;   %  Default
    if (length(ivxgrid) == 2) | (length(ivxgrid) == 3) ;
                                     %  then use input endpoints
      lend = ivxgrid(1) ;
      rend = ivxgrid(2) ;
    end ;
    if length(ivxgrid) == 3 ;          %  then use number of grid points
      nbin = ivxgrid(3) ;
    end ;

    if lend > rend ;    %  Then bad range has been input
      disp('!!!   Error in gpnpr: invalid range chosen  !!!') ;
      bincts = [] ;
    else ;
      bincts = gplbinr(data,[lend,rend,nbin],ieptflag) ;
    end ;

    %  Can do data-based bandwidth selection here, if specified
    if ivh == -1 ;        %  Then use RSW Rule of Thumb
      ivh = bwrswb(data,-1) ;
    elseif min(ivh) <= 0 ;     %  Then be sure to use default,
                               %  RSW Direct Plug In
                               %    (in case an unsupp val was input)
      ivh = bwrswb(data) ;
      if ivh == 0 ;
        disp('!!!   Warning from gpnpr: h_RSW was 0     !!!') ;
        disp('          (going to Rule of Thumb)    ') ;
        ivh = bwrswb(data,-1) ;
      end ;

    end ;

  end ;
  xgrid = linspace(lend,rend,nbin)' ;


  %  Loop through bandwidths
  npr = [] ;
  for ih = 1:length(ivh) ;
    h = ivh(ih) ;

    %  Create vector of kernel values, at equally spaced grid
    delta = (rend - lend) / (nbin - 1) ;    %  binwidth
    k = nbin - 1 ;    %  index of last nonzero entry of kernel vector
    arg = linspace(0,k * delta / h,k + 1)' ;
    kvec0 = exp(-(arg.^2) / 2) / sqrt(2 * pi) ;
    if ipolydeg ~= 0 ;    %  Then need extra stuff for local linear
      arg = arg * h ;
      kvec1 = kvec0 .* arg ;
      kvec2 = kvec1 .* arg ;
    end ;

    %  Do actual kernel smooth
    kvec0 = [flipud(kvec0(2:k+1)); kvec0] ;
          %  construct symmetric kernel
    s0 = conv(bincts(:,1),kvec0) ;
    s0 = s0(k+1:k+nbin) ;
    sy0 = conv(bincts(:,2),kvec0) ;
    sy0 = sy0(k+1:k+nbin) ;
    if ipolydeg ~= 0 ;    %  Then need extra stuff for local linear
      kvec1 = [-flipud(kvec1(2:k+1)); kvec1] ;
          %  skew-symmetric here!
      s1 = conv(bincts(:,1),kvec1) ;
      s1 = s1(k+1:k+nbin) ;
      sy1 = conv(bincts(:,2),kvec1) ;
      sy1 = sy1(k+1:k+nbin) ;

      kvec2 = [flipud(kvec2(2:k+1)); kvec2] ;
          %  construct symmetric kernel
      s2 = conv(bincts(:,1),kvec2) ;
      s2 = s2(k+1:k+nbin) ;
    end ;

    if ipolydeg == 0 ;    %  then do local constant (NW)
      denom = s0 ; 
    else ;                %  then do local linear
      denom = s2 .* s0 - s1 .* s1 ; 
    end ;


    vflag0d = (denom / max([denom; eps])) <= eps ;
          %  ones where denominator is effectively 0

    sflag0d = sum(vflag0d) ;
          %  > 0 when there are some locations with 0 denominator
    if sflag0d > 0 ;
      disp('!!!   Warning from gpnpr:  data too sparse for this h   !!!') ;
      disp('!!!   will interpolate over sparse regions  !!!') ;
      denom(vflag0d) = ones(sflag0d,1) ;
          %  replace 0's by ones (temporarily)
    end ;

    %  main calculation of smooth
    if ipolydeg == 0 ;    %  then do local constant (NW)
      nprh = sy0 ./ denom ;
    else ;                %  then do local linear
      nprh = (s2 .* sy0 - s1 .* sy1) ./ denom ;
    end ;



    if sflag0d > 0 ;     %  then need to come back and fix up 0 denoms
      vflagbf = bincts(:,1) > 0 ;
          %  one where there is some data, in sense that binct > 0

      txbdat = xgrid(vflagbf) ;    
          %  x data, truncated to nonzero bins
      tybdat = bincts(vflagbf,2) ./ bincts(vflagbf,1) ;    
          %  y data as bin averages, truncated to nonzpero bins

      if sflag0d == nbin ;     %  then have denom trouble at each data point,
                               %  so just return interpolant of binned data

        nprh = interp1s(txbdat,tybdat,xgrid) ;

      else ;    %  then have some regions where denom is OK,
                %  so do interpolation only over those regions


        if vflag0d(1) == 1 ;    %  then have data sparsity at left end
                                %  and need to adjust
          [temp,ifge] = max(1 - vflag0d) ;
          %  index of first point with good estimate (where denom is nonzero)
          [temp,ifpbag] = max(vflagbf(ifge:nbin)) ;
          ifpbag = ifpbag + ifge - 1 ;
          %  index of first positive binct after first good estimate

            flag = txbdat < xgrid(ifpbag) ;
            if sum(flag) > 0 ;   %  Then there are some points to int. over
              vx = [txbdat(flag); xgrid(ifpbag)] ;
                  %  x's for interpolation 
              vy = [tybdat(flag); nprh(ifpbag)] ;
                  %  y's for interpolation 
            else ;
              vx = xgrid(ifpbag) ;
                  %  x's for interpolation 
              vy = nprh(ifpbag) ;
                  %  y's for interpolation 
            end ;
          nprh(1:(ifpbag-1)) = interp1s(vx,vy,xgrid(1:(ifpbag-1))) ;
          vflag0d(1:(ifpbag-1)) = zeros(ifpbag-1,1) ;
          %  reset flag, now that values have been fixed
        end ;


        if vflag0d(nbin) == 1 ;    %  then have data sparsity at right end
                                   %  and need to adjust
          [temp,ilge] = max(flipud(1 - vflag0d)) ;
          ilge = nbin + 1 - ilge ;
          %  index of last point with good estimate (where denom is nonzero)
          [temp,ilpbbg] = max(flipud(vflagbf(1:ilge))) ;
          ilpbbg = ilge + 1 - ilpbbg ;
          %  index of last positive binct before last good estimate
          
            flag = txbdat > xgrid(ilpbbg) ;
            if sum(flag) > 0 ;   %  Then there are some points to int. over
              vx = [xgrid(ilpbbg); txbdat(flag)] ;
                  %  x's for interpolation 
              vy = [nprh(ilpbbg); tybdat(flag)] ;
                  %  y's for interpolation 
            else ;
              vx = xgrid(ilpbbg) ;
                  %  x's for interpolation 
              vy = nprh(ilpbbg) ;
                  %  y's for interpolation 
            end ;
          nprh((ilpbbg+1):nbin) = interp1s(txbdat,tybdat, ...
                              xgrid((ilpbbg+1):nbin)) ;
          vflag0d((ilpbbg+1):nbin) = zeros(nbin-ilpbbg,1) ;
          %  reset flag, now that values have been fixed
        end ;


        while sum(vflag0d) > 0 ;   %  loop until all sparsity problems fixed
          [temp,ind0d] = max(vflag0d) ;
          %  an index where still have 0 denom

          [temp,inge] = max(1 - vflag0d(ind0d:nbin)) ;
          inge = inge + ind0d - 1 ;
          %  index of next point with good estimate (where denom is nonzero)
          [temp,inpbag] = max(vflagbf(inge:nbin)) ;
          inpbag = inpbag + inge - 1 ;
          %  index of next positive binct after first good estimate

          [temp,ilge] = max(flipud(1 - vflag0d(1:ind0d))) ;
          ilge = ind0d + 1 - ilge ;
          %  index of last point with good estimate (where denom is nonzero)
          [temp,ilpbbg] = max(flipud(vflagbf(1:ilge))) ;
          ilpbbg = ilge + 1 - ilpbbg ;
          %  index of last positive binct before last good estimate

            flag = (txbdat > xgrid(ilpbbg)) .* (txbdat < xgrid(inpbag)) ;
            if sum(flag) > 0 ;   %  Then there are some points to int. over
              vx = [xgrid(ilpbbg); txbdat(flag); xgrid(inpbag)] ;
                  %  x's for interpolation 
              vy = [nprh(ilpbbg); tybdat(flag); nprh(inpbag)] ;
                  %  y's for interpolation 
            else ;
              vx = [xgrid(ilpbbg); xgrid(inpbag)] ;
                  %  x's for interpolation 
              vy = [nprh(ilpbbg); nprh(inpbag)] ;
                  %  y's for interpolation 
            end ;
          nprh((ilpbbg+1):(inpbag-1)) = interp1(vx,vy, ...
                                      xgrid((ilpbbg+1):(inpbag-1))) ;
          %  linearly interpolate accross gap
          vflag0d((ilpbbg+1):(inpbag-1)) = zeros(inpbag-ilpbbg-1,1) ;
        end ;

      end ;

    end ;

    npr = [npr nprh] ;

  end ;

end ;



%  Create matrix of kernels, if this is needed
%
if nargout == 3 ;
  cent = mean([lend; rend]) ;
          %  centerpoint of evaluation grid
  mih = vec2mat(1 ./ ivh',nbin) ;
  mker = vec2mat((xgrid - cent),length(ivh)) .* mih;
          %  argument of gaussian kernel
  mker = exp(-mker.^2 / 2) .* mih / sqrt(2 * pi) ;
          %  Gaussian kernels with mass 1
  mker = 0.05 * mker ;
          %  Make masses = 0.05
end ;



%  Make plots if no numerical output requested
%
if nargout == 0 ;  %  Then no numerical output, but make a plot
  plot(xgrid,npr) ;
end ;

