function [XGrid, Fitted] = gpanal(data,vrange,nhp) 
% GPANAL, General Purpose data ANALysis
% Inputs:
%     data    - either n x 1 column vector of density estimation data
%                  or n x 2 matrix of regression data:
%                            X's in first column,  Y's in second
%     vrange  - 3 vector:  use minx as vrange(1) and maxx as vrange(2)
%                               and number of grid points as vrange(3)
%     nhp     - value is number in family, and rows in SiZer
%
% Outputs:
%
% Assumes path can find personal functions:
%    bwsjpib.m
%    bwrswb.m
%    gpkde.m
%    gpnpr.m
%    gpsz1.m

%    Copyright (c) J. S. Marron 1997, 1998
%    Modified M Bellis 2001






xdat = data(:,1) ;
ydat = data(:,2) ;
idatyp = 2 ;



mind = vrange(1) ;
maxd = vrange(2) ;
ngrid = vrange(3) ;

Min_pos=find(xdat>=mind);
Min_pos=Min_pos(1);
data=[xdat(Min_pos:end),ydat(Min_pos:end)];

centd = mean([mind;maxd]) ;

%  Set h grid stuff, as in SiZer1
range = maxd - mind ;
binw = range / (ngrid - 1) ;
hmin = 2 * binw ;
hmax = range ;
vh = logspace(log10(hmin),log10(hmax),nhp) ;

[Fitted, XGrid] = gpnpr(data,vh,[mind; maxd; ngrid]) ;

Fitted=[Fitted(1,:);Fitted];
XGrid=[0;XGrid];

