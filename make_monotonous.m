%===========================
% FUNCTION MAKE_MONOTONOUS %
%===========================

% MAKE_MONOTONOUS makes a vector of values monotonous either by keeping
% the last greatest or last smallest value.
% Outliers are eliminated before by a simple procedure (deviation from mean local value).

% INPUT PARAMETERS
% 1      Val : values that must be rendered monotonous
% 2     Type : either inc (monotonous increasing) or dec (monotonous decreasing)
% 3 FlipFlag : if ==1 flip the data before making them monotonous


% OUTPUT PARAMETERS
% 1 Val : values made monotonous


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤%
%                          c) Michel Bellis                                                %
%                          michel.bellis@crbm.cnrs.fr                                      %
%            Affiliation:  CNRS (Centre National de la Recherche Scientifique - France)    %                               
%  Bioinformatic Project:  ARRAYMATIC => http://code.google.com/p/arraymatic               %
%        Code Repository:  GITHUB => http://github.com/mbellis                             %
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤%

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%  THIS CODE IS DISTRIBUTED UNDER THE CeCILL LICENSE, WHICH IS COMPATIBLE WITH       %
%  THE GNU GENERAL PUBLIC LICENCE AND IN ACCORDANCE WITH THE EUROPEAN LEGISLATION.   %
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%


function Val=make_monotonous (Val,Type,varargin)

%verify that Val it is a vector
if nargin==2
    FlipFlag=0;
elseif nargin==3
    FlipFlag=varargin{1};
else
    errordlg('make_monotonous must have 2 or 3 paramters');
    error('process canceled')
end

if ndims(Val)>2
    errordlg('more than two dimensions in Val');
    error('process canceled')
end
if size(Val,1)>1&size(Val,2)>1
    errordlg('Val is not a vector');
    error('process canceled')
end
if FlipFlag==1
    if size(Val,1)>1
        Val=flipud(Val);
    else
        Val=fliplr(Val);
    end
end

%statistics on interval

ValDiff=abs(diff(Val));
%eliminate zero values which can bias the Limit estimation if there are
%numerous
MemValDiff=ValDiff;
ValDiff(ValDiff==0)=[];
Limit=mean(ValDiff)+2*std(ValDiff);
ValDiff=MemValDiff;
%test the first value
if ValDiff(1)>Limit
    Pos=find(ValDiff<=Limit); 
    LastVal=Val(Pos(1));
    Val(1)=LastVal;
else
    LastVal=Val(1);
end
if isequal(Type,'inc')
    for i=2:length(Val)
        if Val(i)<LastVal|ValDiff(i-1)>Limit
            Val(i)=LastVal;
        else
            LastVal=Val(i);
        end
    end
elseif isequal(Type,'dec')
    for i=2:length(Val)
        if Val(i)>LastVal|ValDiff(i-1)>Limit
            Val(i)=LastVal;
        else
            LastVal=Val(i);
        end
    end
else
    errordlg('Type must be inc or dec');
    error('process canceled')
end

if FlipFlag==1
    if size(Val,1)>1
        Val=flipud(Val);
    else
        Val=fliplr(Val);
    end
end

