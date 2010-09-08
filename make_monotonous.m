%&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% FUNCTION make_monotonous
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% make a vector of values monotonous either by keeping
% the last greatest or last smallest value

%INPUT PARAMETERS
%1- Val : values that must be rendered monotonous
%2- Type : either inc (monotonous increasing) or dec (monotonous decreasing)
%3- FlipFlag : if ==1 flip the data before makinf them monotonous


%OUTPUT PARAMETERS
%1- Val : values made monotonous

%VERSIONS
%V01 22-03-2010 Refactorinf of the existing version
%V02 14-04-2010 Add FlipFlag

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

LastVal=Val(1);
if isequal(Type,'inc')
    for i=1:length(Val)
        if Val(i)<LastVal
            Val(i)=LastVal;
        else
            LastVal=Val(i);
        end
    end
elseif isequal(Type,'dec')
    for i=1:length(Val)
        if Val(i)>LastVal
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

