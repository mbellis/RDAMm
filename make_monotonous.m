%&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% FUNCTION make_monotonous
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% make a vector of values monotonous either by keeping
% the last greatest or last smallest value

%INPUT PARAMETERS
%1- Val : values that must be rendered monotonous
%2- Type : either inc (monotonous increasing) or dec (monotonous decreasing)


%OUTPUT PARAMETERS
%1- Val : values made monotonous

%VERSIONS
%V01 22-3-2010 Refactorinf of the existing version

function Val=make_monotonous (Val,Type)

%verify that Val it is a vector
if ndims(Val)>2
    errordlg('more than two dimensions in Val');
    error('process canceled')
end
if size(Val,1)>1&size(Val,2)>1
    errordlg('Val is not a vector');
    error('process canceled')
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


