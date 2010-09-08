% CUMDISTPOS cumulative frequency distribution
function [SortedFValues,FCumulDist,SortedUValues,UCumulDist,MinValue,MaxValue]...
   =cumul_distribution(Values,Limit,FlipFlag,DisplayFlag)

if nargin~=4
    h=errordlg('cumul_distribution needs four parameters !');
    waitfor(h)
end

%if Values are in a matrix format, transform them into a vector format
Values=Values(:);
%clear values <= Limit
if ~isempty(Limit)
   LimitBindex=Values>Limit;
else
   LimitBindex=ones(size(Values));
end
%keep not nan values
ExistBindex=~isnan(Values);
FilterBindex=LimitBindex & ExistBindex;

%filter values
FValues=Values(FilterBindex);
SortedFValues=sort(FValues);

MinValue=min(SortedFValues); 
MaxValue=max(SortedFValues);

% calculate the cumulative distribution
FCumulDist=[1:length(SortedFValues)]'./length(SortedFValues);
if FlipFlag==1
   FCumulDist=flipud(FCumulDist);
end

% Doublons and eventual non monotonic points are removed in order to use
% polyfit function in a later step
DoublonIndex=find(diff(SortedFValues)==0);
SortedUValues=SortedFValues;
SortedUValues(DoublonIndex)=[];
% Do not recalculate the distribution
% because one will get an error since some data have been removed.
% Instead keep the original value of the distribution, simply by deleting
% the distribution points corresponding to the removed SortedFValues
UCumulDist=FCumulDist;
UCumulDist(DoublonIndex)=[];

if DisplayFlag==1
   figure
   hold on
   plot(SortedUValues,UCumulDist,'b+');
   FittedCumulDisp=interp1(SortedUValues,UCumulDist,[MinValue:0.1:MaxValue],'linear');
   FittedCumulDisp(1)=1;
   plot([MinValue:0.1:MaxValue],FittedCumulDisp,'co');
end
