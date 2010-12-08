%=====================
% CUMUL_DISTRIBUTION 
%=====================
% 
% CUMUL_DISTRIBUTION calculates a cumulative frequency distribution.
% 
% INPUT PARAMETERS
%  1      Values: the values to be cumulated
%  2       Limit: values <=Limit are cleared
%  3    FlipFlag: if equal to 1, values are flipped
%  4 DisplayFlag: display cumulative curve
% 
% OUTPUT PARAMETERS
%  1 SortedFValues: sorted values that are filtered (values <= Limit have been cleared)
%  2    FCumulDist: cumulative frequency distribution of filtered values
%  3 SortedUValues: sorted unique values (doublons in filtered values have been cleared)
%  4    UCumulDist: cumulative frequency distribution of unique values
%  5      MinValue: minimal value of filtered values (and of unique values by the way)
%  6      MaxValue: maximal value of filtered values (and of unique values by the way)


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
   h=figure;
   set(gcf,'color',[1,1,1])
   set(h,'name','CUMULATIVE DISTRIBUTION')
   hold on
   plot(SortedUValues,UCumulDist,'b+');
   FittedCumulDisp=interp1(SortedUValues,UCumulDist,[MinValue:0.1:MaxValue],'linear');
   FittedCumulDisp(1)=1;
   plot([MinValue:0.1:MaxValue],FittedCumulDisp,'co');
end
