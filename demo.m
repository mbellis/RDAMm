%=======
% DEMO %
%=======

% DEMO run a demonstration
% demo must be run from inside the main directory:
% 'cd .../RDAMm/'
% 'demo'

% comment/uncomment lines to test different combination of parameters as indicated



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
global DataRanks P S 

RootDir=pwd;
cd(RootDir)
DirContent=dir;
Ok=0;
for i=1:length(DirContent)
    if isequal(DirContent(i).name,'rdam.m')
        Ok=1;
        break
    end
end
if Ok==0
    h=errordlg('you must run demo from inside RDAMm directory');
    waitfor(h)
    error('process canceled')
end
DataDir=fullfile(RootDir,'data');
%load global variable
cd(DataDir)
load global
if exist('CalibSet_01.mat','file')
    load CalibSet_01
else
    S=[];
end
P.dir.data=DataDir;
    
CompScheme={[1,2;1,2];[1,2;2,1]};
TGRankList= [1,2];
CGRankList=[3,4];
LoadDataFlag=1;
RankThreshold=[0,0];
CalibType='idem';
ClearIndex=[];
NormType='quantile';
AnalyseType='transcriptome';
SizerFittingDegrees=7;
SingleCalibPointFlag=0;
SingleCalibCurveFlag=0;
CalibUpdateFlag=0;
CalibSaveFlag=1;
DisplayFlag=1;
ComparisonFlag=1;
ResRank=1;
CalibSchemeFlag=0;

if LoadDataFlag
    cd(DataDir)
    load DataRanks
    cd(RootDir)
end



[ZVar,Pv,Ppv,Fdr,Sensitivity,TotalVar]=rdam(CompScheme,TGRankList,CGRankList,LoadDataFlag,...
    RankThreshold,CalibType,ClearIndex,NormType,AnalyseType,SizerFittingDegrees,SingleCalibPointFlag,...
    SingleCalibCurveFlag,CalibUpdateFlag,CalibSaveFlag,DisplayFlag,ComparisonFlag,ResRank,CalibSchemeFlag);

h=figure;
set(gcf,'color',[1,1,1])
set(h,'name','RESULTS')
subplot(1,2,1)
plot(ZVar,-(-log10(Pv)),'g.')
hold on
plot(ZVar,-(-log10(Ppv)),'y.')
set(gca,'box','on')
title('Pv and Ppv')
xlabel('ZVar')
ylabel('log10 of Pv (green) and Ppv (yellow)')
subplot(1,2,2)
plot(ZVar,Fdr,'r.')
hold on
plot(ZVar,Sensitivity,'m.')
set(gca,'box','on')
title('Fdr and Sensitivity')
xlabel('ZVar')
ylabel('Fdr (red) and Sensitivity (magenta)')

