%================
% FUNCTION DEMO %
%================

% DEMO run a demonstration
% demo must be run from inside the main directory:
% 'cd .../RDAMm/'
% 'demo'
%
% GLOBAL VARIABLES
% DataRanks contains data ranks (if LoadDataFlag==1,  DataRanks is empty)
% P contains metadata (P.point : description of points, P.biol: description of biological conditions ...)
% S contains calibration sets


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
function demo()
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
if exist(fullfile(DataDir,'CalibSet_01.mat'),'file')
    load CalibSet_01
else
    S=[];
end
P.dir.data=DataDir;
    
CompScheme={[1,2;1,2];[1,2;2,1]};
TGRankList= [2,3];
CGRankList=[4,5];
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

if LoadDataFlag==0
    cd(DataDir)
    load DataRanks
    cd(RootDir)
else
    DataRanks=[];
end



[ZVar,Pv,Ppv,Fdr,Sensitivity,TotalVar]=rdam(CompScheme,TGRankList,CGRankList,LoadDataFlag,...
    RankThreshold,CalibType,ClearIndex,NormType,AnalyseType,SizerFittingDegrees,SingleCalibPointFlag,...
    SingleCalibCurveFlag,CalibUpdateFlag,CalibSaveFlag,DisplayFlag,ComparisonFlag,ResRank,CalibSchemeFlag);


