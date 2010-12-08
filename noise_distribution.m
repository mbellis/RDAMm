%==============================
% FUNCTION NOISE_DISTRIBUTION 
%==============================
%
% NOISE_DISTRIBUTION calculates noise distribution:
% The distribution of the variation between two series of rank values (one called a baseline (BL)
% and the other called a highline (HL)) is computed
% The units used for representing the variation  is the difference between the ranks of values 
% (rank diff)
% The script allows to compute the values of the variation unit for which there is a predetermined
% percentage (threshold of 10, 3 and 1) of points with a greater variability : these are the 
% calibration curves of the percentiles (or fractiles)
% These curves which are irregular, because they are obtained by computing the points in a sliding
% window, are smoothed by interpolation (spline)
% 
% INPUT PARAMETERS
%  1             HLValues: the highline values
%  2             BLValues: the baseline values
%  3        RankThreshold: a vector of two values allowing to start the process of constructing 
%                          quantile curves using a defined range of rank values 
%                          (>=RankThreshold(1)&<=RankThreshold(2))
%  4               HLRank: the rank of point used as highline
%  5               BLRank: the rank of point used as baseline
%  6           ClearIndex: allow to clear some points from the points used in the calibration
%                          process
%  7             NormType: either 'standardization' (original method in RDAM = var - mean /std)
%                          or 'quantile' (more general method based on percentile curves)
%  8         AnalyseType : either 'transcriptome' or 'chipchip' (some parameters have to be 
%                          adapted to the type of analysis)
%  9            CalibType: the type of data, either
%                          up or down: distinct biological conditions
%                          up: only the values which are increased in the HL vs BL comparison 
%                              are used
%                          down: only the values which are decreased in the HL vs BL comparison
%                                are used
%                                used for chiphip analysis type and for biological conditions 
%                                without replicates
%                          idem: replicates of the same biological condition
%                                [HL,BL] is compared to [BL,HL] (no distinction between increased
%                                and decreased distribution which should be equal)
%  10 SizerFittingDegrees: used by sizer (the procedure of curve smoothing).Indicates the number of 
%                          increasing fitting degrees used
%  11         DisplayFlag: indicates if figures must be drawn or not
% 
% OUTPUT PARAMETERS
%  1 RankGrid: sampling range of ranks [0.25:0.25:100]
%  2     Grid: smoothed variation curves corresponding to RankGrid
%  3 ZVarGrid: Normalized variations used for 2-D interpolation
%  4     ZVar: normalized variations (ZVar=interp2(RankGrid,VarGrid,ZVarGrid,Rank,Var))
% 
% VARIABLE NAMES
%   Val: values (ranks)
%   Var: variations (rank differences)
%  Perc: percentile
%   Pos: positive variations
%   Neg: negative variations
%    S~: sorted values
%    P~: values that keep the correspondance with sorted value
% 
% SORT USING
%  Y=f(X) => [SX,SortIndex]=sort(X) => PY=Y(SortIndex)
%  the original order can be find with : [temp ReverseIndex]=sort(SortIndex) => X=SX(ReverseIndex)
% 
%  DIFFERENCE BETWEEN INDEX AND BINDEX:
%  a=[0,2,45,0]
%  Index=find(a==0) => Index=[1,4]
%  Bindex=a==0 => Bindex=[1,0,0,1]


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


function [RankGrid,Grid,ZVarGrid,ZVar]=noise_distribution(HLValues,BLValues,RankThreshold,HLRank,BLRank,ClearIndex,NormType,...
    AnalyseType,CalibType,SizerFittingDegree,DisplayFlag)


RankNb=length(BLValues);
% Not measured probe sets have their rank set to -1 and  are eliminated
NegIndex=find(BLValues==-1);
if ~isempty(NegIndex)
    HLValues(NegIndex)=[];
    BLValues(NegIndex)=[];
end
NegIndex=find(HLValues==-1);
if ~isempty(NegIndex)
    HLValues(NegIndex)=[];
    BLValues(NegIndex)=[];
end

if ~isempty(ClearIndex)
    HLValues(ClearIndex)=[];
    BLValues(ClearIndex)=[];
end


if isequal(CalibType,'idem')
    % the two series of values (highline and baseline) are merged to calculate the positive rank differences
    BLInput=[BLValues;HLValues];
    HLInput=[HLValues;BLValues];
elseif isequal(CalibType,'up')
    % highline and baseline are kept distinct (no replicates, or ChipChip data with two channels :
    % test channel in HL and control channel in BL
    BLInput=BLValues;
    HLInput=HLValues;
elseif isequal(CalibType,'down')
    % highline and baseline are kept distinct (no replicates, or ChipChip data with two channels :
    % test channel in HL and control channel in BL
    BLInput=HLValues;
    HLInput=BLValues;    
end


%sort BLInput and keep the corresponance between HL and BL values
[SBLInput,PHLInput,SortIndex,Temp,Temp]=sort2series(BLInput,HLInput);
InputRankNb=length(BLInput);
SHLInput=sort(HLInput);

% selection of positive and negative variation
VarThreshold=SHLInput(1) - SBLInput(1);
if VarThreshold<0
    VarThreshold=0;
end
Var=PHLInput-SBLInput;
PHLThresholdBindex=PHLInput==100/InputRankNb;
SBLThresholdBindex=SBLInput==100/InputRankNb;

PosVarBindex=Var>VarThreshold;
NegVarBindex=Var<VarThreshold;

ThresholdBindex=(PHLThresholdBindex & PosVarBindex) |(SBLThresholdBindex & NegVarBindex);
Var(ThresholdBindex)=VarThreshold;
PosVarIndex=find(Var>VarThreshold);

% adapt percentile range, windows size and windows step to the number of values
StepNb=1500;
WinStep=ceil(length(PosVarIndex)/StepNb);
WinSize=WinStep*10;
if WinSize>=10000
    PercRange=[0.10:0.10:0.90,0.95,0.99,0.999,0.9999,1];
elseif WinSize>=1000
    PercRange=[0.10:0.10:0.90,0.95,0.99,0.999,1];
elseif WinSize>=100
    PercRange=[0.10:0.10:0.90,0.95,0.99,1];
else
    WinSize=100;
    WinStep=10;
    PercRange=[0.10:0.10:0.90,0.95,0.99,1];
end
PercNb=length(PercRange);


[OutputRes,FRank,FVar,RankGrid,Grid]=quantile_curves(SBLInput,Var,NormType,AnalyseType,PosVarIndex,RankThreshold,WinStep,WinSize,PercRange,SizerFittingDegree,DisplayFlag);


if ~isequal(AnalyseType,'chipchip')
    IncDecPos=1:max(SortIndex);
    IncDecPos=IncDecPos(SortIndex);
    IncDecPos=IncDecPos(PosVarIndex);
    IncIndex=find(IncDecPos<=RankNb);
    DecIndex=find(IncDecPos>RankNb);
    if length(IncIndex)>200000
        RandIndex=ceil(rand(200000,1)*length(IncIndex));
        IncIndex=IncIndex(RandIndex);
    end
    if length(DecIndex)>200000
        RandIndex=ceil(rand(200000,1)*length(DecIndex));
        DecIndex=DecIndex(RandIndex);
    end
end

CalibName=sprintf('P%03u(BL) & P%03u(HL)',BLRank,HLRank);


if DisplayFlag==1
    Colors='bgmkcrbgmkcrb';
    %plot Var vs Rank & percentile lines
    h1=figure;
    set(gcf,'color',[1,1,1])
    set(h1,'name','DISTRIBUTION OF RANK DIFF')
    subplot(1,2,1)
    hold on
    if isequal(AnalyseType,'chipchip')
        plot(FRank,FVar,'b.','markersize',3);
    else
        plot(FRank(IncIndex),FVar(IncIndex),'r.','markersize',3);
        plot(FRank(DecIndex),-FVar(DecIndex),'b.','markersize',3);
    end
    if isequal(NormType,'quantile')
        for PercL=1:PercNb
            plot(RankGrid,Grid.perc{PercL},Colors(PercL))
            MaxVal=max(Grid.perc{PercL});
            line([0,100],[MaxVal,MaxVal],'color',Colors(PercL));
        end
    elseif isequal(NormType,'standardization')
        plot(RankGrid,Grid.mean,'g')
        plot(RankGrid,Grid.std,'b')
    end
    set(gca,'box','on')
    xlabel('min(rank(Replicate 1),rank(Replicate 2))')
    ylabel('RankDiff = Rank(Replicate 2)-Rank(Replicate 1)')
    title(sprintf('%s: Distribution of Rank Diff',CalibName))
end

if isequal(NormType,'quantile')
    %fill a (PercNb+1) x length(RankGrid) matrix with  percentile Var values
    %the last line is the limit value (100 for Var==100)
    GridNb=length(RankGrid);
    PercVal=zeros(PercNb+1,GridNb);
    for PercL=1:PercNb
        PercVal(PercL+1,:)=Grid.perc{PercL}';
    end
    PercVal=[PercVal;ones(1,length(RankGrid))*100];


    %Var where  interpolation is to be made
    VarGrid=0:0.5:100;

    %fraction Value to be interpolated
    FractVal=zeros(PercNb+1,1);
    for PercL=1:PercNb
        FractVal(PercL+1)=max(Grid.perc{PercL});
    end
    FractVal=[FractVal;100];

    %Interpolation (column by column, that is rank by rank)of FractVal onto VarGrid
    ZVarGrid=zeros(length(VarGrid),GridNb);
    for GridL=1:GridNb
        try
            ZVarGrid(:,GridL)=interp1(PercVal(:,GridL),FractVal,VarGrid,'linear')';
        catch
            %error if some X values are equal
            IdemRank=find(diff(PercVal(:,GridL))==0);
            while ~isempty(IdemRank)
                PercVal(IdemRank,GridL)=PercVal(IdemRank,GridL)+randperm(length(IdemRank))'*2*eps;              
                try
                    ZVarGrid(:,GridL)=interp1(PercVal(:,GridL),FractVal,VarGrid,'linear')';
                catch
                end
                IdemRank=find(diff(PercVal(:,GridL))==0);
            end
            
        end
    end
    if ~isempty(find(isnan(ZVarGrid)))
        h=errordlg('exist NaN in ZVarGrid ! Process aborted');
        waitfor(h)
    end
    %calculus of equalized variation (all variations representing the same percentile along the
    % Rank axis receive one equalized variation interpolated from the max(Var|percentile) vs percentile graph
    ZVar=interp2(RankGrid,VarGrid,ZVarGrid,FRank,FVar);
elseif isequal(NormType,'standardization')
    ZVarGrid=[];
    ZVar=(FVar-interp1(RankGrid,Grid.mean,FRank))./interp1(RankGrid,Grid.std,FRank);
end

if DisplayFlag==1
    subplot(1,2,2);
    hold on
    if isequal(AnalyseType,'chipchip')|isequal(CalibType,'distinct')
        plot(FRank,ZVar,'r.','markersize',3);
    else
        if isequal(NormType,'quantile')
            plot(FRank(IncIndex),ZVar(IncIndex),'r.','markersize',3);
            plot(FRank(DecIndex),-ZVar(DecIndex),'b.','markersize',3);
        elseif isequal(NormType,'standardization')
            %in standardization procedure ZVar overlap
            plot(FRank(IncIndex),ZVar(IncIndex)+2,'r.','markersize',3);
            plot(FRank(DecIndex),-ZVar(DecIndex)-2,'b.','markersize',3);
        end
    end
    set(gca,'box','on')
    xlabel('min(rank(Replicate 1),rank(Replicate 2))')
    if isequal(NormType,'quantile')
        ylabel('ZVar (quantile normalised Rank Diff)')
        title(sprintf('%s: Distribution of Quantile Normalized RD',CalibName))
    elseif isequal(NormType,'standardization')
        ylabel('ZVar (standardized Rank diff)')
        title(sprintf('%s: Distribution of Standardized RD',CalibName))
    end
    set(h1,'units','normalized')
    set(h1,'Position',[0.05 0.05 0.80 0.60])
    set(h1,'Color',[1 1 1])
end
