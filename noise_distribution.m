%&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% FUNCTION noise_distribution
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% Calculates noise distribution
% The distribution of the variation between two series of rankk values (one called a baseline (BL)
% and the other called a highline (HL)) is computed
% The units used for representing the variation  is the difference between the ranks of values (rank diff)
% The script allows to compute the values of the variation unit for which there is a predetermined percentage
%(threshold of 10%, 3% and 1%) of points with a greater variability : these are the calibration curves of the
% percentiles (or fractiles)
% These curves which are irregular, because they are obtained by computing the points in a sliding window,
% are smoothed by interpolation (spline)

% INPUT PARAMETERS
% HLValues: the highline values
% BLValues: the baseline values
% RankThreshold: a vector of two values allowing to start the process of constructing quantile curves using
%                a defined range of rank values (>=RankThreshold(1)&<=RankThreshold(2))
% HLNname: the highline name
% BLName: the baseline name
% NormalisationType: either 'standardisation' (original method in RDAM = var - mean /std)
%                    or 'quantile' (more general method based on percentile curves)
% AnalyseType : either 'transcriptome' or 'chipchip' (some parameters have to be adapted to the type of analysis)
% DataType: the type of data, either
%   distinct: distinct biological conditions
%              only the values which are increased in the HL vs BL comparison are used
%              used for chiphip analysis type
%              and for biological conditions without replicates
%   replicate: replicates of the same biological condition
%               [HL,BL] is compared to [BL,HL] (no distinction between increased and decreased distribution which should be equal)
%SizerFittingDegrees: used by sizer (the procedure of curve smoothing).Indicates the number of increasing fitting degrees used
%SaveFlag: indicates if the output must be saved or not
%DisplayFlag: indicates if figures must be drawn or not

%OUTPUT PARAMETERS
%RankGrid: sampling range of ranks [0.25:0.25:100]
%Grid: smoothed variation curves corresponding to RankGrid
%ZVarGrid: Normalized variations used for 2-D interpolation 
%ZVar: normalized variations (ZVar=interp2(RankGrid,VarGrid,ZVarGrid,Rank,Var))

% VARIABLE NAME
%Val values (ranks)
%Var variations (rank differences)
%Perc percentile
%Pos positive variations
%Neg negative variations
%S~ sorted values
%P~ values that keep the correspondance with sorted value
%Y=f(X) => [SX,SortIndex]=sort(X) => PY=Y(SortIndex)
%the original order can be find with : [temp ReverseIndex]=sort(SortIndex) => X=SX(ReverseIndex)

%difference between Index and Bindex:
%a=[0,2,45,0]
%Index=find(a==0) => Index=[1,4]
%Bindex=a==0 => Bindex=[1,0,0,1]


function [RankGrid,Grid,ZVarGrid,ZVar]=noise_distribution(HLRanks,BLRanks,RankThreshold,HLName,BLName,AnalyseType,DataType,SizerFittingDegree,SaveFlag,DisplayFlag)
% global variables used to pass parameters (used in the context of an application)
%global F K P S
%[RankGrid,Grid,ZVarGrid,ZVar]=noise_distribution(Data{2}{1}.rank,Data{3}{1}.rank,[0,0],'point2','point3','transcriptome','distinct',7,0,1)

% Not measured probe sets have their rank set to -1 and  are eliminated
RankNb=length(BLRanks);
ClearIndex=find(BLRanks==-1);
if ~isempty(ClearIndex)
    HLRanks(ClearIndex)=[];
    BLRanks(ClearIndex)=[];
end
ClearIndex=find(HLRanks==-1);
if ~isempty(ClearIndex)
    HLRanks(ClearIndex)=[];
    BLRanks(ClearIndex)=[];
end


if isequal(DataType,'replicate')
    % the two series of values (highline and baseline) are merged to calculate the positive rank differences
    BLInput=[BLRanks;HLRanks];
    HLInput=[HLRanks;BLRanks];
else isequal(DataType,'distinct')
    % highline and baseline are kept distinct (no replicates, or ChipChip data with two channels :
    % test channel in HL and control channel in BL
    BLInput=BLRanks;
    HLInput=HLRanks;
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


[OutputRes,FRank,FVar,RankGrid,Grid]=quantile_curves(SBLInput,Var,AnalyseType,PosVarIndex,RankThreshold,WinStep,WinSize,PercRange,SizerFittingDegree,DisplayFlag);


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

CalibName=sprintf('P%s(BL) & P%s(HL)',BLName,HLName);


if DisplayFlag==1
    Colors='bgmkcrbgmkcrb';
    %plot Var vs Rank & percentile lines
    h1=figure;
    subplot(1,2,1)
    hold on
    if isequal(AnalyseType,'chipchip')
        plot(FRank,FVar,'b.','markersize',3);
    else
        plot(FRank(IncIndex),FVar(IncIndex),'r.','markersize',3);
        plot(FRank(DecIndex),-FVar(DecIndex),'b.','markersize',3);
    end
    for PercL=1:PercNb
        plot(RankGrid,Grid.perc{PercL},Colors(PercL))
        MaxVal=max(Grid.perc{PercL});
        line([0,100],[MaxVal,MaxVal],'color',Colors(PercL));
    end
    set(gca,'box','on')
    xlabel('Rank(Test)')
    ylabel('RankDiff = Rank(Control)-Rank(Test)|')
    title(sprintf('%s: Distribution of Rank Diff',CalibName))
end

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
    ZVarGrid(:,GridL)=interp1(PercVal(:,GridL),FractVal,VarGrid,'linear')';
end
if ~isempty(find(isnan(ZVarGrid)))
    h=errordlg('exist NaN in ZVarGrid ! Process aborted');
    waitfor(h)
end
%calculus of equalized variation (all variations representing the same percentile along the
% Rank axis receive one equalized variation interpolated from the max(Var|percentile) vs percentile graph
ZVar=interp2(RankGrid,VarGrid,ZVarGrid,FRank,FVar);
if DisplayFlag==1
    subplot(1,2,2);
    hold on
    if isequal(AnalyseType,'chipchip')|isequal(DataType,'distinct')
        plot(FRank,ZVar,'r.','markersize',3);
    else
        plot(FRank(IncIndex),ZVar(IncIndex),'r.','markersize',3);
        plot(FRank(DecIndex),ZVar(DecIndex),'b.','markersize',3);
    end
    set(gca,'box','on')
    xlabel('Rank(Test)')
    ylabel('quantile normalised Rank Diff')
    title(sprintf('%s: Distribution of Quantile Normalized RD',CalibName))
    set(h1,'units','normalized')
    set(h1,'Position',[0.05 0.05 0.80 0.60])
    set(h1,'Color',[1 1 1])
%     cd(K.dir.rescalib)
%     set(h1,'name',['Quantile Normalisation : ',CalibName])
%     saveas(h1,sprintf('%s_rdqn_%s.bmp',CalibName,P.par.array))
%     close(h1)
end

% RECOVER RESULT IN S.noise structure
% if ~isfield(S,'noise')
%     SPos=1;
%     S.noise.pos{P.par.arrayrank}(BLName,HLName)=SPos;
% else
%     if size(S.noise.pos{P.par.arrayrank},1)>=BLName&size(S.noise.pos{P.par.arrayrank},2)>=HLName
%         SPos=S.noise.pos{P.par.arrayrank}(BLName,HLName);
%         if SPos==0
%             SPos=length(S.noise.mean)+1;
%             S.noise.pos{P.par.arrayrank}(BLName,HLName)=SPos;
%         end
%     else
%         SPos=length(S.noise.mean)+1;
%         S.noise.pos{P.par.arrayrank}(BLName,HLName)=SPos;
%     end
% end
% InterpNb=length(RankGrid);
% 
% if isequal(NormalisationType,'standardisation')
%     NormLine=ones(InterpNb,1);
%     TestLine=(NormLine.*Grid.std)+ Grid.mean;
%     recover current pos in S.noise structure
%     S.noise.testcurve{SPos,1}=TestLine;
%     normalize on the maximum possible (surface under the diagonal = 401 points * 100 (max value possible) /20050) => 0 <= testsurf <= 1
%     S.noise.testsurf{P.par.arrayrank}(BLName,HLName)=sum(TestLine)/20050;
%     S.noise.mean{SPos,1}=Grid.mean;
%     S.noise.std{SPos,1}=Grid.std;
%     GraphicalControl=0;
%     [Temp,Temp,UFiltered_norm_var_calib_BL,UFiltered_plevel_calib_BL,...
%         Min_norm_var_calib_BL,Max_norm_var_calib_BL]...
%         =CUMDISTPOS(Fitted_norm_phased_incdec_var(P.calib.incdecindex),1,GraphicalControl);
%     % %             Graphical_control=0;
%     S.noise.zrd{SPos,1}=UFiltered_norm_var_calib_BL;
%     S.noise.pvalue{SPos,1}=Pv;
%     S.noise.minzrd{P.par.arrayrank}(BLName,HLName)=Min_norm_var_calib_BL;
%     S.noise.maxzrd{P.par.arrayrank}(BLName,HLName)=Max_norm_var_calib_BL;
% elseif isequal(NormalisationType,'quantile')
%     take median for calculating testline
%     TestLine=Grid.perc{5};
%     S.noise.ftestcurve{SPos,1}=TestLine;
%     S.noise.ftestsurf{P.par.arrayrank}(BLName,HLName)=sum(TestLine)/20050;
%     S.noise.mean{SPos,1}=Grid.mean;
%     S.noise.std{SPos,1}=Grid.std;
%     GraphicalControl=0;
%     [Temp,Temp,ZRd,PValue,MinVar,MaxVar]=CUMDISTPOS(ZVar,1,GraphicalControl);
%     if max(ZRd)<100
%         ZRd=[ZRd;100];
%         PValue=[PValue;PValue(end)/10];
%     end
%     RdRange=0:0.25:100;
%     Pv=interp1(ZRd,PValue,RdRange);
%     if isnan(Pv(1))
%         Pv(1)=1;
%     end
%     S.noise.fzrd{SPos,1}=RdRange;
%     S.noise.fpvalue{SPos,1}=Pv;
%     S.noise.fminzrd{P.par.arrayrank}(BLName,HLName)=MinVar;
%     S.noise.fmaxzrd{P.par.arrayrank}(BLName,HLName)=MaxVar;
%     S.noise.zval{SPos,1}=ZVarGrid;
% 
%     if isfield(P.par,'grpname')
%         if ~isempty(P.par.grpname)
%             if isfield(S.noise,'grpname')
%                 GrpRank=strmatch(P.par.grpname,S.noise.grpname,'exact');
%                 if isempty(GrpRank)
%                     GrpRank=length(S.noise.grpname)+1;
%                     S.noise.grpname{GrpRank,1}=P.par.grpname;
%                 else
%                     GrpRank=GrpRank(1);
%                 end
%             else
%                 GrpRank=1;
%                 S.noise.grpname{GrpRank,1}=P.par.grpname;
%             end
%             S.noise.testgrp(GrpRank,P.par.arrayrank)=sum(TestLine)/20050;
%         end
%     end
% 
% end
% if SaveFlag==1
%     cd(K.dir.stat);
%     eval(sprintf('save %s S',P.file.stat))
% end
