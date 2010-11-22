%===========================
% FUNCTION QUANTILE_CURVES %
%===========================

% QUANTILE_CURVES calculates quantile curves:
%Variations (rank differences) are sorted according to the rank to which they correspond.
%Successive, range of ranks are selected to calculate either the mean, the std or several percentiles of all
%the positive variations (indexed by the SelIndex parameter) contained in the current range of ranks
%The percentile that are calculated correspond to the fraction stored in PercRange.
%The ranges are determined by a sliding window of size WinSize, shifted by steps equal to WinStep.
%The current X position corresponds to a rank. The first occurence of this rank
%is the current position.
%The positive variations values used for calculating the fractiles are put in the range of positions going from the
%position used in the previous step (the value of which is stored in MemAllPosition) to the actual current position.
%The window size is constant but the process act on a selection of point (corresponding to positive variations),
%so in the final result a window is represented by rank intervals of different length

%INPUT PARAMETERS
% 1               Rank: ranks
% 2                Var: variations (rank diff)
% 3        AnalyseType: either 'transcriptome' or 'chipchip'
% 4           SelIndex: position of positive variations used
% 5      RankThreshold: range of rank (in general at the beginning of the rank range) to be processed as the first
%                       window position (can be empty)
% 6            WinStep: window step
% 7            WinSize: window size
% 8          PercRange: the fraction to which percentile must be calculated
% 9 SizerFittingDegree: the number of different fitting degree used by SIZER
% 10       DisplayFlag: indicates if figures must be drawn or not

% OUTPUT PARAMETERS
% 1 OutputRes: structure with normalized variations
% 2      Rank: rank of positive variations (indexed by SelIndex)
% 3       Var: positive variations
% 4  RankGrid: sampling range of ranks [0.25:0.25:100]
% 5      Grid: smoothed variation curves corresponding to RankGrid

% OUTPUTRES STRUCTURE
%     OutputRes.rank: rank of the X values
%     OutputRes.perc: series of raw percentile curves in the sliding window
%  OutputRes.fitperc: series of smoothed percentiles curves (SIZER)
%     OutputRes.mean:  mean of the variation in the sliding window
%  OutputRes.fitmean: smoothed mean (SIZER)
%      OutputRes.std: std of the variation in the sliding window
%   OutputRes.fitstd: smoothed std (SIZER)
%  OutputRes.stepvar: standardised variation (original RDAM method) in the sliding window
%   OutputRes.fitvar: smoothed standardised variation  (SIZER)
% OutputRes.fraction: percentage of positive variation in the sliding window
% OutputRes.gridperc: smoothed percentiles curves corresponding to a sampling range of ranks [0.25:0.25:100]

% VARIABLE NAMES
% Pos: position (index)
%   F: filtered (selected with SelIndex)
%   S: sorted


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


function [OutputRes,Rank,Var,RankGrid,Grid]=quantile_curves (Rank,Var,CalibType,AnalyseType,SelIndex,RankThreshold,WinStep,WinSize,PercRange,SizerFittingDegree,DisplayFlag)                                
PercNb=length(PercRange);
RankNb=length(Rank);

%initialise OutputRes
OutputRes.rank=Rank;
if isequal(CalibType,'quantile')
    OutputRes.perc=zeros(RankNb,12);
end
OutputRes.mean=zeros(RankNb,1);
OutputRes.std=zeros(RankNb,1);
OutputRes.fraction=zeros(RankNb,1);

%keep rank and var corresponding to positive variations (indexed by SelIndex)
if ~isempty(SelIndex)
    Rank=Rank(SelIndex);
    Var=Var(SelIndex);
end

%calculate the size of the first sliding window
% which contains all the variation values for which the rank is inside a defined interval (RankTrheshold)
% search selected ranks (by SelIndex)
AllNb=length(find(OutputRes.rank>=RankThreshold(1)&OutputRes.rank<=RankThreshold(2)));
FirstIndex=find(Rank>=RankThreshold(1)&Rank<=RankThreshold(2));
PosNb=length(FirstIndex);
if PosNb>0
    if PosNb<100
        %at least 100 values must be processed
        PosNb=100;
    end
else % assign the first window the default window size
    FirstIndex=1:WinSize;
    %find the rank which marks the end of the window in the whole series of ranks
    AllRankPos=find(OutputRes.rank==Rank(WinSize));
    AllNb=AllRankPos(1);
end
%fraction of positive variations in the first window
Fraction=PosNb/AllNb;


CurrVar=Var(FirstIndex);
CurrMean=mean(CurrVar);
CurrStd=std(CurrVar);
%keep last calculated values into memory variables
MemMean=CurrMean;
MemStd=CurrStd;
MemFraction=Fraction;

if isequal(CalibType,'quantile')
    %sort CurrVar to find percentiles
    CurrVar=sort(CurrVar);
    %find the percentiles corresponding to the range passed in paramater (PercRange)
    MemPerc=zeros(1,PercNb);
    for PercL=1:PercNb
        MemPerc(PercL)=CurrVar(floor(length(CurrVar)*(PercRange(PercL))));
    end


    %all the first PosNb places are filled in this round
    %but in the next round new values will replace them
    %from MemPosRankPos
    for PercL=1:PercNb
        OutputRes.perc(1:PosNb,PercL)=MemPerc(PercL);
    end
end
OutputRes.mean(1:PosNb)=CurrMean;
OutputRes.std(1:PosNb)=CurrStd;
OutputRes.fraction(1:PosNb)=Fraction;

%calculate the first position to be changed in the next round
if isempty(FirstIndex)
    %at the next step fill half of the default window
    MemAllPosition=round(AllNb/2);
    MemPosPosition=round(PosNb/2);
else
    %at the next step fill all interval minus window step
    MemAllPosition=PosNb-WinStep+1;
    if MemAllPosition<WinStep+1
        MemAllPosition=WinStep+1;
    end
    MemPosPosition=PosNb-WinStep+1;
end

if MemAllPosition<WinStep+1
    MemAllPosition=WinStep+1;
end
if MemPosPosition<1
    MemPosPosition=1;
end

%cut the data in several bins
Range=0:5:100;
BinNb=length(Range)-1;
BinContent=zeros(1,BinNb);
for BinL=1:BinNb
    BinContent(BinL)=length(find(Rank>Range(BinL)&Rank<=Range(BinL+1)));
end
%calculate the positions of ranks that delimitate the bins
BinPosition=cumsum(BinContent);

% treat all the bins
for BinL=1:BinNb
    if BinL==BinNb
        BinEndPosition=size(Rank,1);
    else
        BinEndPosition=BinPosition(BinL);
    end
    %adapt window size and window step to the current bin's content
    WinStep=ceil(BinContent(BinL)/100);
    WinSize=WinStep*10;
    if WinSize<100
        WinSize=100;
        WinStep=10;
    end
    for CurrPosition=MemPosPosition+round(WinSize/2):WinStep:BinEndPosition-round(WinSize/2)
        %the window is centered on the current X position => WinSize nb of
        %points are selected
        StartPosition=CurrPosition-round(WinSize/2);
        EndPosition=CurrPosition+round(WinSize/2);
        if EndPosition>BinEndPosition
            EndPosition=BinEndPosition;
        end
        if StartPosition<1
            StartPosition=1;
        end

        PosNb=EndPosition-StartPosition+1;

        CurrVar=(Var(StartPosition:EndPosition));
        CurrMean=mean(CurrVar);
        CurrStd=std(CurrVar);
        CurrVar =sort(CurrVar);
        if isequal(CalibType,'quantile')
            Percentile=zeros(1,PercNb);
            for PercL=1:PercNb
                %Percentile(PercL)=CurrVar(ceil(size(CurrVar,1)*(PercRange(PercL))));
                Percentile(PercL)=CurrVar(floor(size(CurrVar,1)*(PercRange(PercL))));
            end
        end
        CurrentXPosition=find(OutputRes.rank==Rank(CurrPosition));
        CurrentXPosition=CurrentXPosition(1);
        %total nb of selected point
        AllNb=(CurrentXPosition-MemAllPosition+1)*(PosNb/WinStep);
        Fraction=PosNb/AllNb;

        PrevBlocSize=round((CurrentXPosition-MemAllPosition+1)/2);
        PrevBloc=ones(PrevBlocSize,1);
        CurrBlocSize=CurrentXPosition-MemAllPosition-PrevBlocSize+1;
        CurrBloc=ones(CurrBlocSize,1);
        if isequal(CalibType,'quantile')
            for PercL=1:PercNb
                % The forward bin equal to previously calculated value is
                % known only at this step
                OutputRes.perc(MemAllPosition:MemAllPosition+PrevBlocSize-1,PercL)=PrevBloc*MemPerc(PercL);
                % The backward bin equal to the currently calculated value is
                % filled out
                OutputRes.perc(MemAllPosition+PrevBlocSize:CurrentXPosition,PercL)=CurrBloc*Percentile(PercL);
                MemPerc(PercL)=Percentile(PercL);
            end
        end
        OutputRes.mean(MemAllPosition:MemAllPosition+PrevBlocSize-1)=PrevBloc*MemMean;
        OutputRes.mean(MemAllPosition+PrevBlocSize:CurrentXPosition)=CurrBloc*CurrMean;
        OutputRes.std(MemAllPosition:MemAllPosition+PrevBlocSize-1)=PrevBloc*MemStd;
        OutputRes.std(MemAllPosition+PrevBlocSize:CurrentXPosition)=CurrBloc*CurrStd;
        OutputRes.fraction(MemAllPosition:MemAllPosition+PrevBlocSize-1)=PrevBloc*MemFraction;
        OutputRes.fraction(MemAllPosition+PrevBlocSize:CurrentXPosition)=CurrBloc*Fraction;

        MemAllPosition=CurrentXPosition+1;
        MemMean=CurrMean;
        MemStd=CurrStd;
        MemFraction=Fraction;
    end
end


%fill the end of OutputRes with an extrapolation of the best fit line constructed on the
% scatter plot on the last position of the sliding window
% control that extrapolated lines don not cross each other (the last value
% of the  current interpolated line must be greater than the last value of
% the previously interpolated line (LastVar)

LastVar=0;
FitParameter=polyfit(OutputRes.rank(MemAllPosition-1-WinSize:MemAllPosition-1),...
    OutputRes.mean(MemAllPosition-1-WinSize:MemAllPosition-1),1);
Try=polyval(FitParameter,OutputRes.rank(MemAllPosition:end));
if Try(end)<=LastVar
    FitParameter=polyfit([OutputRes.rank(MemAllPosition-1),OutputRes.rank(end)],...
        [OutputRes.mean(MemAllPosition-1),LastVar+eps],1);
end
OutputRes.mean(MemAllPosition:end)=polyval(FitParameter,OutputRes.rank(MemAllPosition:end));

LastVar=0;
FitParameter=polyfit(OutputRes.rank(MemAllPosition-1-WinSize:MemAllPosition-1),...
    OutputRes.std(MemAllPosition-1-WinSize:MemAllPosition-1),1);
Try=polyval(FitParameter,OutputRes.rank(MemAllPosition:end));
if Try(end)<=LastVar
    FitParameter=polyfit([OutputRes.rank(MemAllPosition-1),OutputRes.rank(end)],...
        [OutputRes.std(MemAllPosition-1),LastVar+eps],1);
end
OutputRes.std(MemAllPosition:end)=polyval(FitParameter,OutputRes.rank(MemAllPosition:end));

if isequal(CalibType,'quantile')
    LastVar=0;
    for PercL=1:PercNb
        FitParameter=polyfit(OutputRes.rank(MemAllPosition-1-WinSize:MemAllPosition-1),...
            OutputRes.perc(MemAllPosition-1-WinSize:MemAllPosition-1,PercL),1);
        Try=polyval(FitParameter,OutputRes.rank(MemAllPosition:end));
        if Try(end)<=LastVar
            FitParameter=polyfit([OutputRes.rank(MemAllPosition-1),OutputRes.rank(end)],...
                [OutputRes.perc(MemAllPosition-1,PercL),LastVar+eps],1);
        end
        OutputRes.perc(MemAllPosition:end,PercL)=polyval(FitParameter,OutputRes.rank(MemAllPosition:end));
        LastVar=OutputRes.perc(end,PercL);
    end
end


if DisplayFlag==1
    if isequal(CalibType,'quantile')
        Colors='bgmkcrbgmkcrb';
        h1=figure;
        subplot(2,2,1)
        hold on
        for PercL=1:PercNb
            plot(0:100/(length(OutputRes.perc(:,PercL))-1):100,OutputRes.perc(:,PercL),Colors(PercL))
        end
        title('raw percentile curves')
        xlabel('rank')
        ylabel('percentile of variations')
    end

    h2=figure;
    subplot(2,2,1)
    hold on
    plot(0:100/(length(OutputRes.mean)-1):100,OutputRes.mean,'g')
    plot(0:100/(length(OutputRes.std)-1):100,OutputRes.std,'b')
    title('raw mean (g) and std (b) curves')
    xlabel('rank')
    ylabel('percentile of variations')
end

%filled mean must be positive
% The filled values are not corrected if negative values were already there
% in order to display them
% Also the correction procedure does not function correctly if exists other negative values
% than filled one
NegIndex=find(OutputRes.mean<0);
if ~isempty(NegIndex)
    FirstNegPos=NegIndex(1);
    if FirstNegPos>1
        LastPosVar=OutputRes.mean(FirstNegPos-1);
    else
        LastPosVar=OutputRes.mean(1);
    end
    OutputRes.mean(NegIndex)=LastPosVar;
end
% Step Normalized Var
try
OutputRes.stepvar(SelIndex)=(Var-OutputRes.mean(SelIndex))./OutputRes.std(SelIndex);
catch
    'stop'
end

%Smoothing of percentile curves with sizer
%RankGrid recovers the rank values used for interpolation
%Grid recovers the interpolated and smoothed values (SizerFittingDegree
%series of values with increasing degree of fitting (from loose ,7, to high,0)

[RankGrid,Grid.mean]=gpanal_calib([OutputRes.rank(SelIndex),OutputRes.mean(SelIndex)],[0.25,100,400],SizerFittingDegree);
[RankGrid,Grid.std]=gpanal_calib([OutputRes.rank(SelIndex),OutputRes.std(SelIndex)],[0.25,100,400],SizerFittingDegree);
if isequal(CalibType,'quantile')
    for PercL=1:PercNb
        [RankGrid,Grid.perc{PercL}]=gpanal_calib([OutputRes.rank(SelIndex),OutputRes.perc(SelIndex,PercL)],[0.25,100,400],SizerFittingDegree);
    end
end

%plot percentiles curves at two different degree of fitting
if DisplayFlag==1
    if isequal(CalibType,'quantile')
        figure(h1)
        subplot(2,2,2)
        hold on
        for PercL=1:PercNb
            plot(RankGrid,Grid.perc{PercL}(:,3),[Colors(PercL),':'])
            plot(RankGrid,Grid.perc{PercL}(:,4),[Colors(PercL),'-'])
        end
        title('smoothed percentile curves (3(..)&4(-))')
        xlabel('rank')
        ylabel('percentile of variations')
    end

    figure(h2);
    subplot(2,2,2)
    hold on
    plot(RankGrid,Grid.mean,'g')
    plot(RankGrid,Grid.std,'b')
    title('smoothed mean (g) and std (b) curves')
    xlabel('rank')
    ylabel('percentile of variations')
end


%limit of segments that have to be fitted with different precision for mean
%and std curves
PercRank=[];
PercPos=[];
PercRank(1)=90;
CurrPercPos=find(RankGrid>=PercRank(1));
PercPos(1)=CurrPercPos(1);
PercRank(2)=99.75;
CurrPercPos=find(RankGrid>=PercRank(2));
if ~isempty(CurrPercPos)
    PercRank(2)=99;
    CurrPercPos=find(RankGrid>=PercRank(2));
end
%values of precision used for each segment
Sizer=[];
if isequal(CalibType,'quantile')
    CompositFittedPerc=cell(PercNb,1);
end
if ~isempty(CurrPercPos)
    Sizer(1)=4;
    Sizer(2)=3;
    % Make a composit sizer fit by concatenation of three segments having different precision of fitting
    % Low precision between 1 and percentile PercRank(1), medium between PercRank(1) and PercRank(2) and linear fitting between PercRank(2) and 100.
    PercPos(2)=CurrPercPos(1);
    CompositFittedMean=[Grid.mean(1:PercPos(1)-1,Sizer(1));Grid.mean(PercPos(1):PercPos(2)-1,Sizer(2))];
    CompositFittedStd=[Grid.std(1:PercPos(1)-1,Sizer(1));Grid.std(PercPos(1):PercPos(2)-1,Sizer(2))];
    % correct end of calibration curve to end with a positive value and mean=std
    EndMean=max(Grid.mean(end,1),0);
    EndStd=max(Grid.std(end,1),0);
    EndVal=(EndMean+EndStd)/2;
    if EndVal==0
        EndVal=eps;
    end
    % linearisation of end of CompositFittedMean and CompositFittedStd
    Linear_parameter=polyfit([PercRank(2);100],[Grid.mean(PercPos(2),3);EndVal],1);
    Flip_value=polyval(Linear_parameter,RankGrid(PercPos(2):end));
    CompositFittedMean=[CompositFittedMean;Flip_value];

    Linear_parameter=polyfit([PercRank(2);100],[Grid.std(PercPos(2),3);EndVal],1);
    Flip_value=polyval(Linear_parameter,RankGrid(PercPos(2):end));
    CompositFittedStd=[CompositFittedStd;Flip_value];

    if isequal(CalibType,'quantile')
        if isequal(AnalyseType,'chipchip')
            Sizer(1)=2;
            Sizer(2)=3;
            %percentile curves to be corrected
            PercLimit(1,1)=9;
            PercLimit(1,2)=PercNb;
            %percentile curves used as such
            PercLimit(2,1)=1;
            PercLimit(2,2)=8;

        else
            Sizer(1)=4;
            Sizer(2)=0;
            %percentile curves used as such
            PercLimit(2,1)=1;
            PercLimit(2,2)=PercNb;
        end

        %limit of segments that have to be fitted with different precision for mean and std curves
        % for percentile curves
        PercRank(1)=15;
        CurrPercPos=find(RankGrid>=PercRank(1));
        PercPos(1)=CurrPercPos(1);                    %
        PercRank(2)=90;
        CurrPercPos=find(RankGrid>=PercRank(2));
        PercPos(2)=CurrPercPos(1);                    %


        if Sizer(2)~=0
            %CHIPCHIP PROCESSING
            %ad hoc corrections applied for the particular examples that I had to process (may be not a general solution)
            %the limit of segments is set to the cross between two curves in the first 1=>15 range
            %and last 90=>end range
            %only percentile curves in range (PercLimit(1,:) are corrected
            for PercL=PercLimit(1,1):PercLimit(1,2)
                Pos=[];
                Limit=[];
                %correction of start
                %detect the second cross between 3rd and 2nd sizer line
                StartPerc1=Grid.perc{PercL}(1:PercPos(1),Sizer(1));
                StartPerc2=Grid.perc{PercL}(1:PercPos(1),Sizer(2));
                DiffStart=StartPerc2-StartPerc1;
                Pos(1)=1;
                while DiffStart(Pos(1))==0
                    Pos(1)=Pos(1)+1;
                end
                if DiffStart(Pos(1))<0
                    Sign=1;
                else
                    Sign=-1;
                end
                Continue=1;
                %search the first inversion
                Limit(1)=0;
                Limit(2)=0;
                while Continue==1 & Pos(1)<length(DiffStart)
                    Pos(1)=Pos(1)+1;
                    if DiffStart(Pos(1))*Sign>0
                        Limit(1)=Pos(1);
                        %search the second inversion (pos=>neg)
                        while Continue==1 & Pos(1)<PercPos(1)
                            Pos(1)=Pos(1)+1;
                            if DiffStart(Pos(1))*Sign<length(DiffStart)
                                Limit(2)=Pos(1);
                                Continue=0;
                            end
                        end
                        Continue=0;
                    end
                end
                if Limit(2)>0
                    Pos(1)=Limit(2);
                elseif Limit(1)>0
                    Pos(1)=Limit(1);
                else
                    Pos(1)=1;
                end

                %correction of end
                %detect the second cross between 1rd and 2nd sizer line
                EndPerc1=Grid.perc{PercL}(PercPos(2):end,Sizer(1));
                EndPerc2=Grid.perc{PercL}(PercPos(2):end,Sizer(2));
                DiffEnd=EndPerc2-EndPerc1;
                Pos(2)=length(DiffEnd);
                while DiffEnd(Pos(2))==0
                    Pos(2)=Pos(2)-1;
                end
                if DiffEnd(Pos(2))>0
                    Sign=1;
                else
                    Sign=-1;
                end
                Continue=1;
                %search the first inversion
                Limit(1)=0;
                Limit(2)=0;
                while Continue==1 & Pos(2)>1
                    Pos(2)=Pos(2)-1;
                    if DiffEnd(Pos(2))*Sign<0
                        Limit(1)=Pos(1);
                        %search the second inversion (pos=>neg)
                        while Continue==1 & Pos(2)>1
                            Pos(2)=Pos(2)-1;
                            if DiffEnd(Pos(2))*Sign>0
                                Limit(2)=Pos(1);
                                Continue=0;
                            end
                        end
                        Continue=0;
                    end
                end

                if Limit(2)>0
                    Pos(2)=Limit(2);
                elseif Limit(1)>0
                    Pos(2)=Limit(1);
                else
                    Pos(2)=length(DiffEnd)-1;
                end
                Pos(2)=length(Grid.perc{PercL}(:,2))-(length(DiffEnd)-Pos(2))+1;
                CompositFittedPerc{PercL}=[Grid.perc{PercL}(1:Pos(1)-1,Sizer(1));Grid.perc{PercL}(Pos(1):Pos(2)-1,Sizer(2));Grid.perc{PercL}(Pos(2):end,Sizer(1))];
            end
        end
        for PercL=PercLimit(2,1):PercLimit(2,2)
            CompositFittedPerc{PercL}=Grid.perc{PercL}(:,Sizer(1));
        end
    end
    if DisplayFlag==1
        if isequal(CalibType,'quantile')
            figure(h1)
            subplot(2,2,3)
            hold on
            if Sizer(2)>0
                for PercL=PercLimit(1,1):PercLimit(1,2)
                    plot(RankGrid,Grid.perc{PercL}(:,Sizer(2)),'k')
                end
            end
            for PercL=1:PercNb
                plot(RankGrid,CompositFittedPerc{PercL},Colors(PercL))
            end
            %plot variations of ranks > 90
            a=find(Rank>90);
            plot(Rank(a),Var(a),'k.','markersize',3)
            title('composit percentile curves')
            xlabel('rank')
            ylabel('percentile of variations')
        end
        figure(h2)
        subplot(2,2,3)
        hold on
        plot(RankGrid,CompositFittedMean,'g')
        plot(RankGrid,CompositFittedStd,'b')
        %plot variations of ranks > 90
        a=find(Rank>90);
        plot(Rank(a),Var(a),'k.','markersize',3)
        title('composit mean and std curves')
        xlabel('rank')
        ylabel('percentile of variations')
    end

    if isequal(CalibType,'quantile')
        LastVar=0;
        for PercL=1:PercNb
            RankNb=length(CompositFittedPerc{PercL});
            FitParameter=polyfit(PercPos(2)-1:RankNb,...
                CompositFittedPerc{PercL}(PercPos(2)-1:end)',1);
            Try=polyval(FitParameter,PercPos(2)-1:RankNb);
            if Try(end)<=LastVar
                FitParameter=polyfit([PercPos(2)-1,RankNb],...
                    [CompositFittedPerc{PercL}(PercPos(2)-1),LastVar+eps],1);
            end
            CompositFittedPerc{PercL}(PercPos(2)-1:end)=polyval(FitParameter,PercPos(2)-1:RankNb);
            LastVar=CompositFittedPerc{PercL}(end);
        end
    end

    if DisplayFlag==1
        if isequal(CalibType,'quantile')
            figure(h1)
            subplot(2,2,4)
            hold on
            for PercL=1:PercNb
                plot(RankGrid,CompositFittedPerc{PercL},Colors(PercL))
            end
            a=find(Rank>90);
            plot(Rank(a),Var(a),'k.','markersize',3)
            title('corrrected composit percentile curves')
            xlabel('rank')
            ylabel('percentile of variations')
        end

        figure(h2)
        subplot(2,2,4)
        hold on
        plot(RankGrid,CompositFittedMean,'g')
        plot(RankGrid,CompositFittedStd,'b')
        %plot variations of ranks > 90
        a=find(Rank>90);
        plot(Rank(a),Var(a),'k.','markersize',3)
        title('corrected composit mean and std curves')
        xlabel('rank')
        ylabel('percentile of variations')
    end
else
    h=warndlg('No linear fitting of end of curves. Curves could cross each other !');
    waitfor(h)
    CompositFittedMean=[Grid.mean(1:PercPos(1)-1,4);Grid.mean(PercPos(1):end,3)];
    CompositFittedStd=[Grid.std(1:PercPos(1)-1,4);Grid.std(PercPos(1):end,3)];
    if isequal(CalibType,'quantile')
        for PercL=1:PercNb
            CompositFittedPerc{PercL}=[Grid.perc{PercL}(1:PercPos(1)-1,4);Grid.perc{PercL}(PercPos(1):end,3)];
        end
    end
end

%return results
Grid.mean=CompositFittedMean;
Grid.std=CompositFittedStd;
CompositFittedMean=interp1(RankGrid,CompositFittedMean,Rank);
CompositFittedStd=interp1(RankGrid,CompositFittedStd,Rank);
OutputRes.fitmean(SelIndex)=CompositFittedMean;
OutputRes.fitstd(SelIndex)=CompositFittedStd;
OutputRes.fitvar(SelIndex)=(Var-CompositFittedMean)./CompositFittedStd;
if isequal(CalibType,'quantile')
    for PercL=1:PercNb
        Grid.perc{PercL}=CompositFittedPerc{PercL};
        CompositFittedPerc{PercL}=interp1(RankGrid,CompositFittedPerc{PercL},Rank);
    end
    OutputRes.fitperc=CompositFittedPerc;
    for PercL=1:PercNb
        OutputRes.gridperc{PercL}=Grid.perc{PercL};
    end
end