%&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% FUNCTION rdam
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% Find which variations are statistically significative when comparing two groups
%   of points, each point being a series of rank values. One group is called the test group (TG)
%   and the other is called the control group (CG). The comparison is TG vs CG.  Increased (devreased) probesets are those
%   which are higher (lower) in the test group.
% A group which contains several points is not 'reduced' to a one-point group by using e.g. the mean of all the points.
%    Therefore, there exist multiple way for doing the comparison between two groups which have variable number of replicates.
%    For exemple if two points are available in each condition (TG1,TG2 in the TG group and CG1,CG2 in the CG group,
%    we can do either a 'parallel' comparison (that is a unique one-to-one comparison as with paired sample
%    that uses only TG1vsCG1 and TG2vsCG2). In this case one set of two comparisons is made. In each comparison a p-value is
%    calculated for each probeset, and combined (e.g. by doing their product) to give the finale statistical values (p-values, FDR,
%    sensitivity). If three, four ... replicates were available, we would have done one set of three, four ... comparisons.
%    If the samples are not paired, we can do crossed comparisons that uses all possible comparisons (TG1vsCG1, TG2vsCG2, TG1vsCG2, TG1vsCG2.
%    In this case two sets of two comparisons are made. If the number of replicates is higher, there exist more possibilities and a comparison
%    sheme is necessary in order to indicate how many sets of comparisons are to be made, and which comparisons are made in each set.
% The comparison scheme is coded in a cell array : { 1st set of comparisons, 2nd set of comparisons, ...}
%   a set of comparisons indicates the points that must be used in the TG and CG group to do each comparison
%   {[1,2;1,2]} indicates that the first set of comparisons is TG1vsCG1 and TG2vsCG2
%   {[1,2;1,2],[1,2;2,1]} indicates that the first set of comparisons is TG1 v sCG1 and TG2vsCG2 and the
%   second set of comparisons is TG1vsCG2 and TG2vsCG1
% For each comparison
%    a calibration set
%    (which allows two things:
%       first it is a way of normalizing the variation (Var => zVar,
%           e.g. Rank Difference (RD) => normalized RD (zRD), either by the standardization procedure (orignal method of RDAM), or by a more sophisticated
%           algorithm (surface mapping),
%       and secondly it gives the p-value for any zVar value (PvCurve))
%    must be used (calculated by the noise_distribution script).

%INPUT PARAMETERS
%1- ChipRank indicates the rank of the chip used (e.g. HG U74 has four chips A,B,C,D)
%2- CompScheme : The comparisons to be made
%3- TGRankList and CGRankList: list of data ranks used in the comparisons
%   if TGRankList=[12,15,23], and ChipRank = 1, the Data used for TG point in position 1 is Data{12}{1}.rank
%4&5- LoadDataFlag = 1 indicates that each data point must be loaded when to be used (no permanent storage of Data in the memory)
%6- RankThreshold: a vector of two values allowing to start the process of constructing quantile curves using
%   a defined range of rank values (>=RankThreshold(1)&<=RankThreshold(2))
%7- CalibType indicates the type of couple of two replicates used to calculate the calibration set for the current comparison.
%   The couple noted {G1,G2}, which is passed as parameter, indicates the
%   G1 and G2 will be identified to the HL and BL lines, respectively, in the noise_distribution function.
%   If (TG1,TG2) and (CG1,CG2) are the two couples of replicates for the the TG and CG condition, respectively, we have three possibilities:
%   CalibType='idem' =>  {TG1,TG2} or {CG1,CG2} is used,
%   CalibType='up' =>  {TG1,CG1} or {TG2,CG2} is used,
%   CalibType='down' => {CG1,TG1} or {CG2,TG2} is used.
%   A rule allows to choose the right couple from the two existing ones (FirstCouple, SndCouple).
%8- ClearIndex: allows to clear some points from the points used in the calibration process(e.g. probe sets that are 
%   linked to gender can'b eliminated from calibration process in case where points used as replicates in the
%   calibration process belong to different gender)
%9- NormType: either 'standardization' (original method in RDAM = var - mean /std)
%   or 'quantile' (more general method based on percentile curves)
%10- AnalyseType : either 'transcriptome' or 'chipchip' (some parameters have to be adapted to the type of analysis)
%11- SizerFittingDegrees: used by sizer (the procedure of curve smoothing).Indicates the number of increasing fitting degrees used
%12- SingleCalibPointFlag = 1 indicates that the same TG (CG) point is used to construct the pairs of TG (CG) points that are
%   considered as replicates in the calibration process (noise_distribution script)
%13- SingleCalibCurveFlag = 1 indicates that a single calibration set is used for all the comparisons
%14-  SaveFlag: indicates if the output must be saved or not
%15- DisplayFlag: indicates if figures must be drawn or not
%16- ComparisonFlag: indicates if the comparison must be done (if not the function is simply used to construct calibration sets which are stored in S.calib)

%VARARGIN PARAMETERS
% if SingleCalibPointFlag==1, HLCalibRank = varargin{1} & CGCalibRank = varargin{2} indicates the ranks of points that are
%   systematically used to construct pairs of  calibration points (HLCalibRank (CGCalibRank) must not be equal to any of
%   the points contained in TGRankList (CGRankList))
% if SingleCalibCurveFlag==1, only one calibration set is used =varargin{3}
% if ClearCalibFlag==1, an index of probeset to be cleared from point used in calibration process
%   is passed : ClearCalibIndex=varargin{4}

%GLOBAL VARIABLES
% Data contains data ranks (if LoadDataFlag==1, Data does not exist)
% K contains directorie names
% P contains metadata (P.point : description of points, P.biol: description of biological conditions ...)
% S contains calibration sets




function rdam(ChipRank,CompScheme,TGRankList,CGRankList,LoadDataFlag,RankThreshold,CalibType,ClearIndex,NormType,...
    AnalyseType,SizerFittingDegree,SingleCalibPointFlag,SingleCalibCurveFlag,SaveFlag,DisplayFlag,ComparisonFlag,varargin)
%rdam(1,{[1,2;1,2];[1,2;2,1]},[1,2],[3,4],0,[0,0],'idem',[],'quantile','transcriptome',7,0,0,0,0,0)
%rdam(1,{[1,2;1,2];[1,2;2,1]},[1,2],[3,4],0,[0,0],'up',[],'quantile','transcriptome',7,0,0,0,0,0)
%rdam(1,{[1,2;1,2];[1,2;2,1]},[1,2],[3,4],0,[0,0],'down',[],'standardization','transcriptome',7,0,0,0,0,0)
global Data K P S

%% Verifications
if SingleCalibPointFlag==1
    if nargin~=18
        h=errordlg('must have 18 paramters');
        waitfor(h)
    else
        HLCalibRank=varargin{1};
        BLCalibRank=varargin{2};
        for i=length(CompScheme)
            if ~isempty(find(CompScheme{i}(1,:)==HLCalibRank))
                h=errordlg('HL nb %u is used in Comparison set %u ! Noise distribution can''t be calculated.',HLCalibRank,i);
                waitfor(h)
            end
            if ~isempty(find(CompScheme{i}(2,:)==BLCalibRank))
                h=errordlg('BL nb %u is used in Comparison set %u ! Noise distribution can''t be calculated.',BLCalibRank,i);
                waitfor(h)
            end
        end
    end
elseif SingleCalibCurveFlag==1
    if nargin~=20
        h=errordlg('must have 20 parameters');
        waitfor(h)
    else
        CoupleRank(1,2)=varargin{3};
        CoupleRank(1,1)=varargin{4};
    end
else
     if nargin~=16
        h=errordlg('must have 16 parameters');
        waitfor(h)
     end
end
CoupleRanks=zeros(2,2);

%% Set up LoopNb & CompNb
% For each comparison set the program starts a new run
LoopNb=length(CompScheme)
% A run is made of CompNb comparisons
CompNb=[];
for LoopL=1:LoopNb
    CompNb=[CompNb,size(CompScheme{LoopL},2)];
end

%%

for LoopL=1:LoopNb
    for CompL=1:CompNb(LoopL)
        %====================================
        %{{{ BL,HL selection => Calibration curves => zVar => pvalues
        % one of the CompNb independant comparison to be made
        Comp_rank=sprintf('%.0f',CompL);
        %%%New_comp_flag=1;
        %------------------------------------
        %{{{ BL & HL SELECTION

        % When there are no group comparison and in case of two duplicates comparison,
        % BL1, BL2, HL1, HL2 and calibration sets are calculated once
        % When there is a group comparison, BL1, BL2,HL1 and HL2 are calculated at each CompL value
        % Calibration curves are then calculated once only if BL3 and HL3 exist

        %find the current pairs of BL points and HL points considered as replicate in the current comparison HL{1}vsBL{1}
        BLIndexes=zeros(2,1);
        HLIndexes=zeros(2,1);
        BLRanks=zeros(2,1);
        HLRanks=zeros(2,1);
        %The index of the first member of the current BL and HL pair
        BLIndexes(1)= CompScheme{LoopL}(2,CompL);
        HLIndexes(1)= CompScheme{LoopL}(1,CompL);
        %Recover the corresponding ranks
        BLRanks(1)=CGRankList(BLIndexes(1));
        HLRanks(1)=TGRankList(HLIndexes(1));
        %Search the second member if SingleCalibPointFlag==0
        %otherwise the second member is always the same and is known
        if SingleCalibPointFlag==0
            %The index of the second member of the current BL and HL pair
            if CompL<CompNb(LoopL)
                BLIndexes(2)=CompScheme{LoopL}(2,CompL+1);
            else
                BLIndexes(2)=CompScheme{LoopL}(2,1);
            end
            BLRanks(2)=CGRankList(BLIndexes(2));

            if CompL<CompNb(LoopL)
                HLIndexes(2)=CompScheme{LoopL}(1,CompL+1);
            else
                HLIndexes(2)=CompScheme{LoopL}(1,1);
            end
            HLRanks(2)=TGRankList(HLIndexes(2));
        end

        %Recover the data (ranks) for each BL and HL
        for i=1:3
            HL{i}=[];
            BL{i}=[];
        end
        if SingleCalibPointFlag==1
            StepNb=1;
        else
            StepNb=2;
        end
        for StepL=1:StepNb
            % selection of HL and BL results stocked in Var         
            if LoadDataFlag==1
                cd(K.dir.work)
                eval(sprintf('load Data%04u',BLRanks(StepL)))
                BL{StepL}=CData{ChipRank}.rank;
                eval(sprintf('load Data%04u',HLRanks(StepL)))
                HL{StepL}=CData{ChipRank}.rank;
                % Not measured genes (set to NaN in raw data) are set to -1 in both exprimental points
                % in order to prevent detection of artefactual variations
                clear CData
            else
                BL{StepL}=Data{BLRanks(StepL)}{ChipRank}.rank;
                HL{StepL}=Data{HLRanks(StepL)}{ChipRank}.rank;
            end
            NaNIndex=find(BL{StepL}==-1);
            if ~isempty(NaNIndex)
                HL{StepL}(NaNIndex)=-1;
            end
            NaNIndex=find(HL{StepL}==-1);
            if ~isempty(NaNIndex)
                BL{StepL}(NaNIndex)=-1;
            end
        end


        %Eventually recover data for the single(s) points used to complete the BL and HL pairs or replicates
        if SingleCalibPointFlag==1&LoopL==1&CompL==1
            if LoadDataFlag==1
                cd(K.dir.work)
                eval(sprintf('load Data%04u',HLCalibRank))
                HL{3}=CData{ChipRank}.rank;
                eval(sprintf('load Data%04u',BLCalibRank))
                BL{3}=CData{ChipRank}.rank;
            else
                HL{3}=Data{HLCalibRank}{ChipRank}.rank;
                BL{3}=Data{BLCalibRank}{ChipRank}.rank;
            end
            NaNIndex=find(BL{3}==-1);
            if ~isempty(NaNIndex)
                HL{3}(NaNIndex)=-1;
            end
            NaNIndex=find(HL{3}==-1);
            if ~isempty(NaNIndex)
                BL{3}(NaNIndex)=-1;
            end
        end


        if SingleCalibCurveFlag==0
            %VERIFY THAT NEEDED CALIBRATION SETS EXIST OR CONSTRUCT THEM
            if isequal(CalibType,'idem')
                if SingleCalibPointFlag==1
                    CoupleRanks(1,1)=max(BLCalibRank,BLRanks(1));
                    CoupleRanks(1,2)=min(BLCalibRank,BLRanks(1));
                else
                    CoupleRanks(1,1)=max(BLRanks(1),BLRanks(2));
                    CoupleRanks(1,2)=min(BLRanks(1),BLRanks(2));
                end
                if SingleCalibPointFlag==1
                    CoupleRanks(2,1)=max(HLCalibRank,HLRanks(1));
                    CoupleRanks(2,2)=min(HLCalibRank,HLRanks(1));
                else
                    CoupleRanks(2,1)=max(HLRanks(1),HLRanks(2));
                    CoupleRanks(2,2)=min(HLRanks(1),HLRanks(2));
                end
            elseif isequal(CalibType,'up')
                CoupleRanks(1,1)=HLRanks(1);
                CoupleRanks(1,2)=BLRanks(1);
                if SingleCalibPointFlag==1
                    CoupleRanks(2,1)=HLCalibRank;
                    CoupleRanks(2,2)=BLCalibRank;
                else
                    CoupleRanks(2,1)=HLRanks(2);
                    CoupleRanks(2,2)=BLRanks(2);
                end
            elseif isequal(CalibType,'down')
                CoupleRanks(1,1)=BLRanks(1);
                CoupleRanks(1,2)=HLRanks(1);
                if SingleCalibPointFlag==1
                    CoupleRanks(2,1)=BLCalibRank;
                    CoupleRanks(2,2)=HLCalibRank;
                else
                    CoupleRanks(2,1)=BLRanks(2);
                    CoupleRanks(2,2)=HLRanks(2);
                end
            elseif SingleCalibCurveFlag==0
                h=errordlg('Unknow CalibType %s',CalibType);
                waitfor(h)
            end


            ExistCalib(1)=0;
            ExistCalib(2)=0;
            % if calibration set is constructed on true replicates (same biological conditions) the
            % order of BLRank and HLRank does not matter but there exist only one test curve in
            % S.calib (at S.calib.testcurve{ChipRank}(min(BLRank,HLRank),max(BLRank,HLRank))

            %if calibration set is constructed on points from different biological condition
            %there exist two different test curves
            if isfield(S,'calib')
                if isequal(NormType,'standardization')
                    if isequal(CalibType,'idem')
                        CalibRank=1;
                    elseif isequal(CalibType,'up')
                        CalibRank=2;
                    elseif isequal(CalibType,'down')
                        CalibRank=3;
                    end
                elseif isequal(NormType,'quantile')
                    if isequal(CalibType,'idem')
                        CalibRank=4;
                    elseif isequal(CalibType,'up')
                        CalibRank=5;
                    elseif isequal(CalibType,'down')
                        CalibRank=6;
                    end
                end
                for CoupleL=1:2
                    if P.point.biolrank(CoupleRanks(CoupleL,2))~=P.point.biolrank(CoupleRanks(CoupleL,1))
                        try
                            if S.calib.testsurf{CalibRank}{ChipRank}(CoupleRanks(CoupleL,2),CoupleRanks(CoupleL,1))>0
                                try
                                    if S.calib.testsurf{CalibRank}{ChipRank}(CoupleRanks(CoupleL,1),CoupleRanks(CoupleL,2))>0
                                        ExistCalib(CoupleL)=1;
                                    end
                                catch                                        
                                end
                            end
                        catch
                        end
                    else
                        try
                            if S.calib.testsurf{CalibRank}{ChipRank}(min(CoupleRanks(CoupleL,2),CoupleRanks(CoupleL,1)),max(CoupleRanks(CoupleL,2),CoupleRanks(CoupleL,1)))>0
                                ExistCalib(CoupleL)=1;
                            end
                        catch                            
                        end
                    end
                end
            end
            if ExistCalib(1)==0|ExistCalib(2)==0
                %calculate calibration sets and save results in S.calib
                for CoupleL=1:2
                    if ExistCalib(CoupleL)==0
                        [RankGrid,Grid,ZVarGrid,ZVar]=noise_distribution(Data{CoupleRanks(CoupleL,1)}{ChipRank}.rank,Data{CoupleRanks(CoupleL,2)}{ChipRank}.rank,RankThreshold,CoupleRanks(CoupleL,1),CoupleRanks(CoupleL,2),ClearIndex,NormType,AnalyseType,CalibType,SizerFittingDegree,DisplayFlag);
                    end
                    %Recover results in S.calib structure

                    if isequal(NormType,'standardization')
                        if isequal(CalibType,'idem')
                            CalibRank=1;
                        elseif isequal(CalibType,'up')
                            CalibRank=2;
                        elseif isequal(CalibType,'down')
                            CalibRank=3;
                        end                            
                        %construct the test curve at a distance + one std of the mean
                        %the test curve is used to compare the dispersion of the data in a pair of replicates
                        %it allows to use the pair of replicates with the highest dispersion when a comparison is to be made
                        % between two conditions
                        TestCurve=Grid.std+ Grid.mean;                        
                    elseif isequal(NormType,'quantile')                        
                        if isequal(CalibType,'idem')
                            CalibRank=4;
                        elseif isequal(CalibType,'up')
                            CalibRank=5;
                        elseif isequal(CalibType,'down')
                            CalibRank=6;
                        end
                        %take the median for calculating the test curve
                        TestCurve=Grid.perc{5};
                    end                    
                        S.calib.testcurve{CalibRank}{CoupleRanks(CoupleL,2),CoupleRanks(CoupleL,1)}=TestCurve;
                        %normalize on the maximum possible (surface under the diagonal = 401 points * 100
                        %(max value possible) /20050) => 0 <= testsurf <= 1
                        S.calib.testsurf{CalibRank}{ChipRank}(CoupleRanks(CoupleL,2),CoupleRanks(CoupleL,1))=sum(TestCurve)/20050;
                        S.calib.mean{CalibRank}{CoupleRanks(CoupleL,2),CoupleRanks(CoupleL,1)}=Grid.mean;
                        S.calib.std{CalibRank}{CoupleRanks(CoupleL,2),CoupleRanks(CoupleL,1)}=Grid.std;
                        %recover the p-value corresponding to ZVar
                        [Temp,Temp,ZVar,Pv,MinZVar,MaxZVar]=cumul_distribution(ZVar,[],1,0);
                        if max(ZVar)<100
                            ZVar=[ZVar;100];
                            %add a p-value equal to the 1/10th of the last existing p-value
                            Pv=[Pv;Pv(end)/10];
                        end
                        %Keep p-values corresponding to interpolation values
                        RdRange=0:0.25:100;
                        Pv=interp1(ZVar,Pv,RdRange);
                        if isnan(Pv(1))
                            Pv(1)=1;
                        end                        
                        S.calib.zrd{CalibRank}{CoupleRanks(CoupleL,2),CoupleRanks(CoupleL,1)}=RdRange;
                        S.calib.pvalue{CalibRank}{CoupleRanks(CoupleL,2),CoupleRanks(CoupleL,1)}=Pv;
                        S.calib.minzrd{CalibRank}{ChipRank}(CoupleRanks(CoupleL,2),CoupleRanks(CoupleL,1))=MinZVar;
                        S.calib.maxzrd{CalibRank}{ChipRank}(CoupleRanks(CoupleL,2),CoupleRanks(CoupleL,1))=MaxZVar;
                        S.calib.zval{CalibRank}{CoupleRanks(CoupleL,2),CoupleRanks(CoupleL,1)}=ZVarGrid;                        
                    
                    if SaveFlag==1
                        cd(K.dir.stat);
                        eval(sprintf('save %s S',P.file.stat))
                    end
                end
            end            
            if ComparisonFlag==1
            % select the right calibration set (those which correspond to the pair of replicates having
            % the highest dispersion  (less reproducible pair in order not to overestimate the p-values)            
                if S.calib.stestsurf(CoupleRanks(1,2),CoupleRanks(1,1))<S.calib.stestsurf(CoupleRanks(2,2),CoupleRanks(2,1))
                    CompPv=getpv(NormType,Calib.Type,HLRanks(1),BLRanks(1),CoupleRanks(2,2),CoupleRanks(2,1));
                else
                    CompPv=getpv(NormType,Calib.Type,HLRanks(1),BLRanks(1),CoupleRanks(1,2),CoupleRanks(1,1));
                end
            end
        else %exist a single calibration set used for all comparisons
            CompPv=getpv(HLRanks(1),BLRanks(1),CoupleRanks(1,2),CoupleRanks(1,1));
        end
    end %of Comp l
end %of LoopL


