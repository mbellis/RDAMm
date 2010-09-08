%rdam(1,{[1,2;1,2];[1,2;2,1]},[1,2],[3,4],0,[0,0],'idem',[],'quantile','transcriptome',7,0,0,0,1,1,1);
%rdam(1,{[1,2;1,2];[1,2;2,1]},[1,2],[3,4],0,[0,0],'idem',[],'quantile','transcriptome',7,0,0,0,1,0,1);
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
%3&4- TGRankList and CGRankList: list of data ranks used in the comparisons
%   if TGRankList=[12,15,23], and ChipRank = 1, the Data used for TG point in position 1 is Data{12}{1}.rank
%5- LoadDataFlag = 1 indicates that each data point must be loaded when to be used (no permanent storage of Data in the memory)
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
%14- CalibUpdateFlag: indicates if noise_distribution must be run even if calibration set already exist (should be 0 if clearIndex=[]
%   unless calculus must be redone for a special reason; important if ClearIndex is not empty, because an existing calibration set can be different)
%15- CalibSaveFlag: indicates if the output of noise_distribution must be saved or not
%16- DisplayFlag: indicates if figures must be drawn or not
%17- ComparisonFlag: indicates if the comparison must be done (if not the function is simply used to construct calibration sets which are stored in S.calib)
%18- CalibSchemeFlag: indicates that a particular combination of
%    calibration set is used (e.g. in double chanel technics)


%VARARGIN PARAMETERS
% if SingleCalibPointFlag==1, HLCalibRank = varargin{1} & CGCalibRank = varargin{2} indicates the ranks of points that are
%   systematically used to construct pairs of  calibration points (HLCalibRank (CGCalibRank) must not be equal to any of
%   the points contained in TGRankList (CGRankList))
% if SingleCalibCurveFlag==1, only one calibration set is used zval=varargin{1} pv=varargin{2}

%OUTPUT PARAMETERS
% ZVar : normalized variation (from 0 to 100 for INCREASED and -100 to 0 for DECREASED)
% Ppv : product of p-values (from 0 to 1 for INCREASED and -1 to 0 for DECREASED)
% Pv : p-value of Ppv  (from 0 to 1 for INCREASED and -1 to 0 for DECREASED)
% Sensitivity (from 0 to 1 for INCREASED and -1 to 0 for DECREASED)
% TotalVar : estimate of the total number of increased and decreased probe sets

%GLOBAL VARIABLES
% Data contains data ranks (if LoadDataFlag==1, Data does not exist)
% K contains directory names
% P contains metadata (P.point : description of points, P.biol: description of biological conditions ...)
% S contains calibration sets

%VERSIONS
%V04 14-03-2010 Allows comparison between two points (without replicates)
%               Added CalibScheme to extend calibration possibilities
%V03 19-01-2010 Completely refactored version
%V02 19-01-2010 Completely refactored version but for internal usage (comparison with previous results)
%V01 09-12-2009 - Refactoring of the existing version (not complete)

function [ZVar,Pv,Ppv,Fdr,Sensitivity,TotalVar]=rdam(ChipRank,CompScheme,TGRankList,CGRankList,LoadDataFlag,RankThreshold,CalibType,ClearIndex,NormType,...
    AnalyseType,SizerFittingDegree,SingleCalibPointFlag,SingleCalibCurveFlag,CalibUpdateFlag,CalibSaveFlag,DisplayFlag,ComparisonFlag,CalibSchemeFlag,varargin)

%with SaveFlag=1
%rdam(1,{[1,2;1,2];[1,2;2,1]},[1,2],[3,4],0,[0,0],'idem',[],'quantile','transcriptome',7,0,0,0,1,0,0,0);
%exemple used in processing double chanel chip
%Create calibration sets
%rdam(1,{[1,2;1,2]},[5,6],[1,2],0,[0,0],'idem',[],'quantile','transcriptome',7,0,0,0,1,1,0,0);
%rdam(1,{[1,2;1,2]},[7,8],[3,4],0,[0,0],'idem',[],'quantile','transcriptome',7,0,0,0,1,1,0,0);
%rdam(1,{[1,2;1,2]},[1,5],[2,6],0,[0,0],'idem',[],'quantile','transcriptome',7,0,0,0,1,1,0,0);
%rdam(1,{[1,2;1,2]},[3,7],[4,8],0,[0,0],'idem',[],'quantile','transcriptome',7,0,0,0,1,1,0,0);
%Do analysis by indicating the calibration sets to be used in each comparison
%[ZVar,Pv,Ppv,Fdr,Sensitivity,TotalVar]=rdam(1,{[1,2;1,2]},[5,6],[1,2],0,[0,0],'idem',[],'quantile','transcriptome',7,0,0,0,1,1,1,0);
%[ZVar,Pv,Ppv,Fdr,Sensitivity,TotalVar]=rdam(1,{[1,2;1,2]},[7,8],[3,4],0,[0,0],'idem',[],'quantile','transcriptome',7,0,0,0,1,1,1,0);
%[ZVar,Pv,Ppv,Fdr,Sensitivity,TotalVar]=rdam(1,{[1,2,3,4;1,2,3,4]},[5,6,7,8],[1,2,3,4],0,[0,0],'idem',[],'quantile','transcriptome',7,0,0,0,1,1,1,1,{[1,1,3,3;2,2,4,4;3,3,3,3]});
%[ZVar,Pv,Ppv,Fdr,Sensitivity,TotalVar]=rdam(1,{[1,2,3,4;1,2,3,4]},[5,6,7,8],[1,2,3,4],0,[0,0],'idem',[],'quantile','transcriptome',7,0,0,0,1,1,1,1,{[1,2,3,4;5,6,7,8;3,3,3,3]});

%with CalibUpdateFlag=1 & SaveFlag=0
%rdam(1,{[1,2;1,2];[1,2;2,1]},[1,2],[3,4],0,[0,0],'idem',[],'quantile','transcriptome',7,0,0,1,0,0,0,0);
%with ComparisonFlag=1 & SaveFlag=1 & DisplayFlag=1
%rdam(1,{[1,2;1,2];[1,2;2,1]},[1,2],[3,4],0,[0,0],'idem',[],'quantile','transcriptome',7,0,0,0,1,1,1,0);
%rdam(1,{[1,2;1,2];[1,2;2,1]},[1,2],[3,4],0,[0,0],'up',[],'quantile','transcriptome',7,0,0,0,0,0,0,0);
%rdam(1,{[1,2;1,2];[1,2;2,1]},[1,2],[3,4],0,[0,0],'idem',[],'standardization','transcriptome',7,0,0,0,0,0,0,0);
%rdam(1,{[1,2;1,2];[1,2;2,1]},[1,2],[3,4],0,[0,0],'down',[],'standardization','transcriptome',7,0,0,0,0,0,0,0);

global Data K P S

%% Verifications
if SingleCalibPointFlag==1
    if nargin~=20
        h=errordlg('must have 20 parameters');
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
        CurrZVal=varargin{1};
        CurrZVarCdf=varargin{2};
    end
elseif CalibSchemeFlag==1
    if nargin~=19
        h=errordlg('must have 19 parameters');
        waitfor(h)
    else
        %CalibScheme gives the coordinates of the test curve in S.calib for
        %each comparison (e.g for CompScheme = {[1,2,3,4;1,2,3,4]} & TG=[7,8,9,10] & CG=[11,12,13,14]=>
        %CalibScheme = {[7,7,13,13;8,8,14,14;3,3,3,]} indicates that for
        %comparison between first Test point (rank 7) and the first Control
        %point (rank 11), the test curve used is located at S.calib{3}(7,8)
        CalibScheme=varargin{1};
    end
else
    if nargin~=18
        h=errordlg('must have 18 parameters');
        waitfor(h)
    end
end
if ~isequal(CalibType,'idem')&~isequal(CalibType,'down')&~isequal(CalibType,'up')
    h=errordlg(sprintf('CalibType %s unknown (must be idem, up or down)',CalibType));
    waitfor(h)
end
if ~isequal(NormType,'standardization')&~isequal(NormType,'quantile')
    h=errordlg(sprintf('NormType %s unknown (must be standardization or quantile)',NormType));
    waitfor(h)
end
if ~isequal(AnalyseType,'transcriptome')&~isequal(AnalyseType,'chipchip')
    h=errordlg(sprintf('AnalyzeType %s unknown (must be transcriptome or chipchip)',AnalyseType));
    waitfor(h)
end
if CalibSaveFlag==0
    %CalibSaveFlag=0 indicates that a new calibration set must be calculated, so update must be forced
    CalibUpdateFlag=1;
end



%% Set up RoundNb & CompNb
% For each comparison set the program starts a new round
RoundNb=length(CompScheme);
% A round is made of CompNb comparisons
CompNb=zeros(RoundNb,1);
for RoundL=1:RoundNb
    CompNb(RoundL)=size(CompScheme{RoundL},2);
end

for RoundL=1:RoundNb
    for CompL=1:CompNb(RoundL)
        %find the current pairs of BL points and HL points considered as replicate in the current comparison HL{1}vsBL{1}
        BLIndexes=zeros(2,1);
        HLIndexes=zeros(2,1);
        BLRanks=zeros(2,1);
        HLRanks=zeros(2,1);
        %The index of the first member of the current BL and HL pair
        BLIndexes(1)= CompScheme{RoundL}(2,CompL);
        HLIndexes(1)= CompScheme{RoundL}(1,CompL);
        %Recover the corresponding ranks
        BLRanks(1)=CGRankList(BLIndexes(1));
        HLRanks(1)=TGRankList(HLIndexes(1));
        %Search the second member if SingleCalibPointFlag==0
        %otherwise the second member is always the same and is known
        if SingleCalibPointFlag==0
            %The index of the second member of the current BL and HL pair
            if CompL<CompNb(RoundL)
                BLIndexes(2)=CompScheme{RoundL}(2,CompL+1);
            else
                BLIndexes(2)=CompScheme{RoundL}(2,1);
            end
            BLRanks(2)=CGRankList(BLIndexes(2));

            if CompL<CompNb(RoundL)
                HLIndexes(2)=CompScheme{RoundL}(1,CompL+1);
            else
                HLIndexes(2)=CompScheme{RoundL}(1,1);
            end
            HLRanks(2)=TGRankList(HLIndexes(2));
        end

        if SingleCalibCurveFlag==0
            if isequal(NormType,'standardization')
                if isequal(CalibType,'idem')
                    CalibRank=1;
                else
                    CalibRank=2;
                end
            elseif isequal(NormType,'quantile')
                if isequal(CalibType,'idem')
                    CalibRank=3;
                else
                    CalibRank=4;
                end
            end
            if CalibSchemeFlag==0
                CalibRanks=zeros(2,2);
                %VERIFY THAT NEEDED CALIBRATION SETS EXIST OR CONSTRUCT THEM
                if isequal(CalibType,'idem')
                    if SingleCalibPointFlag==1
                        CalibRanks(1,1)=min(HLCalibRank,HLRanks(1));
                        CalibRanks(1,2)=max(HLCalibRank,HLRanks(1));
                    else
                        CalibRanks(1,1)=min(HLRanks(1),HLRanks(2));
                        CalibRanks(1,2)=max(HLRanks(1),HLRanks(2));
                    end
                    if SingleCalibPointFlag==1
                        CalibRanks(2,1)=min(BLCalibRank,BLRanks(1));
                        CalibRanks(2,2)=max(BLCalibRank,BLRanks(1));
                    else
                        CalibRanks(2,1)=min(BLRanks(1),BLRanks(2));
                        CalibRanks(2,2)=max(BLRanks(1),BLRanks(2));
                    end
                else %up or down case
                    if SingleCalibPointFlag==1
                        CalibRanks(1,1)=HLCalibRank;
                        CalibRanks(1,2)=BLCalibRank;
                    else
                        CalibRanks(1,1)=HLRanks(1);
                        CalibRanks(1,2)=BLRanks(1);
                    end
                end

                % if calibration set is constructed on true replicates (same biological conditions) the
                % order of BLRank and HLRank does not matter but there exist only one test curve in
                % S.calib (at S.calib.testcurve{ChipRank}(min(BLRank,HLRank),max(BLRank,HLRank))

                %if calibration set is constructed on points from different biological condition
                %there exist two different test curves one constructed on increased values in the comparison
                %max(BLRanks,HLRanks) vs min(BLRanks,HLRanks) (at S.calib.testcurve{ChipRank}(min(BLRank,HLRank),max(BLRank,HLRank))
                %and one on decreased values in the comparisons
                %max(BLRanks,HLRanks) vs min(BLRanks,HLRanks) (at S.calib.testcurve{ChipRank}(max(BLRank,HLRank),min(BLRank,HLRank))

                ExistCalib(1)=0;
                ExistCalib(2)=0;                

                if SingleCalibPointFlag==1
                    CalibNb=1;
                else
                    CalibNb=2;
                end

                if isfield(S,'calib')
                    for CalibL=1:CalibNb
                        switch CalibType
                            case 'idem'
                                try
                                    if S.calib.testsurf{ChipRank}{CalibRank}(min(CalibRanks(CalibL,1),CalibRanks(CalibL,2)),max(CalibRanks(CalibL,1),CalibRanks(CalibL,2)))>0
                                        if isempty(ClearIndex)
                                            %the existing calibration set must have its clearflag set to 0
                                            if S.calib.clearflag{ChipRank}{CalibRank}(min(CalibRanks(CalibL,1),CalibRanks(CalibL,2)),max(CalibRanks(CalibL,1),CalibRanks(CalibL,2)))==0
                                                ExistCalib(CalibL)=1;
                                            end
                                        else
                                            %the existing calibration set must have its clearflag set to 1 (it's under the responsability of user to assume that it is the same
                                            %desired calibration set (otherwise CalibUpdateFlag must be set at 1)
                                            if S.calib.clearflag{ChipRank}{CalibRank}(min(CalibRanks(CalibL,1),CalibRanks(CalibL,2)),max(CalibRanks(CalibL,1),CalibRanks(CalibL,2)))==1
                                                ExistCalib(CalibL)=1;
                                            end
                                        end
                                    end
                                catch                                    
                                end
                            case 'up'
                                try
                                    if S.calib.testsurf{ChipRank}{CalibRank}(min(CalibRanks(CalibL,1),CalibRanks(CalibL,2)),max(CalibRanks(CalibL,1),CalibRanks(CalibL,2)))>0
                                        if isempty(ClearIndex)
                                            %the existing calibration set must have its clearflag set to 0
                                            if S.calib.clearindex{ChipRank}{CalibRank}(min(CalibRanks(CalibL,1),CalibRanks(CalibL,2)),max(CalibRanks(CalibL,1),CalibRanks(CalibL,2)))==0
                                                ExistCalib(CalibL)=1;
                                            end
                                        else
                                            %the existing calibration set must have its clearflag set to 1 (it's under the responsability of user to assume that it is the same
                                            %desired calibration set (otherwise CalibUpdateFlag must be set at 1)
                                            if S.calib.clearindex{ChipRank}{CalibRank}(min(CalibRanks(CalibL,1),CalibRanks(CalibL,2)),max(CalibRanks(CalibL,1),CalibRanks(CalibL,2)))==1
                                                ExistCalib(CalibL)=1;
                                            end
                                        end
                                    end
                                catch
                                end
                            case 'down'
                                try
                                    if S.calib.testsurf{ChipRank}{CalibRank}(max(CalibRanks(CalibL,1),CalibRanks(CalibL,2)),min(CalibRanks(CalibL,1),CalibRanks(CalibL,2)))>0
                                        if isempty(ClearIndex)
                                            %the existing calibration set must have its clearflag set to 0
                                            if S.calib.clearflag{ChipRank}{CalibRank}(max(CalibRanks(CalibL,1),CalibRanks(CalibL,2)),min(CalibRanks(CalibL,1),CalibRanks(CalibL,2)))==0
                                                ExistCalib(CalibL)=1;
                                            end
                                        else
                                            %the existing calibration set must have its clearflag set to 1 (it's under the responsability of user to assume that it is the same
                                            %desired calibration set (otherwise CalibUpdateFlag must be set at 1)
                                            if S.calib.clearflag{ChipRank}{CalibRank}(max(CalibRanks(CalibL,1),CalibRanks(CalibL,2)),min(CalibRanks(CalibL,1),CalibRanks(CalibL,2)))==1
                                                ExistCalib(CalibL)=1;
                                            end
                                        end
                                    end
                                catch
                                end
                        end
                    end
                end
                if ExistCalib(1)==0|(ExistCalib(2)==0&CalibNb==2)|CalibUpdateFlag==1
                    %calculate calibration sets and save results in S.calib
                    for CalibL=1:CalibNb
                        if ExistCalib(CalibL)==0|CalibUpdateFlag==1
                            [HL,BL]=getdata(ChipRank,CalibRanks(CalibL,1),CalibRanks(CalibL,2),LoadDataFlag);
                            [RankGrid,Grid,ZVarGrid,ZVar]=noise_distribution(HL,BL,RankThreshold,CalibRanks(CalibL,1),CalibRanks(CalibL,2),ClearIndex,NormType,AnalyseType,CalibType,SizerFittingDegree,DisplayFlag);
                            clear HL BL
                            %Recover results in S.calib structure
                            if isequal(NormType,'standardization')
                                %construct the test curve at a distance + one std of the mean
                                %the test curve is used to compare the dispersion of the data in a pair of replicates
                                %it allows to use the pair of replicates with the highest dispersion when a comparison is to be made
                                % between two conditions
                                TestCurve=Grid.std+ Grid.mean;
                            elseif isequal(NormType,'quantile')
                                %take the median for calculating the test curve
                                TestCurve=Grid.perc{5};
                            end
                            switch CalibType
                                %find the line and column position in S to store the calibration set
                                case 'idem'
                                    SLine=CalibRanks(CalibL,1);
                                    SColumn=CalibRanks(CalibL,2);
                                case 'up'
                                    SLine=min(CalibRanks(CalibL,1),CalibRanks(CalibL,2));
                                    SColumn=max(CalibRanks(CalibL,1),CalibRanks(CalibL,2));
                                case 'down'
                                    SLine=max(CalibRanks(CalibL,1),CalibRanks(CalibL,2));
                                    SColumn=min(CalibRanks(CalibL,1),CalibRanks(CalibL,2));
                            end
                            if isempty(ClearIndex)
                                S.calib.clearflag{ChipRank}{CalibRank}(SLine,SColumn)=uint8(0);
                            else
                                S.calib.clearflag{ChipRank}{CalibRank}(SLine,SColumn)=uint8(1);
                            end
                            S.calib.testcurve{ChipRank}{CalibRank}{SLine,SColumn}=single(TestCurve);
                            %normalize on the maximum possible (surface under the diagonal = 401 points * 100
                            %(max value possible) /20050) => 0 <= testsurf <= 1
                            S.calib.testsurf{ChipRank}{CalibRank}(SLine,SColumn)=single(sum(TestCurve)/20050);
                            S.calib.mean{ChipRank}{CalibRank}{SLine,SColumn}=single(Grid.mean);
                            S.calib.std{ChipRank}{CalibRank}{SLine,SColumn}=single(Grid.std);
                            %recover the p-value corresponding to ZVar
                            [Temp,Temp,ZVar,ZVarCdf,MinZVar,MaxZVar]=cumul_distribution(ZVar,[],1,0);
                            if max(ZVar)<100
                                ZVar=[ZVar;100];
                                %add a p-value equal to the 1/10th of the last existing p-value
                                ZVarCdf=[ZVarCdf;ZVarCdf(end)/10];
                            end
                            %Keep p-values corresponding to interpolation values
                            ZVarRange=0:0.25:100;
                            ZVarCdf=interp1(ZVar,ZVarCdf,ZVarRange);
                            if isnan(ZVarCdf(1))
                                ZVarCdf(1)=1;
                            end
                            S.calib.zvar{ChipRank}{CalibRank}{SLine,SColumn}=single(ZVarRange);
                            S.calib.cdf{ChipRank}{CalibRank}{SLine,SColumn}=single(ZVarCdf);
                            S.calib.minzvar{ChipRank}{CalibRank}(SLine,SColumn)=single(MinZVar);
                            S.calib.maxzvar{ChipRank}{CalibRank}(SLine,SColumn)=single(MaxZVar);
                            S.calib.zval{ChipRank}{CalibRank}{SLine,SColumn}=single(ZVarGrid);

                            if CalibSaveFlag==1
                                cd(K.dir.root);
                                save data Data K P S
                                %in the definitive vertion :
                                %cd(K.dir.stat)
                                %eval(sprintf('save %s S',P.file.stat))
                            else
                                CurrZVal{CalibL}=single(ZVarGrid);
                                CurrZVarCdf{CalibL}=single(ZVarCdf);
                                CurrTestSurf(CalibL)=single(sum(TestCurve)/20050);
                                CurrMinZVar(CalibL)=single(MinZVar);
                                CurrMaxZVar(CalibL)=single(MaxZVar);
                            end
                        end
                    end
                end
            end
        end
        if ComparisonFlag==1
            %Recover the data (ranks) for BL and HL
            [HL,BL]=getdata(ChipRank,HLRanks(1),BLRanks(1),LoadDataFlag);
            PsNb=length(BL);
            if CompL==1
                if RoundL==1
                    ZVar=single(zeros(PsNb,RoundNb));
                    Pv=single(ones(PsNb,RoundNb));
                    Ppv=single(ones(PsNb,RoundNb));
                end
                CurrZVar=single(zeros(PsNb,CompNb(RoundL)));
                CurrPv=single(ones(PsNb,CompNb(RoundL)));
            end

            if SingleCalibCurveFlag==0
                if CalibSchemeFlag==0
                    % select the right calibration set (those which correspond to the pair of replicates having
                    % the highest dispersion  (less reproducible pair in order not to overestimate the p-values)
                    switch CalibType
                        case 'idem'
                            if CalibSaveFlag==1
                                if S.calib.testsurf{ChipRank}{CalibRank}(CalibRanks(1,1),CalibRanks(1,2))>S.calib.testsurf{ChipRank}{CalibRank}(CalibRanks(2,1),CalibRanks(2,2))
                                    CLine=1;
                                else
                                    CLine=2;
                                end
                                CurrZVal=S.calib.zval{ChipRank}{CalibRank}{CalibRanks(CLine,1),CalibRanks(CLine,2)};
                                CurrZVarCdf=S.calib.cdf{ChipRank}{CalibRank}{CalibRanks(CLine,1),CalibRanks(CLine,2)};
                                CurrMinZVar=S.calib.minzvar{ChipRank}{CalibRank}(CalibRanks(CLine,1),CalibRanks(CLine,2));
                                CurrMaxZVar=S.calib.maxzvar{ChipRank}{CalibRank}(CalibRanks(CLine,1),CalibRanks(CLine,2));
                            else
                                if CurrTestSurf(1)>CurrTestSurf(2)
                                    CurrZVal=CurrZVal{1};
                                    CurrZVarCdf=CurrZVarCdf{1};
                                    CurrMinZVar=CurrMinZVar(1);
                                    CurrMaxZVar=CurrMaxZVar(1);
                                else
                                    CurrZVal=CurrZVal{2};
                                    CurrZVarCdf=CurrZVarCdf{2};
                                    CurrMinZVar=CurrMinZVar(2);
                                    CurrMaxZVar=CurrMaxZVar(2);
                                end
                            end
                        case {'up','down'}
                            if CalibSaveFlag==1
                                if isequal(CalibType,'up')
                                    SLine=CalibRanks(1,1);
                                    SColumn=CalibRanks(1,2);
                                else
                                    SLine=CalibRanks(1,2);
                                    SColumn=CalibRanks(1,1);
                                end
                                CurrZVal=S.calib.zval{ChipRank}{CalibRank}{SLine,SColumn};
                                CurrZVarCdf=S.calib.cdf{ChipRank}{CalibRank}{SLine,SColumn};
                                CurrMinZVar=S.calib.minzvar{ChipRank}{CalibRank}(SLine,SColumn);
                                CurrMaxZVar=S.calib.maxzvar{ChipRank}{CalibRank}(SLine,SColumn);
                            else
                                CurrZVal=CurrZVal{1};
                                CurrZVarCdf=CurrZVarCdf{1};
                                CurrMinZVar=CurrMinZVar(1);
                                CurrMaxZVar=CurrMaxZVar(1);
                            end
                    end
                else
                    SLine=CalibScheme{RoundL}(1,CompL);
                    SColumn=CalibScheme{RoundL}(2,CompL);
                    CalibRank=CalibScheme{RoundL}(3,CompL);
                    CurrZVal=S.calib.zval{ChipRank}{CalibRank}{SLine,SColumn};
                    CurrZVarCdf=S.calib.cdf{ChipRank}{CalibRank}{SLine,SColumn};
                    CurrMinZVar=S.calib.minzvar{ChipRank}{CalibRank}(SLine,SColumn);
                    CurrMaxZVar=S.calib.maxzvar{ChipRank}{CalibRank}(SLine,SColumn);
                end
            end
            if length(find(HL==BL))==length(BL)
                h=errodlg('HL (point nb %u) and BL (point nb %u) have identical values !',HLRanks(1),BLRanks(1));
                waitfor(h)
            else
                [CurrZVar(:,CompL),CurrPv(:,CompL)]=getpv(HL,BL,CurrZVal,CurrZVarCdf,CurrMinZVar,CurrMaxZVar,NormType,[],[]);
            end
        end %if ComparisonFlag==1
    end %of CompL
    if ComparisonFlag==1
        if CompNb(RoundL)>1
            [ZVar(:,RoundL),Pv(:,RoundL),Ppv(:,RoundL)]=getppv(CurrZVar,CurrPv,single(P.ppv.val{CompNb(RoundL)}),single(P.ppv.cdf{CompNb(RoundL)}),DisplayFlag);
        else
            ZVar(:,RoundL)=CurrZVar;
            Pv(:,RoundL)=CurrPv;
            Ppv(:,RoundL)=CurrPv;
        end
        if DisplayFlag==1 & RoundL==2

            figure
            IncIndex=ZVar(:,1)>0;
            DecIndex=ZVar(:,1)<0;
            subplot(2,2,1)
            plot(Data{1}{1}.rank-Data{3}{1}.rank,Data{2}{1}.rank-Data{4}{1}.rank,'k.','markersize',3)
            hold on
            plot(Data{1}{1}.rank(IncIndex)-Data{3}{1}.rank(IncIndex),Data{2}{1}.rank(IncIndex)-Data{4}{1}.rank(IncIndex),'r.','markersize',3)
            plot(Data{1}{1}.rank(DecIndex)-Data{3}{1}.rank(DecIndex),Data{2}{1}.rank(DecIndex)-Data{4}{1}.rank(DecIndex),'b.','markersize',3)
            set(gca,'xlim',[-50,50])
            set(gca,'ylim',[-50,50])
            ylabel('rank(HL2)-rank(BL2)')
            xlabel('rank(HL1)-rank(BL1)')
            title('reproducibility of rank difference')

            subplot(2,2,2)
            plot(ZVar(:,1),ZVar(:,2),'k.','markersize',3)
            hold on
            plot(ZVar(IncIndex,1),ZVar(IncIndex,2),'r.','markersize',3)
            plot(ZVar(DecIndex,1),ZVar(DecIndex,2),'b.','markersize',3)
            for i=1:7
                plot(ZVar(P.inczeroindex(P.test(i)),1),ZVar(P.inczeroindex(P.test(i)),2),sprintf('%co',P.color(i)),'markersize',10)
                plot(ZVar(P.deczeroindex(P.test(i)),1),ZVar(P.deczeroindex(P.test(i)),2),sprintf('%c+',P.color(i)),'markersize',10)
            end
            xlabel('ZVar of HL1vsBL1 & HL2vsBL2')
            ylabel('ZVar of HL2vsBL1 & HL1vsBL2')
            set(gca,'xlim',[-50,50])
            set(gca,'ylim',[-50,50])
            title('reproducibility of normalized rank difference (seven corrected points are shown)')

            subplot(2,2,3)
            plot(ZVar(:,1),log10(Pv(:,1)),'k.','markersize',3)
            hold on
            plot(ZVar(IncIndex,1),log10(Pv(IncIndex,1)),'r.','markersize',3)
            plot(ZVar(DecIndex,1),log10(Pv(DecIndex,1)),'b.','markersize',3)
            xlabel('ZVar')
            ylabel('log10(Pv)')
            set(gca,'xlim',[-50,50])
            title('Pv(Ppv)=f(ZVar) is monotonous')

            subplot(2,2,4)
            plot(Data{1}{1}.rank-Data{3}{1}.rank,ZVar(:,1),'k.','markersize',3)
            hold on
            plot(Data{1}{1}.rank(IncIndex)-Data{3}{1}.rank(IncIndex),ZVar(IncIndex(:,1)),'r.','markersize',3)
            plot(Data{1}{1}.rank(DecIndex)-Data{3}{1}.rank(DecIndex),ZVar(DecIndex(:,1)),'b.','markersize',3)
            xlabel('rank(HL1)-rank(BL1)')
            ylabel('ZVar of HL1vsBL1 & HL2vsBL2')
            set(gca,'xlim',[-50,50])
            set(gca,'ylim',[-50,50])
            title('some low rank differences are reassigned (inc<=>dec)')

        end % of if DisplayFlag==1 & RoundL==2
    end %of if ComparisonFlag==1
end %of RoundL
if ComparisonFlag==1
    %normally each round has the same number of comparisons
    %take the mean over all rounds in case it is not true
    MeanCompNb=round(mean(CompNb));
    [ZVar,Pv,Ppv,Fdr,Sensitivity,TotalVar]=result(ZVar,Pv,Ppv,MeanCompNb,DisplayFlag);
end
%% FUNCTION getdata

%INPUT PARAMETERS
%VARARGIN PARAMETERS
%GLOBAL VARIABLES

function [HL,BL]=getdata(ChipRank,HLRank,BLRank,LoadDataFlag)
global Data K
if LoadDataFlag==1
    cd(K.dir.work)
    eval(sprintf('load Data%04u',HLRank))
    HL=CData{ChipRank}.rank;
    eval(sprintf('load Data%04u',BLRank))
    BL=CData{ChipRank}.rank;
    clear CData
else
    HL=Data{HLRank}{ChipRank}.rank;
    BL=Data{BLRank}{ChipRank}.rank;
end
% Not measured genes (set to -1) in one experimental point are set to -1 in both exprimental points
% in order to prevent detection of artefactual variations
HL(BL<0)=-1;
BL(HL<0)=-1;

%% FUNCTION getpv

%INPUT PARAMETERS
%VARARGIN PARAMETERS
%GLOBAL VARIABLES

function [ZVar,Pv]=getpv(HL,BL,ZVal,HoZVarCdf,MinZVar,MaxZVar,NormType,MeanRank,StdRank)

%the calculus is made in two passes : in the first passe Increased
%variations, that is positive variations in the HL vs BL comparison,
%are treated, and in the second passe Decreased variations, that is the
%positive variations in the BL vs HL comparison.

XGrid=0:0.25:100;
YGrid=0:0.5:100;
HoZVar=0:0.25:100;
VarThreshold=0;
ZVar=zeros(length(BL),1);
Pv=ones(length(BL),1);


for TrendL=1:2
    %:::::::::::::::::::::::::::::::::::::::::::
    %CALCULUS OF THE NORMALIZED VARIATION ZVAR
    %:::::::::::::::::::::::::::::::::::::::::::
    %variation is the difference of ranks
    if TrendL==1
        Var=HL-BL;
        Rank=BL;
        Sign=1;
    else
        Var=BL-HL;
        Rank=HL;
        Sign=-1;
    end
    clear Temp
    %find positive variations (missing values (-1) are excluded from the calculus)
    UpBindex=Var>VarThreshold&BL>=0&HL>=0;
    if isequal(NormType,'standardization')
        FittedMean=interp1(XGrid,MeanRank,Rank(UpBindex));
        FittedStd=interp1(XGrid,StdRank,Rank(UpBindex));
        UpZVar=(Var(UpBindex)-FittedMean)./FittedStd;
    elseif isequal(NormType,'quantile')
        UpZVar=interp2(XGrid,YGrid,ZVal,Rank(UpBindex),Var(UpBindex));
    end
    %replace UpZVar at the righ place in ZVar
    ZVar(UpBindex)=UpZVar*Sign;



    %:::::::::::::::::::::
    %CALCULUS of P-VALUES
    %:::::::::::::::::::::

    % Selection of normalized variation which are 'positive'

    % Correct the tails of UpZVar to make them contained in the range of
    % zvar used in the interpolation of p-values
    UpZVar(UpZVar>MaxZVar)=MaxZVar;
    UpZVar(UpZVar<MinZVar)=MinZVar;
    UpPv=interp1(HoZVar,HoZVarCdf,UpZVar);

    % values >1 are replaced
    UpPv(UpPv>1)=1;
    % values <=0 are replaced by eps
    UpPv(UpPv<=0)=eps;
    Pv(UpBindex)=UpPv*Sign;
end

%% FUNCTION getppv
function [ZVar,Pv,Ppv]=getppv(ZVar,Pv,HoPpv,HoPpvCdf,DisplayFlag)
global P

CompNb=size(ZVar,2);

%Unique values for ZVar and Pv are derived
InterpZVar=abs(ZVar(:));
InterpPv=abs(Pv(:));
[InterpPv,SortOrder]=sort(InterpPv);
InterpZVar=InterpZVar(SortOrder);
InterpZVar=make_monotonous(InterpZVar,'dec');
%Reconstruct a complete set of values (Inc and Dec)
InterpPv=[-InterpPv;InterpPv];
InterpZVar=[-InterpZVar;InterpZVar];
[InterpZVar,Index]=unique(InterpZVar);
InterpPv=InterpPv(Index);
%InterpPpv is the product of ppv in case of perfect reproducibility (i.e. with Zvar identical in all the comparisons)
% => InterpPpv=Pv^CompNb
%in order to make below an estimate of Ppv by interpolation of ZVar 
InterpPpv=InterpPv.^CompNb;
[InterpZVar,Index]=unique(InterpZVar);
InterpPpv=InterpPpv(Index);


%:::::::::::::::::::::::::::::::::::::::::::::::::::::
% CALCULATE THE MEAN OF ZVAR AND THE PRODUCT OF PPV
%:::::::::::::::::::::::::::::::::::::::::::::::::::::

%BestPpv=f(ZVar) is identical in each comparison (the two first are plotted)
%BestPpv is the product of ppv in case of perfect reproducibility (Zvar are the same in all the comparisons)
% => BestPpv=Pv^CompNb
if DisplayFlag==1
    figure
    hold on
    for CompL=1:min(CompNb,length(P.color))
    plot(ZVar(:,CompL),log10(Pv(:,CompL).^2),sprintf('%c.',P.color(CompL)),'markersize',3)
    end    
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% combine the different p-values to determine the true direction of variation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[ZVar,Pv]=adjust(ZVar,Pv);
%after adjust Pv=abs(Pv)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% construct product of p-values => ppv and  calculate the corresponding p-values
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ZVar=mean(ZVar,2);
ProdPv=prod(Pv,2);
IncBindex=ZVar>0;
DecBindex=ZVar<0;
%make monotonous the function Ppv=f(ZVar) by interpolating
%Zvar on the curve InterpPpv=f(ZVar) (ZVar corresponds to the maximum
%of reproducibility) in order to get an estimate of Ppv
Ppv=interp1(InterpZVar,InterpPpv,ZVar);

if DisplayFlag==1
    %the real product of ppv
    plot(ZVar,log10(ProdPv),'g.','markersize',3)
    hold on
    %the (ZVar,ppv) of points that have been corrected
    plot(ZVar(P.inczeroindex),log10(ProdPv(P.inczeroindex)),'r+','markersize',3)
    plot(ZVar(P.deczeroindex),log10(ProdPv(P.deczeroindex)),'b+','markersize',3)
    %the seven sampled point with their (ZVar,ppv) before and after correction
    for i=1:7
        plot(P.zvartest(i,1),log10(P.ppvtest(i,1)),sprintf('%co',P.color(i)))
        plot(P.zvartest(i,2),log10(P.ppvtest(i,2)),sprintf('%c+',P.color(i)))
        plot(ZVar(P.inczeroindex(P.test(i))),log10(Ppv(P.inczeroindex(P.test(i)))),sprintf('%co',P.color(i)))
        plot(ZVar(P.deczeroindex(P.test(i))),log10(Ppv(P.deczeroindex(P.test(i)))),sprintf('%c+',P.color(i)))
        line([P.zvartest(i,1),ZVar(P.inczeroindex(P.test(i)))],[log10(P.ppvtest(i,1)),log10(Ppv(P.inczeroindex(P.test(i))))],'color',P.color(i))
        line([P.zvartest(i,2),ZVar(P.deczeroindex(P.test(i)))],[log10(P.ppvtest(i,2)),log10(Ppv(P.deczeroindex(P.test(i))))],'color',P.color(i),'linestyle',':')
    end
    set(gca,'xlim',[-50,50])
    xlabel('mean(ZVar)')
    ylabel('log10(Ppv)')
    title('Correction of Ppv=f(ZVar) and discrepencies between comparisons')
end
%reinitialise Pv to recover pv of ppv
Pv=ones(size(Pv,1),1);
MinPpv=min(HoPpv);
MaxPpv=max(HoPpv);
for TrendL=1:2
    if TrendL==1
        CurrPpv=Ppv(IncBindex);        
    else
        CurrPpv=Ppv(DecBindex);
    end
    %force CurrPpv to be of inside the range of interpolated ppv values
    CurrPpv(CurrPpv<=MinPpv)=MinPpv;
    CurrPpv(CurrPpv>MaxPpv)=MaxPpv;
    CurrPv=interp1(HoPpv,HoPpvCdf,CurrPpv);
    if TrendL==1
        Pv(IncBindex)=CurrPv;
    else
        Pv(DecBindex)=-CurrPv;
    end
end



%% FUNCTION result
function [ZVar,Pv,Ppv,Fdr,Sensitivity,TotalVar]=result(ZVar,Pv,Ppv,CompNb,DisplayFlag)
global P
Fdr=ones(size(Pv,1),1);
Sensitivity=ones(size(Pv,1),1);
RoundNb=size(ZVar,2);

%Unique values for ZVar and Ppv are derived
%in order to make an interpolation of ppv below
[InterpPpv,Index,Temp]=unique(Ppv(:));
InterpZVar=ZVar(:);
InterpZVar=InterpZVar(Index);
[InterpZVar,Index,Temp]=unique(InterpZVar);
InterpPpv=InterpPpv(Index);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% combine the different p-values to determine the true direction of variation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if RoundNb>1
    [ZVar,Pv]=adjust(ZVar,Pv);
end
%Pv = abs(Pv)
%takes ZVar as the best statistics
ZVar=mean(ZVar,2);
%and don't use the calculated ppv
MeanPpv=mean(Ppv,2);
%but the interpolated one instead
ZVar(ZVar<min(InterpZVar))=min(InterpZVar);
ZVar(ZVar>max(InterpZVar))=max(InterpZVar);
Ppv=interp1(InterpZVar,InterpPpv,ZVar);
%don't use the mean of pv(ppv)
MeanPv=mean(Pv,2);
%but the interpolated pv instead
CurrPpv=Ppv;
CurrPpv(Ppv<min(single(P.ppv.val{CompNb})))=min(single(P.ppv.val{CompNb}));
CurrPpv(Ppv>max(single(P.ppv.val{CompNb})))=max(single(P.ppv.val{CompNb}));
Pv=interp1(single(P.ppv.val{CompNb}),single(P.ppv.cdf{CompNb}),CurrPpv);

IncBindex=ZVar>0;
DecBindex=ZVar<0;

if DisplayFlag==1
    figure
    subplot(1,2,1)
    plot(ZVar,log10(MeanPpv),'b+')
    hold on
    plot(ZVar,log10(Ppv),'g.')
    xlabel('final ZVar')
    ylabel('final product of p-values (ppv)')
    title('correction of product of p-values (blue (green) before (after) correction)')
    subplot(1,2,2)
    plot(ZVar,log10(abs(MeanPv)),'b+')
    hold on
    plot(ZVar,log10(abs(Pv)),'g.');
    line([0,0],[0,1],'color','k')
    xlabel('final ZVar')
    ylabel('final p-value of ppv')
    title('correction of p-values (blue (green) before (after) correction)')
end
%forces to display Fdr, Sensitivity curves
DisplayFlag=1;
TotalVar=zeros(1,2);
for TrendL=1:2
    if TrendL==1
        CurrBindex=IncBindex;
        if DisplayFlag==1
            h=figure;
        else
            h=0;
        end
        CurrColor='r';
        Title='INCREASED';
        SubPos=1;
        Sign=1;
    else
        CurrBindex=DecBindex;
        CurrColor='b';
        Title='DECREASED';
        SubPos=2;
        Sign=-1;
    end
    CurrPpv=Ppv(CurrBindex);
    CurrPv=Pv(CurrBindex);
    CurrZVar=ZVar(CurrBindex);

    [SortedPpv DirectSort]=sort(CurrPpv);
    SortedPv=CurrPv(DirectSort);
    [Temp ReverseSort]=sort(DirectSort);
    SortedZVar=CurrZVar(DirectSort);
    

    [Temp,PpvCdf,Temp,Temp,Temp,Temp]=cumul_distribution(CurrPpv,0,0,0);

    [CurrFdr,CurrSensitivity,TruePosNb]= fdr(abs(SortedZVar),SortedPpv,abs(SortedPv),PpvCdf,DisplayFlag,h,CurrColor,Title,SubPos);
    CurrFdr=CurrFdr(ReverseSort);
    Fdr(CurrBindex)=CurrFdr*Sign;
    CurrSensitivity=CurrSensitivity(ReverseSort);
    Sensitivity(CurrBindex)=CurrSensitivity*Sign;
    Pv(CurrBindex)=Pv(CurrBindex)*Sign;
    Ppv(CurrBindex)=Ppv(CurrBindex)*Sign;
    TotalVar(TrendL)=TruePosNb;
end

%% FUNCTION adjust
function [ZVar,Pv]=adjust(ZVar,Pv)
global P

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% determine the main variation trend by using mean of log(p-values)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CompNb=size(ZVar,2);
MemZVar=ZVar;
%use absolute values of p-values
Pv=abs(Pv);
MemPv=Pv;
%as we use mean of log(Pv) to determine the main trend, LogPv must be at zero by default
LogPv=zeros(size(Pv));
for CompL=1:CompNb
    IncBindex=ZVar(:,CompL)>0;
    %use positive values for increased
    LogPv(IncBindex,CompL)=-log10(Pv(IncBindex,CompL));
    DecBindex=ZVar(:,CompL)<0;
    %use negative values for decreased
    LogPv(DecBindex,CompL)=log10(Pv(DecBindex,CompL));
end
MeanLogPv=mean(LogPv,2);
IncBindex=mean(MeanLogPv,2)>0;
DecBindex=mean(MeanLogPv,2)<0;
InvBindex=mean(MeanLogPv,2)==0;
%clear LogPv

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% correct ZVar and Pv for discrepencies
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%recover positions that have to be changed in at least one comparison
IncZeroBindex=zeros(size(Pv,1),1);
DecZeroBindex=zeros(size(Pv,1),1);
%correct all position where the local trend does not corresponds to the actual trend
%as determined with the mean of log(pv)
for CompL=1:CompNb
    CurrBindex=ZVar(:,CompL)>0;
    ZVar(CurrBindex~=IncBindex,CompL)=0;
    IncZeroBindex=IncZeroBindex|CurrBindex~=IncBindex;
    Pv(CurrBindex~=IncBindex,CompL)=1;
    CurrBindex=ZVar(:,CompL)<0;
    ZVar(CurrBindex~=DecBindex,CompL)=0;
    DecZeroBindex=DecZeroBindex|CurrBindex~=DecBindex;
    Pv(CurrBindex~=DecBindex,CompL)=1;
    ZVar(InvBindex,CompL)=0;
    Pv(InvBindex,CompL)=1;
end



P.inczeroindex=find(IncZeroBindex);
P.deczeroindex=find(DecZeroBindex);
%Select seven values to be displayed among the changed positions
ZeroNb=min(length(P.inczeroindex),length(P.deczeroindex));
if ZeroNb>=7
    P.test=round(1:ZeroNb/7:ZeroNb);
    P.zvartest=[mean(MemZVar(P.inczeroindex(P.test),:),2),mean(MemZVar(P.deczeroindex(P.test),:),2)];
    P.ppvtest=[prod(MemPv(P.inczeroindex(P.test),:),2),prod(MemPv(P.deczeroindex(P.test),:),2)];
    %display the values before and after correction
    [LogPv(P.inczeroindex(P.test),:),MeanLogPv(P.inczeroindex(P.test),:),MemZVar(P.inczeroindex(P.test),:),ZVar(P.inczeroindex(P.test),:)]
    [LogPv(P.deczeroindex(P.test),:),MeanLogPv(P.deczeroindex(P.test),:),MemZVar(P.deczeroindex(P.test),:),ZVar(P.deczeroindex(P.test),:)]
end