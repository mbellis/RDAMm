%========
%FILL_S 
%========
%
% FILL_S fills the S strcture
%
% INPUT PARAMETERS
%
% 1  NormType: either 'standardization' (original method in RDAM = var - mean /std)
%              or 'quantile' (more general method based on percentile curves)
% 2 CalibType: indicates the type of couple of two replicates used to calculate
%              the calibration set for the current comparison.
%              The couple noted {G1,G2}, which is passed as parameter, indicates the
%              G1 and G2 will be identified to the HL and BL lines, respectively, in the 
%              noise_distribution function.
%              If (TG1,TG2) and (CG1,CG2) are the two couples of replicates for the the 
%              TG and CG condition, respectively, we have three possibilities:
%                CalibType='idem' =>  {TG1,TG2} or {CG1,CG2} is used,
%                CalibType='up' =>  {TG1,CG1} or {TG2,CG2} is used,
%                CalibType='down' => {CG1,TG1} or {CG2,TG2} is used.
%              A rule allows to choose the right couple from the two existing ones
%             (FirstCouple, SndCouple).
% 3      Grid: A regular grid ([0:0.25:100]) on which Zvar and cdf(ZVar) are interpolated
% 4      ZVar: normalized variation (from 0 to 100 for INCREASED and -100 to 0 for DECREASED)
% 5   ResRank: Calibration set is stored in S and saved as sprintf(CalibSet_02u,ResRank)
%              (allows to disconnect construction of points dendrogram from calculus 
%              on biological conditions which allows to add easily new points).
% 6  SaveFlag: indicates if S must be saved


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

function fill_s(NormType,CalibType,ClearIndexFlag,CalibRanks,Grid,ZVar,ZVarGrid,ResRank,SaveFlag)
global P S

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

%Recover results in S{CalibRank} structure
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
        SLine=CalibRanks(1);
        SColumn=CalibRanks(2);
    case 'up'
        SLine=min(CalibRanks(1),CalibRanks(2));
        SColumn=max(CalibRanks(1),CalibRanks(2));
    case 'down'
        SLine=max(CalibRanks(1),CalibRanks(2));
        SColumn=min(CalibRanks(1),CalibRanks(2));
end
if isempty(S) | length(S)<CalibRank
    S{CalibRank}=[];
    Pos=1;
    S{CalibRank}.position(Pos,:)=[SLine,SColumn,Pos];
else
    if isempty(S{CalibRank})
        Pos=1;
        S{CalibRank}.position(Pos,:)=[SLine,SColumn,Pos];
    else
        Pos=find(S{CalibRank}.position(:,1)==SLine&S{CalibRank}.position(:,2)==SColumn);
    end
end
if isempty(Pos)
    Pos=size(S{CalibRank}.position,1)+1;
    S{CalibRank}.position(Pos,:)=[SLine,SColumn,Pos];
end
if ClearIndexFlag==0
    S{CalibRank}.clearFlag(Pos,1)=uint8(0);
else
    S{CalibRank}.clearFlag(Pos,1)=uint8(1);
end
S{CalibRank}.testCurve{Pos,1}=single(TestCurve);
%normalize on the maximum possible (surface under the diagonal = 401 points * 100
%(max value possible) /20050) => 0 <= testsurf <= 1
S{CalibRank}.testSurf(Pos,1)=single(sum(TestCurve)/20050);
S{CalibRank}.mean{Pos,1}=single(Grid.mean);
S{CalibRank}.std{Pos,1}=single(Grid.std);
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
S{CalibRank}.zvar{Pos,1}=single(ZVarRange);
S{CalibRank}.cdf{Pos,1}=single(ZVarCdf);
S{CalibRank}.minZVar(Pos,1)=single(MinZVar);
S{CalibRank}.maxZVar(Pos,1)=single(MaxZVar);
S{CalibRank}.zval{Pos,1}=single(ZVarGrid);

if SaveFlag==1
    cd(P.dir.data)
    eval(sprintf('save CalibSet_%02u S',ResRank))
end