%&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% FUNCTION getpv
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% Find which variations are statistically significative when comparing two groups
% points, each point being a series of rank values. One group is called the baseline (BL)
% and the other is called the highline (HL). The comparison is HL vs BL.  Increased (decreased) probesets are those
% which are higher (lower) in the highline.

%INPUT PARAMETERS

%VARARGIN PARAMETERS
% if SingleCalibPointFlag==1, HLCalibRank = varargin{1} & BLCalibRank = varargin{2} indicates the ranks of points that are

%GLOBAL VARIABLES


function getpv(HL,BL,ZVal,Pv,MeanRank,StdRank,NormType)

%��������������������������������������������������������������������������������������
% Calculus of variations for the existing comparisons in one of the comparison type (II & X)
%��������������������������������������������������������������������������������������

%calculate sorted, in phase ranks and indexes necessary to make the reverse sorting and
%the sorted to phased sorting
[BLSorted,Phased_HL,BLDirectSort,BLReverseSort,Temp]=Sort2Series(BL,HL);
[HLSorted,Phased_BL,HLDirectSort,HLReverseSort,Temp]=Sort2Series(HL,BL);
clear Temp

XGrid=0:0.25:100;
VarThreshold=0;
BLVar=HL-BL;
HLVar=BL-HL;

BLSortedVar=BLVar(BLDirectSort);
HLSortedKeptZvar=HLVar(HLDirectSort);

IncBindex=BLSortedVar>VarThreshold;
if isequal(NormType,'standardization')
    BLFittedMean=interp1(XGrid,MeanRank,BLSorted);
    BLFittedStd=interp1(XGrid,StdRank,BLSorted);
    BLSortedZVar=(BLSortedVar-BLFittedMean)./BLFittedStd;
elseif isequal(NormType,'quantile')
    BLSortedZVar=interp2(XGrid,YGrid,ZVal,BLSorted,BLSortedVar);
end

BLSortedKeptZvar=zeros(length(BL),1);
BLSortedKeptZvar(IncBindex)=BLSortedZVar(IncBindex);
BLZVar=BLSortedKeptZvar(BLReverseSort);

PlusBindex=BL>=0&HL>=0;
IncBindex=(BLVar>VarThreshold)&PlusBindex;

eval(['VarBL',sprintf('%02.f',CompLoop),'=BLZVar;'])


if P.flag.zscore==1
    BLVar=BLZVar;
end
eval(['BLVar',Comp_rank,IIX_rank,'=BLVar;']);
eval(['BL',Comp_rank,'=BLSorted(BLReverseSort);']);
eval(['BL',Comp_rank,IIX_rank,'=BLSorted(BLReverseSort);']);
eval(['BLZVar',Comp_rank,'=BLZVar;']);
eval(['BLZVar',Comp_rank,IIX_rank,'=BLZVar;']);
eval(['IncBindex',Comp_rank,'=IncBindex;']);
eval(['IncBindex',Comp_rank,IIX_rank,'=IncBindex;']);
eval(['VarThreshold',Comp_rank,IIX_rank,'=VarThreshold;'])
clear BLReverseSort
clear BLVar
clear BLSorted





%P.comp.decbindex=Dec_bindex;
% taking ~IncBindex incorporate zero variation in Dec_bindex
% and produces discontinuity in observed_plevel
%Dec_bindex=~IncBindex;
DecBindex=HLSortedKeptZvar>VarThreshold;
if isequal(NormType,'standardization')
    HLFittedMean=interp1(XGrid,MeanRank,HLSorted);
    HLFittedStd=interp1(XGrid,StdRank,HLSorted);
    HLSortedZVar=(HLSortedKeptZvar-HLFittedMean)./HLFittedStd;
elseif isequal(NormType,'quantile')
    HLSortedZVar=interp2(XGrid,YGrid,ZVal,HLSorted,HLSortedKeptZvar);
end

HLSortedKeptZvar=zeros(size(HL,1),1);
HLSortedKeptZvar(DecBindex)=HLSortedZVar(DecBindex);
HLZVar=HLSortedKeptZvar(HLReverseSort);

DecBindex=(HLVar>VarThreshold)&PlusBindex;
eval(['VarHL',sprintf('%02.f',CompLoop),'=HLZVar;'])


if P.flag.zscore==1
    HLVar=HLZVar;
end
eval(['HLVar',Comp_rank,IIX_rank,'=HLVar;']);
eval(['HLZVar',Comp_rank,'=HLZVar;']);
eval(['HLZVar',Comp_rank,IIX_rank,'=HLZVar;']);
eval(['HL',Comp_rank,'=HLSorted(HLReverseSort);']);
eval(['HL',Comp_rank,IIX_rank,'=HLSorted(HLReverseSort);']);
eval(['DecBindex',Comp_rank,'=DecBindex;']);
eval(['DecBindex',Comp_rank,IIX_rank,'=DecBindex;']);
eval(['VarThreshold',Comp_rank,IIX_rank,'=VarThreshold;'])
clear HLReverseSort
clear HLVar
clear HLSorted



%}}}
%------------------------------------
%------------------------------------
%{{{  CALCULUS of PVALUES
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{{{ Increased genes
%++++++++++++++++++++++++++++++++++++++
%{{{ Interpolation of pv(zVar) , computation of cdf(zVar)

colori=['rm'];
% Selection of normalized variation which are 'positive'
eval(['Selection_bindex_BL',Comp_rank,'=IncBindex',Comp_rank,IIX_rank,';']);
Increased_number=size(find(IncBindex),1);

FNormalized_act_var_BL=BLZVar(IncBindex);
Interpolated_variation=FNormalized_act_var_BL;
Max_norm_var_bindex=Interpolated_variation>Max_norm_var_calib_BL;
Min_norm_var_bindex=Interpolated_variation<Min_norm_var_calib_BL;
Interpolated_variation(Max_norm_var_bindex)=Max_norm_var_calib_BL;
Interpolated_variation(Min_norm_var_bindex)=Min_norm_var_calib_BL;
FExpected_plevel=interp1(UFiltered_norm_var_calib_BL,UFiltered_plevel_calib_BL,...
    Interpolated_variation);

% values >1 are replaced
Sup_one_index=find(FExpected_plevel>1);
FExpected_plevel(Sup_one_index)=1;
% values <=0 are replaced by eps
Inf_zero_index=find(FExpected_plevel<=0);
FExpected_plevel(Inf_zero_index)=eps;


% don't pass the P.par.varthresh parameter at this step
% because we are interesting to recover a not filtered Cum_dist with the
% same size than Expected_plevel
%'Observed_plevel_BL'
Graphical_control=0;
[B,FSorted_observed_plevel,Temp,Temp,Temp_min,Temp_max]...
    =CUMDISTPOS(FNormalized_act_var_BL,0,Graphical_control);
clear temp
%}}}
%++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++
%{{{ Graphical evaluation of total change
%[X ,SortIndex]=sort(Interpolated_variation);
%Cdf0=FExpected_plevel(SortIndex);
%FigTitle='Individual comparisons';
%IIXRank=sprintf('%.f',IIXLoop);
%CompRank=sprintf('%.f',CompLoop);
%SubTitle=['Inc ',IIXRank,'-',CompRank];
%LineNb=IIXLoopLimit;
%ColNb=CompLoopLimit*2;
%SubRank=(IIXLoop-1)*ColNb+(CompLoop-1)*2+1;

%}}}
%++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++
%{{{ Variable assignation (FObserved_plevelBL...)

% values >1 are replaced
Sup_one_index=find(FSorted_observed_plevel>1);
FSorted_observed_plevel(Sup_one_index)=1;
% replace zeros by eps
Inf_zero_index=find(FSorted_observed_plevel<=0);
FSorted_observed_plevel(Inf_zero_index)=eps;

[FSorted_norm_act_var_BL  Direct_sort_index]=sort(FNormalized_act_var_BL);
Indirect_sort_index=flipud(Direct_sort_index);
[temp Inverse_sort_index]=sort(Indirect_sort_index);
FSorted_expected_plevel=FExpected_plevel(Indirect_sort_index);

FObserved_plevel_BL{CompLoop}=FSorted_observed_plevel(Inverse_sort_index);
FExpected_plevel_BL{CompLoop}=FSorted_expected_plevel(Inverse_sort_index);
[FSorted_norm_act_variation Temp_index]=sort(BLZVar(IncBindex)); %}}}
%clear BLZVar
clear FSorted_norm_act_var_BL
clear FObserved_plevel_BL
clear FSorted_norm_act_variation
clear temp

%++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++
%{{{ Graphics
jump_display=1;
if jump_display==0
    Fig_title='Increased (red)';
    FigH=figure;
    [Error1,Error2,Sensitivity,Specificity,Chi2_probability,TP_Tnb,Max_useful_FRP]=STAT_CALIB_PREV(Increased_number,...
        FNormalized_act_var_BL(Indirect_sort_index),FSorted_expected_plevel,...
        FSorted_observed_plevel,1,FigH,'r',Fig_title,FSorted_observed_plevel,1);
end
%}}}
%++++++++++++++++++++++++++++++++++++++
%}}}
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{{{ Decreased genes
%++++++++++++++++++++++++++++++++++++++
%{{{  Interpolation of pv(zVar) , computation of cdf(zVar)
colord=['bc'];
% Selection of normalized variation which are 'positive'
eval(['Selection_bindex_HL',Comp_rank,'=DecBindex',Comp_rank,IIX_rank,';']);
FNormalized_act_var_HL=HLZVar(DecBindex);

Decreased_number=size(find(DecBindex),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation of plevel
%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%Interpolated_variation=HLZVar(Selection_bindex);
Interpolated_variation=HLZVar(DecBindex);
Max_norm_var_bindex=Interpolated_variation>Max_norm_var_calib_HL;
Min_norm_var_bindex=Interpolated_variation<Min_norm_var_calib_HL;
Interpolated_variation(Max_norm_var_bindex)=Max_norm_var_calib_HL;
Interpolated_variation(Min_norm_var_bindex)=Min_norm_var_calib_HL;
clear Max_norm_var_bindex
clear Min_norm_var_bindex

FExpected_plevel=interp1(UFiltered_norm_var_calib_HL,UFiltered_plevel_calib_HL,...
    Interpolated_variation);
clear Interpolated_variation

% values >1 are replaced
Sup_one_index=find(FExpected_plevel>1);
FExpected_plevel(Sup_one_index)=1;

% values <=0 are replaced by eps
Inf_zero_index=find(FExpected_plevel<=0);
size(Inf_zero_index);
FExpected_plevel(Inf_zero_index)=eps;

Graphical_control=0;
[Temp,FSorted_observed_plevel,Temp,Temp,Temp_min,Temp_max]...
    =CUMDISTPOS(FNormalized_act_var_HL,0,Graphical_control);
clear Temp
Graphical_control=0;
%}}}
%++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++
%{{{ Graphical evaluation of total change
%[X ,SortIndex]=sort(Interpolated_variation);
%Cdf0=FExpected_plevel(SortIndex);
%FigTitle='Individual comparisons';
%IIXRank=sprintf('%.f',IIXLoop);
%CompRank=sprintf('%.f',CompLoop);
%SubTitle=['Dec ',IIXRank,'-',CompRank];
%LineNb=IIXLoopLimit;
%ColNb=CompLoopLimit*2;
%SubRank=(IIXLoop-1)*ColNb+(CompLoop-1)*2+2;
%TOTALCHANGE(X,Cdf0,1-FSorted_observed_plevel,FigTitle,SubTitle,LineNb,ColNb,SubRank)
%}}}
%++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++
%{{{ Variable assignation (FObserved_plevelBL...)
% values >1 are replaced
Sup_one_index=find(FSorted_observed_plevel>1);
FSorted_observed_plevel(Sup_one_index)=1;

% replace zeros by eps
Inf_zero_index=find(FSorted_observed_plevel<=0);
size(Inf_zero_index);
FSorted_observed_plevel(Inf_zero_index)=eps;


[FSorted_norm_act_var_HL  Direct_sort_index]=sort(FNormalized_act_var_HL);
Indirect_sort_index=flipud(Direct_sort_index);
[temp Inverse_sort_index]=sort(Indirect_sort_index);
FSorted_expected_plevel=FExpected_plevel(Indirect_sort_index);
FObserved_plevel_HL{CompLoop}=FSorted_observed_plevel(Inverse_sort_index);
FExpected_plevel_HL{CompLoop}=FSorted_expected_plevel(Inverse_sort_index);
[FSorted_norm_act_variation Temp_index]=sort(HLZVar(DecBindex));
%clear HLZVar
clear FSorted_norm_act_var_HL
clear FObserved_plevel_HL
clear FSorted_norm_act_variation
clear temp
%}}}
%++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++
%{{{ Graphics
jump_display=1;
if jump_display==0
    if Super_loop==1
        Fig_title='Comparison II : increased (red and magenta) and decreased (blue and cyan)';
    else
        Fig_title='Comparison X : increased (red and magenta) and decreased (blue and cyan)';
    end

    [Error1,Error2,Sensitivity,Specificity,Chi2_probability,TP_Tnb,Max_useful_FRP]=STAT_CALIB_PREV(Decreased_number,...
        FNormalized_act_var_HL(Indirect_sort_index),FSorted_expected_plevel,...
        FSorted_observed_plevel,1,F.gh.prev.stat1,colord(Comparison_loop),Fig_title,2);
end
%}}}
%++++++++++++++++++++++++++++++++++++++
%}}}
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{{{ MULTI (fill P.comp.variation, .plevel, .bindex_inc, .bindex_dec)
if CompLoop==1 &  IIXLoop==1
    %++++++++++++++++++++++++++++++++++++++
    %{{{ Creation of P.comp variable
    % which are permanent inside CompLoop
    % and volatile outside
    % allow to register all signals
    Structure_size=size(BL,1);
    First_dim=Structure_size;
    Second_dim=max(P.par.blgrpnb,P.par.hlgrpnb);
    %%%Third_dim=IIXLoopLimit;
    if IIXLoopLimit==1
        IIX_size=2;
    else
        IIX_size=IIXLoopLimit;
    end
    if P.par.big==0
        P.comp.blgroup.signal=zeros(First_dim,Second_dim);
        P.comp.blgroup.cdfpos=P.comp.blgroup.signal;
        P.comp.hlgroup.signal=zeros(First_dim,Second_dim);
        P.comp.hlgroup.cdfpos=P.comp.hlgroup.signal;
        P.comp.variation=zeros(Structure_size,CompLoopLimit,IIX_size);
        P.comp.plevel=ones(Structure_size,CompLoopLimit,IIX_size);
        P.comp.var=zeros(Structure_size,IIX_size);
        P.comp.corr_variation=ones(Structure_size,CompLoopLimit,IIX_size)*VarThreshold;
        P.comp.corr_plevel=ones(Structure_size,CompLoopLimit,IIX_size);
        P.comp.product=ones(Structure_size,IIX_size);
        P.comp.corr_nb=ones(Structure_size,IIX_size);
        P.comp.pvalue=ones(Structure_size,IIX_size);
        P.comp.pvalue_cdf=zeros(Structure_size,IIX_size);
        P.comp.single_outlier=zeros(Structure_size,CompLoopLimit,IIX_size);
        P.comp.double_outlier=zeros(Structure_size,CompLoopLimit,IIX_size);
    else
        P.comp.var=zeros(Structure_size,IIX_size);
        P.comp.plevel=ones(Structure_size,CompLoopLimit,IIX_size);
        P.comp.corr_plevel=ones(Structure_size,CompLoopLimit,IIX_size);
    end
    %}}}
    %++++++++++++++++++++++++++++++++++++++
end
if IIXLoop==1
    %++++++++++++++++++++++++++++++++++++++
    %{{{ Fill signal & cdfpos in BL(HL)_group
    %P.comp.blgroup.signal(:,CompLoop)=BL;
    if P.par.big==0
        P.comp.blgroup.signal(:,CompLoop)=interp1(K.chip.ref.rank,K.chip.ref.signal,BL);
        P.comp.blgroup.cdfpos(:,CompLoop)=BL;
        %P.comp.hlgroup.signal(:,CompLoop)=HL;
        P.comp.hlgroup.signal(:,CompLoop)=interp1(K.chip.ref.rank,K.chip.ref.signal,HL);
        P.comp.hlgroup.cdfpos(:,CompLoop)=HL;
    end
    %}}}
    %++++++++++++++++++++++++++++++++++++++
end
%++++++++++++++++++++++++++++++++++++++
%{{{ Comp assignation: .variation, .plevel, .bindex_inc, .bindex_dec
if P.par.big==0
    P.comp.variation(IncBindex,CompLoop,IIXLoop)=BLZVar(IncBindex);
    P.comp.variation(DecBindex,CompLoop,IIXLoop)=HLZVar(DecBindex);
    P.comp.plevel(IncBindex,CompLoop,IIXLoop)=FExpected_plevel_BL{CompLoop};
    clear FExpected_plevel_BL
    P.comp.plevel(DecBindex,CompLoop,IIXLoop)=FExpected_plevel_HL{CompLoop};
    clear FExpected_plevel_HL
    P.comp.bindex_inc{CompLoop,IIXLoop}=IncBindex;
    P.comp.bindex_dec{CompLoop,IIXLoop}=DecBindex;
else
    P.comp.var(IncBindex,CompLoop,IIXLoop)=BLZVar(IncBindex);
    P.comp.var(DecBindex,CompLoop,IIXLoop)=HLZVar(DecBindex);
    clear BLZVar
    clear HLZVar
    P.comp.bindex_inc{CompLoop,IIXLoop}=IncBindex;
    clear IncBindex
    P.comp.bindex_dec{CompLoop,IIXLoop}=DecBindex;
    clear DecBindex
    P.comp.plevel(P.comp.bindex_inc{CompLoop,IIXLoop},CompLoop,IIXLoop)=FExpected_plevel_BL{CompLoop};
    FExpected_plevel_BL{CompLoop}=[];
    clear FExpected_plevel_BL
    P.comp.plevel(P.comp.bindex_dec{CompLoop,IIXLoop},CompLoop,IIXLoop)=FExpected_plevel_HL{CompLoop};
    FExpected_plevel_HL{CompLoop}=[];
    clear FExpected_plevel_HL
end

�������        end % of CompLoop
        
        
        if ~isequal(P.par.comptype,'clustering')&Continue==1
            
            %====================================
            %{{{ Treatment of MULTI
            if IIXLoop==1
                %------------------------------------
                %{{{ Creation of P.comp.BL(HL)_mean_signal, std_signal, v_signal, mean_cdfpos, std_cdf, cv_cdf
                %cdfpos = rank
                if P.par.big==0
                    Structure_size=size(P.comp.blgroup.signal,1);
                    P.comp.BL_mean_signal=zeros(Structure_size,IIXLoopLimit);
                    P.comp.BL_std_signal=P.comp.BL_mean_signal;
                    P.comp.BL_cv_signal=P.comp.BL_mean_signal;
                    P.comp.BL_mean_cdfpos=P.comp.BL_mean_signal;
                    P.comp.BL_std_cdfpos=P.comp.BL_mean_signal;
                    P.comp.BL_cv_cdfpos=P.comp.BL_mean_signal;
                    
                    Structure_size=size(P.comp.hlgroup.signal,1);
                    P.comp.HL_mean_signal=zeros(Structure_size,IIXLoopLimit);
                    P.comp.HL_std_signal=P.comp.HL_mean_signal;
                    P.comp.HL_cv_signal=P.comp.HL_mean_signal;
                    P.comp.HL_mean_cdfpos=P.comp.HL_mean_signal;
                    P.comp.HL_std_cdfpos=P.comp.HL_mean_signal;
                    P.comp.HL_cv_cdfpos=P.comp.HL_mean_signal;                
                end
                %}}}
                %------------------------------------
            end
            %------------------------------------
            %{{{ Fill Comp
            %cdfpos = rank
            if P.par.big==0
                P.comp.BL_mean_signal(:,IIXLoop)=mean(P.comp.blgroup.signal(:,1:P.par.blgrpnb),2);
                P.comp.BL_std_signal(:,IIXLoop)=std(P.comp.blgroup.signal(:,1:P.par.blgrpnb),1,2);
                P.comp.BL_cv_signal(:,IIXLoop)=(P.comp.BL_std_signal(:,IIXLoop)*100)./P.comp.BL_mean_signal(:,IIXLoop);
                P.comp.BL_mean_cdfpos(:,IIXLoop)=mean(P.comp.blgroup.cdfpos(:,1:P.par.blgrpnb),2);
                P.comp.BL_std_cdfpos(:,IIXLoop)=std(P.comp.blgroup.cdfpos(:,1:P.par.blgrpnb),1,2);
                P.comp.BL_cv_cdfpos(:,IIXLoop)=(P.comp.BL_std_cdfpos(:,IIXLoop)*100)./P.comp.BL_mean_cdfpos(:,IIXLoop);
                
                P.comp.HL_mean_signal(:,IIXLoop)=mean(P.comp.hlgroup.signal(:,1:P.par.hlgrpnb),2);
                P.comp.HL_std_signal(:,IIXLoop)=std(P.comp.hlgroup.signal(:,1:P.par.hlgrpnb),1,2);
                P.comp.HL_cv_signal(:,IIXLoop)=(P.comp.HL_std_signal(:,IIXLoop)*100)./P.comp.HL_mean_signal(:,IIXLoop);
                P.comp.HL_mean_cdfpos(:,IIXLoop)=mean(P.comp.hlgroup.cdfpos(:,1:P.par.hlgrpnb),2);
                P.comp.HL_std_cdfpos(:,IIXLoop)=std(P.comp.hlgroup.cdfpos(:,1:P.par.hlgrpnb),1,2);
                P.comp.HL_cv_cdfpos(:,IIXLoop)=(P.comp.HL_std_cdfpos(:,IIXLoop)*100)./P.comp.HL_mean_cdfpos(:,IIXLoop);
            end
            %}}}
            %------------------------------------
            %------------------------------------
            %{{{ treatment of discrepencies
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %{{{ Direction rule (mean(+log(1/pv),-log(1/pv))
            
            Plevel=P.comp.plevel(:,:,IIXLoop);
            for Plevel_loop=1:CompLoopLimit
                Log_plevel=log(1./Plevel(P.comp.bindex_inc{Plevel_loop,IIXLoop},Plevel_loop));
                Infinite_bindex=~isfinite(Log_plevel);
                Log_plevel(Infinite_bindex)=0;
                Plevel(P.comp.bindex_inc{Plevel_loop,IIXLoop},Plevel_loop)=Log_plevel;
                
                Log_plevel=-log(1./Plevel(P.comp.bindex_dec{Plevel_loop,IIXLoop},Plevel_loop));
                Infinite_bindex=~isfinite(Log_plevel);
                Log_plevel(Infinite_bindex)=0;
                Plevel(P.comp.bindex_dec{Plevel_loop,IIXLoop},Plevel_loop)=Log_plevel;
                
                Plevel(~P.comp.bindex_inc{Plevel_loop,IIXLoop}&~P.comp.bindex_dec{Plevel_loop,IIXLoop},Plevel_loop)=0;
            end
            
            Pvalue_log=mean(Plevel,2);
            IncBindex=Pvalue_log>0;      
            DecBindex=Pvalue_log<0;
            P.comp.inc_bindex{IIXLoop}=IncBindex;
            P.comp.dec_bindex{IIXLoop}=DecBindex;
            %}}}
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %{{{ Correction of corr_pvalue & variation according to Inc_ and DecBindex
            % 2007 07 13 !!!! IN FACT CORRECTION WAS NOT MADE !!!! wrong position of end corresponding to if mod(CompLoopLimit,2)==0&CompLoopLimit>1 (L4957)            
            P.comp.corr_plevel(:,:,IIXLoop)=P.comp.plevel(:,:,IIXLoop);
            if P.par.big==0
                P.comp.corr_variation(:,:,IIXLoop)=P.comp.variation(:,:,IIXLoop);
            end
            
            %}}}
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %{{{ Treat indeterminated cases
            % in case of even number of comparison it may happens that there is an equal number of increased and decreased
            % assignation for a gene and that all the plevel are equal. It's then impossible to determine the change status
            % for these genes
            if mod(CompLoopLimit,2)==0&CompLoopLimit>1
                % Search for genes having the same pvalue for all the comparisons
                Equal_bindex=P.comp.plevel(:,1,IIXLoop);
                for Corr_loop=[2:CompLoopLimit]
                    Equal_bindex=Equal_bindex&P.comp.plevel(:,Corr_loop,IIXLoop)==P.comp.plevel(:,1,IIXLoop);
                end
                
                if size(find(Equal_bindex),1)>0
                    % Compute independently the number of increased and decreased assignation
                    Discrepancy_inc_bindex=P.comp.bindex_inc{1,IIXLoop};
                    for Corr_loop=2:CompLoopLimit
                        Discrepancy_inc_bindex=Discrepancy_inc_bindex+P.comp.bindex_inc{Corr_loop,IIXLoop};
                    end
                    Discrepancy_bindex=Discrepancy_inc_bindex==CompLoopLimit/2;
                    Discrepancy_equal_bindex=Discrepancy_bindex&Equal_bindex;
                    if size(find(Discrepancy_equal_bindex),1)>0
                        for Corr_loop=[1:CompLoopLimit]
                            P.comp.corr_variation(Discrepancy_equal_bindex,Corr_loop,IIXLoop)=VarThreshold;
                            P.comp.corr_plevel(Discrepancy_equal_bindex,Corr_loop,IIXLoop)=1;
                        end
                    end
                end
                %% 2007 07 13 !!!! IN FACT CORRECTION WAS NOT MADE !!!! correction of end corresponding to if mod(CompLoopLimit,2)==0&CompLoopLimit>1 (L4957)                
            end
            
            %}}}
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %{{{ % Compute the nb of discrepancy, correct variation, pv & calculate ppv
            Discrepancy_nb=zeros(size(IncBindex));
            for Corr_loop=1:CompLoopLimit
                Discrepancy_bindex=(IncBindex&P.comp.bindex_dec{Corr_loop,IIXLoop})|...
                    (DecBindex&P.comp.bindex_inc{Corr_loop,IIXLoop});
                P.comp.corr_variation(Discrepancy_bindex,Corr_loop,IIXLoop)=VarThreshold;
                P.comp.corr_plevel(Discrepancy_bindex,Corr_loop,IIXLoop)=1;
                Discrepancy_nb(Discrepancy_bindex)=Discrepancy_nb(Discrepancy_bindex)+1;
            end
            P.comp.corr_nb(:,IIXLoop)=Discrepancy_nb;
            %% 2007 07 13 !!!! IN FACT CORRECTION WAS NOT MADE !!!! wrong position of end corresponding to if mod(CompLoopLimit,2)==0&CompLoopLimit>1 (L4957)                
            %%end
            % Calculus of product of plevel from corrected values of plevel
            P.comp.product(:,IIXLoop)=prod(P.comp.corr_plevel(:,:,IIXLoop),2);
            if P.par.big==0
                % Calculus of var from corrected values of variation
                P.comp.var(:,IIXLoop)=mean(P.comp.corr_variation(:,:,IIXLoop),2);
                
            end
            %}}}
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %}}}
            %------------------------------------
            %------------------------------------
            %{{{ Calculus of pvalue(ppv)|Ho for increased genes and decreased genes
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %{{{ Increased genes
            %++++++++++++++++++++++++++++++++++++++
            %{{{ Selection of treated ppv; compute cdf(ppv)
            P.comp.maxpvbindex=max(P.comp.corr_plevel(:,:,IIXLoop),[],2)<=P.par.maxpvalue;
            P.comp.minpvbindex=~P.comp.maxpvbindex;
            Plevel_product=P.comp.product(:,IIXLoop);
            Plevel_product(P.comp.minpvbindex)=1;
            Plevel_product=Plevel_product(P.comp.inc_bindex{IIXLoop});
            [Sorted_plevel_product Direct_sort_index]=sort(Plevel_product);
            [temp Reverse_sort_index]=sort(Direct_sort_index);
            [temp Made_obs_dist temp temp temp]=CUMDISTPOS(Sorted_plevel_product,0);
            clear temp
            Made_obs_dist=Made_obs_dist(Reverse_sort_index);
            P.comp.pvalue_cdf(P.comp.inc_bindex{IIXLoop},IIXLoop)=Made_obs_dist;
            %}}}
            %++++++++++++++++++++++++++++++++++++++
            %++++++++++++++++++++++++++++++++++++++
            %{{{ Calculus of expected product of plevel cdf(ppv)|Ho
            if CompLoopLimit>length(P.calib.ppv)
                [P.calib.ppv{CompLoopLimit},P.calib.cdf{CompLoopLimit}]=MAKE_CDF_PPV(1,CompLoopLimit);
            else
                if isempty(P.calib.ppv{CompLoopLimit})
                    [P.calib.ppv{CompLoopLimit},P.calib.cdf{CompLoopLimit}]=MAKE_CDF_PPV(1,CompLoopLimit);
                end
            end
            
            
            Infered_plevel_BL=P.calib.ppv{CompLoopLimit};
            Infered_plevel_dist_BL=P.calib.cdf{CompLoopLimit};
            Min_plevel_BL=min(Infered_plevel_BL);
            Min_plevel_BL=Min_plevel_BL(1);
            Max_plevel_BL=max(Infered_plevel_BL);
            Max_plevel_BL=Max_plevel_BL(1);
            Inferior_index=find(Sorted_plevel_product<Min_plevel_BL);
            Sorted_plevel_product(Inferior_index)=Min_plevel_BL;
            Superior_index=find(Sorted_plevel_product>Max_plevel_BL);
            Sorted_plevel_product(Superior_index)=Max_plevel_BL;
            Diff_index=find(diff(Infered_plevel_BL)==0);
            while ~isempty(Diff_index)
                Infered_plevel_BL(Diff_index)=[];
                Infered_plevel_dist_BL(Diff_index)=[];
                Diff_index=find(diff(Infered_plevel_BL)==0);
            end
            %}}}
            %++++++++++++++++++++++++++++++++++++++
            %++++++++++++++++++++++++++++++++++++++
            %{{{ Abnormal situation for very low probability we get
            % Infered_plevel_BL =[1.0000    1.0000    1.0000]
            % and Infered_plevel_dist_BL =[0.0000    1.0000    1.0000]
            % and the identity of Infered_plevel_BL is not detected, provoking an error of interp1
            if size(Infered_plevel_BL,2)==3
                warndlg('abnormal pv(ppv|Ho)','COM_PREV')
                Infered_plevel_BL =[0   eps    1.0000];
                Infered_plevel_dist_BL =[0.0000    1.0000    1.0000];
            end
            %}}}
            %++++++++++++++++++++++++++++++++++++++
            %++++++++++++++++++++++++++++++++++++++
            %{{{ Calculate pv(ppv); Correct for abnormal values of pv(ppv|Ho); assignation of P.comp.pvalue
            Expected_plevel=interp1(Infered_plevel_BL,Infered_plevel_dist_BL,Sorted_plevel_product);
            Expected_plevel=Expected_plevel(Reverse_sort_index);
            Zero_index=find(Expected_plevel<=0);
            One_index=find(Expected_plevel>=1);
            Nan_bindex=isnan(Expected_plevel);
            Expected_plevel(Zero_index)=eps;
            Expected_plevel(One_index)=1;
            Expected_plevel(Nan_bindex)=1;
            P.comp.pvalue(P.comp.inc_bindex{IIXLoop},IIXLoop)=Expected_plevel;
            %}}}
            %++++++++++++++++++++++++++++++++++++++
            %}}}
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %{{{ Decreased genes
            %++++++++++++++++++++++++++++++++++++++
            %{{{ Selection of treated ppv; compute cdf(ppv)
            Plevel_product=P.comp.product(:,IIXLoop);
            Plevel_product(P.comp.minpvbindex)=1;
            Plevel_product=Plevel_product(P.comp.dec_bindex{IIXLoop});
            [Sorted_plevel_product Direct_sort_index]=sort(Plevel_product);
            [temp Reverse_sort_index]=sort(Direct_sort_index);
            [temp Made_obs_dist temp temp temp]=CUMDISTPOS(Sorted_plevel_product,0);
            clear temp
            Made_obs_dist=Made_obs_dist(Reverse_sort_index);
            P.comp.pvalue_cdf(P.comp.dec_bindex{IIXLoop},IIXLoop)=Made_obs_dist;
            %}}}
            %++++++++++++++++++++++++++++++++++++++
            %++++++++++++++++++++++++++++++++++++++
            %{{{ Calculus of expected product of plevel cdf(ppv)|Ho
            if CompLoopLimit>length(P.calib.ppv)
                [P.calib.ppv{CompLoopLimit},P.calib.cdf{CompLoopLimit}]=MAKE_CDF_PPV(1,CompLoopLimit);
            else
                if isempty(P.calib.ppv{CompLoopLimit})
                    [P.calib.ppv{CompLoopLimit},P.calib.cdf{CompLoopLimit}]=MAKE_CDF_PPV(1,CompLoopLimit);
                end
            end
            
            
            Infered_plevel_HL=P.calib.ppv{CompLoopLimit};
            Infered_plevel_dist_HL=P.calib.cdf{CompLoopLimit};
            
            Min_plevel_HL=min(Infered_plevel_HL);
            Min_plevel_HL=Min_plevel_HL(1);
            Max_plevel_HL=max(Infered_plevel_HL);
            Max_plevel_HL=Max_plevel_HL(1);
            Inferior_index=find(Sorted_plevel_product<Min_plevel_HL);
            Sorted_plevel_product(Inferior_index)=Min_plevel_HL;
            Superior_index=find(Sorted_plevel_product>Max_plevel_HL);
            Sorted_plevel_product(Superior_index)=Max_plevel_HL;
            
            Diff_index=find(diff(Infered_plevel_HL)==0);
            while ~isempty(Diff_index)
                Infered_plevel_HL(Diff_index)=[];
                Infered_plevel_dist_HL(Diff_index)=[];
                Diff_index=find(diff(Infered_plevel_HL)==0);
            end
            %}}}
            %++++++++++++++++++++++++++++++++++++++
            %++++++++++++++++++++++++++++++++++++++
            %{{{ Abnormal situation for very low probability we get
            % for very low probability we get
            % Infered_plevel_HL =[1.0000    1.0000    1.0000]
            % and Infered_plevel_dist_HL =[0.0000    1.0000    1.0000]
            % and the identity of Infered_plevel_HL is not detected, provoking an error of interp1
            if size(Infered_plevel_HL,2)==3
                Infered_plevel_HL =[0   eps    1.0000];
                Infered_plevel_dist_HL =[0.0000    1.0000    1.0000];
            end
            %}}}
            %++++++++++++++++++++++++++++++++++++++
            %++++++++++++++++++++++++++++++++++++++
            %{{{ Calculate pv(ppv); Correct for abnormal values of pv(ppv|Ho); assignation of P.comp.pvalue
            Expected_plevel=interp1(Infered_plevel_HL,Infered_plevel_dist_HL,Sorted_plevel_product);
            Expected_plevel=Expected_plevel(Reverse_sort_index);
            Zero_index=find(Expected_plevel<=0);
            One_index=find(Expected_plevel>=1);
            Nan_bindex=isnan(Expected_plevel);
            Expected_plevel(Zero_index)=eps;
            Expected_plevel(One_index)=1;
            Expected_plevel(Nan_bindex)=1;
            P.comp.pvalue(P.comp.dec_bindex{IIXLoop},IIXLoop)=Expected_plevel;
            P.comp.inc_ppv{IIXLoop}=Infered_plevel_BL;
            P.comp.inc_ppv_pv{IIXLoop}=Infered_plevel_dist_BL;
            P.comp.dec_ppv{IIXLoop}=Infered_plevel_HL;
            P.comp.dec_ppv_pv{IIXLoop}=Infered_plevel_dist_HL;
            %}}}
            %++++++++++++++++++++++++++++++++++++++
            %}}}
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %}}}
            %------------------------------------
            %}}}
            %====================================
            
            if Continue==0
                return
            end
        end
        %}}}
        %������������������������������������
    end % IIX loop
    
    %clear unused variables
    clear ActVarBL
    clear ActVarHL
    clear   BL 
    clear   Actual_rel_pos_BL1
    clear   Actual_rel_pos_BL1_1
    clear   HL 
    clear   Actual_rel_pos_HL1
    clear   Actual_rel_pos_HL1_1
    clear Actual_var_BL1_1
    clear Actual_var_HL1_1
    clear B
    clear BL
    clear BL1
    clear BL2
    clear CurrChip
    clear  DecBindex
    clear  IncBindex
    if P.par.big==1
        clear IncBindex
        clear DecBindex
    end
    clear Dec_bindex1
    clear Dec_bindex1_1
    clear Direct_sort_index
    clear Discrepancy_bindex
    clear Discrepancy_nb
    clear Expected_plevel
    clear FExpected_plevel
    clear FNormalized_act_var_BL
    clear FNormalized_act_var_HL
    clear FSorted_expected_plevel
    clear FSorted_observed_plevel
    clear HL
    clear HL1
    clear HL2
    clear Inc_bindex1
    clear Inc_bindex1_1
    clear Indirect_sort_index
    clear Infinite_bindex
    clear Inverse_sort_index
    clear Log_plevel
    clear Mad_obs_dist
    clear Nan_bindex
    clear Made_obs_dist
    clear Normalized_act_var_BL1
    clear Normalized_act_var_BL1_1        
    clear Normalized_act_var_HL1
    clear Normalized_act_var_HL1_1        
    clear PhasedRelPosBL
    clear PhasedRelPosHL        
    clear Phased_BL
    clear Phased_HL        
    clear BLSortedKeptZvar
    clear HLSortedKeptZvar        
    clear BLSortedZVar
    clear HLSortedZVar        
    clear Plevel
    clear Plevel_product
    clear PlusBindex
    clear Pvalue_log
    clear Reverse_sort_index
    clear Selection_bindex_BL1
    clear Selection_bindex_HL1        
    clear SortIndex
    clear SortedRelPosBL
    clear SortedRelPosHL
    clear BLSorted
    clear HLSorted
    clear Sorted_plevel_product
    clear Sorted_to_phased_index_BL
    clear Sorted_to_phased_index_HL
    clear Temp_index
    clear VarBL
    clear VarBL01
    clear VarHL
    clear VarHL01
    clear X
    clear Zero_index
    clear ZVal
    
    if ~isequal(P.par.comptype,'clustering')&Continue==1
        %������������������������������������
        %{{{ STATCALIB
        
        %====================================
        %{{{ FINAL COMPARISON SYNTHESIS FOR MULTI
        %------------------------------------
        %{{{ Description of content
        % first column : mean of rank
        % second column : mean of normalized variation
        % third column : interpolated variation from plevel values (not used)
        % fourth column : ppv
        % fifth : cdf(ppv)
        % sixth column : pv(ppv)
        %seventh column : FDR
        %eighth column : mean(max(pv) in each IIX comparison)
        %}}}
        %------------------------------------
        if isequal(P.par.comptype,'II & X')
            %------------------------------------
            %{{{ Several IIX
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %{{{ Determination of variation (inc or dec) with log(pvalue)
            Plevel=P.comp.product;
            for Plevel_loop=[1:IIXLoopLimit]
                Log_plevel=log(1./Plevel(P.comp.inc_bindex{Plevel_loop},Plevel_loop));
                Infinite_bindex=~isfinite(Log_plevel);
                Log_plevel(Infinite_bindex)=0;
                Plevel(P.comp.inc_bindex{Plevel_loop},Plevel_loop)=Log_plevel;
                Log_plevel=-log(1./Plevel(P.comp.dec_bindex{Plevel_loop},Plevel_loop));
                Infinite_bindex=~isfinite(Log_plevel);
                Log_plevel(Infinite_bindex)=0;
                Plevel(P.comp.dec_bindex{Plevel_loop},Plevel_loop)=Log_plevel;
                Plevel(~P.comp.inc_bindex{Plevel_loop}&~P.comp.dec_bindex{Plevel_loop},Plevel_loop)=0;
            end
            Mean_plevel=mean(Plevel,2);
            IncBindex=Mean_plevel>0;
            DecBindex=Mean_plevel<0;
            P.comp.incbindex=IncBindex;
            P.comp.decbindex=DecBindex;
            %}}}
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %{{{ Correction of product
            P.comp.corr_product=P.comp.product;
            P.comp.corr_var=P.comp.var;
            P.comp.corr_pvalue=P.comp.pvalue;
            for Corr_loop=[1:IIXLoopLimit]
                Discrepancy_bindex=(IncBindex&P.comp.dec_bindex{Corr_loop})|(DecBindex&P.comp.inc_bindex{Corr_loop});
                P.comp.corr_product(Discrepancy_bindex,Corr_loop)=1;
                P.comp.corr_var(Discrepancy_bindex,Corr_loop)=VarThreshold;
                P.comp.corr_pvalue(Discrepancy_bindex,Corr_loop)=1;
            end
            %}}}
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %{{{ Fill Comparison
            Comp_size=length(P.comp.product(:,1));
            if P.par.big==1
                P.comparison=[mean(P.comp.corr_product,2),...
                        ones(Comp_size,1)];
            else
                P.comparison=[mean(P.comp.BL_mean_cdfpos,2),...
                        mean(P.comp.corr_var,2),...
                        zeros(Comp_size,1),...
                        mean(P.comp.corr_product,2),...
                        ones(Comp_size,1),...
                        mean(P.comp.corr_pvalue,2),...
                        ones(Comp_size,1),...
                        ones(Comp_size,1),...
                        ones(Comp_size,1)];
                P.comparison(:,8)=mean(mean(max(P.comp.corr_plevel,[],2),2),3);
            end
            %}}}
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %}}}
            %------------------------------------
        elseif isequal(P.par.comptype,'II')|isequal(P.par.comptype,'II all')
            %------------------------------------
            %{{{ One IIX
            P.comp.incbindex=P.comp.inc_bindex{1};
            P.comp.decbindex=P.comp.dec_bindex{1};
            if P.par.big==1
                P.comp.inc_bindex{1}=[];
                P.comp.dec_bindex{1}=[];
            end
            Comp_size=size(P.comp.product(:,1),1);
            if P.par.big==1
                P.comp.bindex_inc=[];
                P.comp.bindex_dec=[];
                P.comp.inc_bindex=[];
                P.comp.dec_bindex=[];
                P.comp.corr_nb=[];
                P.comp.maxpvbindex=[];
                P.comp.minpvbindex=[];
                P.comp.pvalue_cdf=[];
                P.comp.pvalue=[];
                P.comparison=zeros(Comp_size,5);
                P.comparison(:,1)=P.comp.product(:,1);
                P.comparison(:,5)=P.comp.var(:,1);
            else
                P.comparison=[mean(P.comp.BL_mean_cdfpos,2),...
                        P.comp.var(:,1),...
                        zeros(Comp_size,1),...
                        P.comp.product(:,1),...
                        ones(Comp_size,1),...
                        P.comp.pvalue(:,1),...
                        ones(Comp_size,1),...
                        ones(Comp_size,1),...
                        ones(Comp_size,1)];
                P.comparison(:,8)=max(P.comp.corr_plevel(:,:,1),[],2);
            end
            %}}}
            %------------------------------------
        end
        %------------------------------------     
        %{{{ Correct abnormal values for P.comparison(:,4) (mean(ppv))
        % % % %             Zero_index=find(P.comparison(:,4)<=0);
        % % % %             One_index=find(P.comparison(:,4)>1);
        % % % %             Nan_bindex=isnan(P.comparison(:,4));
        % % % %             Nan_index=find(Nan_bindex);
        % % % %             if ~isempty(Zero_index)
        % % % %                 %'zero in comparison(:,4)'
        % % % %                 Expected_plevel(Zero_index)=eps;
        % % % %             end
        % % % %             if ~isempty(One_index)
        % % % %                 %'one in comparison(:,4)'
        % % % %                 Expected_plevel(One_index)=1;
        % % % %             end
        % % % %             if ~isempty(Nan_index)
        % % % %                 %'nan in comparison(:,4)'
        % % % %                 Expected_plevel(Nan_index)=1;
        % % % %             end
        % Store original values of ppp and ppv_pv
        if P.par.big==0
            P.comparison(:,9)=P.comparison(:,4);
            P.comparison(:,12)=P.comparison(:,6);
            P.comp.maxpvbindex=P.comparison(:,8)<=P.par.maxpvalue;
            
            %}}}
            %------------------------------------
            %------------------------------------
            %{{{ Correction of pv(ppv) and ppv if MaxPv<1
            
            % due to the procedure used for the calculus of max(plevel)
            % a point could have mean(max(plevel))< Maxplevel_limit_PREP
            % although have a max(plevel) >Maxplevel_limit_PREP in
            % one of the II,X combinations. The corresponding pv(ppv) will be
            % one and the final mean(pv(ppv)) = ~0.5 in double comparison
            % So the pv(ppv) of these points must be corrected
            
            
            % Bindex of points with a variation pvalue inferior to P.par.maxpvalue
            % in all comparisons in all IIX combinations
            Included_bindex=ones(Comp_size,1)==1;
            
            for IIXLoop=[1:IIXLoopLimit]
                for CompLoop=[1:CompLoopLimit]
                    Included_bindex=Included_bindex&P.comp.plevel(:,CompLoop,IIXLoop)<=P.par.maxpvalue;
                end
            end
            
            % Bindex of points with have at least one variation pvalue superior to P.par.maxpvalue
            % Some of these points could ultimately have a mean(variation pvalue)<=P.par.maxpvalue
            % and then be incorporated in P.comp.maxpvbindex
            Excluded_bindex=~Included_bindex;
            Ppv_bindex=zeros(Comp_size,1);
            Ppv_bindex=Ppv_bindex==1;
            Interpolation_bindex=Ppv_bindex;
            for Incdec_loop=[1:2]
                if Incdec_loop==1
                    Incdec_rank='inc';
                else
                    Incdec_rank='dec';
                end
                
                for IIXLoop=1:IIXLoopLimit
                    % Points which have at least one variation pvalue > P.par.maxpvalue in the IIXLoop currently treated
                    % but which have been selected because mean(pvalue) is <= P.par.maxpvalue
                    eval(['Ppv_bindex=(Ppv_bindex)|(P.comp.maxpvbindex&Excluded_bindex&P.comp.',Incdec_rank,'_bindex{IIXLoop});']);
                    % bindex of data used to do interpolation of pv(ppv) on ppv
                    eval(['Interpolation_bindex=(Interpolation_bindex)|(P.comp.maxpvbindex&Included_bindex&P.comp.',Incdec_rank,'_bindex{IIXLoop});']);
                    % selection of data used to do interpolation
                end % of for IIXLoop=[1:IIXLoopLimit]
                
                if ~isempty(find(Interpolation_bindex))&~isempty(find(Ppv_bindex))
                    Ppv=P.comparison(Interpolation_bindex,4);
                    Ppv_pv=P.comparison(Interpolation_bindex,6);
                    Corr_ppv=P.comparison(Ppv_bindex,4);
                    Undermax_bindex=Corr_ppv>P.par.maxpvalue^CompLoopLimit;
                    Corr_ppv(Undermax_bindex)=(P.par.maxpvalue^CompLoopLimit)-eps;
                    
                    
                    % assure that values to be corrected (Corr_ppv)
                    % can be interpolated
                    Min_ppv=min(Ppv);
                    Min_ppv=Min_ppv(1);
                    Max_ppv=max(Ppv);
                    Max_ppv=Max_ppv(1);
                    Inferior_index=find(Corr_ppv<Min_ppv);
                    Corr_ppv(Inferior_index)=Min_ppv;
                    Superior_index=find(Corr_ppv>Max_ppv);
                    Corr_ppv(Superior_index)=Max_ppv;
                    
                    % keep Ppv_pv strictly increasing by propagating local minimum
                    [Ppv Direct_sort_index]=sort(Ppv);
                    Ppv_pv=Ppv_pv(Direct_sort_index);
                    Ppv_pv=MINSTEP(Ppv_pv,1);
                    % remove duplicated values of Ppv
                    Diff_index=find(diff(Ppv)==0);
                    while ~isempty(Diff_index)
                        Ppv(Diff_index)=[];
                        Ppv_pv(Diff_index)=[];
                        Diff_index=find(diff(Ppv)==0);
                    end
                    
                    % rmove duplicated values of Ppv_pv
                    [Ppv_pv Direct_sort_index]=sort(Ppv_pv);
                    Ppv=Ppv(Direct_sort_index);
                    Diff_index=find(diff(Ppv_pv)==0);
                    while ~isempty(Diff_index)
                        Ppv(Diff_index)=[];
                        Ppv_pv(Diff_index)=[];
                        Diff_index=find(diff(Ppv_pv)==0);
                    end
                    % interpolation of ppv_pv for ppv values to be corrected
                    Corr_ppv_pv=interp1(Ppv,Ppv_pv,Corr_ppv);
                    P.comparison(Ppv_bindex,6)=Corr_ppv_pv;
                end
            end % of for Incdec_loop=[1:2]
            
            
            % finish by keeping  Ppv_pv strictly increasing
            % important to filter with P.comp.maxpvbindex
            if size(find(P.comp.incbindex&P.comp.maxpvbindex),1)>0
                Ppv=P.comparison(P.comp.incbindex&P.comp.maxpvbindex,9);
                Ppv_pv=P.comparison(P.comp.incbindex&P.comp.maxpvbindex,6);
                [Ppv Direct_sort_index]=sort(Ppv);
                [Temp Reverse_sort_index]=sort(Direct_sort_index);
                Ppv_pv=Ppv_pv(Direct_sort_index);
                Ppv_pv=MAXSTEP(Ppv_pv,0);
                Ppv_pv=Ppv_pv(Reverse_sort_index);
                P.comparison(P.comp.incbindex&P.comp.maxpvbindex,6)=Ppv_pv;
            end
            
            if size(find(P.comp.decbindex&P.comp.maxpvbindex),1)>0
                Ppv=P.comparison(P.comp.decbindex&P.comp.maxpvbindex,9);
                Ppv_pv=P.comparison(P.comp.decbindex&P.comp.maxpvbindex,6);
                [Ppv Direct_sort_index]=sort(Ppv);
                [Temp Reverse_sort_index]=sort(Direct_sort_index);
                Ppv_pv=Ppv_pv(Direct_sort_index);
                Ppv_pv=MAXSTEP(Ppv_pv,0);
                Ppv_pv=Ppv_pv(Reverse_sort_index);
                P.comparison(P.comp.decbindex&P.comp.maxpvbindex,6)=Ppv_pv;
            end
            
            %}}}
            %------------------------------------
            %}}}
            %====================================
            
            %====================================
            %{{{ Correction Norm_var_s{1,2,3} & DISPLAY of result
            P.flag.calib=0;
            P.par.displaytype='Multi';
            Display_selection_flag=0;
            if P.flag.display==1
                DISPLAY_MULTI_PREV
            end
        end
        %}}}
        %====================================
        %====================================
        %{{{ FINAL STATISTICS
        Skip_statistics=0;
        if Skip_statistics==0
            %------------------------------------
            %{{{ Figure selection
            if P.flag.display==1
                if F.h.stat(3)==0
                    h=warndlg('pause; OK and >dbcont');
                    waitfor(h)
                    F.h.stat(3)=figure;
                    set(F.h.stat(3),'tag','STAT3_fiPREP')
                    set(F.h.stat(3),'Color',[1 1 1])
                else
                    try
                        figure(F.h.stat(3))
                    catch
                        F.h.stat(3)=figure;
                        set(F.h.stat(3),'tag','STAT3_fiPREP')
                    end
                end
                F.gh.prev.delete=get(F.h.stat(3),'children');
                delete(F.gh.prev.delete);
            end
            %}}}
            %------------------------------------
            %------------------------------------
            %{{{  ppv & cdf(ppv) for increased
            % First render monotonous ppv vs pv(ppv)
            if P.par.big==0
                Corr_ppv=P.comparison(:,4);
                Corr_ppv(P.comp.minpvbindex)=1;
                P.comparison(P.comp.incbindex,4)=Corr_ppv(P.comp.incbindex);
                P.comparison(P.comp.decbindex,4)=Corr_ppv(P.comp.decbindex);
                Ppv=P.comparison(:,4);
                
                
                if CompLoopLimit>length(P.calib.ppv)
                    [P.calib.ppv{CompLoopLimit},P.calib.cdf{CompLoopLimit}]=MAKE_CDF_PPV(1,CompLoopLimit);
                else
                    if isempty(P.calib.ppv{CompLoopLimit})
                        [P.calib.ppv{CompLoopLimit},P.calib.cdf{CompLoopLimit}]=MAKE_CDF_PPV(1,CompLoopLimit);
                    end
                end
                
                
                PpvHo=P.calib.ppv{CompLoopLimit};
                PvHo=P.calib.cdf{CompLoopLimit};
                MinPpvHo=min(PpvHo);
                MinPpvHo=MinPpvHo(1);
                MinIndex=find(Ppv<MinPpvHo);
                if ~isempty(MinIndex)
                    Ppv(MinIndex)=MinPpvHo;
                end
                Pv=interp1(PpvHo,PvHo,Ppv);
                P.comparison(:,6)=Pv;
                
                % Correct pv(ppv) to allow display
                Corr_ppv_pv=P.comparison(:,6);
                Corr_ppv_pv(P.comp.minpvbindex)=1;
                %Corr_ppv_pv=Corr_ppv_pv(P.comp.incbindex);
                P.comparison(P.comp.incbindex,6)=Corr_ppv_pv(P.comp.incbindex);
                P.comparison(P.comp.decbindex,6)=Corr_ppv_pv(P.comp.decbindex);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % cdf(ppv) for Increased
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % Calculate cdf of ppv
            %[Sorted_corr_ppv Direct_sort_index]=sort(Corr_ppv);
            if P.par.big==0
                [Sorted_corr_ppv Direct_sort_index]=sort(P.comparison(P.comp.incbindex,4));
            else
                [Sorted_corr_ppv Direct_sort_index]=sort(P.comparison(P.comp.incbindex,1));
            end
            [temp Reverse_sort_index]=sort(Direct_sort_index);
            [temp ppv_cdf temp temp temp]=CUMDISTPOS(Sorted_corr_ppv,0);
            ppv_cdf=ppv_cdf(Reverse_sort_index);
            %P.comparison(P.comp.incbindex & P.comp.maxpvbindex,5)=ppv_cdf;
            if P.par.big==0
                P.comparison(P.comp.incbindex,5)=ppv_cdf;
            else
                P.comparison(P.comp.incbindex,2)=ppv_cdf;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % cdf(ppv) for Decreased
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %[Sorted_corr_ppv Direct_sort_index]=sort(Corr_ppv);
            if P.par.big==0
                [Sorted_corr_ppv Direct_sort_index]=sort(P.comparison(P.comp.decbindex,4));
            else
                [Sorted_corr_ppv Direct_sort_index]=sort(P.comparison(P.comp.decbindex,1));
            end
            [temp Reverse_sort_index]=sort(Direct_sort_index);
            [temp ppv_cdf temp temp temp]=CUMDISTPOS(Sorted_corr_ppv,0);
            ppv_cdf=ppv_cdf(Reverse_sort_index);
            if P.par.big==0
                P.comparison(P.comp.decbindex,5)=ppv_cdf;
            else
                P.comparison(P.comp.decbindex,2)=ppv_cdf;
            end
            %}}}
            %------------------------------------
            %------------------------------------
            %{{{ Calculus of FDR for increased genes
            %Calculus of FP rate for increased genes
            Variation_number=size(find(P.comp.incbindex),1);
            if P.par.big==0
                XUnit=P.comparison(P.comp.incbindex,4);
            else
                XUnit=P.comparison(P.comp.incbindex,1);
            end
            %%%Max_Xunit=max(XUnit);
            Exp_plevel=interp1(P.calib.ppv{CompLoopLimit},P.calib.cdf{CompLoopLimit},XUnit);
            [Sorted_Xunit Direct_sort_index]=sort(XUnit);
            clear XUnit
            %%%Max_exp=max(Exp_plevel);
            Phased_exp_plevel=Exp_plevel(Direct_sort_index);
            clear Exp_plevel
            [Temp Reverse_sort_index]=sort(Direct_sort_index);
            clear Temp
            if P.par.big==0
                Obs_plevel=P.comparison(P.comp.incbindex,5);
            else
                Obs_plevel=P.comparison(P.comp.incbindex,2);
            end
            %Max_obs=max(Obs_plevel);
            Phased_obs_plevel=Obs_plevel(Direct_sort_index);
            clear Direct_sort_index
            Phased_smooth_obs_plevel=MAXSTEP(Phased_obs_plevel,0);
            clear Phased_obs_plevel
            % Phased_variation not used in STAT_CALIB_PREV
            %%%Phased_variation=1;
            
            %%%X=log(1./Sorted_Xunit);
            %%%Cdf0=Phased_exp_plevel;
            %%%Cdf=Phased_smooth_obs_plevel;
            %%%[X ,SortIndex]=sort(X);
            %%%Cdf0=Cdf0(SortIndex);
            %%%Cdf=Cdf(SortIndex);
            %%%FigTitle='Product of pv';
            %%%FigTag='I_TOTALCH_FIG';
            %%%FigTitle='Increased';
            %%%LineNb=1;
            %%%ColNb=2;
            %%%SubRank=1;
            %TOTALCHANGE(X,Cdf0,Cdf,FigTag,FigTitle,SubTitle,LineNb,ColNb,SubRank)
            
            
            %[Inc_error1,Inc_error2,Inc_sensitivity,Inc_specificity,Inc_chi2_probability,P.comp.inctotalnb,Inc_max_useful_FPR]=STAT_CALIB_PREV(Variation_number,...
            %    Phased_variation,Phased_exp_plevel,Phased_smooth_obs_plevel,1,F.h.stat(3),'r','',Sorted_Xunit,0);
            if P.flag.display==1
                HFig=F.h.stat(3);
            else
                HFig=0;
            end
            [Inc_error1,Inc_error2,Inc_sensitivity,Inc_specificity,P.comp.inctotalnb,Inc_max_useful_FPR]=STAT_CALIB_PREV(Variation_number,...
                Phased_exp_plevel,Phased_smooth_obs_plevel,P.flag.display,HFig,'r','',Sorted_Xunit,0,1);
            Inc_error1=Inc_error1(Reverse_sort_index);
            %%%Inc_error2=Inc_error2(Reverse_sort_index);
            Inc_sensitivity=Inc_sensitivity(Reverse_sort_index);
            %%%Inc_specificity=Inc_specificity(Reverse_sort_index);
            %%%Inc_chi2_probability=Inc_chi2_probability(Reverse_sort_index);
            %}}}
            %------------------------------------
            %------------------------------------
            %{{{ Calculus of FDR for decreased genes
            
            Variation_number=size(find(P.comp.decbindex),1);
            if P.par.big==0
                XUnit=P.comparison(P.comp.decbindex,4);
            else
                XUnit=P.comparison(P.comp.decbindex,1);
            end
            Exp_plevel=interp1(P.calib.ppv{CompLoopLimit},P.calib.cdf{CompLoopLimit},XUnit);
            [Sorted_Xunit Direct_sort_index]=sort(XUnit);
            clear XUnit
            Phased_exp_plevel=Exp_plevel(Direct_sort_index);
            clear Exp_plevel
            [Temp Reverse_sort_index]=sort(Direct_sort_index);
            clear Temp
            if P.par.big==0
                Obs_plevel=P.comparison(P.comp.decbindex,5);
            else
                Obs_plevel=P.comparison(P.comp.decbindex,2);
            end
            Phased_obs_plevel=Obs_plevel(Direct_sort_index);
            clear Obs_plevel
            clear Direct_sort_index
            Phased_smooth_obs_plevel=MAXSTEP(Phased_obs_plevel,0);
            clear Phased_obs_plevel
            
            %X=log(1./Sorted_Xunit);
            %Cdf0=Phased_exp_plevel;
            %Cdf=Phased_smooth_obs_plevel;
            %[X ,SortIndex]=sort(X);
            %Cdf0=Cdf0(SortIndex);
            %Cdf=Cdf(SortIndex);
            %FigTitle='Product of pv';
            %FigTitle='Decreased';
            %FigTag='D_TOTALCH_FIG';
            %LineNb=1;
            %ColNb=2;
            %SubRank=2;
            %TOTALCHANGE(X,Cdf0,Cdf,FigTag,FigTitle,SubTitle,LineNb,ColNb,SubRank)
            
            
            %                 [Dec_error1,Dec_error2,Dec_sensitivity,Dec_specificity,Dec_chi2_probability,P.comp.dectotalnb,Dec_max_useful_FPR]=...
            %                     STAT_CALIB_PREV(Variation_number,Phased_variation,...
            %                     Phased_exp_plevel,Phased_smooth_obs_plevel,1,...
            %                     F.h.stat(3),'b','Mean of Comparison II and Comparison X (Increased in red, Decreased in blue)',...
            %                     Sorted_Xunit,P.comp.inctotalnb);
            [Dec_error1,Dec_error2,Dec_sensitivity,Dec_specificity,P.comp.dectotalnb,Dec_max_useful_FPR]=...
                STAT_CALIB_PREV(Variation_number,...
                Phased_exp_plevel,Phased_smooth_obs_plevel,P.flag.display,...
                HFig,'b','Mean of Comparison II and Comparison X (Increased in red, Decreased in blue)',...
                Sorted_Xunit,P.comp.inctotalnb,2);
            if isfield(P.flag,'abacus')
                if P.flag.abacus==1
                    cd(K.dir.comp)
                    ResName=strrep(sprintf('%s_vs_%s',P.par.hlgrpname{1},P.par.blgrpname{1}),' ','_');
                    set(F.h.stat(3),'name',ResName)
                    set(F.h.stat(3),'units','normalized')
                    set(F.h.stat(3),'Position',[0.05 0.05 0.60 0.40])
                    ResName=strrep(ResName,'/','-');
                    saveas(F.h.stat(3),sprintf('%s_abaque_chip%c',ResName,P.par.array),'png')
                end
            end
            
            
            
            Dec_error1=Dec_error1(Reverse_sort_index);
            %Dec_error2=Dec_error2(Reverse_sort_index);
            Dec_sensitivity=Dec_sensitivity(Reverse_sort_index);
            %Dec_specificity=Dec_specificity(Reverse_sort_index);
            %Dec_chi2_probability=Dec_chi2_probability(Reverse_sort_index);
            
            %P.comparison(P.comp.incbindex & P.comp.maxpvbindex,7)=Inc_error1;
            if P.par.big==0
                P.comparison(P.comp.incbindex,7)=Inc_error1;
                P.comparison(P.comp.incbindex,13)=Inc_sensitivity;
                %P.comparison(P.comp.decbindex & P.comp.maxpvbindex,7)=Dec_error1;
                P.comparison(P.comp.decbindex,7)=Dec_error1;
                P.comparison(P.comp.decbindex,13)=Dec_sensitivity;
            else
                P.comparison(P.comp.incbindex,3)=Inc_error1;
                P.comparison(P.comp.incbindex,4)=Inc_sensitivity;
                %P.comparison(P.comp.decbindex & P.comp.maxpvbindex,7)=Dec_error1;
                P.comparison(P.comp.decbindex,3)=Dec_error1;
                P.comparison(P.comp.decbindex,4)=Dec_sensitivity;
            end
            %}}}
            %------------------------------------
        end
        %------------------------------------
        %{{{ Determination of number of outliers
        if P.par.big==0
            if isequal(P.par.comptype,'II & X')
                P.comparison(:,10)=mean(sum(P.comp.single_outlier,2),3);
                P.comparison(:,11)=mean(sum(P.comp.double_outlier,2),3);
            else
                P.comparison(:,10)=sum(P.comp.single_outlier(:,:,1),2);
                P.comparison(:,11)=sum(P.comp.double_outlier(:,:,1),2);
            end
        end
        %}}}
        %------------------------------------
        %}}}
        %====================================
        %====================================
        %{{{ Variables necessary to display points selected by MAS5 comparison
        if P.par.big==0
            if P.flag.display==1
                figure(F.gh.prev.PREV_fi)
            end
            
            
            % max of FPR for increased and decreased
            Max_FPR_inc=(size(find(P.comp.incbindex),1)-P.comp.inctotalnb)/size(find(P.comp.incbindex),1);
            Max_FPR_dec=(size(find(P.comp.decbindex),1)-P.comp.dectotalnb)/size(find(P.comp.decbindex),1);
            % max of useful FPR
            % it is not advisable to select FPR greater than the one where all the TP are selected (TP_Tnb)
            % and the minimum of FP
            % These limit values are in global K variables Inc_max_useful_FPR and Dec_max_useful_FPR
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Increased variations
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Recover Ppv_pv
            Ppv_pv_inc=P.comparison(P.comp.incbindex,6);
            FPR_inc=P.comparison(P.comp.incbindex,7);
            % FPR vs Ppv_pv is monotonic (increasing)
            
            % sort on Ppv_pv which is the most important criterium
            % and the most rigourously calculated
            % FPR is not strictily monotonic in the reality
            % It is corrected to allow sensible selection on this criterium
            [P.comp.incpv4fdr Direct_sort_index]=sort(Ppv_pv_inc);
            P.comp.incfdr4pv=FPR_inc(Direct_sort_index);
            
            % data are rendered strictly increasing
            Zero_index=find(diff(P.comp.incfdr4pv)==0);
            P.comp.incpv4fdr(Zero_index)=[];
            P.comp.incfdr4pv(Zero_index)=[];
            
            Zero_index=find(diff(P.comp.incpv4fdr)==0);
            P.comp.incpv4fdr(Zero_index)=[];
            P.comp.incfdr4pv(Zero_index)=[];
            
            % an extra data is added to be sure that the range of real FPR data will be
            % entirely covered
            Min_sorted_ppv_pv_inc=min(P.comp.incpv4fdr);
            Min_phased_FPR_inc=min(P.comp.incpv4fdr);
            if Min_sorted_ppv_pv_inc>eps &  Min_phased_FPR_inc>eps
                % sometimes exists values < eps
                % in tis case don't add eps to keep FPR vs Ppv_pv strictly increasing
                P.comp.incpv4fdr=[eps;P.comp.incpv4fdr];
                P.comp.incfdr4pv=[eps;P.comp.incfdr4pv];
            end
            
            % if does not exist ppv_pv of value one add it
            One_index=find(P.comp.incpv4fdr==1);
            if isempty(One_index)
                P.comp.incpv4fdr=[P.comp.incpv4fdr;1];
                if Max_FPR_inc>max(P.comp.incfdr4pv)
                    P.comp.incfdr4pv=[P.comp.incfdr4pv;Max_FPR_inc];
                elseif Max_FPR_inc==max(P.comp.incfdr4pv)
                    P.comp.incfdr4pv=[P.comp.incfdr4pv;max(P.comp.incfdr4pv)+eps];
                elseif Max_FPR_inc<max(P.comp.incfdr4pv)
                    P.comp.incfdr4pv=[P.comp.incfdr4pv;max(P.comp.incfdr4pv)+eps];
                end
            end
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Decreased variations
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            Ppv_pv_dec=P.comparison(P.comp.decbindex,6);
            FPR_dec=P.comparison(P.comp.decbindex,7);
            
            [P.comp.decpv4fdr Direct_sort_index]=sort(Ppv_pv_dec);
            P.comp.decfdr4pv =FPR_dec(Direct_sort_index);
            
            
            Zero_index=find(diff(P.comp.decfdr4pv)==0);
            P.comp.decpv4fdr(Zero_index)=[];
            P.comp.decfdr4pv(Zero_index)=[];
            
            Zero_index=find(diff(P.comp.decpv4fdr)==0);
            P.comp.decpv4fdr(Zero_index)=[];
            P.comp.decfdr4pv(Zero_index)=[];
            
            Min_sorted_ppv_pv_dec=min(P.comp.decpv4fdr);
            Min_phased_FPR_dec=min(P.comp.decpv4fdr);
            if Min_sorted_ppv_pv_dec>eps &  Min_phased_FPR_dec>eps
                P.comp.decpv4fdr=[eps;P.comp.decpv4fdr];
                P.comp.decfdr4pv=[eps;P.comp.decfdr4pv];
            end
            % if does not exist ppv_pv of value one add it
            One_index=find(P.comp.decpv4fdr==1);
            if isempty(One_index)
                P.comp.decpv4fdr=[P.comp.decpv4fdr;1];
                if Max_FPR_dec>max(P.comp.decfdr4pv)
                    P.comp.decfdr4pv=[P.comp.decfdr4pv;Max_FPR_dec];
                elseif Max_FPR_dec==max(P.comp.decfdr4pv)
                    P.comp.decfdr4pv=[P.comp.decfdr4pv;max(P.comp.decfdr4pv)+eps];
                elseif Max_FPR_dec<max(P.comp.decfdr4pv)
                    % should not happen as STAT_CALIB_PREV copes with this problem
                    P.comp.decfdr4pv=[P.comp.decfdr4pv;max(P.comp.decfdr4pv)+eps];
                end
            end
        end
        
        %}}}
        %====================================
        %}}}
        %������������������������������������
    end
    
    %}}}
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %{{{ construct the systematic cross BLxHL
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    Skip_dispersion=0;
    if Skip_dispersion==0&Continue==1
        if isequal(P.par.comptype,'II & X')&isequal(P.par.calibtype,'Single comparison')&P.flag.multi==1
            BL_nb=length(P.par.blgrprank{1});
            HL_nb=length(P.par.hlgrprank{1});
            BL_mean=mean(P.comp.BL_mean_signal(:,:,1),2);
            HL_mean=mean(P.comp.HL_mean_signal(:,:,1),2);
            Com.Xvariation=zeros(size(P.comparison(:,1),1),BL_nb,HL_nb);
            Com.Yvariation=Com.Xvariation;
            
            
            for BL_loop=[1:BL_nb]
                for HL_loop=[1:HL_nb]
                    for XY_loop=[1:2]
                        if XY_loop==1
                            HL=P.comp.blgroup.signal(:,BL_loop);
                            BL=BL_mean;
                        else
                            HL=P.comp.hlgroup.signal(:,HL_loop);
                            BL=HL_mean;
                        end
                        
                        % calculate variation of X position
                        [BLSorted,Phased_HL,Direct_sort_index_BL,BLReverseSort,...
                                Sorted_to_phased_index_HL]...
                            =Sort2Series(BL,HL);
                        [HLSorted,Phased_BL,Direct_sort_index_HL,HLReverseSort,...
                                Sorted_to_phased_index_BL]...
                            =Sort2Series(HL,BL);                        
                        BLSorted=RANKED(BL,Sort_flag,P.par.signalthresh);
                        HLSorted=RANKED(HL,Sort_flag,P.par.signalthresh);
                        
                        BL=BLSorted(BLReverseSort);
                        HL=HLSorted(HLReverseSort);
                        
                        VarThreshold=HLSorted(1)-BLSorted(1);
                        if VarThreshold<0
                            VarThreshold=0;
                        end
                        BLVar=HL-BL;
                        
                        VarThreshold=BLSorted(1)-HLSorted(1);
                        if VarThreshold<0
                            VarThreshold=0;
                        end
                        HLVar=BL-HL;
                        
                        BLSortedKeptZvar=HLSorted(Sorted_to_phased_index_HL)...
                            -BLSorted;
                        HLSortedKeptZvar=BLSorted(Sorted_to_phased_index_BL)...
                            - HLSorted;
                        
                        
                        Inc_bind=BLVar>VarThreshold;
                        Cut_phased_inc_var_bind=BLSortedKeptZvar>VarThreshold;
                        BLFittedMean=interp1(XGrid,MeanRank,BLSorted);
                        BLFittedStd=interp1(XGrid,StdRank,BLSorted);
                        BLSortedZVar=(BLSortedKeptZvar-BLFittedMean)./BLFittedStd;
                        
                        BLSortedKeptZvar=zeros(size(BL,1),1);
                        BLSortedKeptZvar(Cut_phased_inc_var_bind)...
                            =BLSortedZVar(Cut_phased_inc_var_bind);
                        BLZVar=BLSortedKeptZvar(BLReverseSort);
                        
                        Dec_bind=HLVar>VarThreshold;
                        % taking ~Inc_bind incorporate zero variation in Dec_bind
                        % and produces discontinuity in observed_plevel
                        %Dec_bind=~Inc_bind;
                        Cut_phased_dec_var_bind=HLSortedKeptZvar>VarThreshold;
                        HLFittedMean=interp1(XGrid,MeanRank,HLSorted);
                        HLFittedStd=interp1(XGrid,StdRank,HLSorted);
                        HLSortedZVar=(HLSortedKeptZvar-HLFittedMean)./HLFittedStd;
                        
                        HLSortedKeptZvar=zeros(size(HL,1),1);
                        HLSortedKeptZvar(Cut_phased_dec_var_bind)...
                            =HLSortedZVar(Cut_phased_dec_var_bind);
                        HLZVar=HLSortedKeptZvar(HLReverseSort);
                        
                        if XY_loop==1
                            P.comp.Xvariation(Inc_bind,BL_loop,HL_loop)=BLZVar(Inc_bind);
                            P.comp.Xvariation(Dec_bind,BL_loop,HL_loop)=-HLZVar(Dec_bind);
                        else
                            P.comp.Yvariation(Inc_bind,BL_loop,HL_loop)=BLZVar(Inc_bind);
                            P.comp.Yvariation(Dec_bind,BL_loop,HL_loop)=-HLZVar(Dec_bind);
                        end
                    end
                end
            end
            if P.flag.display==1
                h=figure;
                set(h,'name','detection of abnormal experiment');
                hold on
                for HL_loop=[1:HL_nb]
                    for BL_loop=[1:BL_nb]
                        subplot(BL_nb,HL_nb,(HL_loop-1)*HL_nb+BL_loop)
                        xlabel([P.par.blgrprank{BL_loop},': ',P.point.name{str2num(P.par.blgrprank{BL_loop})}])
                        ylabel([P.par.hlgrprank{HL_loop},': ',P.point.name{str2num(P.par.hlgrprank{HL_loop})}])
                        hold on
                        plot(P.comp.Xvariation(P.comp.incbindex,BL_loop,HL_loop),P.comp.Yvariation(P.comp.incbindex,BL_loop,HL_loop),'r.','markersize',3)
                        plot(P.comp.Xvariation(P.comp.decbindex,BL_loop,HL_loop),P.comp.Yvariation(P.comp.decbindex,BL_loop,HL_loop),'b.','markersize',3)
                        Max_X=max(P.comp.Xvariation(:,BL_loop,HL_loop));
                        Min_X=min(P.comp.Xvariation(:,BL_loop,HL_loop));
                        Max_Y=max(P.comp.Yvariation(:,BL_loop,HL_loop));
                        Min_Y=min(P.comp.Yvariation(:,BL_loop,HL_loop));
                        h=line([Max_X Min_X],[0 0]);
                        set(h,'color',[0 1 0],'linewidth',2)
                        h=line([0 0],[Max_Y Min_Y]);
                        set(h,'color',[1 1 0],'linewidth',2)
                        h=get(h,'parent');
                    end
                end
            end
        end
    end
    
    
    if ~isequal(P.par.comptype,'clustering')&Continue==1&P.flag.display==1
        if P.par.big==0
            FPR1_string=get(F.gh.prev.fdr_pm,'string');
            FPR1=findcellb(FPR1_string,'0.010');
            FPR1=find(FPR1);
            set(F.gh.prev.fdr_pm,'value',FPR1)
            COM_PREV('display selection')
        end
    end
    EndTime=clock;
    etime(EndTime,StartTime)
    %}}}
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    %}}}
    %))))))))))))))))))))))))))))))))))))
    
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %{{{test classecomp