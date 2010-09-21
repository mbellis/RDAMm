%&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% FUNCTION fdr
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% calculate fdr, sensitivity

%INPUT PARAMETERS
%1- Ppv: the product of p-values
%2- PpvPv: the p-value of Ppv (calculated in case of null hypothesis)
%3- PpvCdf: the observed cumulative distribution frequency of Ppv
%4- DisplayFlag: indicates if figures must be drawn or not
%5- FigH: the figure handle
%7- SubPos: subplot position

%OUTPUT PARAMETERS
%1- Fdr: False Discovery Rate of Ppv
%2- Fnr: False Undiscovery Rage of Ppv
%3- Sensitivity: sensitivity of Ppv
%4- Specificity: specificity of Ppv
%5- TruePosNb: estimate of the number of true variation

%VERSIONS
%V02 31-8-2010 Detect abnormal position of total variation estimation
%V01 22-3-2010 Refactoring of the existing version


function [Fdr,Sensitivity,TruePosNb]= fdr(ZVar,Ppv,PpvPv,PpvCdf,DisplayFlag,FigH,SubPos)

DataNb=length(PpvPv);

% These verifications should not be necessary (they were used before
% the refactoring)
% Verify that PpvPv is correct
if PpvPv(end)~=1
    PpvPv(end)=1;
end
%Verify that PpvPv and PpvCdf are  monotonous
% correct entry data which are flawed (some data
% with very low ppv could have their PpvCdf or PpvPv equal to one)
PpvPvMin=min(PpvPv);
PpvPvMinPos=find(PpvPv==PpvPvMin);
PpvPvMinPos=PpvPvMinPos(1);
if PpvPvMinPos>1 %abnormal position of minimum
    h=warndlg(sprintf('fdr: position of first min(PpvPv) >1 (%u)',PpvPvMinPos));         
    waitfor(h)
    PpvPv(1:PpvPvMinPos)=PpvPvMin;
end

PpvCdfMin=min(PpvCdf);
PpvCdfMinPos=find(PpvCdf==PpvCdfMin);
PpvCdfMinPos=PpvCdfMinPos(1);

if PpvCdfMinPos>1
    h=warndlg(sprintf('fdr: position of first min(PpvCdf) >1 (%u)',PpvCdfMinPos));
    waitfor(h)
    PpvCdf(1:PpvCdfMinPos)=PpvCdfMin;
end

%Conformation of PpvPv to prevent division by zero
ZeroIndex=find(PpvPv==0);
% to prevent from division by zero
PpvPv(ZeroIndex)=eps;

%ESTIMATION OF TRUE POSITIVE


%If ppvCdf>PpvPv in the last 5% portion of ppvCdf change the ppvCdf values because it is not possible that all
%probe sets have significant variation
PlusPos=find(PpvCdf(end-round(DataNb*0.05):end)>PpvPv(end-round(DataNb*0.05):end));
if ~isempty(PlusPos);
    % find the first inversion
    PlusPos=DataNb-(round(DataNb*0.05)-min(PlusPos))-2;   
    PpvCdf(PlusPos:DataNb)=PpvPv(PlusPos:DataNb);
end

% True positive are the difference between observed and expected
TpFreq=PpvCdf-PpvPv;

% figure
% plot(1./Ppv,(1-PpvCdf)./(1-PpvPv),'b')
% hold on
% plot(1./Ppv,TpFreq,'r')
% plot(1./Ppv, PpvCdf - ((1-PpvCdf)./(1-PpvPv)).*PpvPv,'g')
% plot(1./Ppv, (PpvCdf - ((1-PpvCdf)./(1-PpvPv)).*PpvPv)./TpFreq,'m')
% line([1,0.001],[0,0],'color','k')
% line([1,0.001],[1,1],'color','k')

% difference between Obs and Ho are highly significant for low value of Ho
% when computing TpFreq, starting from low value of Ho, TpFreq can't decrease
MTpFreq=make_monotonous(TpFreq,'inc');
% Number of true
TpNb=MTpFreq*DataNb;

% True positive can't be negative
NegIndex=find(TpNb<0);
 if ~isempty(NegIndex)
%     warndlg(['Exist ',sprintf('%.f',size(NegIndex,1)),' negative nb of Tp'],'STAT_CALIB_PREV')
     %%%MTpFreq(NegIndex)=0;
     TpNb(NegIndex)=0;
 end

%Find the  
StartPos=1;
EndPos=DataNb; 
FoundMax=0;
FirstRound=1;
while FoundMax==0    
    % Restrict search for Max nb of TP to positions before InvPos : pb if inversion at the beginning : find nothing !
    TruePosNb=max(TpNb(StartPos:EndPos));
    %TruePosNb=max(TpNb);
    if ~isempty(TruePosNb)
        %TruePosNb=TruePosNb(1);    
        if TruePosNb>0
            TpMaxNbPos=find(TpNb>=TruePosNb);
            FirstTpMaxNbPos=TpMaxNbPos(1);
            LastTpMaxNbPos=TpMaxNbPos(end);
        else
            FoundMax=1;
            FirstTpMaxNbPos=1;
            LastTpMaxNbPos=1;
        end
    else
        FoundMax=1;
        FirstTpMaxNbPos=1;
        LastTpMaxNbPos=1;
    end
    
    % refine the estimation of TP by convergence procedure 
    % the first estimate is the upper limit (all genes follow the null hypothesis)
    Fo=PpvPv(FirstTpMaxNbPos);
    F=PpvCdf(FirstTpMaxNbPos);
    if Fo==1
        Alpha=1;
    else
        Alpha=(1-F)/(1-Fo);
        if Alpha>1
            Alpha=1;
        end
    end
    
    TpFreq=PpvCdf-Alpha.*PpvPv;
    MTpFreq=make_monotonous(TpFreq,'inc');
    %NewTpNb=TpFreq(FirstTpMaxNbPos)*DataNb;
    %NewTpNb=TpFreq*DataNb;
    NewTpNb=MTpFreq*DataNb;
    if TpNb(FirstTpMaxNbPos)>0
        %TpRatio=NewTpNb/TpNb(FirstTpMaxNbPos);
        TpRatio=NewTpNb(FirstTpMaxNbPos)/TpNb(FirstTpMaxNbPos);
    else
        TpRatio=1;
    end
    
    if TpRatio<2
        FoundMax=1;
    else
        % eliminate intervale containing this current max by searching the next inversion
        if FirstRound==1 % construct smoothed F-Fo
            [XGrid,GridTpFreq]=gpanal_calib([log10(1./Ppv),(PpvCdf-PpvPv)],[0,ceil(max(log10(1./Ppv))),501],10);
            InterpXData=log(1./Ppv);
            XDataDiff=diff(InterpXData);
            ZeroIndex=find(XDataDiff==0);
            while ~isempty(ZeroIndex)
                InterpXData(ZeroIndex)=[];
                XDataDiff=diff(InterpXData);
                ZeroIndex=find(XDataDiff==0);
            end
            
            SmoothDegree=2;
            GridTpFreqDiff=diff(GridTpFreq(:,SmoothDegree));
            ZeroIndex=find(GridTpFreqDiff==0);
            while ~isempty(ZeroIndex)
                GridTpFreq(ZeroIndex,:)=[];
                XGrid(ZeroIndex)=[];
                GridTpFreqDiff=diff(GridTpFreq(:,SmoothDegree));
                ZeroIndex=find(GridTpFreqDiff==0);
            end
            
            SmoothTpFreq=interp1(XGrid,GridTpFreq(:,SmoothDegree),InterpXData);                           
            NegSlopeIndex=find([diff(SmoothTpFreq)<0;0]&SmoothTpFreq>=0);
            FirstRound=0;            
        end
        % search next interval starting with an increasing TpFreq
        NextPos=find(NegSlopeIndex<FirstTpMaxNbPos);
        if ~isempty(NextPos)
            EndPos=NegSlopeIndex(NextPos(end));
%             plot(10.^InterpXData(EndPos),SmoothTpFreq(EndPos),'mo')
%             plot(10.^InterpXData(EndPos),SmoothTpFreq(EndPos),'m+')
        else
            FoundMax=1;
            FirstTpMaxNbPos=1;
            LastTpMaxNbPos=1;
        end    
    end
end % of while FoundMax==0 

PpvPvMem=PpvPv;
PpvPv=PpvPv*Alpha;

%RECALCULATE TP AFTER CORRECTION OF FO
% True positive are the difference between observed and expected
TpFreq=PpvCdf-PpvPv;
MTpFreq=make_monotonous(TpFreq,'inc');
%take floor for TruePosNb and ceil for TpNb to make sure that the first
%position is not misplaced at the end
TruePosNb=floor(max(MTpFreq*DataNb));
TpNb=ceil(MTpFreq*DataNb);
NegIndex=find(TpNb<0);
if ~isempty(NegIndex)
    TpNb(NegIndex)=0;
end


%TruePosNb=max(TpNb);
if ~isempty(TruePosNb)
    %TruePosNb=TruePosNb(1);
    if TruePosNb>0
        TpMaxNbPos=find(TpNb>=TruePosNb);
        FirstTpMaxNbPos=TpMaxNbPos(1);
        LastTpMaxNbPos=TpMaxNbPos(end);
    else
        FirstTpMaxNbPos=1;
        LastTpMaxNbPos=1;
    end
else
    FirstTpMaxNbPos=1;
    LastTpMaxNbPos=1;
end


TruePosNb=TpNb(FirstTpMaxNbPos);
TpNb(FirstTpMaxNbPos:end)=TruePosNb;
%}}}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{{{ Estimation of false positive FpNb

% estimated False Positive Number
FpNb=(PpvCdf*DataNb)-TpNb;

%FpNb(LastTpMaxNbPos:end)=(PpvCdf(LastTpMaxNbPos:end)*DataNb)-TruePosNb;

% should not exist
%Correction_index=find(FpNb>TnMaxNb);
%size(Correction_index)
%FpNb(Correction_index)=TnMaxNb;
%}}}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{{{ Estimation of false negative FnNb
FnNb=TruePosNb-TpNb;
% should not exist
%Correction_index=find(FnNb>TruePosNb);
%size(Correction_index)
%FnNb(Correction_index)=TruePosNb;
%}}}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{{{ Estimation of true negative TnNb
TnNb=DataNb-(TpNb+FnNb+FpNb);
%Correction_index=find(TnNb>TnMaxNb);
%TnNb(Correction_index)=TnMaxNb;
%}}}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{{{ Correction of limit value (0,1)
CorrFpNb=FpNb;
CorrTpNb=TpNb;
CorrFnNb=FnNb;
CorrTnNb=TnNb;

ZeroIndex=find(FpNb==0);
CorrFpNb(ZeroIndex)=eps;
ZeroIndex=find(TpNb==0);
CorrTpNb(ZeroIndex)=eps;
ZeroIndex=find(FnNb==0);
CorrFnNb(ZeroIndex)=eps;
ZeroIndex=find(TnNb==0);
CorrTnNb(ZeroIndex)=eps;

%}}}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{{{ Estimation of statistical parameters (FDR, FNR, Sensitivity, Specificity, Chi2)

FDR=FpNb./(CorrFpNb+CorrTpNb);            
FDR(FDR>1)=1;
FDR(FDR<0)=0;

FNR=FnNb./(CorrFnNb+CorrTnNb);
FNR(FNR>1)=1;
FNR(FNR<0)=0;
%Fdr=make_monotonous(FDR,'inc');
Fdr=make_monotonous(FDR,'dec',1);
Fnr=make_monotonous(FNR,'dec');

Fdr10pc=Fdr>=0.10;
if ~isempty(find(Fdr10pc))
    Fdr10pc=find(Fdr10pc);
    Fdr10pc=Fdr10pc(1);
else
    if Fdr(1)>0.10
        Fdr10pc=1;
    else
        Fdr10pc=max(size(FDR));
    end
end



PositiveNb=TpNb+FnNb;
ZeroIndex=find(PositiveNb==0);
PositiveNb(ZeroIndex)=eps;
SENS=TpNb./(PositiveNb);
% When True positive is small, Sensitivity could be superior to one for small value of false negative
SENS(SENS>1)=1;
Sensitivity=make_monotonous(SENS,'inc');

NegativeNb=TnNb+FpNb;
ZeroIndex=find(NegativeNb==0);
NegativeNb(ZeroIndex)=eps;
SPEC=TnNb./(NegativeNb);            
% Same correction but very much more improbable (TnNb = zero)
SPEC(SPEC>1)=1;
Specificity=make_monotonous(SPEC,'dec');


%GRAPHICS
if DisplayFlag==1
%     figure(FigH)   
%     set(FigH,'name',Title)    
%     subplot(1,2,SubPos)    
%     if SubPos==1       
%         title('FDR(-), Sensitivity(..), cdf((ppv)|Ho)(black -) and cdf(ppv)(--)')    
%         xlabel(sprintf('Inc: %.f    log(1/product of the plevels)',TruePosNb))
%     else
%         xlabel(sprintf('Dec: %.f    log(1/product of the plevels)',TruePosNb))
%     end
%     set(gca,'xscale','log')
%     set(gca,'box','on')    
%     hold on   
%     %  
%     % expected cdf of Ppv under Ho hypothesis
%     plot(1./Ppv,PpvPv,'k-');
%     plot(1./Ppv,PpvPvMem,'k:');
%     
%     % observed cdf of Ppv
%     eval(['plot(1./Ppv,PpvCdf,''',CurrColor,'--'');']);   
%     
%     % indicates the position of TP max nb    
%     plot(1./Ppv(FirstTpMaxNbPos),PpvCdf(FirstTpMaxNbPos),sprintf('%co',CurrColor))
%     plot(1./Ppv(FirstTpMaxNbPos),PpvCdf(FirstTpMaxNbPos),sprintf('%c+',CurrColor))
%     line([1./Ppv(FirstTpMaxNbPos),1./Ppv(FirstTpMaxNbPos)],[0,Fdr(FirstTpMaxNbPos)],'color',CurrColor)
%     line([1,1./Ppv(FirstTpMaxNbPos)],[Fdr(FirstTpMaxNbPos),Fdr(FirstTpMaxNbPos)],'color',CurrColor)    
%    line([1./Ppv(Fdr10pc),1./Ppv(Fdr10pc)],[0,Sensitivity(Fdr10pc)],'color',CurrColor,'linestyle',':')
%    line([1,1./Ppv(Fdr10pc)],[Sensitivity(Fdr10pc),Sensitivity(Fdr10pc)],'color',CurrColor,'linestyle',':')
%           
%     % Fdr & Sensibility
%     plot(1./Ppv,Fdr,sprintf('%c-',CurrColor))
%     plot(1./Ppv,Sensitivity,sprintf('%c:',CurrColor))
%     
    figure(FigH)   
    subplot(1,2,SubPos)
    hold on   
    if SubPos==1                     
        Label=sprintf('zVar of Inc (%.f)',TruePosNb);
        CurrColor='r';
        title('                                                             FDR(-), Sensitivity(..), cdf((ppv)|Ho)(black -) and cdf(ppv)(--)')    
    else
        Label=sprintf('zVar of Dec (%.f)',TruePosNb);
        CurrColor='b';
    end    
    %  
    % expected cdf of Ppv under Ho hypothesis
    plot(ZVar,PpvPv,'k-');
    plot(ZVar,PpvPvMem,'k:');
    
    % observed cdf of Ppv
    eval(['plot(ZVar,PpvCdf,''',CurrColor,'--'');']);   
    
    % indicates the position of TP max nb    
    plot(ZVar(FirstTpMaxNbPos),PpvCdf(FirstTpMaxNbPos),sprintf('%co',CurrColor))
    plot(ZVar(FirstTpMaxNbPos),PpvCdf(FirstTpMaxNbPos),sprintf('%c+',CurrColor))
    line([ZVar(FirstTpMaxNbPos),ZVar(FirstTpMaxNbPos)],[0,Fdr(FirstTpMaxNbPos)],'color',CurrColor)
    line([1,ZVar(FirstTpMaxNbPos)],[Fdr(FirstTpMaxNbPos),Fdr(FirstTpMaxNbPos)],'color',CurrColor)    
   line([ZVar(Fdr10pc),ZVar(Fdr10pc)],[0,Sensitivity(Fdr10pc)],'color',CurrColor,'linestyle',':')
   line([1,ZVar(Fdr10pc)],[Sensitivity(Fdr10pc),Sensitivity(Fdr10pc)],'color',CurrColor,'linestyle',':')
          
    % Fdr & Sensibility
    plot(ZVar,Fdr,sprintf('%c-',CurrColor))
    plot(ZVar,Sensitivity,sprintf('%c:',CurrColor))
    
    set(gca,'box','on')    
    xlabel(Label)
    set(gcf,'color',[1,1,1])
end
