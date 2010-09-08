% calculate and display some statitistics from cdf(obs) & cdf(theorical)
% type I error, type II error, sensibility, specificity, Chi2
% !!! INACTIVATE line  188

function [MFDR,MFNR,MSensitivity,MSpecificity,TpMaxNb,Max_useful_FPR]= fdr(DataNb,...
    HoCdf,ObsCdf,GraphFlag,FigH,Color,Title,XData,SubPos)
%GraphFlag=0;
% Value is inverse_sorted
% ObsCdf is ObsCdf of Value and phased with them (greatest value in first position, assigned with the lowest ObsCdf value)
% HoCdf of Value idem
%((((((((((((((((((((((((((((((((((((
%{{{ STATISTICS
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{{{ Verify that HoCdf is correct
if HoCdf(end)~=1
    HoCdf(end)=1;
end
%}}}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{{{ Verify that HoCdf and ObsCdf are  monotonous
% exist there a patch to correct entry data which are flawed (some data with very low product of variation pvalues and
% observed cdf of ppv or pvalue of ppv|Ho equal to one !

HoCdfMin=min(HoCdf);
HoCdfMin=HoCdfMin(1);
HoCdfMinPos=find(HoCdf==HoCdfMin);
HoCdfMinPos=HoCdfMinPos(1);
if HoCdfMinPos>1 %abnormal position of minimum
    %warndlg(['position of first min(HoCdf) >1 in STAT_CALIB : ',sprintf('%.f',HoCdfMinPos)],'STAT_CALIB_PREV');         
    HoCdf(1:HoCdfMinPos)=HoCdfMin;
end


ObsCdfMin=min(ObsCdf);
ObsCdfMin=ObsCdfMin(1);
ObsCdfMinPos=find(ObsCdf==ObsCdfMin);
ObsCdfMinPos=ObsCdfMinPos(1);
if ObsCdfMinPos>1
    %warndlg(['position of first min(ObsCdf) >1 in STAT_CALIB : ',sprintf('%.f',ObsCdfMinPos)],'STAT_CALIB_PREV');
    ObsCdf(1:ObsCdfMinPos)=ObsCdfMin;
end

%}}}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{{{ Conformation of HoCdf to prevent division by zero
ZeroIndex=find(HoCdf==0);
%OneIndex=find(HoCdf==1);
% to prevent from division by zero
HoCdf(ZeroIndex)=eps;
%HoCdf(OneIndex)=1-eps; 
%}}}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{{{ Estimation of true positive Tp

% True positive are the difference between observed and expected
TpFreq=ObsCdf-HoCdf;
% figure
% plot(1./XData,(1-ObsCdf)./(1-HoCdf),'b')
% hold on
% plot(1./XData,TpFreq,'r')
% plot(1./XData, ObsCdf - ((1-ObsCdf)./(1-HoCdf)).*HoCdf,'g')
% plot(1./XData, (ObsCdf - ((1-ObsCdf)./(1-HoCdf)).*HoCdf)./TpFreq,'m')
% line([1,0.001],[0,0],'color','k')
% line([1,0.001],[1,1],'color','k')

% difference between Obs and Ho are highly significant for low value of Ho
% when computing TpFreq, starting from low value of Ho, TpFreq can't decrease
MTpFreq=MAXSTEP_PREV(TpFreq,0);
% Number of true
TpNb=MTpFreq*DataNb;

% True positive can't be negative
NegIndex=find(TpNb<0);
 if ~isempty(NegIndex)
%     warndlg(['Exist ',sprintf('%.f',size(NegIndex,1)),' negative nb of Tp'],'STAT_CALIB_PREV')
     %%%MTpFreq(NegIndex)=0;
     TpNb(NegIndex)=0;
 end

StartPos=1;
EndPos=DataNb; 
FoundMax=0;
FirstRound=1;
while FoundMax==0
    %InvPos=PpvLimit;
    
    % Restrict InvPos to low values of XData
    %if InvPos>round(DataNb/10)
    %    InvPos=round(DataNb/10);
    %end
    % Restrict search for Max nb of TP to positions before InvPos : pb if inversion at the beginning : find nothing !
    TpMaxNb=max(TpNb(StartPos:EndPos));
    %TpMaxNb=max(TpNb);
    if ~isempty(TpMaxNb)
        TpMaxNb=TpMaxNb(1);    
        if TpMaxNb>0
            TpMaxNbPos=find(TpNb>=TpMaxNb);
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
    Fo=HoCdf(FirstTpMaxNbPos);
    F=ObsCdf(FirstTpMaxNbPos);
    if Fo==1
        Alpha=1;
    else
        Alpha=(1-F)/(1-Fo);
        if Alpha>1
            Alpha=1;
        end
    end
    
    TpFreq=ObsCdf-Alpha.*HoCdf;
    NewTpNb=TpFreq(FirstTpMaxNbPos)*DataNb;
    if TpNb(FirstTpMaxNbPos)>0
        TpRatio=NewTpNb/TpNb(FirstTpMaxNbPos);
    else
        TpRatio=1;
    end
    
    if TpRatio<2
        FoundMax=1;
    else
        % eliminate intervale containing this current max by searching the next inversion
        if FirstRound==1 % construct smoothed F-Fo
            [XGrid,GridTpFreq]=gpanal_calib([log10(1./XData),(ObsCdf-HoCdf)],[0,ceil(max(log10(1./XData))),501],10);
            InterpXData=log(1./XData);
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
            
            
%             figure
%             plot(1./XData,(1-ObsCdf)./(1-HoCdf),'b')
%             hold on
%             plot(1./XData,ObsCdf-HoCdf,'r')
%             plot(1./XData, ObsCdf - ((1-ObsCdf)./(1-HoCdf)).*HoCdf,'g')
%             plot(1./XData, (ObsCdf - ((1-ObsCdf)./(1-HoCdf)).*HoCdf)./(ObsCdf-HoCdf),'m')                        
%             plot(10.^InterpXData,SmoothTpFreq,'k')
%             line([1,0.001],[0,0],'color','k')
%             line([1,0.001],[1,1],'color','k')

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
end

% !!!!!THE NEXT LINE MUST BE INACTIVATED !!!!
% FOR testing FDR estimator quality one impose the true number of variant genes
% in case of artificial data
%Alpha=(DataNb-500)/DataNb;


HoCdfMem=HoCdf;
HoCdf=HoCdf*Alpha;

%RECALCULATE TP AFTER CORRECTION OF FO
% True positive are the difference between observed and expected
TpFreq=ObsCdf-HoCdf;
MTpFreq=MAXSTEP_PREV(TpFreq,0);
TpNb=MTpFreq*DataNb;
NegIndex=find(TpNb<0);
 if ~isempty(NegIndex)
     MTpFreq(NegIndex)=0;
     TpNb(NegIndex)=0;
 end
TpMaxNb=TpNb(FirstTpMaxNbPos);

% calculate Standard error to bracket the estimation of True Positive
TpFreqSE=sqrt((HoCdf.*(1-HoCdf)+ObsCdf.*(1-ObsCdf))/DataNb);
TpFreqMin=MTpFreq-TpFreqSE*1.96;
TpFreqMax=MTpFreq+TpFreqSE*1.96;
TpNbSE=TpFreqSE*DataNb;
TpNbSup=TpNb+TpNbSE*1.96;
TpNbInf=TpNb-TpNbSE*1.96;

TpNb(FirstTpMaxNbPos:end)=TpMaxNb;
TpNbSE(FirstTpMaxNbPos:end)=TpNbSE(LastTpMaxNbPos);
TpMaxNbSE=TpNbSE(LastTpMaxNbPos);
TpMaxNbSup=TpMaxNb+TpMaxNbSE*1.96;
TpMaxNbInf=TpMaxNb-TpMaxNbSE*1.96;
%}}}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{{{ Estimation of maximum nb of true negative TnMaxNb
TnMaxNb=DataNb-TpMaxNb;
TnMaxNbSup=DataNb-TpMaxNbInf;
TnMaxNbInf=DataNb-TpMaxNbSup;
%}}}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{{{ Estimation of false positive FpNb

% Expected False Positive Number
ExpFpNb=HoCdf*DataNb;
% estimated False Positive Number
FpNb=(ObsCdf*DataNb)-TpNb;

%FpNb(LastTpMaxNbPos:end)=(ObsCdf(LastTpMaxNbPos:end)*DataNb)-TpMaxNb;
FpNbSup=FpNb+TpNbSE*1.96;
FpNbInf=FpNb-TpNbSE*1.96;

% should not exist
%Correction_index=find(FpNb>TnMaxNb);
%size(Correction_index)
%FpNb(Correction_index)=TnMaxNb;
%}}}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{{{ Estimation of false negative FnNb
FnNb=TpMaxNb-TpNb;
FnNbSup=max(TpMaxNbSup-TpNbInf);
FnNbInf=max(TpMaxNbInf-TpNbSup);
% should not exist
%Correction_index=find(FnNb>TpMaxNb);
%size(Correction_index)
%FnNb(Correction_index)=TpMaxNb;
%}}}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{{{ Estimation of true negative TnNb
TnNb=DataNb-(TpNb+FnNb+FpNb);
TnNbSup=DataNb-(TpNbInf+FnNbInf+FpNbInf);
TnNbInf=DataNb-(TpNbSup+FnNbSup+FpNbSup);
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

CorrFpNbSup=FpNbSup;
CorrTnNbSup=TnNbSup;
ZeroIndex=find(FpNbSup==0);
CorrFpNbSup(ZeroIndex)=eps;
ZeroIndex=find(TnNbSup==0);
CorrTnNbSup(ZeroIndex)=eps;

CorrFpNbInf=FpNbInf;
CorrTnNbInf=TnNbInf;
ZeroIndex=find(FpNbInf<=0);
CorrFpNbInf(ZeroIndex)=eps;
ZeroIndex=find(TnNbInf<=0);
CorrTnNbInf(ZeroIndex)=eps;
%}}}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{{{ Estimation of statistical parameters (FDR, FNR, Sensitivity, Specificity, Chi2)

FDR=FpNb./(CorrFpNb+CorrTpNb);            
FDR(FDR>1)=1;
FDRSup=FpNbSup./(CorrFpNb+CorrTpNb);
FDRInf=FpNbInf./(CorrFpNb+CorrTpNb);

FNR=FnNb./(CorrFnNb+CorrTnNb);
FNR(FNR>1)=1;
MFDR=MINSTEP(FDR,1);
MFDRSup=MINSTEP(FDRSup,1);
MFDRInf=MINSTEP(FDRInf,1);
MFNR=MINSTEP(FNR,0);

FDR10pc=MFDR>=0.10;
if ~isempty(find(FDR10pc))
    FDR10pc=find(FDR10pc);
    FDR10pc=FDR10pc(1);
else
    if MFDR(1)>0.10
        FDR10pc=1;
    else
        FDR10pc=max(size(FDR));
    end
end

%  h=figure;
%  set(h,'name','effect of monotonic transformation on FDR(r), and FNR(b)')
%  hold on
%  plot(XData,FDR,'r')
%  plot(XData,MFDR,'b')
% plot(XData,FNR,'b+')
% plot(XData,MFNR,'bo')
% plot(FDR,MFDR,'m.')
% plot(FNR,MFNR,'c.')


PositiveNb=TpNb+FnNb;
ZeroIndex=find(PositiveNb==0);
PositiveNb(ZeroIndex)=eps;
Sensitivity=TpNb./(PositiveNb);
% When True positive is small, Sensitivity could be superior to one for small value of false negative
Sensitivity(Sensitivity>1)=1;
MSensitivity=MAXSTEP(Sensitivity,0);

NegativeNb=TnNb+FpNb;
ZeroIndex=find(NegativeNb==0);
NegativeNb(ZeroIndex)=eps;
Specificity=TnNb./(NegativeNb);            
% Same correction but very much more improbable (TnNb = zero)
Specificity(Specificity>1)=1;
MaxSpecificity=TnNbSup./(NegativeNb);
MinSpecificity=TnNbInf./(NegativeNb);
MSpecificity=MINSTEP(Specificity,0);
MMinSpecificity=MINSTEP(MinSpecificity,0);
MMaxSpecificity=MINSTEP(MaxSpecificity,0);
% used only for subplot not displayed presently
% % % Chi2=DataNb*((TpFreq.^2./HoCdf)+(TpFreq.^2./(1-HoCdf)));
% % % Chi2Stat=1-chi2Cdf(Chi2,1);
%}}}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%}}}
%))))))))))))))))))))))))))))))))))))
Max_useful_FPR=MFDR(FirstTpMaxNbPos);
%((((((((((((((((((((((((((((((((((((
%{{{ GRAPHICS
if GraphFlag==1
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %{{{ Subplot 1
    figure(FigH)   
    set(FigH,'name',Title)    
    subplot(1,2,SubPos)
    if SubPos==1
        title('FDR(-), Sensitivity(..), cdf((ppv)|Ho)(black -) and cdf(ppv)(--)')    
        xlabel(sprintf('Inc: %.f    log(1/product of the plevels)',TpMaxNb))
    else
        xlabel(sprintf('Dec: %.f    log(1/product of the plevels)',TpMaxNb))
    end
    set(gca,'xscale','log')
    hold on   
    %  
    % expected cdf of XData under Ho hypothesis
    plot(1./XData,HoCdf,'k-');
    plot(1./XData,HoCdfMem,'k:');
    
    % observed cdf of XData
    eval(['plot(1./XData,ObsCdf,''',Color,'--'');']);   
    
    % indicates the position of TP max nb    
    plot(1./XData(FirstTpMaxNbPos),ObsCdf(FirstTpMaxNbPos),sprintf('%co',Color))
    plot(1./XData(FirstTpMaxNbPos),ObsCdf(FirstTpMaxNbPos),sprintf('%c+',Color))
    line([1./XData(FirstTpMaxNbPos),1./XData(FirstTpMaxNbPos)],[0,MFDR(FirstTpMaxNbPos)],'color',Color)
    line([1,1./XData(FirstTpMaxNbPos)],[MFDR(FirstTpMaxNbPos),MFDR(FirstTpMaxNbPos)],'color',Color)    
   line([1./XData(FDR10pc),1./XData(FDR10pc)],[0,MSensitivity(FDR10pc)],'color',Color,'linestyle',':')
   line([1,1./XData(FDR10pc)],[MSensitivity(FDR10pc),MSensitivity(FDR10pc)],'color',Color,'linestyle',':')
          
    % MFDR & Sensibility
    plot(1./XData,MFDR,sprintf('%c-',Color))
    plot(1./XData,MSensitivity,sprintf('%c:',Color))

    %}}}
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    Skip_following=1;
    if Skip_following==0
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        %{{{ Subplot 2
        subplot(2,2,2)
        title('nb of TP(-), TN(..), FP(-.) and FN(--)')
        xlabel('log (1/product of pvalues))')
        set(gca,'xscale','log')
        hold on
        % "cumulative" nb of TP
        eval(['plot(1./XData,TpNb,''',Color,'-'');']);        
        % raw values of TP nb
        % eval(['plot(1./XData,TpFreq*DataNb,''',Color,'.'',''markersize'',2);']);
        plot(1./XData,TnNb,sprintf('%c:',Color))
        plot(1./XData,FpNb,sprintf('%c-.',Color))
        plot(1./XData,FnNb,sprintf('%c--',Color))
        %}}}
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        %{{{ Subplot 3
        subplot(2,2,3)
        title('Sensitivity(-), Specificity(..), FDR(--), FRR(:)')
        xlabel('Inverse of the product of the plevels')
        set(gca,'xscale','log')
        hold on
        eval(['plot(1./XData,MSensitivity,''',Color,'-'');']);
        %hold on
        plot(1./XData,MSpecificity,sprintf('%c:',Color))
        plot(1./XData,MFDR,sprintf('%c--',Color))
        plot(1./XData,MFNR,sprintf('%c:',Color))
        %eval(['plot(1./XData,Chi2Stat,''',Color,'--'');']);
        %}}} 
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        %{{{ Subplot 4   
        HoCdfDist=RANKED(HoCdf,0);    
        HoCdfDist=HoCdfDist*DataNb;
        MFDRDist=RANKED(MFDR,0);
        MFDRDist=MFDRDist*DataNb;
        MFDRInfDist=RANKED(MFDRInf,0);
        MFDRInfDist=MFDRInfDist*DataNb;
        MFDRSupDist=RANKED(MFDRSup,0);
        MFDRSupDist=MFDRSupDist*DataNb;   
        MSpecificityDist=RANKED(MSpecificity,0);
        MSpecificityDist=MSpecificityDist*DataNb;
        MMinSpecificityDist=RANKED(MMinSpecificity,0);
        MMinSpecificityDist=MMinSpecificityDist*DataNb;
        MMaxSpecificityDist=RANKED(MMaxSpecificity,0);
        MMaxSpecificityDist=MMaxSpecificityDist*DataNb;
        
        subplot(2,2,4)
        title('cdf of p-value of ppv(-), FP rate(--), Specificity(-.)')
        xlabel('log (1/value)')
        set(gca,'xscale','log')
        hold on
        MFDRIndex=find(FDR>0);
        MSpecificityIndex=find(MSpecificity>0);
        plot((1./HoCdf),HoCdfDist,sprintf('%c-',Color))
        plot(1./MFDR(MFDRIndex),MFDRDist(MFDRIndex),sprintf('%c--',Color))
        plot(1./MSpecificity(MSpecificityIndex),MSpecificityDist(MSpecificityIndex),sprintf('%c-.',Color))
        
        X1_max=nanmax(1./HoCdf);
        X1_max=X1_max(1);
        X2_max=nanmax(1./MFDR(MFDRIndex));
        X2_max=X1_max(1);
        X3_max=nanmax(1./MSpecificity(MSpecificityIndex));
        X3_max=X1_max(1);
        X_max=max([X1_max,X2_max,X3_max]);
        X_max=X_max(1);
        PointPos=10.^[0:log10(X_max)];
        plot(PointPos,ones(size(PointPos))*TpMaxNb,sprintf('%c-',Color))
        plot(PointPos,ones(size(PointPos))*TpMaxNbSup,sprintf('%c:',Color))
        plot(PointPos,ones(size(PointPos))*TpMaxNbInf,sprintf('%c:',Color))
    end   
    
    %TpMaxNb=round(TpMaxNb);
    Max_useful_FPR=MFDR(FirstTpMaxNbPos);
    %}}} 
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
end
%}}}
%)))))))))))))))))))))))))))))))))))) 
