PsNb=length(Data{1}{1}.rank);
Frequency=cumsum(ones(1,PsNb))/PsNb;
h=figure;
set(h,'color',[1,1,1])
set(h,'name','DIFFERENT PLOT TYPE')
subplot(2,2,1)
hold on
XData=min(Data{1}{1}.signal,Data{2}{1}.signal);
plot(XData,abs(log2(Data{2}{1}.signal-log2(Data{1}{1}.signal))),'c.','markersize',3)
plot(sort(XData),Frequency*14,'k');
set(gca,'box','on')
set(gca,'ylim',[0,15])
title('log2(fold change) vs signal')
xlabel('min(signal(R1),signal(R2))')
ylabel('abs(log2(signal(R1))-log2(signal(R2)))')
subplot(2,2,2)
hold on
XData=log2(min(Data{1}{1}.signal,Data{2}{1}.signal));
plot(XData,abs(log2(Data{2}{1}.signal-log2(Data{1}{1}.signal))),'b.','markersize',3)
plot(sort(XData),Frequency*14,'k');
set(gca,'box','on')
set(gca,'ylim',[0,15])
title('log2(fold change) vs log2(signal)')
xlabel('log2(min(signal(R1),signal(R2)))')
ylabel('abs(log2(signal(R1))-log2(signal(R2)))')
subplot(2,2,3)
hold on
XData=min(Data{1}{1}.rank,Data{2}{1}.rank);
plot(XData,abs(log2(Data{2}{1}.signal-log2(Data{1}{1}.signal))),'g.','markersize',3)
plot(sort(XData),Frequency*14,'k');
set(gca,'box','on')
set(gca,'ylim',[0,15])
title('log2(fold change) vs rank')
xlabel('min(rank(R1),rank(R2))')
ylabel('abs(log2(signal(R1))-log2(signal(R2)))')
subplot(2,2,4)
hold on
XData=min(Data{1}{1}.rank,Data{2}{1}.rank);
plot(XData,abs(Data{2}{1}.rank-Data{1}{1}.rank),'r.','markersize',3)
plot(sort(XData),Frequency*50,'k');
set(gca,'box','on')
set(gca,'ylim',[0,60])
title('difference of ranks vs rank')
xlabel('min(rank(R1),rank(R2))')
ylabel('abs(rank(R1))-rank(R2))')