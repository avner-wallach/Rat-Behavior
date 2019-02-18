function plot_pauses_bouts(hf_struct,comp_struct,ctl_struct,wdata,all_sigs)
% figures for overall pauses and bouts
F=figure;
COL=colormap('lines');
close(F);
COL2=[218 124 48;62 150 81;146 36 40;107 76 154]/255;
minamp=2;
fsize=22;
%%  example of dynamics, envelope
dynamics=0;
if(dynamics)
wsig=wdata{54};
decomposed_struct=decompose_whisking_new(wsig,[],500);
window=[7.5 12.5];
ind=find(decomposed_struct.t>=window(1) & decomposed_struct.t<=window(2));
ywindow=[min(wsig(ind))*0.85 max(wsig(ind))*1.25];
t0=decomposed_struct.t(ind(1));
indw_start=find(diff(decomposed_struct.wamp(ind)<minamp)==1);
indw_end=find(diff(decomposed_struct.wamp(ind)<minamp)==-1);
if(numel(indw_start)<numel(indw_end))
    indw_start=[1 indw_start];
end
if(numel(indw_end)<numel(indw_start))
    indw_end=[indw_end numel(ind)];
end
if(indw_end(1)<indw_start(1))
    indw_end=[indw_end numel(ind)];
    indw_start=[1 indw_start];
end

A=axes;
for i=1:numel(indw_start)
    H=area([decomposed_struct.t(ind(indw_start(i))) decomposed_struct.t(ind(indw_end(i)))]-t0,[ywindow(2) ywindow(2)],ywindow(1));
    set(H,'FaceColor',[0 0 0],'FaceAlpha',0.05,'LineStyle','none');
    hold on;
end
x=decomposed_struct.t(ind)-t0;
H=plot(x,decomposed_struct.wsig(ind));
set(H,'Color',0.5*[1 1 1],'LineWidth',2);
H=plot(x,decomposed_struct.woff(ind));
set(H,'Color',[0 0 0],'LineWidth',3);
H=plot(x,decomposed_struct.woff(ind)+decomposed_struct.wamp(ind));
set(H,'Color',COL(4,:),'LineWidth',3);
H=plot(x,decomposed_struct.woff(ind)-decomposed_struct.wamp(ind));
set(H,'Color',COL(4,:),'LineWidth',3);

%amplitude
y0=ywindow(2);
H=plot(x,y0*ones(size(x)));
set(H,'Color',[0 0 0]);
y=decomposed_struct.wamp(ind)/max(decomposed_struct.wamp(ind))*15+y0;
my=max(y);
H=plot(x,y);
set(H,'Color',COL(4,:),'LineWidth',2);
ytick=get(A,'YTick');
yticklab=get(A,'YTickLabel');
indt=find(ytick<130);
ytick=[ytick(indt) y0 my];
yticklab={yticklab{indt} '0' num2str(round(max(decomposed_struct.wamp(ind))))};
set(A,'YTick',ytick,'YTickLabel',yticklab);
set(A,'Xlim',window-window(1));
set(A,'Ylim',[ywindow(1) my*1.01]);
set(A,'FontSize',fsize);
hold off;
end
%%  amplitude histograms hf vs. freely moving
idletime=1;
if(idletime)
free_wamp=[comp_struct.wamp;ctl_struct.wamp];
hf_whisk=mean(hf_struct.wamp>minamp);
comp_whisk=mean(comp_struct.wamp>minamp);
ctl_whisk=mean(ctl_struct.wamp>minamp);
free_whisk=mean(free_wamp>minamp);
amp_edges=[0:0.5:40];
amp_bins=mean([amp_edges(1:end-1);amp_edges(2:end)],1)
F=figure;
A=axes;
h=histcounts(hf_struct.wamp,amp_edges,'normalization','pdf');
ind=find(amp_bins<=2);
H=bar(amp_bins(ind),h(ind));
set(H,'FaceColor',0.65*COL(1,:),'BarWidth',1);
hold on;
ind=find(amp_bins>2);
H=bar(amp_bins(ind),h(ind));
set(H,'FaceColor',COL(1,:),'BarWidth',1);
% H=histogram(free_wamp,[0:1:40],'normalization','pdf');
% set(H,'FaceColor',0.5*[1 1 1]);
% H=histogram(comp_struct.wamp,[0:0.5:40],'normalization','pdf');
h=histcounts(ctl_struct.wamp,amp_edges,'normalization','pdf');
ind=find(amp_bins<=2);
H=bar(amp_bins(ind),h(ind));
set(H,'FaceColor',0.65*COL(2,:),'BarWidth',1);
hold on;
ind=find(amp_bins>2);
H=bar(amp_bins(ind),h(ind));
set(H,'FaceColor',COL(2,:),'BarWidth',1);

set(A,'FontSize',fsize,'XTick',[0 2 10:10:40]);
%inset
F=figure;
A=axes;
B=bar([1],100*[1-hf_whisk],'stacked');
B.FaceColor=COL(1,:);
% B.FaceAlpha=0.5;
hold on;
B=bar([2],100*[1-comp_whisk],'stacked');
B.FaceColor=COL(2,:);
% B=bar([2],100*[1-free_whisk],'stacked');
% B.FaceColor=0.5*[1 1 1];
% B.FaceAlpha=0.5;
set(A,'FontSize',round(fsize));
set(F,'Position',[680 558 150 150]);
set(A,'XTick',[],'XLim',[0 3],'Ylim',[0 30]);
1-hf_whisk
1-free_whisk
end
%% bout analysis
bouts=0;
if(bouts)
Mdur=max(hf_struct.all_bout_dur);
Mwhisk=max(hf_struct.all_bout_whisks);
whisk_lables=[0:10:Mwhisk];
F=figure;
A1=axes;
H=cdfplot(hf_struct.all_bout_dur);
set(H,'Color',COL2(1,:));
set(A1,'Xlim',[0 Mdur]);
set(A1,'FontSize',fsize);
dur_lables=get(A1,'XTickLabel');
dur_ticks=get(A1,'XTick');
H=cdfplot(hf_struct.all_bout_whisks);
set(H,'Color',COL2(2,:));
set(A1,'FontSize',fsize);
set(A1,'Xlim',[0 Mwhisk],'Box','off');
set(A1,'YTick',[],'XLabel',[],'YLabel',[],'Title',[]);
set(get(A1,'XAxis'),'Color',COL2(1,:),'TickDirection','out')
apos=get(A1,'Position');

A2=axes;
H1=cdfplot(hf_struct.all_bout_whisks/Mwhisk);
set(H1,'LineWidth',3,'Color',COL2(1,:));
hold on;
H2=cdfplot(hf_struct.all_bout_dur/Mdur);
set(H2,'LineWidth',3,'Color',COL2(2,:));
set(A2,'Xlim',[0 1]);
set(A2,'XAxisLocation','top','XTick',dur_ticks/Mdur,'XTickLabel',dur_lables);
set(get(A2,'XAxis'),'Color',COL2(2,:),'TickDirection','out')
set(A2,'Position',apos,'XLabel',[],'YLabel',[],'Title',[]);
set(A2,'XGrid','off','YGrid','off');
set(A2,'FontSize',fsize,'Box','off');
%inset markers
inset_d=[0.44 7.3 1.042];
y=get(H2,'YData');
x=get(H2,'XData');
for i=1:3;
    ind=find(x>=inset_d(i)/Mdur);
    H=plot(x(ind(1)),y(ind(1)),'.');
    set(H,'Color',0.5*[1 1 1],'MarkerSize',35);
end
%insets
ind=[5;6;33;36];
windows=[.2 1;.95 8;8.2 9.6;3.87 5.25]; 
for i=1:4;
    plot_inset(wdata{ind(i)},windows(i,:));
end

%inset: PDF in loglog to show heavy tail
binedges=1:50;
binctrs=1.5:1:49.5;
F=figure;
A=axes;
h=histogram(hf_struct.all_bout_whisks,binedges,'Normalization','pdf');
y=get(h,'Values');
f1=fit(binctrs(:),y(:),'Power1');
f2=fit(binctrs(:),y(:),'Exp1');
H=plot(binctrs,y,'.');
set(H,'MarkerSize',20,'Color',COL2(1,:));
hold on;
H=plot(binctrs,feval(f1,binctrs));
set(H,'LineWidth',2,'Color',0*[1 1 1])
H=plot(binctrs,feval(f2,binctrs),':');
set(H,'LineWidth',2,'Color',0.5*[1 1 1])
set(A,'XScale','log','YScale','log','FontSize',fsize,'YLim',[1e-3 1]);
set(F,'Position',[680 558 200 200]);
end
%%  example of dynamics, envelope
dynamics=0;
if(dynamics)
    F=figure;
    wsig=wdata{16};
    decomposed_struct=decompose_whisking_new(wsig,[],500);
    window=[0.07 0.75];
    ind=find(decomposed_struct.t>=window(1) & decomposed_struct.t<=window(2));
    ind_peaks=decomposed_struct.segmented_struct.ind_peak(decomposed_struct.segmented_struct.ind_peak>=ind(1) & ...
        decomposed_struct.segmented_struct.ind_peak<=ind(end));
    ind_troughs=decomposed_struct.segmented_struct.ind_trough(decomposed_struct.segmented_struct.ind_trough>=ind(1) & ...
        decomposed_struct.segmented_struct.ind_trough<=ind(end));
    ind_pumps=decomposed_struct.ind_pumps(decomposed_struct.ind_pumps>=ind(1) & ...
        decomposed_struct.ind_pumps<=ind(end));
    
    if(ind_troughs(1)<ind_peaks(1))
        ind_peaks=[1 ind_peaks];
    end
    if(ind_peaks(end)>ind_troughs(end))
        ind_troughs=[ind_troughs ind(end)];
    end
    ywindow=[min(wsig(ind))*0.9 max(wsig(ind))*2.1];
    t0=decomposed_struct.t(ind(1));
    dwsig=diff(decomposed_struct.wsig)*500/1000;
    A=axes;
    for i=1:numel(ind_peaks)
        H=area([decomposed_struct.t(ind_peaks(i)) decomposed_struct.t(ind_troughs(i))]-t0,[ywindow(2) ywindow(2)],ywindow(1));
        set(H,'FaceColor',[0 0 0],'FaceAlpha',0.05,'LineStyle','none');
        hold on;
    end

    H=plot(decomposed_struct.t(ind)-t0,decomposed_struct.wsig(ind));
    set(H,'LineWidth',3,'Color',[0 0 0]);    
    H=plot(decomposed_struct.t(ind_pumps)-t0,decomposed_struct.wsig(ind_pumps),'.');
    set(H,'MarkerSize',30,'Color',COL2(3,:));
    M=max(decomposed_struct.wsig(ind));
    m=min(decomposed_struct.wsig(ind));
    set(A,'FontSize',fsize,'YLim',[m M]);
    T1=[60 80 100];
    set(A,'YTick',T1);
    TL1=get(A,'YTickLabel');
    
%     hold off;

%     my=min(dwsig(ind));
%     ywindow=[my-.1*abs(my) max(dwsig(ind))*1.1];    
%     F=figure;
%     A=axes;
%     for i=1:numel(ind_peaks)
%         H=area([decomposed_struct.t(ind_peaks(i)) decomposed_struct.t(ind_troughs(i))]-t0,[ywindow(2) ywindow(2)],ywindow(1));
%         set(H,'FaceColor',[0 0 0],'FaceAlpha',0.05,'LineStyle','none');
%         hold on;
%     end
    y=dwsig(ind);
    %plot to get labels;
    F2=figure;
    A2=axes;
    H=plot(decomposed_struct.t(ind)-t0,y);
    set(H,'Color',[0 0 0],'LineWidth',3);
    set(A2,'FontSize',fsize,'YLim',[min(y) max(y)]);
    T2=get(A2,'YTick');
    TL2=get(A2,'YTickLabel');
    close(F2);

    MM=max(abs(y));
    y=(y/MM/3+1.4)*M;
    Y=dwsig(ind_pumps);
    Y=(Y/MM/3+1.4)*M;
    T=[T1 (T2/MM/3+1.4)*M];
    TL={TL1{:} TL2{:}};
    
    H=plot(decomposed_struct.t(ind)-t0,y);
    set(H,'Color',[0 0 0],'LineWidth',3);
    H=plot(decomposed_struct.t(ind_pumps)-t0,Y,'.');
    set(H,'MarkerSize',30,'Color',COL2(3,:));
    H=plot(decomposed_struct.t(ind)-t0,ones(size(ind))*1.4*M);
    set(H,'LineWidth',1,'Color',0.5*[1 1 1]);    
%     line(window-window(1),[0 0],'Color',0.5*[1 1 1]);
%     set(A,'XAxisLocation','origin')
%     set(A,'Xlim',window-window(1));
%     set(A,'Ylim',ywindow);
%     set(A,'FontSize',fsize);
    hold off;
    set(A,'Xlim',window-window(1));
    set(A,'Ylim',[m*.9 1.6*M]);
    set(A,'FontSize',fsize);
    set(A,'YTick',T,'YTickLabel',TL);
    set(A,'Box','off');
end


end

function plot_inset(wsig,window)
%insets
fsize=18;
COL=colormap('lines');
% minamp=2;
decomposed_struct=decompose_whisking_new(wsig,[],500);
ind=find(decomposed_struct.t>=window(1) & decomposed_struct.t<=window(2));
ywindow=[min(wsig(ind))*0.9 max(wsig(ind))*1.1];
t0=decomposed_struct.t(ind(1));

F=figure;
A=axes;
H=plot(decomposed_struct.t(ind)-t0,decomposed_struct.wsig(ind));
set(H,'Color',0.5*[1 1 1],'LineWidth',2);
set(A,'Xlim',window-window(1));
set(A,'Ylim',ywindow);
set(A,'FontSize',fsize,'XTick',[],'YTick',[]);
hold off;
set(A,'FontSize',round(fsize*.75));
set(F,'Position',[680 558 100*diff(window) 150]);

end
