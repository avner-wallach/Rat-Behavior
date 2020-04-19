function plot_figure7()
% figures for tirmming effect
global m M midoff minamp maxamp k fsize;
load('../Data/all_analyzed_structs_head.mat');
load('../Data/headfixed_c2h_data.mat');
m=round(min(cellfun(@nanmin,wdata)));
M=round(max(cellfun(@nanmax,wdata)));
midoff=(M+m)/2;
maxamp=(M-m)/2;
minamp=2.5;
fsize=22;
F=figure;
COL=colormap('lines');
close(F);
head=2;
%% trimming effect: offset, amplitude (figure 7)

A=plotCDFComp(ctl_struct.all_diluted_off,comp_struct.all_diluted_off,COL(2:3,:));
set(A(1),'XLim',[40 120]);
A=plotCDFComp(ctl_struct.all_diluted_amp,comp_struct.all_diluted_amp,COL(2:3,:));
set(A(1),'XLim',[2 30]);
end

function A=plotCDFComp(data1,data2,COL)
global fsize
mean1=mean(data1);
mean2=mean(data2);
var1=std(data1);
var2=std(data2);
pval_mean=boottest(data1,data2,1e4,@mean,'double')
pval_var=boottest(data1,data2,1e4,@var,'double')
F=figure;
A(1)=axes;
H=cdfplot(data1);
set(H,'Color',COL(1,:),'LineWidth',3);
hold on;
H=cdfplot(data2);
set(H,'Color',COL(2,:),'LineWidth',3);
set(A(1),'FontSize',fsize);
set(A(1),'XGrid','off');
set(A(1),'YGrid','off');
set(A(1),'Title',[],'XLabel',[],'YLabel',[]);

%insets
F=figure;
A(2)=axes;
B=bar([1],mean1);
B.FaceColor=COL(1,:);
hold on;
B=bar([2],mean2);
B.FaceColor=COL(2,:);
set(A(2),'FontSize',16);
set(F,'Position',[680 558 150 150]);
set(A(2),'XTick',[],'Xlim',[0 3]);
set(A(2),'YLim',[0 max(mean1,mean2)*1.25]);

F=figure;
A(3)=axes;
B=bar([1],var1);
B.FaceColor=COL(1,:);

hold on;
B=bar([2],var2);
B.FaceColor=COL(2,:);
set(A(3),'FontSize',16);
set(F,'Position',[680 558 150 150]);
set(A(3),'XTick',[],'Xlim',[0 3]);
set(A(3),'YLim',[0 max(var1,var2)*1.25]);
end
