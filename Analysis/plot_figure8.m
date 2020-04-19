function plot_figure8()
%% plots for headfixing effects: bouts and pauses
global fsize;
load('../Data/all_analyzed_structs_head.mat');
fsize=22;
minamp=2.5;
figure;
COL=colormap('lines');
close(gcf);
%% velcity /rotation for whisk / non whisk
% amp_off_v_da=[ctl_struct.all_diluted_amp_off];
amp_off_v_da=[ctl_struct.all_diluted_amp_off;comp_struct.all_diluted_amp_off];
hf_amp_off=hf_struct.all_diluted_amp_off;

ind_nw=find(amp_off_v_da(:,1)<=minamp);
ind_w=find(amp_off_v_da(:,1)>minamp);
hf_ind_nw=find(hf_amp_off(:,1)<=minamp);
hf_ind_w=find(hf_amp_off(:,1)>minamp);

A=plotCDFComp(amp_off_v_da(ind_nw,3),amp_off_v_da(ind_w,3),[0.5*COL(2,:);COL(2,:)]);
set(A(1),'XLim',[0 10]);

A=plotCDFComp(abs(amp_off_v_da(ind_nw,4)),abs(amp_off_v_da(ind_w,4)),[0.5*COL(2,:);COL(2,:)]);
set(A(1),'XLim',[0 120]);

A=plotCDFComp(amp_off_v_da(ind_nw,2),hf_amp_off(hf_ind_nw,2),[0.5*COL(2,:);0.5*COL(1,:)]);
set(A(1),'XLim',[30 120]);

A=plotCDFComp(amp_off_v_da(ind_w,2),hf_amp_off(hf_ind_w,2),[COL(2,:);COL(1,:)]);
set(A(1),'XLim',[30 120]);
end

function A=plotCDFComp(data1,data2,COL)
global fsize
% COL=colormap('lines');
mean1=mean(data1);
mean2=mean(data2);
% var1=var(data1);
% var2=var(data2);
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
