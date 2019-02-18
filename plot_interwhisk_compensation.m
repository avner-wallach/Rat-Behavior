function plot_interwhisk_compensation(hf_struct,comp_struct,ctl_struct)
% figures for spatial compensation
global m M midoff minamp maxamp k;
m=round(min(hf_struct.amp_off(:,2)));
M=round(max(hf_struct.amp_off(:,2)));
midoff=(M+m)/2;
maxamp=(M-m)/2;
minamp=2;
F=figure;
COL=colormap('lines');
close(F);
%%  xcov plots
xcov_plot=0;
varnames={'amp','off','dur'};
if(xcov_plot)
    for i=1:numel(varnames)
        plotXcor({hf_struct,comp_struct,ctl_struct},varnames{i});
    end
end

%% cdf plot
cdfs=0;
if(cdfs)
A=plotCDFComp(hf_struct.all_diluted_prot*1e3,comp_struct.all_diluted_prot*1e3,ctl_struct.all_diluted_prot*1e3);
set(A(1),'XLim',[25 150]);
A=plotCDFComp(hf_struct.all_diluted_ret*1e3,comp_struct.all_diluted_ret*1e3,ctl_struct.all_diluted_ret*1e3);
set(A(1),'XLim',[25 150]);
end

%% duration vs. lag
dur_corr=1;
varnames={'amp','off','dur'};
if(dur_corr)
    for i=1:numel(varnames)
        plotDurCorr({hf_struct,comp_struct,ctl_struct},varnames{i});
    end
end

end
function plotXcor(strcts,varname)
COL=colormap('lines');
x_field=['all_X_',varname];
l_field=['all_lag_',varname];
for i=1:numel(strcts)
    z=getfield(strcts{i},x_field);
    l=getfield(strcts{i},l_field);
    s=strcts{i}.segmented_diluted_length;
    x=[1:size(z,2)];
    y=[1:size(z,1)];
    [lags{i},I]=sort(l);
    X=linspace(1,x(end),1e3);
    Y=linspace(1,y(end),1e3+1);
    Z=interp2(x,y,z(I,:),X,Y','linear');
    [XX,YY]=meshgrid(X,Y);

    F=figure;
    A=axes;
%     imagesc(x(I,:));
    surf(XX,YY,zeros(size(Z)),'CData',Z,'LineStyle','none');
    view(0,90);
    hold on;
    H=plot(lags{i},y,'--w');
    set(H,'LineWidth',3);
    set(A,'FontSize',14);
    set(A,'Xlim',[1,20]);
    set(A,'Ylim',[1,y(end)]);
    set(A,'YTick',[]);
    colormap('hot');
    set(A,'XGrid','off','YGrid','off');
    set(F,'Position',[600,486,300,420])
end
% 
% figure;
% plot(s(I)
F=figure;
A=axes;
for i=1:numel(strcts)
    H=cdfplot(lags{i});
    mean_l(i)=mean(lags{i});
    median_l(i)=median(lags{i});
    set(H,'Color',COL(i,:),'LineWidth',2);
    hold on;
end
set(A,'FontSize',14);
set(A,'XGrid','off','Xlim',[0,20]);
set(A,'YGrid','off');
set(A,'Title',[],'XLabel',[],'YLabel',[]);

pval_mean=[boottest(lags{1},lags{2},5e3,@mean,'double'),boottest(lags{1},lags{3},5e3,@mean,'double'),boottest(lags{2},lags{3},5e3,@mean,'double')];
pval_median=[boottest(lags{1},lags{2},5e3,@median,'double'),boottest(lags{1},lags{3},5e3,@median,'double'),boottest(lags{2},lags{3},5e3,@median,'double')];
eval(['pval_mean_',varname,'=pval_mean']);
eval(['pval_median_',varname,'=pval_median']);

%inset
F=figure;
A=axes;
for i=1:numel(strcts)
    B=bar([i],mean_l(i));
    B.FaceColor=COL(i,:);
    hold on;
end
set(A,'FontSize',10);
set(F,'Position',[680 558 150 150]);
set(A,'XTick',[],'Xlim',[0 4]);
set(A,'YLim',[0 max(mean_l)*1.5]);

end
function [pvalue]=xcov_bootstrap(strct,varname)
N=1e2;
mode='median';
l_field=['all_lag_',varname];
for i=1:3
    l{i}=getfield(strct{i},l_field);
end
pvalue=nan(2,3);
for i=1:2
    for j=(i+1):3
        pvalue(i,j)=boottest(l{i},l{j},N,mode);
    end
end
end

function plotDurCorr(strcts,varname)
COL=colormap('lines');
F=figure;
A=axes;
a=0.01;

l_field=['all_lag_',varname];
for i=1:numel(strcts)
    l=getfield(strcts{i},l_field);
    l=l+a*randn(size(l));
    s=strcts{i}.segment_dur;
    H=plot(s,l,'.');
    set(H,'MarkerSize',18);
    hold on;
    if(i==1)
        f=fit(s(:),l(:),'poly1');
        H=plot(f);
        set(H,'LineWidth',2,'Color',0.5*[1 1 1],'LineStyle',':');
        legend('off');
        [c,p]=corr(s(:),l(:))
    end
end
set(A,'XLabel',[],'YLabel',[],'FontSize',14);
end

function A=plotCDFComp(data1,data2,data3)
COL=colormap('lines');
median1=median(data1);
median2=median(data2);
median3=median(data3);
mean1=mean(data1);
mean2=mean(data2);
mean3=mean(data3);
var1=var(data1);
var2=var(data2);
var3=var(data3);
pval_median=[boottest(data1,data2,5e3,@median,'double'),boottest(data1,data3,5e3,@median,'double'),boottest(data2,data3,5e3,@median,'double')]
pval_mean=[boottest(data1,data2,5e3,@mean,'double'),boottest(data1,data3,5e3,@mean,'double'),boottest(data2,data3,5e3,@mean,'double')]
pval_var=[boottest(data1,data2,5e3,@var,'double'),boottest(data1,data3,5e3,@var,'double'),boottest(data2,data3,5e3,@var,'double')]
F=figure;
A(1)=axes;
H=cdfplot(data1);
set(H,'Color',COL(1,:),'LineWidth',2);
hold on;
H=cdfplot(data2);
set(H,'Color',COL(2,:),'LineWidth',2);
H=cdfplot(data3);
set(H,'Color',COL(3,:),'LineWidth',2);
set(A(1),'FontSize',14);
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
B=bar([3],mean3);
B.FaceColor=COL(3,:);
set(A(2),'FontSize',10);
set(F,'Position',[680 558 150 150]);
set(A(2),'XTick',[],'Xlim',[0 4]);
set(A(2),'YLim',[0 max([mean1,mean2,mean3])*1.5]);

% F=figure;
% A(3)=axes;
% B=bar([1],var1);
% B.FaceColor=COL(1,:);
% hold on;
% B=bar([2],var2);
% B.FaceColor=COL(2,:);
% B=bar([3],var3);
% B.FaceColor=COL(3,:);
% set(A(3),'FontSize',10);
% set(F,'Position',[680 558 150 150]);
% set(A(3),'XTick',[],'Xlim',[0 4]);
% set(A(3),'YLim',[0 max([var1,var2,var3])*1.5]);
% 
% F=figure;
% A(4)=axes;
% B=bar([1],median1);
% B.FaceColor=COL(1,:);
% hold on;
% B=bar([2],median2);
% B.FaceColor=COL(2,:);
% B=bar([3],median3);
% B.FaceColor=COL(3,:);
% set(A(4),'FontSize',10);
% set(F,'Position',[680 558 150 150]);
% set(A(4),'XTick',[],'Xlim',[0 4]);
% set(A(4),'YLim',[25 max([median1,median2,median3])*1.25]);
end
