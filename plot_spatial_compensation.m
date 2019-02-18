function plot_spatial_compensation(hf_struct,comp_struct,ctl_struct,all_sigs)
% figures for spatial compensation
global m M midoff minamp maxamp k fsize;
% m=round(min(hf_struct.amp_off(:,2)));
% M=round(max(hf_struct.amp_off(:,2)));
m=round(min(all_sigs{1}));
M=round(max(all_sigs{1}));
midoff=(M+m)/2;
maxamp=(M-m)/2;
minamp=2;
fsize=22;
F=figure;
COL=colormap('lines');
close(F);
head=2;
%%  envelope distributions
env_dist=1;
if(env_dist)
% plotAmpOff(hf_struct,0);
plotAmpOff(comp_struct,0);
plotAmpOff(ctl_struct,0);
plotAmpOff(comp_struct,head);
plotAmpOff(ctl_struct,head);
end
%%  entropy comparison
entropy_comp=1;
if(entropy_comp)
N=75;
P=0.3;
% [h_vals,H_vals,pvals]=compare_entropy_all({hf_struct,comp_struct,ctl_struct},N,P,0);
% pvals
% plot_metrics(H_vals,h_vals);
% 
[h_vals,H_vals,pvals]=compare_entropy_all({comp_struct,ctl_struct},N,P,0);
pvals
plot_metrics(H_vals,h_vals);

[h_vals,H_vals,pvals]=compare_entropy_all({comp_struct,ctl_struct},N,P,head);
pvals
plot_metrics(H_vals,h_vals);

% [h_vals,H_vals,pvals]=compare_entropy_head({hf_struct,ctl_struct},N,P,head);
% pvals
% plot_metrics(H_vals,h_vals);

% [d_vals,D_vals,pvals]=compare_deviation_all({hf_struct,comp_struct,ctl_struct},N);
% pvals
% plot_metrics(D_vals,d_vals);
end
%% trimming effect: offset, amplitude
trim_effect=0;
if(trim_effect)
A=plotCDFComp(ctl_struct.all_diluted_off,comp_struct.all_diluted_off,COL(2:3,:));
set(A(1),'XLim',[40 120]);
A=plotCDFComp(ctl_struct.all_diluted_amp,comp_struct.all_diluted_amp,COL(2:3,:));
set(A(1),'XLim',[2 30]);
% A=plotCDFComp(all_sigs{2},all_sigs{3});
A=plotCDFComp(hf_struct.all_diluted_off,[ctl_struct.all_diluted_off],[COL(1,:);COL(2,:)]);
set(A(1),'XLim',[40 120]);

% %amp_off=comp_struct.all_diluted_amp_off;
% amp_off=ctl_struct.all_diluted_amp_off;
% bottoms=amp_off(:,2)-amp_off(:,1);
% tops=amp_off(:,2)+amp_off(:,1);
% F=figure;
% A=axes;
% h_bot=cdfplot(bottoms);
% x_bot=h_bot.XData;
% y_bot=h_bot.YData;
% hold on;
% h_top=cdfplot(tops);
% x_top=h_top.XData;
% y_top=h_top.YData;
% hold off;
% [X,I]=sort([x_bot(1:end-1),x_top(1:end-1)]);
% y=[diff(y_bot),-diff(y_top)];
% Y=y(I);
% Z=cumsum(Y);
% plot(X,Z);
end
end

function plotAmpOff(strct,head)
global m M midoff maxamp fsize;
minamp=2;
all_amp=strct.amp_off(1:end,1);
if(~head)
    all_offset=strct.amp_off(1:end,2);
elseif(head==1)
    all_offset=strct.amp_hoff(1:end,2);
else
    all_offset=strct.amp_hoff_ctr(1:end,2);
end
x=linspace(minamp,maxamp,25);
y=linspace(m,M,26);
d=mean(diff(x))*mean(diff(y));
ind=find(all_amp>2 & all_offset<112);
Nx=hist3([all_amp(ind) all_offset(ind)],{x,y});
Nx=Nx/sum(Nx(:))/d;
X=linspace(minamp,maxamp,1e3);
Y=linspace(m,M,1e3+1);
Z=interp2(x,y,Nx',X,Y','linear')*100;
[XX,YY]=meshgrid(X,Y);
a=(M-midoff)/(minamp-maxamp);
b=M-a*minamp;
Z(YY>a*XX+b)=NaN;
a=(m-midoff)/(minamp-maxamp);
b=m-a*minamp;
Z(YY<a*XX+b)=NaN;

F=figure;
A=axes;
surf(XX,YY,zeros(size(Z)),'CData',Z,'FaceColor','Interp','LineStyle','none');
% surf(x,y,zeros(size(Nx')),'CData',Nx','FaceColor','Interp','LineStyle','none');
view(0,90);
% set(A,'FontSize',16);
set(A,'Clim',[0 0.3]);
% L=line([minamp maxamp],[midoff midoff]);
% L=line([minamp minamp],[m M]);
set(A,'Xlim',[minamp,maxamp]);
set(A,'Ylim',[m,M]);
colormap('hot');
set(A,'XGrid','off','YGrid','off');
% set(A,'XTick',[],'Ytick',[],'XColor','none','YColor','none');
set(A,'FontSize',fsize);
set(F,'Position',[680 558 420 420]);
end

function [h_free,H,p]=compare_entropy(strct_hf,strct_free,N)
global m M midoff maxamp fsize;
minamp=2;
%maxoff=120;
ctrs{1}=linspace(minamp,maxamp,25);
ctrs{2}=linspace(m,M,25);
hf_amp_off=strct_hf.amp_off;
free_amp_off=strct_free.amp_off;

S=numel(strct_hf.segment_length);
inds=[1 cumsum(strct_hf.segment_length)];
for n=1:N   
    I=randperm(S);
    a=[];
    i=1;
    while(size(a,1)<size(free_amp_off,1))
        a=[a;hf_amp_off(inds(I(i)):inds(I(i)+1),:)];
        i=i+1;
    end
    ind=find(a(:,1)>minamp & a(:,2)<M);
    Nx=hist3([a(ind,1) a(ind,2)],ctrs);
    Nx=Nx/sum(Nx(:));    
    H(n)=ent(Nx(:));
    if(H(n)==max(H))
        Nx_top=Nx;
    end
end

ind_free=find(free_amp_off(:,1)>minamp & free_amp_off(:,2)<M);
N_free=hist3([free_amp_off(ind_free,1) free_amp_off(ind_free,2)],ctrs);
N_free=N_free/sum(N_free(:));
h_free=ent(N_free(:));

p=sum(H<=h_free)/N;
end

function [h_vals,H_vals,pvals]=compare_entropy_all(strcts,N,P,head)
global m M midoff minamp maxamp fsize;
%maxoff=120;
ctrs{1}=linspace(minamp,maxamp,25);
ctrs{2}=linspace(m,M,25);
for i=1:numel(strcts)
    if(~head)
        amp_off{i}=strcts{i}.amp_off;
    elseif (head==1 | i==1)
        amp_off{i}=strcts{i}.amp_hoff;
    else
        amp_off{i}=strcts{i}.amp_hoff_ctr;
    end
    S(i)=numel(strcts{i}.segment_length);
    inds{i}=[1 cumsum(strcts{i}.segment_length)];
    L(i)=numel(amp_off{i});
end

%get enthropies of entire datasets
for i=1:numel(strcts)
    ind=find(amp_off{i}(:,1)>minamp & amp_off{i}(:,2)<M);
    Nx=hist3([amp_off{i}(ind,1) amp_off{i}(ind,2)],ctrs);
    Nx=Nx/sum(Nx(:));    
    H_vals(i)=ent(Nx(:));    
end

L_boot=min(L*P);
for i=1:numel(strcts)
    for n=1:N   
%         I=randperm(S(i));
        I=randi(S(i),1e3,1);
        a=[];
        k=1;
        while(size(a,1)<L_boot)
            a=[a;amp_off{i}(inds{i}(I(k)):inds{i}(I(k)+1),:)];
            k=k+1;
        end
        ind=find(a(:,1)>minamp & a(:,2)<M);
        Nx=hist3([a(ind,1) a(ind,2)],ctrs);
        Nx=Nx/sum(Nx(:));    
        h_vals{i}(n)=ent(Nx(:));
    end
end

for i=1:(numel(strcts)-1)
    for j=(i+1):numel(strcts)
        pvals(i,j-1)=boottest(h_vals{i},h_vals{j},1e4,@mean,'double');
    end
end
end

function [h_vals,H_vals,pvals]=compare_entropy_head(strcts,N,P,head)
global m M midoff minamp maxamp fsize;
%maxoff=120;
ctrs{1}=linspace(minamp,maxamp,25);
ctrs{2}=linspace(m,M,25);
j=1;
for i=1:numel(strcts)
    if(i==1) %headfixed
        j=2; %middle, 3 for all conditions
        amp_off{j}=strcts{i}.amp_off;
        S(j)=numel(strcts{i}.segment_length);
        inds{j}=[1 cumsum(strcts{i}.segment_length)];
        L(j)=numel(amp_off{j});
        j=1;
    else
        %whisker only
        amp_off{j}=strcts{i}.amp_off;
        S(j)=numel(strcts{i}.segment_length);
        inds{j}=[1 cumsum(strcts{i}.segment_length)];
        L(j)=numel(amp_off{j});
        
        % whisker + head
        j=j+2; %+3 for all conditions
        if(head==1)
            amp_off{j}=strcts{i}.amp_hoff;
        else
            amp_off{j}=strcts{i}.amp_hoff_ctr;
        end
        S(j)=numel(strcts{i}.segment_length);
        inds{j}=[1 cumsum(strcts{i}.segment_length)];
        L(j)=numel(amp_off{j});
        j=j-2;        
    end
end
J=numel(amp_off);

%get enthropies of entire datasets
for i=1:J
    ind=find(amp_off{i}(:,1)>minamp & amp_off{i}(:,2)<M);
    Nx=hist3([amp_off{i}(ind,1) amp_off{i}(ind,2)],ctrs);
    Nx=Nx/sum(Nx(:));    
    H_vals(i)=ent(Nx(:));    
end

L_boot=min(L*P);
for i=1:J
    for n=1:N   
%         I=randperm(S(i));
        I=randi(S(i),1e3,1);
        a=[];
        k=1;
        while(size(a,1)<L_boot)
            a=[a;amp_off{i}(inds{i}(I(k)):inds{i}(I(k)+1),:)];
            k=k+1;
        end
        ind=find(a(:,1)>minamp & a(:,2)<M);
        Nx=hist3([a(ind,1) a(ind,2)],ctrs);
        Nx=Nx/sum(Nx(:));    
        h_vals{i}(n)=ent(Nx(:));
    end
end

for i=1:(J-1)
    for j=(i+1):J
        pvals(i,j-1)=boottest(h_vals{i},h_vals{j},1e4,@mean,'double');
    end
end
end

function H=ent(p)
    H=nansum(-log2(p).*p);
end

function [d_vals,D_vals,pvals]=compare_deviation_all(strcts,N)
global m M midoff minamp maxamp;
minamp=2;
for i=1:numel(strcts)
    amp_off{i}=strcts{i}.amp_off;
    S(i)=numel(strcts{i}.segment_length);
    inds{i}=[1 cumsum(strcts{i}.segment_length)];
    L(i)=numel(amp_off{i});
end

%get enthropies of entire datasets
for i=1:numel(strcts)
    ind=find(amp_off{i}(:,1)>minamp);
    D_vals(i)=sqrt(det(cov(amp_off{i}(ind,:))));    
end

L_boot=min(L*0.3);
for i=1:numel(strcts)
    for n=1:N   
%         I=randperm(S(i));
        I=randi(S(i),1e3,1);
        a=[];
        k=1;
        while(size(a,1)<L_boot)
            a=[a;amp_off{i}(inds{i}(I(k)):inds{i}(I(k)+1),:)];
            k=k+1;
        end
        ind=find(a(:,1)>minamp);
        d_vals{i}(n)=sqrt(det(cov(a(ind,:))));
    end
end

for i=1:(numel(strcts)-1)
    for j=(i+1):numel(strcts)
        pvals(i,j-1)=boottest(d_vals{i},d_vals{j},5e3,@mean,'double');
    end
end
end

function plot_metrics(metric,mk)
global fsize
F=figure;
AA=axes;
X=[];
G=[];
COL=colormap('lines');
% COL1=COL([2 3 1 2 3],:);
COL1=COL([2 1 2],:);
for i=1:numel(mk)
    X=[X;mk{i}'];
    G=[G;ones(numel(mk{i}),1)*i];
end
% positions=[1 1.1];
positions=1:numel(metric);
boxplot(X,G,'plotstyle','compact','symbol','','colorgroup',G,'median','line','positions',positions,...
    'colors',1-(1-COL1)/2,'Widths',.5,'FactorGap',[],'FactorSeparator',[],...
    'orientation','vertical','LabelOrientation','inline');
A=findobj('Tag','Box');
for i=1:numel(A)
    if(A(i).Parent.Parent==AA)
        A(i).LineWidth=40;
    end
end
A=findobj('Tag','Whisker');
for i=1:numel(A)
    if(A(i).Parent.Parent==AA)    
        A(i).LineWidth=2;
    end
end
A=findobj('Tag','Median');
for i=1:numel(A)
    if(A(i).Parent.Parent==AA)
        A(i).Color=1-(1-A(i).Color)*2;
    end
end
hold on;
for i=1:numel(mk)
    H=plot(positions(i),metric(i),'.');
    hold on;
    set(H,'MarkerSize',25,'Color',COL1(i,:));
    m(i)=min(mk{i});
    M(i)=max(mk{i});
end
axis([0.5 numel(metric)+.5 min(m)*.9 max(M)*1.1])
box('off');
set(AA,'FontSize',fsize);
set(F,'Position',[680 558 420 420]);

end

function [r_free,R,p]=compare_rms(strct_hf,strct_free,N)
global m M midoff maxamp;
minamp=2;
hf_amp_off=strct_hf.amp_off;
free_amp_off=strct_free.amp_off;

S=numel(strct_hf.segment_length);
inds=[1 cumsum(strct_hf.segment_length)];
for n=1:N   
    I=randperm(S);
    a=[];
    i=1;
    while(size(a,1)<size(free_amp_off,1))
        a=[a;hf_amp_off(inds(I(i)):inds(I(i)+1),:)];
        i=i+1;
    end
    ind=find(a(:,1)>minamp & a(:,2)<M);
    a=a(ind,:)-ones(numel(ind),1)*mean(a(ind,:));
    R(n)=sqrt(mean((a.^2*[1;1])));
end

ind_free=find(free_amp_off(:,1)>minamp & free_amp_off(:,2)<M);
a=free_amp_off(ind_free,:);
a=a-ones(size(a,1),1)*mean(a);
R_free=a.^2*[1;1];
r_free=sqrt(mean(a.^2*[1;1]));

ind_hf=find(hf_amp_off(:,1)>minamp & hf_amp_off(:,2)<M);
a=hf_amp_off(ind_hf,:);
a=a-ones(size(a,1),1)*mean(a);
R_hf=a.^2*[1;1];
r_hf=sqrt(mean(a.^2*[1;1]));


p=sum(R<=r_free)/N;
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
