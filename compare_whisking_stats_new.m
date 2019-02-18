function compare_whisking_stats_new(strct_comp,strct_ctl,strct_hf)
global m M midoff maxamp k;
%% load data
L={'free, 1 row','free, 3 whiskers','headfixed, 3 whiskers'};
m=min(strct_hf.amp_off(:,2));
M=max(strct_hf.amp_off(:,2));
midoff=(M+m)/2;
maxamp=(M-m)/2;
k=1;

COL=colormap('lines');
%% amplitude distributions
plotAmp(strct_comp,COL(1,:));
plotAmp(strct_ctl,COL(2,:));
plotAmp(strct_hf,COL(3,:));
legend(L);

N=1e3;
mode='mean';
p_amp_ctl_comp=boottest(strct_comp.all_diluted_amp,strct_ctl.all_diluted_amp,N,mode)
p_amp_ctl_hf=boottest(strct_hf.all_diluted_amp,strct_ctl.all_diluted_amp,N,mode)
p_amp_comp_hf=boottest(strct_comp.all_diluted_amp,strct_hf.all_diluted_amp,N,mode)

%% offset distributions
plotOff(strct_comp,COL(1,:));
plotOff(strct_ctl,COL(2,:));
plotOff(strct_hf,COL(3,:));
legend(L);

N=1e3;
mode='mean';
p_off_ctl_comp=boottest(strct_comp.all_diluted_off,strct_ctl.all_diluted_off,N,mode)
p_off_ctl_hf=boottest(strct_hf.all_diluted_off,strct_ctl.all_diluted_off,N,mode)
p_off_comp_hf=boottest(strct_comp.all_diluted_off,strct_hf.all_diluted_off,N,mode)

%% durations, ratios, protraction, retraction distributions
strct_comp=plotRatDur(strct_comp,COL(1,:));
strct_ctl=plotRatDur(strct_ctl,COL(2,:));
strct_hf=plotRatDur(strct_hf,COL(3,:));


%% amplitude-offset coverage

plotLowOff(strct_comp,1);
plotLowOff(strct_ctl,2);
plotLowOff(strct_hf,3);

%randomize hf collections  
N=5e3;
p_H_comp=compare_entropy(strct_hf,strct_comp,N)
p_H_ctl=compare_entropy(strct_hf,strct_ctl,N)

%% pumps
[v_comp_prot,v_comp_ret]=plotPump(strct_comp,COL(1,:));
[v_ctl_prot,v_ctl_ret]=plotPump(strct_ctl,COL(2,:));
[v_hf_prot,v_hf_ret]=plotPump(strct_hf,COL(3,:));
legend(L);

figure(50);
subplot(1,2,1);
bar([v_comp_prot;v_ctl_prot;v_hf_prot],'stacked');
subplot(1,2,2);
bar([v_comp_ret;v_ctl_ret;v_hf_ret],'stacked');

[pval_pumps]=pump_bootstrap({strct_comp,strct_ctl,strct_hf})

%% xcor
varnames={'amp','off','dur','pr'};
for i=1:numel(varnames)
    F1=figure;
    F2=figure;
    varname=varnames{i};
    plotXcor(strct_comp,COL(1,:),1,varname,F1,F2);
    plotXcor(strct_ctl,COL(2,:),2,varname,F1,F2);
    plotXcor(strct_hf,COL(3,:),3,varname,F1,F2);
    legend(L);
    pval_xcov(:,:,i)=xcov_bootstrap({strct_comp,strct_ctl,strct_hf},varname)
end


end

function plotAmp(strct,c)
ctrs{1}=linspace(0,50,24);
all_amp=strct.all_diluted_amp;
Namp=hist(all_amp',ctrs{1});
Namp=Namp/sum(Namp);
figure(10);
H=plot(ctrs{1},Namp);
set(H,'Color',c,'LineWidth',2);
hold on;
H=xlabel('Amplitude (deg)');
set(H,'FontSize',12);
end

function plotOff(strct,c)
ctrs{1}=linspace(0,120,24);
all_off=strct.all_diluted_off;
Noff=hist(all_off',ctrs{1});
Noff=Noff/sum(Noff);
figure(20);
H=plot(ctrs{1},Noff);
set(H,'Color',c,'LineWidth',2);
hold on;
H=xlabel('Offset (deg)');
set(H,'FontSize',12);
end

function plotLowOff(strct,n)
global m M midoff maxamp;
minamp=2;
all_amp=strct.amp_off(1:end,1);
all_offset=strct.amp_off(1:end,2);
ctrs{1}=linspace(minamp,45,25);
ctrs{2}=linspace(30,112,25);
ind=find(all_amp>2 & all_offset<112);
Nx=hist3([all_amp(ind) all_offset(ind)],ctrs);
Nx=Nx/sum(Nx(:));
figure(11)
H=subplot(1,3,n);
image(ctrs{1},ctrs{2},Nx'*5e3);
set(H,'YDir','normal');
hold on;
H=ylabel('Offset (deg)');
set(H,'FontSize',12);
H=xlabel('Amplitude (deg)');
set(H,'FontSize',12);
line([2 maxamp],[m midoff]);
line([2 maxamp],[M midoff]);
[X Y]=meshgrid(ctrs{1},ctrs{2});
grid_area=sum(withinT(X(:),Y(:),2,maxamp,m,M));
support=sum(Nx(:)>0)/grid_area;
title(['support=',num2str(support)]);
end

function v=withinT(X,Y,x1,x2,y1,y2)
%checks whether points in p matrx are in triangle
X=X-x1;
Y=Y-y1;
y2=y2-y1;
x2=x2-x1;
a=y2/2/x2;
v=(Y>=a*X & Y<=y2-a*X);
end

function plotXcor(strct,c,n,varname,F1,F2)
x_field=['all_X_',varname];
l_field=['all_lag_',varname];
x=getfield(strct,x_field);
l=getfield(strct,l_field);
[lags,I]=sort(l);
figure(F1)
subplot(1,3,n);
imagesc(x(I,:));

figure(F2)
H=cdfplot(lags);
set(H,'Color',c,'LineWidth',2);
hold on;
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


function strct=plotRatDur(strct,c)
% dur=strct.dur_rat(:,1);
% pr_ratio=strct.dur_rat(:,2);
dur=strct.all_diluted_dur;
pr_ratio=strct.all_diluted_pr;
prot=strct.all_diluted_prot;
ret=strct.all_diluted_ret;
ctrs=linspace(3,18,50);
Nd=hist(1./dur,ctrs);
Nd=Nd/sum(Nd);
figure(30);
subplot(2,2,1);
H=plot(ctrs,Nd);
set(H,'Color',c,'LineWidth',2);
xlabel('Frequency (Hz)');
hold on;
ctrs=linspace(0,5,50);
Nr=hist(pr_ratio,ctrs);
Nr=Nr/sum(Nr);
subplot(2,2,2);
H=plot(ctrs,Nr);
set(H,'Color',c,'LineWidth',2);
hold on;
xlabel('Protraction-Retraction Ratio');
ctrs=linspace(0,200,50);
Nr=hist(prot*1e3,ctrs);
Nr=Nr/sum(Nr);
subplot(2,2,3);
H=plot(ctrs,Nr);
set(H,'Color',c,'LineWidth',2);
hold on;
xlabel('Protraction Duration (ms)');
ctrs=linspace(0,200,50);
Nr=hist(ret*1e3,ctrs);
Nr=Nr/sum(Nr);
subplot(2,2,4);
H=plot(ctrs,Nr);
set(H,'Color',c,'LineWidth',2);
hold on;
xlabel('Retraction Duration (ms)');
strct.prot=prot;
strct.ret=ret;
end

function H=ent(p)
    H=nansum(-log2(p).*p);
end

function [v_prot,v_ret]=plotPump(strct,c)
prot_pump_str=strct.all_pump_stregth(:,1);
ret_pump_str=strct.all_pump_stregth(:,2);

figure(40)
A=subplot(1,2,1)
H=cdfplot(prot_pump_str);
set(H,'Color',c,'LineWidth',2);
hold on;
set(A,'Xlim',[0 2.5]);
A=subplot(1,2,2)
H=cdfplot(ret_pump_str);
set(H,'Color',c,'LineWidth',2);
hold on;
set(A,'Xlim',[0 2.5]);

v_prot=get_pumptypes(prot_pump_str);
v_ret=get_pumptypes(ret_pump_str);
end

function [pvalue]=pump_bootstrap(strct)
N=1e4;
for i=1:3
    pump_str{i,1}=strct{i}.all_pump_stregth(:,1);    
    pump_str{i,2}=strct{i}.all_pump_stregth(:,2);
end

pvalue=nan(2,3,2);
for k=1:2
    for i=1:2
        for j=(i+1):3
            x=pump_str{i,k};
            y=pump_str{j,k};
            Nx=numel(x);
            Ny=numel(y);
            Px=sum(x>0)/Nx;
            Py=sum(y>0)/Ny;
            D=Px-Py;
            V=[x(:);y(:)];
            for l=1:N               
%                 VV=V(randperm(numel(V))); %no replacments
                VV=V(randi(numel(V),1,numel(V))); %with replacments

                xx=VV(1:Nx);
                yy=VV((Nx+1):end);
                p_xx=sum(xx>0)/Nx;
                p_yy=sum(yy>0)/Ny;
                dd(l)=p_xx-p_yy;
            end
            pvalue(i,j,k)=sum(abs(dd)>=abs(D))/N;
%             if(D>0)
%                 pvalue(i,j,k)=sum(dd>=D)/N;
%             else
%                 pvalue(i,j,k)=sum(dd<=D)/N;
%             end


        end
    end
end
end

function v=get_pumptypes(pump_vec)
%v=[sum(pump_vec==0) sum(pump_vec>0 & pump_vec<=1) sum(pump_vec>1)]/numel(pump_vec);
%v=[sum(pump_vec==0) sum(pump_vec>0)]/numel(pump_vec);
v=[sum(pump_vec>0)]/numel(pump_vec);
end

function p=compare_entropy(strct_hf,strct_free,N)
minamp=2;
maxoff=120;
ctrs{1}=linspace(minamp,45,25);
ctrs{2}=linspace(30,maxoff,25);
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
    ind=find(a(:,1)>minamp & a(:,2)<maxoff);
    Nx=hist3([a(ind,1) a(ind,2)],ctrs);
    Nx=Nx/sum(Nx(:));    
    H(n)=ent(Nx(:));
    if(H(n)==max(H))
        Nx_top=Nx;
    end
end

ind_free=find(free_amp_off(:,1)>minamp & free_amp_off(:,2)<maxoff);
N_free=hist3([free_amp_off(ind_free,1) free_amp_off(ind_free,2)],ctrs);
N_free=N_free/sum(N_free(:));
h_free=ent(N_free(:));

p=sum(H<=h_free)/N;
end