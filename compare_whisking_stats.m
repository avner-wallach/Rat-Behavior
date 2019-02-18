function compare_whisking_stats(strct_comp,strct_ctl,strct_hf)
global m M midoff maxamp k;
%% load data
% load('comparative_c2_data');
% wdata_comp=wdata;
% 
% load('control_c2_data');
% wdata_ctl=wdata;
% 
% load('headfixed_c2_data');
% wdata_hf=wdata;
% 
% %% analyze each
% strct_comp=get_whisking_stats(wdata_comp,0);
% strct_ctl=get_whisking_stats(wdata_ctl,0);
% strct_hf=get_whisking_stats(wdata_hf,0);
L={'free, 1 row','free, 3 whiskers','headfixed, 3 whiskers'};
m=min(strct_hf.amp_off(:,2));
M=max(strct_hf.amp_off(:,2));
midoff=(M+m)/2;
maxamp=(M-m)/2;
k=1;
%% vizualize
plotAmp(strct_comp,'b');
plotAmp(strct_ctl,'m');
plotAmp(strct_hf,'r');
legend(L);

plotLowOff(strct_comp,1);
plotLowOff(strct_ctl,2);
plotLowOff(strct_hf,3);

plotXcor(strct_comp,'b');
plotXcor(strct_ctl,'m');
plotXcor(strct_hf,'r');
% 
% plotPump(strct_comp,'b');
% plotPump(strct_ctl,'m');
% plotPump(strct_hf,'r');

strct_comp=plotRatDur(strct_comp,'b');
strct_ctl=plotRatDur(strct_ctl,'m');
strct_hf=plotRatDur(strct_hf,'r');
figure(31);
p_prot=my_boxplot({strct_comp.prot',strct_ctl.prot',strct_hf.prot'},{'Free,1 row','Free, 3 whiskers','Headfixed'})
p_ret=my_boxplot({strct_comp.ret',strct_ctl.ret',strct_hf.ret'},{'Free,1 row','Free, 3 whiskers','Headfixed'})
%% randomize hf collections  
N=10;
maxlag=500;

ctrs{1}=linspace(0,45,25);
ctrs{2}=linspace(30,112,25);
amp_off=strct_hf.amp_off;
comp_amp_off=strct_comp.amp_off;
ctl_amp_off=strct_ctl.amp_off;
inds=find(isnan(amp_off(1:end-1,1)) & ~isnan(amp_off(2:end,1)))+1;
inde=find(isnan(amp_off(2:end,1)) & ~isnan(amp_off(1:end-1,1)));
S=numel(inde);
dr=size(strct_ctl.dur_rat,1);
for n=1:N   
    n
    I=randperm(S);
    a=[];
    i=1;
    while(numel(a)<numel(find(~isnan(ctl_amp_off))))
        a=[a;amp_off(inds(I(i)):inde(I(i)),:)];
        [x_a(i,:),lags]=xcov(amp_off(inds(I(i)):inde(I(i)),1),maxlag,'coeff');
        [x_o(i,:),lags]=xcov(amp_off(inds(I(i)):inde(I(i)),2),maxlag,'coeff');        
        i=i+1;
    end
    ind=find(a(:,1)>2 & a(:,2)<112);
    Nx=hist3([a(ind,1) a(ind,2)],ctrs);
    Nx=Nx/sum(Nx(:));    
    H(n)=ent(Nx(:));
    if(H(n)==max(H))
        Nx_top=Nx;
    end
    
    X_a(n,:)=nanmean(x_a,1);        
    X_a(n,:)=X_a(n,:)/max(X_a(n,:));
    X_o(n,:)=nanmean(x_o,1);        
    X_o(n,:)=X_o(n,:)/max(X_o(n,:));
end

ind_comp=find(comp_amp_off(:,1)>2 & comp_amp_off(:,2)<112);
N_comp=hist3([comp_amp_off(ind_comp,1) comp_amp_off(ind_comp,2)],ctrs);
N_comp=N_comp/sum(N_comp(:));
h_comp=ent(N_comp(:));
ind_ctl=find(ctl_amp_off(:,1)>2 & ctl_amp_off(:,2)<112);
N_ctl=hist3([ctl_amp_off(ind_ctl,1) ctl_amp_off(ind_ctl,2)],ctrs);
N_ctl=N_ctl/sum(N_ctl(:));
h_ctl=ent(N_ctl(:));

figure(42)
subplot(1,2,1)
plot(ctrs{1},sum(N_comp,2),'b');
hold on;
plot(ctrs{1},sum(N_ctl,2),'m');
% p_amp=boottest(comp_amp_off(ind_comp,1),ctl_amp_off(ind_ctl,1),1e4)
subplot(1,2,2)
plot(ctrs{2},sum(N_comp,1),'b');
hold on;
plot(ctrs{2},sum(N_ctl,1),'m');
% p_off=boottest(comp_amp_off(ind_comp,2),ctl_amp_off(ind_ctl,2),1e4,'std')


figure(40), hist(H,100);
hold on;
plot(h_comp,0,'v');
plot(h_ctl,0,'vm');

figure(41)
h1=subplot(1,3,1);
image(ctrs{1},ctrs{2},N_comp'*5e3);
set(h1,'YDir','normal');
h1=ylabel('Offset (deg)');
set(h1,'FontSize',12);
h1=xlabel('Amplitude (deg)');
set(h1,'FontSize',12);
h1=subplot(1,3,2);
image(ctrs{1},ctrs{2},N_ctl'*5e3);
set(h1,'YDir','normal');
h1=ylabel('Offset (deg)');
set(h1,'FontSize',12);
h1=xlabel('Amplitude (deg)');
set(h1,'FontSize',12);
h1=subplot(1,3,3);
image(ctrs{1},ctrs{2},Nx_top'*5e3);
set(h1,'YDir','normal');
h1=ylabel('Offset (deg)');
set(h1,'FontSize',12);
h1=xlabel('Amplitude (deg)');
set(h1,'FontSize',12);

figure(50);
% subplot(1,2,1);
% my_plotWithConf(lags/500,nanmean(X_a,1),nanstd(X_a,[],1),[1 0 0]);
% hold on;
% plot(strct_ctl.x_lags,strct_ctl.x_amp,'m');
% plot(strct_comp.x_lags,strct_comp.x_amp,'b');
hist(sum(abs(X_a),2),100);
hold on;
plot(sum(abs(strct_ctl.x_amp_m)),0,'vm');
plot(sum(abs(strct_comp.x_amp_m)),0,'vb');

% subplot(1,2,2);
% my_plotWithConf(lags/500,nanmean(X_o,1),nanstd(X_o,[],1),[1 0 0]);
% hold on;
% plot(strct_ctl.x_lags,strct_ctl.x_offset,'m');
% plot(strct_comp.x_lags,strct_comp.x_offset,'b');

%% randomize hf collections - durations and pr
N=100;
maxlag=20;
dur_rat=strct_hf.dur_rat;
comp_dur_rat=strct_comp.dur_rat;
ctl_dur_rat=strct_ctl.dur_rat;
inds=find(isnan(dur_rat(1:end-1,1)) & ~isnan(dur_rat(2:end,1)))+1;
inde=find(isnan(dur_rat(2:end,1)) & ~isnan(dur_rat(1:end-1,1)));
S=numel(inde);
for n=1:N    
    I=randperm(S);
    dp=[];
    i=1;
    while(numel(dp)<numel(find(~isnan(ctl_dur_rat))))
        dp=[dp;dur_rat(inds(I(i)):inde(I(i)),:)];
        [x_d(i,:),lags]=xcov(dur_rat(inds(I(i)):inde(I(i)),1),maxlag,'coeff');
        [x_p(i,:),lags]=xcov(dur_rat(inds(I(i)):inde(I(i)),2),maxlag,'coeff');        
        i=i+1;
    end
    
    X_d(n,:)=nanmean(x_d,1);        
    X_d(n,:)=X_d(n,:)/max(X_d(n,:));
    X_p(n,:)=nanmean(x_p,1);        
    X_p(n,:)=X_p(n,:)/max(X_p(n,:));    
end

clear x_d x_p;
inds=find(isnan(comp_dur_rat(1:end-1,1)) & ~isnan(comp_dur_rat(2:end,1)))+1;
inde=find(isnan(comp_dur_rat(2:end,1)) & ~isnan(comp_dur_rat(1:end-1,1)));
S=numel(inde);
for s=1:S
    [x_d(i,:),lags]=xcov(comp_dur_rat(inds(s):inde(s),1),maxlag,'none');
    [x_p(i,:),lags]=xcov(comp_dur_rat(inds(s):inde(s),2),maxlag,'none');        
end
X_comp_d=nanmean(x_d,1);        
X_comp_d=X_comp_d/max(X_comp_d);
X_comp_p=nanmean(x_p,1);        
X_comp_p=X_comp_p/max(X_comp_p);    

clear x_d x_p;
inds=find(isnan(ctl_dur_rat(1:end-1,1)) & ~isnan(ctl_dur_rat(2:end,1)))+1;
inde=find(isnan(ctl_dur_rat(2:end,1)) & ~isnan(ctl_dur_rat(1:end-1,1)));
S=numel(inde);
for s=1:S
    [x_d(i,:),lags]=xcov(ctl_dur_rat(inds(s):inde(s),1),maxlag,'none');
    [x_p(i,:),lags]=xcov(ctl_dur_rat(inds(s):inde(s),2),maxlag,'none');        
end
X_ctl_d=nanmean(x_d,1);        
X_ctl_d=X_ctl_d/max(X_ctl_d);
X_ctl_p=nanmean(x_p,1);        
X_ctl_p=X_ctl_p/max(X_ctl_p);    

figure(60);
subplot(1,2,1);
my_plotWithConf(lags,nanmean(X_d,1),nanstd(X_d,[],1),[1 0 0]);
hold on;
plot(lags,X_comp_d,'b');
plot(lags,X_ctl_d,'m');

subplot(1,2,2);
my_plotWithConf(lags,nanmean(X_p,1),nanstd(X_p,[],1),[1 0 0]);
hold on;
plot(lags,X_comp_p,'b');
plot(lags,X_ctl_p,'m');

end

function plotAmp(strct,c)
all_amp=strct.amp_off(1:end,1);
all_offset=strct.amp_off(1:end,2);
ctrs{1}=linspace(0,50,50);
ctrs{2}=linspace(20,120,50);
Nx=hist3([all_amp all_offset],ctrs);
Nx=Nx/sum(Nx(:));
Namp=hist(all_amp',ctrs{1});
Namp=Namp/sum(Namp);
Noff=hist(all_offset',ctrs{2});
Noff=Noff/sum(Noff);
figure(10);
% H=axes('position',[0.35 0.35 0.6 0.6]);
% imagesc(ctrs{1},ctrs{2},Nx');
% set(H,'YDir','normal');
% set(H,'YTick',[],'YTickLabel',[]);
% set(H,'XTick',[],'XTickLabel',[]);
% H=axes('position',[0.1 0.35 0.2 0.6]);
% barh(ctrs{2},Noff);
% set(H,'XDir','reverse');
% set(H,'YLim',[20,120]);
% H=ylabel('Offset (deg)');
% set(H,'FontSize',12);
% H=axes('position',[0.35 0.1 0.6 0.2]);
plot(ctrs{1},Namp,c);
% set(H,'XLim',[-1.5 50]);
hold on;
H=xlabel('Amplitude (deg)');
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

function plotXcor(strct,c)
global k;
lags=strct.x_lags;
x_amp=strct.x_amp_m;
x_dur=strct.x_duration;
x_off=strct.x_offset;
figure(20)
subplot(1,3,1);
plot(lags,x_amp,c);
hold on;
subplot(1,3,2);
plot(lags,x_off,c);
hold on;
subplot(1,3,3);
plot(lags,x_dur,c);
hold on;

% [y,I]=sort(strct.x_zc_amp);
% x=(1:size(strct.x_amp(1,500:end),2))*2;
% y=(1:size(strct.x_amp(1,500:end),1));
% figure(21)
% subplot(1,3,k);
% imagesc(x,y,abs(strct.x_amp(I,500:end)));
% 
% figure(22);
% H=cdfplot(strct.x_zc_amp);
% set(H,'Color',c);
% hold on;
% 
% [y,I]=sort(strct.x_zc_dur);
% x=(1:size(strct.x_duration(1,500:end),2))*2;
% y=(1:size(strct.x_duration(1,500:end),1));
% figure(23)
% subplot(1,3,k);
% imagesc(x,y,abs(strct.x_duration(I,500:end)));
% xlabel('Time lag (ms)');
% ylabel('Segment');
% 
% figure(24);
% H=cdfplot(strct.x_zc_dur);
% set(H,'Color',c);
% hold on;
% k=k+1;
end

function strct=plotRatDur(strct,c)
dur=strct.dur_rat(:,1);
pr_ratio=strct.dur_rat(:,2);
prot=dur.*pr_ratio./(1+pr_ratio);
ret=dur./(1+pr_ratio);
ctrs=linspace(3,18,50);
Nd=hist(1./dur,ctrs);
Nd=Nd/sum(Nd);
figure(30);
subplot(2,2,1);
plot(ctrs,Nd,c);
xlabel('Frequency (Hz)');
hold on;
ctrs=linspace(0,5,50);
Nr=hist(pr_ratio,ctrs);
Nr=Nr/sum(Nr);
subplot(2,2,2);
plot(ctrs,Nr,c);
hold on;
xlabel('Protraction-Retraction Ratio');
ctrs=linspace(0,200,50);
Nr=hist(prot*1e3,ctrs);
Nr=Nr/sum(Nr);
subplot(2,2,3);
plot(ctrs,Nr,c);
hold on;
xlabel('Protraction Duration (ms)');
ctrs=linspace(0,200,50);
Nr=hist(ret*1e3,ctrs);
Nr=Nr/sum(Nr);
subplot(2,2,4);
plot(ctrs,Nr,c);
hold on;
xlabel('Retraction Duration (ms)');
strct.prot=prot;
strct.ret=ret;
end

function H=ent(p)
    H=nansum(-log2(p).*p);
end

function plotPump(strct,c)
th=linspace(0,0.5,200);
pp=fliplr(strct.prot_pump);
pr=fliplr(strct.ret_pump);
figure(70);
subplot(1,2,1)
plot(th,pp,c);
xlabel('pump threshold')
ylabel('precentage of  protraction pumps')
hold on;
subplot(1,2,2)
plot(th,pr,c);
xlabel('pump threshold')
ylabel('precentage of  retraction pumps')
hold on;
end
