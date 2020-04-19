function plot_figure3_S2()
global COL scale fsize
fsize=22;
F=figure;
g=0;
COL=[zeros(1,3);g*ones(1,3)];

close(F);
framerate=500;

% correlation of pumps with head motion
all_gamma=cell(2,1);
pump_gamma=cell(2,1);
rpump_gamma=cell(2,1);
ppump_gamma=cell(2,1);

all_V=cell(2,1);
pump_V=cell(2,1);
rpump_V=cell(2,1);
ppump_V=cell(2,1);

all_Vl=cell(2,1);
pump_Vl=cell(2,1);
rpump_Vl=cell(2,1);
ppump_Vl=cell(2,1);

all_Vt=cell(2,1);
pump_Vt=cell(2,1);
rpump_Vt=cell(2,1);
ppump_Vt=cell(2,1);

all_da=cell(2,1);
pump_da=cell(2,1);
ppump_da=cell(2,1);
rpump_da=cell(2,1);

pump_str=cell(2,1);
ppump_str=cell(2,1);
rpump_str=cell(2,1);
%% load data
hD=cell(2,1);
pD=cell(2,1);
heads=cell(2,1);

side='r_'; %right
load(['../Data/comparative_c2h_',side,'data.mat']);
comp_hdata=hdata;
comp_pdata=pdata;
comp_headsize=headsize;
load(['../Data/control_c2h_',side,'data.mat']);
ctl_hdata=hdata;
ctl_pdata=pdata;
ctl_headsize=headsize;
hD{1}={comp_hdata{:} ctl_hdata{:}};
pD{1}={comp_pdata{:} ctl_pdata{:}};
heads{1}=[comp_headsize ctl_headsize];

side=''; %left
load(['../Data/comparative_c2h_',side,'data.mat']);
comp_hdata=hdata;
comp_pdata=pdata;
comp_headsize=headsize;
load(['../Data/control_c2h_',side,'data.mat']);
ctl_hdata=hdata;
ctl_pdata=pdata;
ctl_headsize=headsize;
hD{2}={comp_hdata{:} ctl_hdata{:}};
pD{2}={comp_pdata{:} ctl_pdata{:}};
heads{2}=[comp_headsize ctl_headsize];

%% get stats
M=11;
p=linspace(0,1,M);
c=mean([p(1:end-1);p(2:end)]);
scale=100; %100- percent, 1- probability

for s=1:2
    for i=1:numel(hD{s})
        i
        hdat=hD{s}{i};
        pdat=pD{s}{i};
        beta=hdat(:,1);
        X=hdat(:,2);
        Y=hdat(:,3);
        Vx=diff(X)*framerate;
        Vy=diff(Y)*framerate;
        Vx=[Vx(1);Vx]; Vy=[Vy(1);Vy];
        V=sqrt(Vx.^2+Vy.^2)/heads{s}(i);
        alpha=(atan2(Vy,Vx))/pi*180;
        da=diff(beta)*framerate;
        da=[da(1);da(:)];
        gamma=alpha-beta(1:end);
        Vl=V.*cosd(gamma);%longitudinal motion
        Vt=V.*sind(gamma);%transverse motion

        ind=pdat(:,1);
        str=pdat(:,2);
        phs=pdat(:,3);
        msp=pdat(:,4);
        rf=pdat(:,5);
        
        %remove small amplitude pumps
        minamp=3;
        ind=ind(rf>minamp);
        str=str(rf>minamp);
        phs=phs(rf>minamp);
        msp=msp(rf>minamp);        

        all_gamma{s}=[all_gamma{s};gamma];
        pump_gamma{s}=[pump_gamma{s};gamma(ind)];
        ppump_gamma{s}=[ppump_gamma{s};gamma(ind(phs==0))];
        rpump_gamma{s}=[rpump_gamma{s};gamma(ind(phs==1))];

        all_V{s}=[all_V{s};V];
        pump_V{s}=[pump_V{s};V(ind)];
        ppump_V{s}=[ppump_V{s};V(ind(phs==0))];
        rpump_V{s}=[rpump_V{s};V(ind(phs==1))];

        all_da{s}=[all_da{s};da];
        pump_da{s}=[pump_da{s};da(ind)]; 
        ppump_da{s}=[ppump_da{s};da(ind(phs==0))];
        rpump_da{s}=[rpump_da{s};da(ind(phs==1))];

        pump_str{s}=[pump_str{s};str];    
        ppump_str{s}=[ppump_str{s};str(phs==0)];
        rpump_str{s}=[rpump_str{s};str(phs==1)];

        all_Vl{s}=[all_Vl{s};Vl];
        pump_Vl{s}=[pump_Vl{s};Vl(ind)];
        ppump_Vl{s}=[ppump_Vl{s};Vl(ind(phs==0))];
        rpump_Vl{s}=[rpump_Vl{s};Vl(ind(phs==1))];

        all_Vt{s}=[all_Vt{s};Vt];
        pump_Vt{s}=[pump_Vt{s};Vt(ind)];
        ppump_Vt{s}=[ppump_Vt{s};Vt(ind(phs==0))];
        rpump_Vt{s}=[rpump_Vt{s};Vt(ind(phs==1))];
        
    end
    p_ppump(s)=numel(ppump_V{s})/numel(all_V{s}); %prior protraction pump
    p_rpump(s)=numel(rpump_V{s})/numel(all_V{s}); %prior retraction pump
end

%% get histograms

edges_V=quantile(flatten(all_V),p);
ctrs_V=quantile(flatten(all_V),c);
edges_g=quantile(flatten(all_gamma),p);
ctrs_g=quantile(flatten(all_gamma),c);
edges_da=quantile(flatten_mirror(all_da),p);
ctrs_da=quantile(flatten_mirror(all_da),c);
edges_Vl=quantile(flatten(all_Vl),p);
ctrs_Vl=quantile(flatten(all_Vl),c);
edges_Vt=quantile(flatten_mirror(all_Vt),p);
ctrs_Vt=quantile(flatten_mirror(all_Vt),c);

%% get 2D histograms and correlations (S2 Fig)
K=100;
x_vt=linspace(edges_Vt(2),edges_Vt(end-1),K);
x_da=linspace(edges_da(2),edges_da(end-1),K);
x_vl=linspace(edges_Vl(2),edges_Vl(end-3),K);
[X,Y]=meshgrid(x_vt,x_da);
f_da_vt=reshape(ksdensity([all_Vt{1},all_da{1}],[X(:),Y(:)]),K,K);
plot_2d(X,Y,f_da_vt);
f=fit(all_Vt{1},all_da{1},'poly1')
H=plot(f);
legend('off');
set(H,'LineWidth',3,'Color',[1 1 1],'LineStyle','--');
xlabel('');
ylabel('');

[X,Y]=meshgrid(x_vl,x_da);
f_da_vl=reshape(ksdensity([all_Vl{1},all_da{1}],[X(:),Y(:)]),K,K);
plot_2d(X,Y,f_da_vl);
f=fit(all_Vl{1},all_da{1},'poly1')
H=plot(f);
legend('off');
set(H,'LineWidth',3,'Color',[1 1 1],'LineStyle','--');
xlabel('');
ylabel('');

%% pooling both sides
display('rotation:')
plot_tuning(-ctrs_da,edges_da,flatten_mirror(all_da),flatten_mirror(ppump_da),flatten_mirror(rpump_da));    
display('long. vel:')
plot_tuning(ctrs_Vl,edges_Vl,all_Vl{s},ppump_Vl{s},rpump_Vl{s});    
display('trans. vel:')
plot_tuning(-ctrs_Vt,edges_Vt,flatten_mirror(all_Vt),flatten_mirror(ppump_Vt),flatten_mirror(rpump_Vt));    

end

function vec_out=flatten(cell_in)
vec_out=[];
for c=1:numel(cell_in)
    vec_out=[vec_out;cell_in{c}];
end
end

function vec_out=flatten_mirror(cell_in)
vec_out=[cell_in{1};-cell_in{2}];
end

% function plot_da_tuning(ctrs,hppump,hrpump,hda)
function plot_tuning(ctrs,edges_da,da,ppump,rpump)
global COL scale fsize
scale=1;
N=1e4;
p_prior=numel(ppump)/numel(da);
r_prior=numel(rpump)/numel(da);
[hda,edges_da]=histcounts(da,edges_da);
hppump=histcounts(ppump,edges_da)./hda/p_prior;
hrpump=histcounts(rpump,edges_da)./hda/r_prior;
display('prot-');
[p_pval,p_Q]=gen_monte(edges_da,da,numel(ppump),N,hppump);
display('ret-');
[r_pval,r_Q]=gen_monte(edges_da,da,numel(rpump),N,hrpump);

M=max([hppump(:);hrpump(:)])+0.1;
m=min([hppump(:);hrpump(:)])-0.1;
F=figure;
A=axes;
hold on;
% my_plotWithConfQ(ctrs,p_Q,COL(1,:));
% my_plotWithConfQ(ctrs,r_Q,COL(2,:));
H=plot(ctrs,hppump*scale,'-');
set(H,'LineWidth',3,'Color',COL(1,:),'MarkerSize',40,'Marker','.');
H=plot(ctrs,hrpump*scale,'o:');
set(H,'LineWidth',3,'Color',COL(2,:),'MarkerSize',10,'Marker','o');
H=plot([ctrs(1)*2 ctrs(end)*2],[1 1],'--');
set(H,'Color',[0 0 0]);
H=plot([0 0],[m/2 M*2],':');
set(H,'Color',[0 0 0]);
% ind=find(p_pval<0.05);
% H=plot(ctrs(ind),(m.*(hppump(ind)<1)+M.*(hppump(ind)>1)),'*');
% set(H,'Color',COL(1,:));
% ind=find(r_pval<0.05);
% H=plot(ctrs(ind),(m.*(hrpump(ind)<1)+M.*(hrpump(ind)>1)),'*');
% set(H,'Color',COL(2,:));

set(A,'FontSize',fsize,'XLim',[min(ctrs)*1.1 max(ctrs)*1.1],'Ylim',[m-0.15,M+0.15]);
end

function [pval,Q]=gen_monte(edges,data,pnum,N,hp)
q=[0.05 0.5 0.95];
prior=pnum/numel(data);
pumps=data(randi(numel(data),N,pnum));
monte=zeros(N,numel(edges)-1);
[h,edges]=histcounts(data,edges);
for n=1:N
    h_pump(n,:)=histcounts(pumps(n,:),edges)./h/prior;
end
pval_up=sum(h_pump>=ones(N,1)*hp)/N;
pval_down=sum(h_pump<=ones(N,1)*hp)/N;
pval=pval_up.*(hp>mean(h_pump))+pval_down.*(hp<=mean(h_pump))
Q=quantile(h_pump,q,1);
end

function plot_2d(X,Y,Z)
    global COL scale fsize
    F=figure;
    A=axes;
    hold on;
%     S=surf(X,Y,zeros(size(Z)),Z,'FaceColor','interp','LineStyle','none');
    imagesc(X(1,:),Y(:,1)',Z);
%     view(0,90);
    colormap('parula');
    set(A,'FontSize',fsize,'XLim',[min(X(:)) max(X(:))],'Ylim',[min(Y(:)) max(Y(:))]);
    set(F,'Position',[680 558 420 420]);
    set(A,'CLim',[0 0.03])
end
