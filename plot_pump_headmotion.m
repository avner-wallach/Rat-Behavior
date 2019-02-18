function plot_pump_headmotion()
global COL scale fsize
fsize=22;
F=figure;
% COL1=colormap('parula');
% COL2=colormap('copper');
% COL=colormap('lines');
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
load(['comparative_c2h_',side,'data.mat']);
comp_hdata=hdata;
comp_pdata=pdata;
comp_headsize=headsize;
load(['control_c2h_',side,'data.mat']);
ctl_hdata=hdata;
ctl_pdata=pdata;
ctl_headsize=headsize;
hD{1}={comp_hdata{:} ctl_hdata{:}};
pD{1}={comp_pdata{:} ctl_pdata{:}};
heads{1}=[comp_headsize ctl_headsize];

side=''; %left
load(['comparative_c2h_',side,'data.mat']);
comp_hdata=hdata;
comp_pdata=pdata;
comp_headsize=headsize;
load(['control_c2h_',side,'data.mat']);
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
%         plot_trajectories(hdat,pdat);
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
        %    k=get_curvature(X,Y);

        ind=pdat(:,1);
        str=pdat(:,2);
        phs=pdat(:,3);
        msp=pdat(:,4);
        rf=pdat(:,5);
        
        %remove small amplitude pumps
        minamp=2;
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

%% by side
byside=0;
if(byside)
    for s=1:2

        [h_V,edges_V]=histcounts(all_V{s},edges_V);
        h_pump_V=histcounts(pump_V{s},edges_V);
        h_ppump_V=histcounts(ppump_V{s},edges_V);
        h_rpump_V=histcounts(rpump_V{s},edges_V);

        [h_g,edges_g]=histcounts(all_gamma{s},edges_g);
        h_pump_g=histcounts(pump_gamma{s},edges_g);
        h_ppump_g=histcounts(ppump_gamma{s},edges_g);%/p_ppump;
        h_rpump_g=histcounts(rpump_gamma{s},edges_g);%/p_rpump;


    %     h_pump_da=histcounts(pump_da{s},edges_da);
    %     [h_da,edges_da]=histcounts(all_da{s},edges_da);
    %     h_ppump_da=histcounts(ppump_da{s},edges_da);
    %     h_rpump_da=histcounts(rpump_da{s},edges_da);
    %     plot_da_tuning(ctrs_da,h_ppump_da,h_rpump_da,h_da);
        plot_tuning(ctrs_da,edges_da,all_da{s},ppump_da{s},rpump_da{s});    
        plot_tuning(ctrs_Vl,edges_Vl,all_Vl{s},ppump_Vl{s},rpump_Vl{s});    
        plot_tuning(ctrs_Vt,edges_Vt,all_Vt{s},ppump_Vt{s},rpump_Vt{s});    

    %     plot_V_tuning(ctrs_V,h_ppump_V,h_rpump_V,h_V);    
    %     plot_V_tuning(ctrs_Vl,h_ppump_V,h_rpump_V,h_V);    
    %     plot_gamma_tuning(ctrs_g,h_ppump_g,h_rpump_g,h_g);

    end
end
%% pooling both sides
display('rotation:')
plot_tuning(-ctrs_da,edges_da,flatten_mirror(all_da),flatten_mirror(ppump_da),flatten_mirror(rpump_da));    
display('long. vel:')
plot_tuning(ctrs_Vl,edges_Vl,all_Vl{s},ppump_Vl{s},rpump_Vl{s});    
display('trans. vel:')
plot_tuning(-ctrs_Vt,edges_Vt,flatten_mirror(all_Vt),flatten_mirror(ppump_Vt),flatten_mirror(rpump_Vt));    

%% by prot/ret
% [h_V,edges_V]=histcounts(flatten(all_V),edges_V);
% h_pump_V=histcounts(flatten(pump_V),edges_V);
% h_ppump_V=histcounts(flatten(ppump_V),edges_V);
% h_rpump_V=histcounts(flatten(rpump_V),edges_V);
% plot_V_tuning(ctrs_V,h_ppump_V,h_rpump_V,h_V);
% 
% [h_g,edges_g]=histcounts(flatten(all_gamma),edges_g);
% h_pump_g=histcounts(flatten(pump_gamma),edges_g);
% h_ppump_g=histcounts(flatten(ppump_gamma),edges_g);%/p_ppump;
% h_rpump_g=histcounts(flatten(rpump_gamma),edges_g);%/p_rpump;
% plot_gamma_tuning(ctrs_g,h_ppump_g,h_rpump_g,h_g);
% 
% h_right_pump_g=histcounts([ppump_gamma{1};rpump_gamma{1}],edges_g);
% h_left_pump_g=histcounts([ppump_gamma{2};rpump_gamma{2}],edges_g);
% plot_gamma_tuning(ctrs_g,h_right_pump_g,h_left_pump_g,h_g);
% 
% 
end

function plot_trajectories(hdat,pdat)
%%  plot head trajectory + direction
F=figure;
COL1=colormap('parula');
COL2=colormap('copper');
COL=colormap('lines');
K=20;
R=10;
plot(hdat(:,2),hdat(:,3));
hold on;
ind=[1:K:size(hdat,1)];
X=hdat(:,2);
Y=hdat(:,3);
U=R*cos(hdat(:,1)/180*pi);
V=R*sin(hdat(:,1)/180*pi);
for i=1:numel(ind)
    H=quiver(X(ind(i)),Y(ind(i)),U(ind(i)),V(ind(i)));
    set(H,'Color',COL1(ceil((i)/numel(ind)*64),:));
end

%% plot pumps
ind=pdat(:,1);
for i=1:numel(ind)
    if(pdat(i,3)) %retraction
        H=plot(X(ind(i)),Y(ind(i)),'v');
    else %protaction
        H=plot(X(ind(i)),Y(ind(i)),'*');
    end
%     H=quiver(X(ind(i)),Y(ind(i)),U(ind(i)),V(ind(i)));
    set(H,'Color',COL2(ceil(min(pdat(i,2),1)*64),:));
end
%% get velocity + heading changes
% Vx=diff(X);
% Vy=diff(Y);
% V=sqrt(Vx.^2+Vy.^2);
% alpha=(atan2(Vx,Vy))/pi*180;
% 
% figure;
% plot(V);
% hold on;
% plot(ind,zeros(size(ind)),'*k');
% figure;
% plot(alpha);
% hold on;
% plot(hdat(:,1));
% plot(ind,zeros(size(ind)),'*k');
% figure;
% plot(hdat(2:end,1)-alpha(:));
% hold on;
% plot(ind,zeros(size(ind)),'*k');

end

function k=get_curvature(X,Y)
J=1; %jump factor
for i=1:numel(X)-2*J
    x1=X(i); x2=X(i+J); x3=X(i+2*J);
    y1=Y(i); y2=Y(i+J); y3=Y(i+2*J);
    k(i) = 2*((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)) ./ ...
        sqrt(((x2-x1).^2+(y2-y1).^2)*((x3-x1).^2+(y3-y1).^2)*((x3-x2).^2+(y3-y2).^2));
end
k=[k(1)*ones(J,1);k(:);k(end)*ones(J,1)];

% p=linspace(0,1,64);
% q=quantile(k,p);
% C=discretize(k,q);
% figure;
% COL=colormap('parula');
% for i=1:numel(X)-2
% H=plot(X(i:i+1),Y(i:i+1));
% set(H,'Color',COL(C(i),:));
% hold on;
% end
end

function [U,V,Th,Rho]=meanvec(theta,rho)
Vec=sum(rho.*exp(i*theta));
U=real(Vec);
V=imag(Vec);
Th=angle(Vec);
Rho=abs(Vec);
end

function vec_out=flatten(cell_in)
vec_out=[];
for c=1:numel(cell_in)
    vec_out=[vec_out;cell_in{c}];
end
end

function vec_out=flatten_mirror(cell_in)
% vec_out=[];
% for c=1:numel(cell_in)
%     vec_out=[vec_out;cell_in{c}];
% end
vec_out=[cell_in{1};-cell_in{2}];
end

function plot_V_tuning(ctrs,hppump,hrpump,hv)
global COL scale
F=figure;
A=axes;
hold on;
H=plot(ctrs,hppump./hv*scale,'^:');
set(H,'LineWidth',2,'Color',COL(1,:));
H=plot(ctrs,hrpump./hv*scale,'v:');
set(H,'LineWidth',2,'Color',COL(2,:));
set(A,'XScale','log','FontSize',14,'XLim',[ctrs(1) ctrs(end)],'Ylim',[0 1.2]);
set(A,'XTick',[0.1 0.5 1 2 4],'XTickLabel',{'0.1','0.5','1','2','4'});
end    

function plot_gamma_tuning(ctrs,hppump,hrpump,hg)
global COL scale
F=figure;
theta=ctrs/180*pi;
theta=[theta theta(1)];
rho=hppump./hg*scale;
rho=[rho rho(1)];
H=polarplot(theta,rho);
set(H,'Color',COL(1,:),'LineWidth',2);
hold on;
[U,V,Th,Rho]=meanvec(theta,rho);
R=3; Ar_angle=pi/20; Ar_size=0.8;
H=polarplot([Th Th],[0 Rho/R]);
set(H,'Color',COL(1,:),'LineWidth',3);
H=polarplot([Th Th+Ar_angle],[Rho/R Rho/R*Ar_size]);
set(H,'Color',COL(1,:),'LineWidth',3);
H=polarplot([Th Th-Ar_angle],[Rho/R Rho/R*Ar_size]);
set(H,'Color',COL(1,:),'LineWidth',3);

rho=hrpump./hg*scale;
rho=[rho rho(1)];
H=polarplot(theta,rho);
set(H,'Color',COL(2,:),'LineWidth',2);
[U,V,Th,Rho]=meanvec(theta,rho);
H=polarplot([Th Th],[0 Rho/R]);
set(H,'Color',COL(2,:),'LineWidth',3);
H=polarplot([Th Th+Ar_angle],[Rho/R Rho/R*Ar_size]);
set(H,'Color',COL(2,:),'LineWidth',3);
H=polarplot([Th Th-Ar_angle],[Rho/R Rho/R*Ar_size]);
set(H,'Color',COL(2,:),'LineWidth',3);
A=gca;
A.ThetaZeroLocation='top';
set(A,'FontSize',14);
% A.ThetaTickLabel=[];
A.ThetaDir='clockwise';
A.RAxisLocation=210;
A.LineWidth=1
A.RLim=[0 1.2];
A.RTick=[0:0.2:1.2];

end

% function plot_da_tuning(ctrs,hppump,hrpump,hda)
function plot_tuning(ctrs,edges_da,da,ppump,rpump)
global COL scale fsize
scale=1;
N=1e4;
p_prior=numel(ppump)/numel(da);
r_prior=numel(rpump)/numel(da);
[hda,edges_da]=histcounts(da,edges_da);
hppump=histcounts(ppump,edges_da)./hda;%/p_prior;
hrpump=histcounts(rpump,edges_da)./hda;%/r_prior;
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