function plot_figure4EF_S3CD()
global COL scale fsize
fsize=22;
msize=12;
F=figure;
g=0; gg=0.75;
COL=[zeros(1,3);g*ones(1,3);gg*ones(1,3)];

close(F);
framerate=500;

%% load data
hD=cell(2,1);
pD=cell(2,1);
heads=cell(2,1);

side='r_'; %right
Str=load(['../Data/comparative_c2h_',side,'data.mat']);
comp_hdata=Str.hdata;
comp_pdata=Str.pdata;
comp_headsize=Str.headsize;
Str=load(['../Data/control_c2h_',side,'data.mat']);
ctl_hdata=Str.hdata;
ctl_pdata=Str.pdata;
ctl_headsize=Str.headsize;
hD{1}={comp_hdata{:} ctl_hdata{:}};
pD{1}={comp_pdata{:} ctl_pdata{:}};
heads{1}=[comp_headsize ctl_headsize];

%head turn params
hminamp=20;
hmaxint=20000;

side=''; %left
Str=load(['../Data/comparative_c2h_',side,'data.mat']);
comp_hdata=Str.hdata;
comp_pdata=Str.pdata;
comp_headsize=Str.headsize;
Str=load(['../Data/control_c2h_',side,'data.mat']);
ctl_hdata=Str.hdata;
ctl_pdata=Str.pdata;
ctl_headsize=Str.headsize;
hD{2}={comp_hdata{:} ctl_hdata{:}};
pD{2}={comp_pdata{:} ctl_pdata{:}};
heads{2}=[comp_headsize ctl_headsize];

%% get stats
tpre=0.25; %pre window s
tpost=0.4; %post window s
intpre=tpre*framerate;
intpost=tpost*framerate;
ips_turns=cell(0);
cont_turns=cell(0);
head_avel=cell(0);
pind=cell(0);
rind=cell(0);
get_head_stats;

%% get time histograms
KC=5000; %number of control replicas
e=50;
edges=[-intpre:e:intpost];
bins=mean([edges(1:end-1);edges(2:end)])/framerate;
qh=[0.25 0.5 0.75];
p=[0.05 0.5 0.95];
ctl=0; 
[ppump_ipsi,ppump_contra,rpump_ipsi,rpump_contra,Qturn_p,Qturn_r]=get_turn_dynamics(ctl);
ctl=1;
if(ctl)
    tic;
    for k=1:KC
        k
        [ppump_ipsi_ctl(k,:),ppump_contra_ctl(k,:),rpump_ipsi_ctl(k,:),rpump_contra_ctl(k,:),Qp,Qr]=get_turn_dynamics(ctl);
        Qturn_p_ctl(k,:)=Qp(2,:);
        Qturn_r_ctl(k,:)=Qr(2,:);
        toc
        tic;
    end
    save('to_turns_ctl','ppump_ipsi_ctl','ppump_contra_ctl','rpump_ipsi_ctl','rpump_contra_ctl','Qturn_p_ctl','Qturn_r_ctl');
else
    S=load('to_turns_ctl.mat')
    ppump_ipsi_ctl=S.ppump_ipsi_ctl;
    ppump_contra_ctl=S.ppump_contra_ctl;
    rpump_ipsi_ctl=S.rpump_ipsi_ctl;
    rpump_contra_ctl=S.rpump_contra_ctl;
    Qturn_p_ctl=S.Qturn_p_ctl;
    Qturn_r_ctl=S.Qturn_r_ctl;
end

pval_p_ipsi=get_pval(ppump_ipsi,ppump_ipsi_ctl) 
pval_r_ipsi=get_pval(rpump_ipsi,rpump_ipsi_ctl) 
pval_p_contra=get_pval(ppump_contra,ppump_contra_ctl) 
pval_r_contra=get_pval(rpump_contra,rpump_contra_ctl) 
pvalQp=get_pval(Qturn_p(2,:),Qturn_p_ctl);
pvalQr=get_pval(Qturn_r(2,:),Qturn_r_ctl);

Qpctl=quantile([ppump_ipsi_ctl;ppump_contra_ctl],p,1);
Qrctl=quantile([rpump_ipsi_ctl;rpump_contra_ctl],p,1);
% normalize to mean ctl
ppump_ipsi=ppump_ipsi./Qpctl(2,:);
ppump_contra=ppump_contra./Qpctl(2,:);
Qpctl=Qpctl./([1;1;1]*Qpctl(2,:));
rpump_ipsi=rpump_ipsi./Qrctl(2,:);
rpump_contra=rpump_contra./Qrctl(2,:);
Qrctl=Qrctl./([1;1;1]*Qrctl(2,:));

t_turn=[-intpre:1:intpost-1];

F=figure;
A=axes;
yyaxis right;
my_plotWithConfQ(t_turn/framerate,Qturn_p,COL(3,:));
hold on;
Mq=max(Qturn_p(:))*1.1*ones(1,length(t_turn));
Mq(pvalQp>=0.025)=nan;
H=plot(t_turn/framerate,Mq,'-');
set(H,'LineWidth',4,'Color',COL(3,:));
set(A,'YColor',.5*COL(3,:),'YLim',1.2*max(abs(Qturn_p(:)))*[-1 1]);
set(A,'FontSize',fsize,'YLabel',[],'XLabel',[]);
yyaxis left;
H=plot(bins,ppump_ipsi,'^-');
set(H,'LineWidth',3,'MarkerSize',msize,'Marker','^','Color',COL(1,:),'MarkerFaceColor',COL(1,:));
hold on;
H=plot(bins,ppump_contra,'v--');
set(H,'LineWidth',3,'MarkerSize',msize,'Color',COL(1,:),'MarkerFaceColor',COL(1,:));
H=plot([0 0],[0 2],'-');
set(H,'Color',COL(2,:))
set(A,'YLim',0.45*[-1 1]+1,'YColor',COL(1,:),'Ylabel',[],'Xlabel',[]);
set(A,'YTick',[0.5:0.25:1.5]);
set(A,'XLim',[-0.25 0.35])

F=figure;
A=axes;
yyaxis right;
my_plotWithConfQ(t_turn/framerate,Qturn_r,COL(3,:));
hold on;
Mq=max(Qturn_r(:))*1.1*ones(1,length(t_turn));
Mq(pvalQr>=0.025)=nan;
H=plot(t_turn/framerate,Mq);
set(H,'LineWidth',4,'Color',COL(3,:));

set(A,'YColor',.5*COL(3,:),'YLim',1.2*max(abs(Qturn_p(:)))*[-1 1]);
set(A,'FontSize',fsize,'YLabel',[],'XLabel',[]);
yyaxis left;
H=plot(bins,rpump_ipsi,'^-');
set(H,'LineWidth',3,'MarkerSize',msize,'Marker','^','Color',COL(1,:));
hold on;
H=plot(bins,rpump_contra,'--v');
set(H,'LineWidth',3,'MarkerSize',msize,'Color',COL(1,:));
H=plot([0 0],[0 2],'-');
set(H,'Color',COL(2,:))
set(A,'YLim',0.45*[-1 1]+1,'YColor',COL(1,:));
set(A,'YTick',[0.5:0.25:1.5]);
set(A,'XLim',[-0.25 0.35])
drawnow;
%% threshold dependence
KC=5000;
minamp_v=linspace(5,50,10);
p5=[0.05 0.5 0.95];
p1=[0.01 0.5 0.99];
for a=1:numel(minamp_v)
    hminamp=minamp_v(a)
    get_head_stats;
    ctl=0;
    [post_p_ips_turns(a),post_p_cont_turns(a),post_r_ips_turns(a),post_r_cont_turns(a)]=get_turn_post(ctl);       
    ctl=1;
    if(ctl)
        for k=1:KC
                [post_p_ips_turns_ctl(k,a),post_p_cont_turns_ctl(k,a),...
                    post_r_ips_turns_ctl(k,a),post_r_cont_turns_ctl(k,a)]=get_turn_post(ctl);       
        end
        save('toturns_th_ctl.mat','post_p_ips_turns_ctl','post_p_cont_turns_ctl',...
            'post_r_ips_turns_ctl','post_r_cont_turns_ctl');    
    else
        S=load('toturns_th_ctl.mat')
        post_p_ips_turns_ctl=S.post_p_ips_turns_ctl;
        post_p_cont_turns_ctl=S.post_p_cont_turns_ctl;
        post_r_ips_turns_ctl=S.post_r_ips_turns_ctl;
        post_r_cont_turns_ctl=S.post_r_cont_turns_ctl;
    end
    
    Q5_p_ips_pre(a,:)=quantile(post_p_ips_turns_ctl(:,a),p5);
    Q5_p_cont_pre(a,:)=quantile(post_p_cont_turns_ctl(:,a),p5);
    Q5_r_ips_pre(a,:)=quantile(post_r_ips_turns_ctl(:,a),p5);
    Q5_r_cont_pre(a,:)=quantile(post_r_cont_turns_ctl(:,a),p5);
    Q1_p_ips_pre(a,:)=quantile(post_p_ips_turns_ctl(:,a),p1);
    Q1_p_cont_pre(a,:)=quantile(post_p_cont_turns_ctl(:,a),p1);
    Q1_r_ips_pre(a,:)=quantile(post_r_ips_turns_ctl(:,a),p1);
    Q1_r_cont_pre(a,:)=quantile(post_r_cont_turns_ctl(:,a),p1);    
end

%normalize
post_p_ips_turns=post_p_ips_turns./Q5_p_ips_pre(:,2)';
Q5_p_ips_pre=Q5_p_ips_pre./(Q5_p_ips_pre(:,2)*[1 1 1]);
post_p_cont_turns=post_p_cont_turns./Q5_p_cont_pre(:,2)';
Q5_p_cont_pre=Q5_p_cont_pre./(Q5_p_cont_pre(:,2)*[1 1 1]);
post_r_ips_turns=post_r_ips_turns./Q5_r_ips_pre(:,2)';
Q5_r_ips_pre=Q5_r_ips_pre./(Q5_r_ips_pre(:,2)*[1 1 1]);
post_r_cont_turns=post_r_cont_turns./Q5_r_cont_pre(:,2)';
Q5_r_cont_pre=Q5_r_cont_pre./(Q5_r_cont_pre(:,2)*[1 1 1]);

Q1_p_ips_pre=Q1_p_ips_pre./(Q1_p_ips_pre(:,2)*[1 1 1]);
Q1_p_cont_pre=Q1_p_cont_pre./(Q1_p_cont_pre(:,2)*[1 1 1]);
Q1_r_ips_pre=Q1_r_ips_pre./(Q1_r_ips_pre(:,2)*[1 1 1]);
Q1_r_cont_pre=Q1_r_cont_pre./(Q1_r_cont_pre(:,2)*[1 1 1]);

F=figure;
A=axes;
my_plotWithConfQ(minamp_v,Q5_p_ips_pre',COL(2,:));
hold on;
my_plotWithConfQ(minamp_v,Q1_p_ips_pre',COL(3,:));
H=plot(minamp_v,post_p_ips_turns,'-');
set(H,'LineWidth',3,'MarkerSize',msize,'Marker','^','Color',COL(1,:),'MarkerFaceColor',COL(1,:));
H=plot(minamp_v,post_p_cont_turns,'--');
set(H,'LineWidth',3,'MarkerSize',msize,'Marker','v','Color',COL(1,:),'MarkerFaceColor',COL(1,:));
set(A,'Ylim',[0.5 1.5],'Xlim',[minamp_v(1) minamp_v(end)],'FontSize',fsize,'Ylabel',[],'Xlabel',[]);

F=figure;
A=axes;
my_plotWithConfQ(minamp_v,Q5_r_ips_pre',COL(2,:));
hold on;
my_plotWithConfQ(minamp_v,Q1_r_ips_pre',COL(3,:));
H=plot(minamp_v,post_r_ips_turns,'-');
set(H,'LineWidth',3,'MarkerSize',msize,'Marker','^','Color',COL(1,:));
H=plot(minamp_v,post_r_cont_turns,'--');
set(H,'LineWidth',3,'MarkerSize',msize,'Marker','v','Color',COL(1,:));
set(A,'Ylim',[0.3 2],'Xlim',[minamp_v(1) minamp_v(end)],'FontSize',fsize,'Ylabel',[],'Xlabel',[]);


%% functions

    function [hist_p_ipsi,hist_p_contra,hist_r_ipsi,hist_r_contra,Qhead_p,Qhead_r]=get_turn_dynamics(ctl)
        
        total_ppumps=0;
        total_rpumps=0;

        ppump_ipsi_raster=[];
        rpump_ipsi_raster=[];
        ppump_cont_raster=[];
        rpump_cont_raster=[];
        
        head_dynamics_p=[];
        head_dynamics_r=[];
        
        total_frames=0;
        total_ips_turns=0;
        total_cont_turns=0;
        for j=1:numel(ips_turns)
                    
            %accumulate events for normalization
            total_frames=total_frames+numel(head_avel{j});
            iturns=ips_turns{j};            
            cturns=cont_turns{j};
            total_ips_turns=total_ips_turns+numel(iturns);
            total_cont_turns=total_cont_turns+numel(cturns);

            if(ctl==0)
                ppumps=pind{j};
                rpumps=rind{j};
            else
                ppumps=randi(numel(head_avel{j}),numel(pind{j}),1);
                rpumps=randi(numel(head_avel{j}),numel(rind{j}),1);
            end
            total_ppumps=total_ppumps+numel(pind{j});
            total_rpumps=total_rpumps+numel(rind{j});     
            
            %ppump
            if(numel(ppumps))
                I=find(ppumps>intpre & ppumps<(numel(head_avel{j})-intpost)); %find pumps with window around them
                %go over all pumps
                for n=1:numel(I)
                    T=ppumps(I(n));
                    Ii=find((iturns-T)>(-intpre) & (iturns-T)<intpost);
                    Ic=find((cturns-T)>(-intpre) & (cturns-T)<intpost);
                    ppump_ipsi_raster=[ppump_ipsi_raster   ; iturns(Ii)-T];
                    ppump_cont_raster=[ppump_cont_raster ; cturns(Ic)-T];
                    head_dynamics_p=[head_dynamics_p;head_avel{j}(T-intpre:T+intpost-1)'];
                end
            end

            %rpump
            if(numel(rpumps))
                I=find(rpumps>intpre & rpumps<(numel(head_avel{j})-intpost)); %find rpumps with window around them
                for n=1:numel(I)
                    T=rpumps(I(n));
                    Ii=find((iturns-T)>(-intpre) & (iturns-T)<intpost);
                    Ic=find((cturns-T)>(-intpre) & (cturns-T)<intpost);
                    rpump_ipsi_raster=[rpump_ipsi_raster   ; iturns(Ii)-T];
                    rpump_cont_raster=[rpump_cont_raster ; cturns(Ic)-T];
                    head_dynamics_r=[head_dynamics_r;head_avel{j}(T-intpre:T+intpost-1)'];
                end
            end
            
        end
        
        p_iturns=total_ips_turns/total_frames; %prior ipsi turns
        p_cturns=total_cont_turns/total_frames; %prior contra turns

        hist_p_ipsi=histcounts(ppump_ipsi_raster,edges)/total_ppumps/e/p_iturns;
        hist_p_contra=histcounts(ppump_cont_raster,edges)/total_ppumps/e/p_cturns;
        hist_r_ipsi=histcounts(rpump_ipsi_raster,edges)/total_rpumps/e/p_iturns;
        hist_r_contra=histcounts(rpump_cont_raster,edges)/total_rpumps/e/p_cturns;
        
        Qhead_p=[1;1;1]*nanmean(head_dynamics_p)+[-1;0;1]*nanstd(head_dynamics_p)/sqrt(size(head_dynamics_p,2));
        Qhead_r=[1;1;1]*nanmean(head_dynamics_r)+[-1;0;1]*nanstd(head_dynamics_r)/sqrt(size(head_dynamics_r,2));
    end

    function [post_p_ips_turns,post_p_cont_turns,post_r_ips_turns,post_r_cont_turns]=get_turn_post(ctl)       
                
        total_ppumps=0;
        total_rpumps=0;
        total_frames=0;

        total_ips_turns=0;
        total_cont_turns=0;

        post_p_ips_turns=0;
        post_p_cont_turns=0;
        post_r_ips_turns=0;
        post_r_cont_turns=0;
        
        for j=1:numel(ips_turns)
            
            %accumulate events for normalization
            total_frames=total_frames+numel(head_avel{j});
            iturns=ips_turns{j};            
            cturns=cont_turns{j};
            total_ips_turns=total_ips_turns+numel(iturns);
            total_cont_turns=total_cont_turns+numel(cturns);
                    
            if(ctl==0)
                ppumps=pind{j};
                rpumps=rind{j};
            else
                ppumps=randi(numel(head_avel{j}),numel(pind{j}),1);
                rpumps=randi(numel(head_avel{j}),numel(rind{j}),1);
            end
            total_ppumps=total_ppumps+numel(pind{j});
            total_rpumps=total_rpumps+numel(rind{j});     
            
            %ppump
            if(numel(ppumps))
                I=find(ppumps>intpre & ppumps<(numel(head_avel{j})-intpost)); %find pumps with window around them
                %go over all pumps
                for n=1:numel(I)
                    T=ppumps(I(n));
                    Ii=find((iturns-T)>(e/2) & (iturns-T)<(3*e/2));
                    Ic=find((cturns-T)>(e/2) & (cturns-T)<(3*e/2));
                    post_p_ips_turns=post_p_ips_turns+numel(Ii);
                    post_p_cont_turns=post_p_cont_turns+numel(Ic);
                end
            end

            %rpump
            if(numel(rpumps))
                I=find(rpumps>intpre & rpumps<(numel(head_avel{j})-intpost)); %find rpumps with window around them
                for n=1:numel(I)
                    T=rpumps(I(n));
                    Ii=find((iturns-T)>(e/2) & (iturns-T)<(3*e/2));
                    Ic=find((cturns-T)>(e/2) & (cturns-T)<(3*e/2));
                    post_r_ips_turns=post_r_ips_turns+numel(Ii);
                    post_r_cont_turns=post_r_cont_turns+numel(Ic);
                end
            end
            
        end

        p_iturns=total_ips_turns/total_frames; %prior ipsi turns
        p_cturns=total_cont_turns/total_frames; %prior contra turns

        post_p_ips_turns=post_p_ips_turns/e/total_ppumps/p_iturns;
        post_p_cont_turns=post_p_cont_turns/e/total_ppumps/p_cturns;
        post_r_ips_turns=post_r_ips_turns/e/total_rpumps/p_iturns;
        post_r_cont_turns=post_r_cont_turns/e/total_rpumps/p_cturns;
                
    end

    function get_head_stats()
        allamps=[];
        allints=[];
        j=1;
        for s=1:2 %both sides
            for i=1:numel(hD{s}) %go over all sections
                hdat=hD{s}{i};
                pdat=pD{s}{i};
                beta=hdat(:,1); %head direction
                da=diff(beta)*framerate; %head rotation: negative=leftward
                da=[da(1);da(:)];
                X=hdat(:,2);
                Y=hdat(:,3);
                Vx=diff(X)*framerate;
                Vy=diff(Y)*framerate;
                V=sqrt(Vx.^2+Vy.^2)/heads{s}(i);
                alpha=(atan2(Vy,Vx))/pi*180; %head motion direction
                gamma=alpha-beta(2:end); %head translational egocentric motion direciton, negative=leftward
                Vl=V.*cosd(gamma);%longitudinal motion
                Vt=V.*sind(gamma);%transverse motion

                %get turns
                [lturns,rturns,amps,ints]=find_head_turns(da,hminamp,hmaxint);
                allamps=[allamps;amps];
                allints=[allints;ints];

                ind=pdat(:,1);
                phs=pdat(:,3);
                rf=pdat(:,5); %pump am

                %remove small amplitude pumps
                minamp=2.5;
                ind=ind(rf>minamp);
                phs=phs(rf>minamp);
                pind{j}=ind(phs==0);
                rind{j}=ind(phs==1);

                dda=diff(da)*framerate; %head accel.
                dda=[dda(1);dda(:)];
                if s==1
                    ips_turns{j}=lturns;
                    cont_turns{j}=rturns;
                    head_avel{j}=dda;
                else
                    ips_turns{j}=rturns;
                    cont_turns{j}=lturns;
                    head_avel{j}=-dda;
                end
                j=j+1;

            end
        end
    end

    function pval=get_pval(p,pctl)        
        pval_up=sum(pctl>=ones(KC,1)*p)/KC;
        pval_down=sum(pctl<=ones(KC,1)*p)/KC;
        pval=pval_up.*(p>mean(pctl))+pval_down.*(p<=mean(pctl));
    end

end