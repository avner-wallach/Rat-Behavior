function plot_pump_headttrans()
global COL scale fsize
fsize=22;
F=figure;
% COL1=colormap('parula');
% COL2=colormap('copper');
% COL=colormap('lines');
g=0; gg=0.75;
COL=[zeros(1,3);g*ones(1,3);gg*ones(1,3)];

close(F);
framerate=500;

%% load data
hD=cell(2,1);
pD=cell(2,1);
heads=cell(2,1);

side='r_'; %right
Str=load(['comparative_c2h_',side,'data.mat']);
comp_hdata=Str.hdata;
comp_pdata=Str.pdata;
comp_headsize=Str.headsize;
Str=load(['control_c2h_',side,'data.mat']);
ctl_hdata=Str.hdata;
ctl_pdata=Str.pdata;
ctl_headsize=Str.headsize;
hD{1}={comp_hdata{:} ctl_hdata{:}};
pD{1}={comp_pdata{:} ctl_pdata{:}};
heads{1}=[comp_headsize ctl_headsize];

%head turn params
hminamp=.1;
hmaxint=20000;

side=''; %left
Str=load(['comparative_c2h_',side,'data.mat']);
comp_hdata=Str.hdata;
comp_pdata=Str.pdata;
comp_headsize=Str.headsize;
Str=load(['control_c2h_',side,'data.mat']);
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

%% plot turn amps and ints
F=figure;
A=subplot(1,3,1);
H=histogram(abs(allamps),0.1:.1:10);
set(A,'FontSize',fsize);
A=subplot(1,3,2);
H=histogram(allints/500,0:.04:1);
set(A,'FontSize',fsize);
A=subplot(1,3,3);
H=histogram(abs(allamps)./allints*500,0:2:30);
set(A,'FontSize',fsize);

%% get time histograms
e=50;
edges=[-intpre:e:intpost];
bins=mean([edges(1:end-1);edges(2:end)])/framerate;
qh=[0.25 0.5 0.75];
ctl=0;

[ppump_ipsi,ppump_contra,rpump_ipsi,rpump_contra,Qturn]=get_pump_dynamics(ctl);
ctl=1;
KC=1000; %number of control replicas
tic
for k=1:KC
    k
     [ppump_ipsi_ctl(k,:),ppump_contra_ctl(k,:),rpump_ipsi_ctl(k,:),rpump_contra_ctl(k,:)]=get_pump_dynamics(ctl);
     toc
     tic
end
save('headttrans_ctl.mat','ppump_ipsi_ctl','ppump_contra_ctl','rpump_ipsi_ctl','rpump_contra_ctl');

p=[0.05 0.5 0.95];

pval_p_ipsi=get_pval(ppump_ipsi,ppump_ipsi_ctl) 
pval_r_ipsi=get_pval(rpump_ipsi,rpump_ipsi_ctl) 
pval_p_contra=get_pval(ppump_contra,ppump_contra_ctl) 
pval_r_contra=get_pval(rpump_contra,rpump_contra_ctl) 

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
my_plotWithConfQ(t_turn/framerate,Qturn,COL(3,:));
set(A,'YColor',.5*COL(3,:),'YLim',1.2*max(abs(Qturn(:)))*[-1 1]);
set(A,'FontSize',fsize,'YLabel',[],'XLabel',[]);
yyaxis left;
H=plot(bins,ppump_ipsi,'-');
set(H,'LineWidth',3,'MarkerSize',40,'Marker','.','Color',COL(1,:));
hold on;
H=plot(bins,rpump_ipsi,'o:');
set(H,'LineWidth',3,'MarkerSize',10,'Color',COL(2,:));
H=plot([0 0],[0 2],'-');
set(H,'Color',COL(2,:))
% H=plot(bins,Qpctl,':k');
set(A,'YLim',0.5*[-1 1]+1,'YColor',COL(1,:));
set(A,'YTick',[0.5:0.25:1.5]);
set(A,'XLim',[-0.25 0.35])

F=figure;
A=axes;
yyaxis right;
my_plotWithConfQ(t_turn/framerate,Qturn,COL(3,:));
set(A,'YColor',.5*COL(3,:),'YLim',1.2*max(abs(Qturn(:)))*[-1 1]);
set(A,'FontSize',fsize,'YLabel',[],'XLabel',[]);
yyaxis left;
H=plot(bins,ppump_contra,'-');
set(H,'LineWidth',3,'MarkerSize',40,'Marker','.','Color',COL(1,:));
hold on;
H=plot(bins,rpump_contra,':o');
set(H,'LineWidth',3,'MarkerSize',10,'Color',COL(2,:));
H=plot([0 0],[0 2],'-');
set(H,'Color',COL(2,:))
set(A,'YLim',0.5*[-1 1]+1,'YColor',COL(1,:));
set(A,'YTick',[0.5:0.25:1.5]);
set(A,'XLim',[-0.25 0.35])
% H=plot(bins,Qrctl,':k');
return;
%% threshold dependence
minamp_v=linspace(0.1,2,10);
for a=1:numel(minamp_v)
%     [lturns,rturns,amps,ints]=find_head_turns(beta,minamp_v(i),hmaxint);
    hminamp=minamp_v(a);
    get_head_stats;
    ctl=0;
    [pre_ips_turns(a),pre_cont_turns(a)]=get_pump_pre(ctl);   
    ctl=1;
    for k=1:KC
        [pre_ips_turns_ctl(k),pre_cont_turns_ctl(k)]=get_pump_pre(ctl);   
    end
    
    Q_ips_pre(a,:)=quantile(pre_ips_turns_ctl,p);
    Q_cont_pre(a,:)=quantile(pre_cont_turns_ctl,p);
end

%normalize
pre_ips_turns=pre_ips_turns./Q_ips_pre(:,2)';
Q_ips_pre=Q_ips_pre./(Q_ips_pre(:,2)*[1 1 1]);
pre_cont_turns=pre_cont_turns./Q_cont_pre(:,2)';
Q_cont_pre=Q_cont_pre./(Q_cont_pre(:,2)*[1 1 1]);


F=figure;
A=axes;
H=plot(minamp_v,pre_ips_turns,'--*');
hold on;
H=plot(minamp_v,pre_cont_turns,'+:');
H=plot(minamp_v,Q_ips_pre(:,1),':k');
H=plot(minamp_v,Q_ips_pre(:,3),':k');

%% functions

    function [hist_p_ipsi,hist_p_contra,hist_r_ipsi,hist_r_contra,Qhead]=get_pump_dynamics(ctl)
        
        total_frames=0;
        total_ppumps=0;
        total_rpumps=0;

        ppump_ipsi_raster=[];
        rpump_ipsi_raster=[];
        ppump_cont_raster=[];
        rpump_cont_raster=[];
        
        head_dynamics=[];
        
        total_frames=0;
        total_ips_turns=0;
        total_cont_turns=0;
        for j=1:numel(ips_turns)
                    
            total_frames=total_frames+numel(head_avel{j});
            total_ppumps=total_ppumps+numel(pind{j});
            total_rpumps=total_rpumps+numel(rind{j});     
            if(ctl==0)
                iturns=ips_turns{j};
                cturns=cont_turns{j};
            else
                iturns=randi(numel(head_avel{j}),numel(ips_turns{j}));
                cturns=randi(numel(head_avel{j}),numel(cont_turns{j}));
            end
            total_ips_turns=total_ips_turns+numel(iturns);
            total_cont_turns=total_cont_turns+numel(cturns);
            
            %ipsi head turns
            if(numel(iturns))
                I=find(iturns>intpre & iturns<(numel(head_avel{j})-intpost)); %find turns with window around them
                for n=1:numel(I)
                    T=iturns(I(n));
                    Ip=find((pind{j}-T)>(-intpre) & (pind{j}-T)<intpost);
                    Ir=find((rind{j}-T)>(-intpre) & (rind{j}-T)<intpost);
                    ppump_ipsi_raster=[ppump_ipsi_raster   ; pind{j}(Ip)-T];
                    rpump_ipsi_raster=[rpump_ipsi_raster ; rind{j}(Ir)-T];
                    head_dynamics=[head_dynamics;head_avel{j}(T-intpre:T+intpost-1)'];
                end
            end

            %contra head turns
            if(numel(cturns))
                I=find(cturns>intpre & cturns<(numel(head_avel{j})-intpost)); %find turns with window around them
                for n=1:numel(I)
                    T=cturns(I(n));
                    Ip=find((pind{j}-T)>(-intpre) & (pind{j}-T)<intpost);
                    Ir=find((rind{j}-T)>(-intpre) & (rind{j}-T)<intpost);
                    ppump_cont_raster=[ppump_cont_raster   ; pind{j}(Ip)-T];
                    rpump_cont_raster=[rpump_cont_raster ; rind{j}(Ir)-T];
%                     if(s==1)
%                         lturn_da=[lturn_da;da(T-intpre:T+intpost-1)'];
%                     end
                end
            end
            
        end
        
        p_ppump=total_ppumps/total_frames; %prior protraction pump
        p_rpump=total_rpumps/total_frames; %prior retraction pump

        hist_p_ipsi=histcounts(ppump_ipsi_raster,edges)/total_ips_turns/e/p_ppump;
        hist_p_contra=histcounts(ppump_cont_raster,edges)/total_cont_turns/e/p_ppump;
        hist_r_ipsi=histcounts(rpump_ipsi_raster,edges)/total_ips_turns/e/p_rpump;
        hist_r_contra=histcounts(rpump_cont_raster,edges)/total_cont_turns/e/p_rpump;
        
        if(ctl==0)
            Qhead=quantile(head_dynamics,qh);
        else
            Qhead=[];
        end
    end

    function [pre_ips_turns,pre_cont_turns]=get_pump_pre(ctl)       
                
        total_ips_turns=0;
        total_cont_turns=0;
        pre_ips_turns=0;
        pre_cont_turns=0;
        for j=1:numel(ips_turns)
                    
            if(ctl==0)
                iturns=ips_turns{j};
                cturns=cont_turns{j};
            else
                iturns=randi(numel(head_avel{j}),numel(ips_turns{j}));
                cturns=randi(numel(head_avel{j}),numel(cont_turns{j}));
            end
            total_ips_turns=total_ips_turns+numel(iturns);
            total_cont_turns=total_cont_turns+numel(cturns);
            
            %ipsi head turns
            if(numel(iturns))
                I=find(iturns>intpre & iturns<(numel(head_avel{j})-intpost)); %find turns with window around them
                for n=1:numel(I)
                    T=iturns(I(n));
                    Ip=find((pind{j}-T)>(-e/2) & (pind{j}-T)<e/2);
                    Ir=find((rind{j}-T)>(-e/2) & (rind{j}-T)<e/2);
                    pre_ips_turns=pre_ips_turns+numel(Ip);
                end
            end

            %contra head turns
            if(numel(cturns))
                I=find(cturns>intpre & cturns<(numel(head_avel{j})-intpost)); %find turns with window around them
                for n=1:numel(I)
                    T=cturns(I(n));
                    Ip=find((pind{j}-T)>(-e/2) & (pind{j}-T)<e/2);
                    Ir=find((rind{j}-T)>(-e/2) & (rind{j}-T)<e/2);
                    pre_cont_turns=pre_cont_turns+numel(Ip);
                end
            end
            
        end
        
        pre_ips_turns=pre_ips_turns/e/total_ips_turns;
        pre_cont_turns=pre_cont_turns/e/total_cont_turns;
    end

    function get_head_stats()
        allamps=[];
        allints=[];
        j=1;
        for s=1:2 %both sides
            for i=1:numel(hD{s}) %go over all sections
                hdat=hD{s}{i};
                pdat=pD{s}{i};
                beta=hdat(:,1);
                da=diff(beta)*framerate;
                da=[da(1);da(:)];
                X=hdat(:,2);
                Y=hdat(:,3);
                Vx=diff(X)*framerate;
                Vy=diff(Y)*framerate;
                V=sqrt(Vx.^2+Vy.^2)/heads{s}(i);
                alpha=(atan2(Vy,Vx))/pi*180;
                gamma=alpha-beta(2:end);
                Vl=V.*cosd(gamma);%longitudinal motion
                Vt=V.*sind(gamma);%transverse motion

                %get turns
                [lturns,rturns,amps,ints]=find_head_turns(Vt,hminamp,hmaxint);
                allamps=[allamps;amps];
                allints=[allints;ints];

                ind=pdat(:,1);
                phs=pdat(:,3);
                rf=pdat(:,5); %pump am

                %remove small amplitude pumps
                minamp=2;
                ind=ind(rf>minamp);
                phs=phs(rf>minamp);
                pind{j}=ind(phs==0);
                rind{j}=ind(phs==1);

                if s==1
                    ips_turns{j}=lturns;
                    cont_turns{j}=rturns;
                    head_avel{j}=da;
                else
                    ips_turns{j}=rturns;
                    cont_turns{j}=lturns;
                    head_avel{j}=-da;
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