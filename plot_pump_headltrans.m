function plot_pump_headltrans()
% pump prob in time, rel to onset of lognitudinal direction change
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
hminamp=0.1;
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
for_turns=cell(0);
back_turns=cell(0);
head_avel=cell(0);
pind=cell(0);
rind=cell(0);
get_head_stats;

%% plot turn amps and ints
F=figure;
A=subplot(1,3,1);
H=histogram(abs(allamps),0:.4:10);
set(A,'FontSize',fsize);
A=subplot(1,3,2);
H=histogram(allints/500,0:.04:1);
set(A,'FontSize',fsize);
A=subplot(1,3,3);
H=histogram(abs(allamps)./allints*500,0:2:50);
set(A,'FontSize',fsize);

%% get time histograms
e=50;
edges=[-intpre:e:intpost];
bins=mean([edges(1:end-1);edges(2:end)])/framerate;
qh=[0.25 0.5 0.75];
p=[0.05 0.5 0.95];
ctl=0;
t_turn=[-intpre:1:intpost-1];

[ppump_for,ppump_back,rpump_for,rpump_back,Qturnf,Qturnb]=get_pump_dynamics(ctl);
ctl=1;
KC=1000; %number of control replicas
% for k=1:KC
%     k
%     [ppump_for_ctl(k,:),ppump_back_ctl(k,:),rpump_for_ctl(k,:),rpump_back_ctl(k,:)]=get_pump_dynamics(ctl);
%     toc;
%     tic;
% end
% save('headltrans_ctl','ppump_for_ctl','ppump_back_ctl','rpump_for_ctl','rpump_back_ctl');
S=load('headltrans_ctl','ppump_for_ctl','ppump_back_ctl','rpump_for_ctl','rpump_back_ctl');
ppump_for_ctl=S.ppump_for_ctl;
ppump_back_ctl=S.ppump_back_ctl;
rpump_for_ctl=S.rpump_for_ctl;
rpump_back_ctl=S.rpump_back_ctl;

pval_p_for=get_pval(ppump_for,ppump_for_ctl) 
pval_r_for=get_pval(rpump_for,rpump_for_ctl) 
pval_p_back=get_pval(ppump_back,ppump_back_ctl) 
pval_r_back=get_pval(rpump_back,rpump_back_ctl) 

%forward motion- ret and prot pumps
Qpctl=quantile(ppump_for_ctl,p,1);
Qrctl=quantile(rpump_for_ctl,p,1);
% normalize to mean ctl
ppump_for=ppump_for./Qpctl(2,:);
Qpctl=Qpctl./([1;1;1]*Qpctl(2,:));
rpump_for=rpump_for./Qrctl(2,:);
Qrctl=Qrctl./([1;1;1]*Qrctl(2,:));

F=figure;
A=axes;
yyaxis right;
my_plotWithConfQ(t_turn/framerate,Qturnf,COL(3,:));
set(A,'FontSize',fsize,'YColor',COL(2,:),'YLabel',[],'XLabel',[],'YColor',0.5*COL(3,:));
set(A,'YLim',max(abs(Qturnf(:)))*1.2*[-1 1]);
yyaxis left;
H=plot(bins,ppump_for,'-');
set(H,'LineWidth',3,'MarkerSize',40,'Marker','.','Color',COL(1,:));
hold on;
H=plot(bins,rpump_for,'o:');
set(H,'LineWidth',3,'MarkerSize',10,'Color',COL(2,:));
% H=plot(bins,Qpctl,':k');
% H=plot(bins,Qrctl,':r');
set(A,'YLim',1+0.35*[-1 1],'YColor',COL(1,:));
set(A,'Xlim',[-0.25 0.35]);

%backward motion- prot and ret pumps
Qpctl=quantile(ppump_back_ctl,p,1);
Qrctl=quantile(rpump_back_ctl,p,1);
% normalize to mean ctl
ppump_back=ppump_back./Qpctl(2,:);
Qpctl=Qpctl./([1;1;1]*Qpctl(2,:));
rpump_back=rpump_back./Qrctl(2,:);
Qrctl=Qrctl./([1;1;1]*Qrctl(2,:));

F=figure;
A=axes;
yyaxis right;
my_plotWithConfQ(t_turn/framerate,Qturnb,COL(3,:));
set(A,'FontSize',fsize,'YColor',COL(2,:),'YLabel',[],'XLabel',[],'YColor',0.5*COL(3,:));
set(A,'YLim',max(abs(Qturnb(:)))*1.2*[-1 1]);
yyaxis left;
H=plot(bins,ppump_back,'-');
set(H,'LineWidth',3,'MarkerSize',40,'Marker','.','Color',COL(1,:));
hold on;
H=plot(bins,rpump_back,'o:');
set(H,'LineWidth',3,'MarkerSize',10,'Color',COL(2,:));
set(A,'YLim',1+0.35*[-1 1],'YColor',COL(1,:));
set(A,'Xlim',[-0.25 0.35]);
% H=plot(bins,Qpctl,':k');
% H=plot(bins,Qrctl,':r');
% set(A,'YLim',[0.7 1.3]);

drawnow;
%% threshold dependence
KC=1000;
minamp_v=linspace(0.1,1,10);
p5=[0.05 0.5 0.95];
p1=[0.01 0.5 0.99];
for a=1:numel(minamp_v)
%     [lturns,rturns,amps,ints]=find_head_turns(beta,minamp_v(i),hmaxint);
    hminamp=minamp_v(a)
    get_head_stats;
    ctl=0;
    [pre_for_turns(a),pre_back_turns(a)]=get_pump_pre(ctl);   
    ctl=1;
    for k=1:KC
        [pre_for_turns_ctl(k,a),pre_back_turns_ctl(k,a)]=get_pump_pre(ctl);   
    end
    save('headltrans_th_ctl.mat','pre_for_turns_ctl','pre_back_turns_ctl');
    Q5_for_pre(a,:)=quantile(pre_for_turns_ctl(:,a),p5);
    Q5_back_pre(a,:)=quantile(pre_back_turns_ctl(:,a),p5);
    Q1_for_pre(a,:)=quantile(pre_for_turns_ctl(:,a),p1);
    Q1_back_pre(a,:)=quantile(pre_back_turns_ctl(:,a),p1);

end

%normalize
pre_for_turns=pre_for_turns./Q5_for_pre(:,2)';
Q5_for_pre=Q5_for_pre./(Q5_for_pre(:,2)*[1 1 1]);
Q1_for_pre=Q1_for_pre./(Q1_for_pre(:,2)*[1 1 1]);

pre_back_turns=pre_back_turns./Q5_back_pre(:,2)';
Q5_back_pre=Q5_back_pre./(Q5_back_pre(:,2)*[1 1 1]);
Q1_back_pre=Q1_back_pre./(Q1_back_pre(:,2)*[1 1 1]);


F=figure;
A=axes;
my_plotWithConfQ(minamp_v,Q5_for_pre',COL(2,:));
hold on;
my_plotWithConfQ(minamp_v,Q1_for_pre',COL(3,:));
H=plot(minamp_v,pre_for_turns,'-');
set(H,'LineWidth',3,'MarkerSize',40,'Marker','.','Color',COL(1,:));
% set(A,'YLim',[0.8 1.5]);

% F=figure;
% A=axes;
% my_plotWithConfQ(minamp_v,Q5_back_pre',COL(2,:));
% hold on;
% my_plotWithConfQ(minamp_v,Q1_back_pre',COL(3,:));
H=plot(minamp_v,pre_back_turns,'--');
set(H,'LineWidth',3,'MarkerSize',22,'Marker','*','Color',COL(1,:));
set(A,'YLim',[0.8 1.5],'Xlim',[minamp_v(1) minamp_v(end)],'Font',fsize,'Ylabel',[],'Xlabel',[]);

%% functions

    function [hist_p_for,hist_p_back,hist_r_for,hist_r_back,Qhead_for,Qhead_back]=get_pump_dynamics(ctl)
        
        total_frames=0;
        total_ppumps=0;
        total_rpumps=0;

        ppump_for_raster=[];
        rpump_for_raster=[];
        ppump_back_raster=[];
        rpump_back_raster=[];
        
        head_dynamics_for=[];
        head_dynamics_back=[];
        
        total_frames=0;
        total_for_turns=0;
        total_back_turns=0;
        for j=1:numel(for_turns)
                    
            total_frames=total_frames+numel(head_avel{j});
            total_ppumps=total_ppumps+numel(pind{j});
            total_rpumps=total_rpumps+numel(rind{j});     
            if(ctl==0)
                iturns=for_turns{j};
                cturns=back_turns{j};
            else
                iturns=randi(numel(head_avel{j}),numel(for_turns{j}));
                cturns=randi(numel(head_avel{j}),numel(back_turns{j}));
            end
            total_for_turns=total_for_turns+numel(iturns);
            total_back_turns=total_back_turns+numel(cturns);
            
            %forward head turns
            if(numel(iturns))
                I=find(iturns>intpre & iturns<(numel(head_avel{j})-intpost)); %find turns with window around them
                for n=1:numel(I)
                    T=iturns(I(n));
                    Ip=find((pind{j}-T)>(-intpre) & (pind{j}-T)<intpost);
                    Ir=find((rind{j}-T)>(-intpre) & (rind{j}-T)<intpost);
                    ppump_for_raster=[ppump_for_raster   ; pind{j}(Ip)-T];
                    rpump_for_raster=[rpump_back_raster ; rind{j}(Ir)-T];
                    head_dynamics_for=[head_dynamics_for;head_avel{j}(T-intpre:T+intpost-1)'];
                end
            end

            %back head turns
            if(numel(cturns))
                I=find(cturns>intpre & cturns<(numel(head_avel{j})-intpost)); %find turns with window around them
                for n=1:numel(I)
                    T=cturns(I(n));
                    Ip=find((pind{j}-T)>(-intpre) & (pind{j}-T)<intpost);
                    Ir=find((rind{j}-T)>(-intpre) & (rind{j}-T)<intpost);
                    ppump_back_raster=[ppump_back_raster   ; pind{j}(Ip)-T];
                    rpump_back_raster=[rpump_back_raster ; rind{j}(Ir)-T];
                    head_dynamics_back=[head_dynamics_back;head_avel{j}(T-intpre:T+intpost-1)'];
                end
            end
            
        end
        
        p_ppump=total_ppumps/total_frames; %prior protraction pump
        p_rpump=total_rpumps/total_frames; %prior retraction pump

        hist_p_for=histcounts(ppump_for_raster,edges)/total_for_turns/e/p_ppump;
        hist_p_back=histcounts(ppump_back_raster,edges)/total_back_turns/e/p_ppump;
        hist_r_for=histcounts(rpump_for_raster,edges)/total_for_turns/e/p_rpump;
        hist_r_back=histcounts(rpump_back_raster,edges)/total_back_turns/e/p_rpump;
        
        if(ctl==0)
            Qhead_for=quantile(head_dynamics_for,qh);
            Qhead_back=quantile(head_dynamics_back,qh);
        else
            Qhead=[];
        end
    end

    function [pre_for_turns,pre_back_turns]=get_pump_pre(ctl)       
                
        total_for_turns=0;
        total_back_turns=0;
        pre_for_turns=0;
        pre_back_turns=0;
        for j=1:numel(for_turns)
                    
            if(ctl==0)
                iturns=for_turns{j};
                cturns=back_turns{j};
            else
                iturns=randi(numel(head_avel{j}),numel(for_turns{j}));
                cturns=randi(numel(head_avel{j}),numel(back_turns{j}));
            end
            total_for_turns=total_for_turns+numel(iturns);
            total_back_turns=total_back_turns+numel(cturns);
            
            %ipsi head turns
            if(numel(iturns))
                I=find(iturns>intpre & iturns<(numel(head_avel{j})-intpost)); %find turns with window around them
                for n=1:numel(I)
                    T=iturns(I(n));
                    Ip=find((pind{j}-T)>(-3*e/2) & (pind{j}-T)<-e/2);
                    Ir=find((rind{j}-T)>(-3*e/2) & (rind{j}-T)<-e/2);
                    pre_for_turns=pre_for_turns+numel(Ip);
                end
            end

            %contra head turns
            if(numel(cturns))
                I=find(cturns>intpre & cturns<(numel(head_avel{j})-intpost)); %find turns with window around them
                for n=1:numel(I)
                    T=cturns(I(n));
                    Ip=find((pind{j}-T)>(-3*e/2) & (pind{j}-T)<-e/2);
                    Ir=find((rind{j}-T)>(-3*e/2) & (rind{j}-T)<-e/2);
                    pre_back_turns=pre_back_turns+numel(Ip);
                end
            end
            
        end
        
        pre_for_turns=pre_for_turns/e/total_for_turns;
        pre_back_turns=pre_back_turns/e/total_back_turns;
    end

    function get_head_stats()
        allamps=[];
        allints=[];
        j=1;
        for s=2:2 % ONLY LEFT SIDE both sides
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
                [fturns,bturns,amps,ints]=find_head_turns(Vl,hminamp,hmaxint);
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

%                 if s==1
                    for_turns{j}=fturns;
                    back_turns{j}=bturns;
                    head_avel{j}=Vl;
%                 else
%                     ips_turns{j}=rturns;
%                     cont_turns{j}=lturns;
%                     head_avel{j}=-da;
%                 end
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