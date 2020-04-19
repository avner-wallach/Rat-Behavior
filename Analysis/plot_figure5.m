function plot_figure5()
% figures for spatial compensation
global minamp;
minamp=2.5;
load('../Data/all_analyzed_structs_head.mat');
F=figure;
COL=colormap('lines');
COL2=[linspace(1,COL(2,1),64)',linspace(1,COL(2,2),64)',linspace(1,COL(2,3),64)'];
close(F);
fsize=22;

%% pump prob. and ratio by velocity
prob_rat_vel=1;
if(prob_rat_vel)
    vel=ctl_struct.all_pump_hspeed;
    pstr=ctl_struct.all_pump_stregth;
    vbins=[(1/3):(2/3):4 max(vel(:))];
    vinf=4.33;
    for i=1:numel(vbins)
        for j=1:2
            ind=find(vel(:,j)<vbins(i));
            pump_prob(i,j)=sum(pstr(ind,j)>0)/numel(ind);  
            [pval_pumps(i,j),n]=pump_bootstrap_oneside({hf_struct,ctl_struct},{[1:numel(hf_struct.all_pump_stregth(:,1))],ind},j);        
            N_pumps(i,j)=numel(ind);
        end
        ind=find(vel(:,1)<vbins(i) & vel(:,2)<vbins(i));    
        pump_rat(i)=sum(pstr(ind,1)>0)/sum(pstr(ind,2)>0);
    end
    hf_prob(1)=sum(hf_struct.all_pump_stregth(:,1)>0)/size(hf_struct.all_pump_stregth(:,1),1);
    hf_prob(2)=sum(hf_struct.all_pump_stregth(:,2)>0)/size(hf_struct.all_pump_stregth(:,2),1);
    figure;
    H=plot([0],hf_prob(1),'.');
    hold on;
    set(H,'Color',COL(1,:),'Marker','.','MarkerSize',40,'LineStyle','-','LineWidth',3);
    H=plot([0],hf_prob(2));
    set(H,'Color',COL(1,:),'Marker','o','MarkerSize',10,'LineStyle','-','LineWidth',3);
    H=plot(vbins(1:end-1),pump_prob(1:end-1,:));
    set(H(1),'Color',COL(2,:),'Marker','.','MarkerSize',40,'LineStyle','-','LineWidth',3);
    set(H(2),'Color',COL(2,:),'Marker','o','MarkerSize',10,'LineStyle','--','LineWidth',3);
    H=plot([vinf],pump_prob(end,1));
    set(H,'Color',COL(2,:),'Marker','.','MarkerSize',40,'LineStyle','-','LineWidth',3);
    H=plot([vinf],pump_prob(end,2));
    set(H,'Color',COL(2,:),'Marker','o','MarkerSize',10,'LineStyle','-','LineWidth',3);
    set(gca,'XLim',[-0.5 5],'YLim',[0 0.55],'FontSize',fsize);
    set(gca,'XTick',[0 1:4 vinf],'XTickLabel',{'0','1','2','3','4','\infty'})    

    figure;
    H=plot([0],hf_prob(1)/hf_prob(2),'*');
    hold on;
    set(H,'Color',COL(1,:),'Marker','*','MarkerSize',10,'LineStyle','-','LineWidth',3);
    H=plot(vbins(1:end-1),pump_prob(1:end-1,1)./pump_prob(1:end-1,2));
    set(H,'Color',COL(2,:),'Marker','*','MarkerSize',10,'LineStyle',':','LineWidth',3);
    H=plot([vinf],pump_prob(end,1)./pump_prob(end,2));
    set(H,'Color',COL(2,:),'Marker','*','MarkerSize',10,'LineStyle','-','LineWidth',3);
    set(gca,'XLim',[-0.5 5],'YLim',[.8 5],'FontSize',fsize);
    set(gca,'XTick',[0 1:4 vinf],'XTickLabel',{'0','1','2','3','4','\infty'})    
    set(gca,'FontSize',16);
    set(gcf,'Position',[680 558 150 150]);

    %phase duration
    dur=ctl_struct.dur_rat(:,1);
    rat=ctl_struct.dur_rat(:,2);
    P_R(:,1)=dur./(1+1./rat);
    P_R(:,2)=P_R(:,1)./rat;
    p=[0.25 0.5 0.75];
    rpts=linspace(0.5,2.5,1e2);
    pts=linspace(1,200,1e2)*1e-3;
    for i=1:numel(vbins)
        for j=1:2
            ind=find(vel(:,j)<vbins(i));
            q_ctl(i,:,j)=quantile(P_R(ind,j),p);
            p_ctl(i,:,j)=ksdensity(P_R(ind,j),pts);
        end
        ind=find(vel(:,1)<vbins(i) & vel(:,2)<vbins(i));
        q_ctl_rat(i,:)=quantile(rat(ind),p);
        p_ctl_rat(i,:)=ksdensity(rat(ind),rpts);
    end

    hf_dur=hf_struct.dur_rat(:,1);
    hf_rat=hf_struct.dur_rat(:,2);
    clear P_R;
    P_R(:,1)=hf_dur./(1+1./hf_rat);
    P_R(:,2)=P_R(:,1)./hf_rat;
    q_hf(1,:)=quantile(P_R(:,1),p);
    q_hf(2,:)=quantile(P_R(:,2),p);
    q_hf_rat=quantile(hf_rat,p);

    x_hf=[-0.1 0.1];
    x_inf=x_hf+vinf;

    figure;
    my_plotWithConfQ(vbins(1:end-1),q_ctl_rat(1:end-1,:)',COL(2,:));
    hold on;
    H=plot(vbins(1:end-1),q_ctl_rat(1:end-1,2)','.');
    set(H,'Color',COL(2,:),'Marker','*','MarkerSize',10,'LineStyle','none','LineWidth',3);
    my_plotWithConfQ(x_inf,([1;1]*q_ctl_rat(end,:))',COL(2,:));
    H=plot(vinf,q_ctl_rat(end,2)','.');
    set(H,'Color',COL(2,:),'Marker','*','MarkerSize',10,'LineStyle','none','LineWidth',3);
    my_plotWithConfQ(x_hf,([1;1]*q_hf_rat)',COL(1,:));
    H=plot(0,q_hf_rat(2)','.');
    set(H,'Color',COL(1,:),'Marker','*','MarkerSize',10,'LineStyle','none','LineWidth',3);
    set(gca,'XLim',[-0.5 5],'YLim',[.8 1.8],'FontSize',fsize);
    set(gca,'XTick',[0 1:4 vinf],'XTickLabel',{'0','1','2','3','4','\infty'},'XLabel',[],'YLabel',[])    

    figure;
    my_plotWithConfQ(vbins(1:end-1),q_ctl(1:end-1,:,1)'*1e3,COL(2,:));
    hold on;
    H=plot(vbins(1:end-1),q_ctl(1:end-1,2,1)'*1e3,'.');
    set(H,'Color',COL(2,:),'Marker','.','MarkerSize',40,'LineStyle','none','LineWidth',3);
    my_plotWithConfQ(vbins(1:end-1),q_ctl(1:end-1,:,2)'*1e3,COL(2,:)/2);
    H=plot(vbins(1:end-1),q_ctl(1:end-1,2,2)'*1e3,'.');
    set(H,'Color',COL(2,:)/2,'Marker','o','MarkerSize',10,'LineStyle','none','LineWidth',3);
    my_plotWithConfQ(x_inf,([1;1]*q_ctl(end,:,1))'*1e3,COL(2,:));
    H=plot(vinf,q_ctl(end,2,1)'*1e3,'.');
    set(H,'Color',COL(2,:),'Marker','.','MarkerSize',40,'LineStyle','none','LineWidth',3);
    my_plotWithConfQ(x_inf,([1;1]*q_ctl(end,:,2))'*1e3,COL(2,:)/2);
    H=plot(vinf,q_ctl(end,2,2)'*1e3,'.');
    set(H,'Color',COL(2,:)/2,'Marker','o','MarkerSize',10,'LineStyle','none','LineWidth',3);

    my_plotWithConfQ(x_hf,([1;1]*q_hf(1,:))'*1e3,COL(1,:));
    H=plot(0,q_hf(1,2)'*1e3,'.');
    set(H,'Color',COL(1,:),'Marker','.','MarkerSize',40,'LineStyle','none','LineWidth',3);
    my_plotWithConfQ(x_hf,([1;1]*q_hf(2,:))'*1e3,COL(1,:)/2);
    set(gca,'XLim',[-0.5 5],'YLim',[45 100],'FontSize',fsize);
    H=plot(0,q_hf(2,2)'*1e3,'.');
    set(H,'Color',COL(1,:)/2,'Marker','o','MarkerSize',10,'LineStyle','none','LineWidth',3);
    set(gca,'XTick',[0 1:4 vinf],'XTickLabel',{'0','1','2','3','4','\infty'},'XLabel',[],'YLabel',[])    
end
end

function [pvalue,N_pumps]=pump_bootstrap_oneside(strct,inds,side)
global minamp;
S=numel(strct);
if(nargin<2)
    for i=1:S
        inds{i}=1:numel(strct{i}.rise_fall(:,side));
    end
end
N=5e4;
for i=1:S
    ind=inds{i}(find(strct{i}.rise_fall(inds{i},side)>minamp));
    N_pumps(i)=numel(ind);
    pump_str{i}=strct{i}.all_pump_stregth(ind,side);    
end

pvalue=nan(S-1,S-1);
% for k=1:2
    for i=1:S-1
        for j=(i+1):S
            x=pump_str{i};
            y=pump_str{j};
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
            pvalue(i,j-1)=sum(abs(dd)>=abs(D))/N;
%             if(D>0)
%                 pvalue(i,j,k)=sum(dd>=D)/N;
%             else
%                 pvalue(i,j,k)=sum(dd<=D)/N;
%             end


        end
    end
end