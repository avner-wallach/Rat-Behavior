function plot_intrawhisker_compensation(hf_struct,comp_struct,ctl_struct)
% figures for spatial compensation
global minamp;
minamp=2;
F=figure;
COL=colormap('lines');
close(F);
fsize=22;
%% retraction/protraction durations
prot_ret=1;
if(prot_ret)
A=plotCDFComp({hf_struct.all_diluted_prot*1e3,ctl_struct.all_diluted_prot*1e3});
set(A(1),'XLim',[25 150]);
A=plotCDFComp({hf_struct.all_diluted_ret*1e3,ctl_struct.all_diluted_ret*1e3});
set(A(1),'XLim',[25 150]);
end

%% pumps prob.
pumps=1;
if(pumps)
ind_hf=find(hf_struct.rise_fall(:,1)>minamp);
ind_comp=find(comp_struct.rise_fall(:,1)>minamp);
ind_ctl=find(ctl_struct.rise_fall(:,1)>minamp);
[prot_percent]=plotPump({hf_struct.all_pump_stregth(ind_hf,1),ctl_struct.all_pump_stregth(ind_ctl,1)})
ind_hf=find(hf_struct.rise_fall(:,2)>minamp);
ind_comp=find(comp_struct.rise_fall(:,2)>minamp);
ind_ctl=find(ctl_struct.rise_fall(:,2)>minamp);
[ret_percent]=plotPump({hf_struct.all_pump_stregth(ind_hf,2),ctl_struct.all_pump_stregth(ind_ctl,2)})

[pval_pumps,N_pumps]=pump_bootstrap({hf_struct,ctl_struct})
end
%% jerk
jerk=0;
if(jerk)
ind_hf=find(hf_struct.rise_fall(:,1)>minamp);
ind_comp=find(comp_struct.rise_fall(:,1)>minamp);
ind_ctl=find(ctl_struct.rise_fall(:,1)>minamp);
[prot_percent]=plotPump({hf_struct.all_jerk_costs(ind_hf,1),comp_struct.all_jerk_costs(ind_comp,1),ctl_struct.all_jerk_costs(ind_ctl,1)});
ind_hf=find(hf_struct.rise_fall(:,2)>minamp);
ind_comp=find(comp_struct.rise_fall(:,2)>minamp);
ind_ctl=find(ctl_struct.rise_fall(:,2)>minamp);
[ret_percent]=plotPump({hf_struct.all_jerk_costs(ind_hf,2),comp_struct.all_jerk_costs(ind_comp,2),ctl_struct.all_jerk_costs(ind_ctl,2)});
end
%% pump clustering
pump_clust=0;
if(pump_clust)
% plot_pump_xcor(hf_struct.all_pump_stregth>0,hf_struct.segment_whisks,COL(1,:));
% plot_pump_xcor(comp_struct.all_pump_stregth>0,comp_struct.segment_whisks,COL(2,:));
% plot_pump_xcor(ctl_struct.all_pump_stregth>0,ctl_struct.segment_whisks,COL(3,:));

%all data
all_pumps=[hf_struct.all_pump_stregth>0;comp_struct.all_pump_stregth>0;ctl_struct.all_pump_stregth>0];
all_whisks=[hf_struct.segment_whisks comp_struct.segment_whisks ctl_struct.segment_whisks];
plot_pump_xcor(all_pumps,all_whisks,0*[1 1 1]);

end
%% pump correlations with duration and offset
pumpef=0;
if(pumpef)
struct.all_pump_stregth=[hf_struct.all_pump_stregth;comp_struct.all_pump_stregth;ctl_struct.all_pump_stregth];
struct.rise_fall=[hf_struct.rise_fall; comp_struct.rise_fall; ctl_struct.rise_fall];
struct.dur_rat=[hf_struct.dur_rat; comp_struct.dur_rat; ctl_struct.dur_rat];
struct.midpoints=[hf_struct.midpoints; comp_struct.midpoints; ctl_struct.midpoints];
struct.idx=cumsum([size(hf_struct.all_pump_stregth,1) size(comp_struct.all_pump_stregth,1) size(ctl_struct.all_pump_stregth,1)]);
% plot_pump_effects(hf_struct);
% plot_pump_effects(comp_struct);
% plot_pump_effects(ctl_struct);
plot_pump_effects(struct);

end
%% pump residency
pres=0;
if(pres)
struct.all_pump_stregth=[hf_struct.all_pump_stregth;comp_struct.all_pump_stregth;ctl_struct.all_pump_stregth];
struct.segment_whisks=[hf_struct.segment_whisks comp_struct.segment_whisks ctl_struct.segment_whisks];
plot_pump_runs(struct);
end
pcor=0;
if(pcor)
M=5;
m=1:M;
[hf_corrs,hf_pvals,hf_all_c,hf_all_p]=get_pump_smoothing(hf_struct.all_pump_stregth>0,hf_struct.segment_whisks,M);
[comp_corrs,comp_pvals,comp_all_c,comp_all_p]=get_pump_smoothing(comp_struct.all_pump_stregth>0,comp_struct.segment_whisks,M);
[ctl_corrs,ctl_pvals,ctl_all_c,ctl_all_p]=get_pump_smoothing(ctl_struct.all_pump_stregth>0,ctl_struct.segment_whisks,M);
F=figure;
A=axes;
plot_pump_cors(m,hf_all_c,hf_all_p,COL(1,:))
plot_pump_cors(m,comp_all_c,comp_all_p,COL(2,:))
plot_pump_cors(m,ctl_all_c,ctl_all_p,COL(3,:))
set(A,'FontSize',fsize,'Xlim',[0.5 M+0.5]);
end
end

function plot_pump_cors(m,all_c,all_p,C)
H=plot(m,all_c);
set(H,'LineWidth',2,'Color',C);
ind=find(all_p>=0.05);
H=plot(m(ind),all_c(ind));
set(H,'LineWidth',2,'Color',C,'Marker','o','MarkerSize',8);
hold on;
ind=find(all_p<0.05);
H=plot(m(ind),all_c(ind));
set(H,'LineWidth',2,'Color',C,'Marker','*','MarkerSize',10);
end

function [percent]=plotPump(pump_str)
fsize=22;
N=numel(pump_str);
COL=colormap('lines');
F=figure;
A=axes;
%area([0.001,1],[.99,.99],'FaceColor',0.9*[1 1 1],'LineStyle','none','BaseValue',0.001,'ShowBaseLine','off');
hold on;
for i=1:N
    H=cdfplot(pump_str{i});
    set(H,'Color',COL(i,:),'LineWidth',3);
    Y=get(H,'YData');
    Y(2)=nan;
    set(H,'YData',Y);
    nopump(i)=sum(pump_str{i}<=0)/numel(pump_str{i});
    delayedp(i)=sum(pump_str{i}>0 & pump_str{i}<=1)/numel(pump_str{i});
    doublep(i)=sum(pump_str{i}>1)/numel(pump_str{i});
    H=plot(0,nopump(i),'.');
    set(H,'Color',COL(i,:),'MarkerSize',24);
end

set(A,'Xlim',[0 2],'Ylim',[0 1],'XLabel',[],'Title',[],'Ylabel',[],'Xgrid','off','Ygrid','off','Box','on');
set(A,'FontSize',fsize);

percent=delayedp + doublep;

%inset
F=figure;
A=axes;
for i=1:N
    B=bar([i],[100*(percent(i))]);
    B.FaceColor=COL(i,:);
    hold on;
end
% for i=1:N
%     B=bar([i+N],[100*delayedp(i)]);
%     B.FaceColor=COL(i,:);
%     hold on;
% end
% for i=1:N
%     B=bar([i+2*N],[100*doublep(i)]);
%     B.FaceColor=COL(i,:);
%     hold on;
% end
set(A,'FontSize',16);
set(F,'Position',[680 558 150 150]);
set(A,'XTick',[]);
set(A,'Ylim',[0 max(100*(percent))*1.25],'Xlim',[0 N+1]);
end

function [pvalue,N_pumps]=pump_bootstrap(strct)
global minamp;
N=5e4;
for i=1:numel(strct)
    ind=find(strct{i}.rise_fall(:,1)>minamp);
    N_pumps(i,1)=numel(ind);
    pump_str{i,1}=strct{i}.all_pump_stregth(ind,1);    
    ind=find(strct{i}.rise_fall(:,2)>minamp);
    N_pumps(i,2)=numel(ind);
    pump_str{i,2}=strct{i}.all_pump_stregth(ind,2);    
end

pvalue=nan(2,3,2);
for k=1:2
    for i=1:numel(strct)-1
        for j=(i+1):numel(strct)
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


function A=plotCDFComp(data)
COL=colormap('lines');
fsize=22;
for i=1:numel(data)
    med(i)=median(data{i});
    men(i)=mean(data{i})
    vr(i)=var(data{i});
end
for i=1:numel(data)-1
    for j=i+1:numel(data)        
        pval_median(i,j)=boottest(data{i},data{j},5e3,@median,'double');
        pval_mean(i,j)=boottest(data{i},data{j},5e3,@mean,'double');
        pval_var(i,j)=boottest(data{i},data{j},5e3,@var,'double');        
    end
end
pval_median;
pval_mean;
pval_var;

F=figure;
A(1)=axes;
wn=1;
for i=1:numel(data)
    H=cdfplot(data{i}+randn(size(data{i}))*wn);
    set(H,'Color',COL(i,:),'LineWidth',3);
    hold on;
end
set(A(1),'FontSize',fsize);
set(A(1),'XGrid','off');
set(A(1),'YGrid','off');
set(A(1),'Title',[],'XLabel',[],'YLabel',[]);

%insets
F=figure;
A(2)=axes;
for i=1:numel(data)
    B=bar([i],men(i));
    B.FaceColor=COL(i,:);
    hold on;
end
set(A(2),'FontSize',16);
set(F,'Position',[680 558 150 150]);
set(A(2),'XTick',[],'Xlim',[0 numel(data)+1]);
set(A(2),'YLim',[0 max(men)*1.25]);

% F=figure;
% A(4)=axes;
% for i=1:numel(data)
%     B=bar([i],med(i));
%     B.FaceColor=COL(i,:);
%     hold on;
% end
% set(A(4),'FontSize',10);
% set(F,'Position',[680 558 150 150]);
% set(A(4),'XTick',[],'Xlim',[0 4]);
% set(A(4),'YLim',[25 max(med)*1.25]);
end

% function [matrices,ctl_matrices]=get_pump_martices(S)
% pumps=S.all_pump_stregth>0;
% matrices=zeros(2,2,2,2);
% k=1;
% for i=1:numel(S.segment_whisks)
%     p=pumps(k:(k+S.segment_whisks(i)-1),:);
%     p1=p(1:end-1,:);
%     p2=p(2:end,:);
%     for a=1:2
%         for b=1:2
%             m(1,1,a,b)=sum(~p1(:,a) & ~p2(:,b));
%             m(1,2,a,b)=sum(~p1(:,a) & p2(:,b));
%             m(2,1,a,b)=sum(p1(:,a) & ~p2(:,b));
%             m(2,2,a,b)=sum(p1(:,a) & p2(:,b));
%         end
%     end
%     matrices=matrices+m;
%     k=k+S.segment_whisks(i);
% end
% p=mean(pumps);
% for a=1:2
%     for b=1:2
%         matrices(:,:,a,b)=matrices(:,:,a,b)/sum(sum(matrices(:,:,a,b)));
%         ctl_matrices(:,:,a,b)=[1-p(a);p(a)]*[1-p(b) p(b)];
%     end
% end
% 
% 
% end

function [pump_xcor]=get_pump_xcor(pumps,segment_whisks,K);
pump_xcor=zeros(4,4*K); %1st row: prot pump, 2nd row: prot nupump, 3rd & 4th:ret pump/nopump; columns: number of cycles*2 phases*past and future 
counters=zeros(4,4*K);
k=1;
for i=1:numel(segment_whisks)
    p=pumps(k:(k+segment_whisks(i)-1),:)';
%     t=ctl(k:(k+S.segment_whisks(i)-1),:)';
    P=[zeros(1,2*K) p(:)' zeros(1,2*K)];
%     T=[zeros(1,2*K) t(:)' zeros(1,2*K)];
    C=[zeros(1,2*K) ones(size(p(:)')) zeros(1,2*K)];
    for j=1:numel(p(:))
        r=(~mod(j,2)*2+~p(j)+1); %row number
        ind=[(j-2*K):(j-1) (j+1):(j+2*K)]+2*K;
        pump_xcor(r,:)=pump_xcor(r,:)+P(ind);
%         ctl_xcor(r,:)=ctl_xcor(r,:)+T(ind);
        counters(r,:)=counters(r,:)+C(ind);
    end
    k=k+segment_whisks(i);
end                 
            
pump_xcor=pump_xcor./counters;

end

function plot_pump_xcor(pumps,segments_whisks,COL)
K=10;
fsize=22;
[pump_xcor]=get_pump_xcor(pumps,segments_whisks,K);
N=2000;
for n=1:N
    ctl=[pumps(randperm(size(pumps,1)),1) pumps(randperm(size(pumps,1)),2)];
    ctl_xcor(:,:,n)=get_pump_xcor(ctl,segments_whisks,K);
end
m_ctl_xcor=mean(ctl_xcor,3);
p=[0.01 0.05 0.1 0.5 0.9 0.95 0.99];
q_ctl_xcor=quantile(ctl_xcor,p,3);
pumpx_norm=(pump_xcor-m_ctl_xcor);
ind=find(pumpx_norm<0);
q=q_ctl_xcor(:,:,2)-m_ctl_xcor;
pumpx_norm(ind)=pumpx_norm(ind)./abs(q(ind));
ind=find(pumpx_norm>0);
q=q_ctl_xcor(:,:,6)-m_ctl_xcor;
pumpx_norm(ind)=pumpx_norm(ind)./abs(q(ind));
% K=10;
x=[-2*K:-1 1:2*K];
xx=x/2;
ind_even=find(~mod(x,2));
ind_odd=find(mod(x,2));
ind_even_pos=ind_even(x(ind_even)>0);
ind_even_neg=ind_even(x(ind_even)<0);
ind_odd_pos=ind_odd(x(ind_odd)>0);
ind_odd_neg=ind_odd(x(ind_odd)<0);
F=figure;
A=axes;
%ywindow=[min(min(pumpx_norm(1:2,:)))*1.2 max(max(pumpx_norm(1:2,:)))*1.2];
ywindow=[-6 8];
fa=0.05;
for i=1:numel(ind_odd)
    X=[xx(ind_odd(i))-0.25 xx(ind_odd(i))+0.25];    
    H=area(X,[ywindow(2) ywindow(2)],ywindow(1)-0.5);
    set(H,'FaceColor',[0 0 0],'FaceAlpha',fa,'LineStyle','none');
    hold on;
end

% H=plot(xx(ind_even),pumpx_norm(1:2,ind_even)');
on_p_marker='.'; on_p_size=40; on_p_style='-';
on_r_marker='o'; on_r_size=10; on_r_style=':';
off_p_marker='none'; off_p_size=32; off_p_style='none';
off_r_marker='none'; off_r_size=8; off_r_style='none';
lw=3;
%same phase
H=plot(xx(ind_even_neg),pumpx_norm(1:2,ind_even_neg)');
set(H(1),'Color',COL,'Marker',on_p_marker,'MarkerSize',on_p_size,'LineWidth',lw,'LineStyle',on_p_style);
set(H(2),'Color',COL,'Marker',off_p_marker,'MarkerSize',off_p_size,'LineWidth',lw,'LineStyle',off_p_style);
hold on;
H=plot(xx(ind_even_pos),pumpx_norm(1:2,ind_even_pos)');
set(H(1),'Color',COL,'Marker',on_p_marker,'MarkerSize',on_p_size,'LineWidth',lw,'LineStyle',on_p_style);
set(H(2),'Color',COL,'Marker',off_p_marker,'MarkerSize',off_p_size,'LineWidth',lw,'LineStyle',off_p_style);

%opposite phase
H=plot(xx(ind_odd_neg),pumpx_norm(1:2,ind_odd_neg)');
set(H(1),'Color',COL,'Marker',on_r_marker,'MarkerSize',on_r_size,'LineWidth',lw,'LineStyle',on_r_style);
set(H(2),'Color',COL,'Marker',off_r_marker,'MarkerSize',off_r_size,'LineWidth',lw,'LineStyle',off_r_style);
H=plot(xx(ind_odd_pos),pumpx_norm(1:2,ind_odd_pos)');
set(H(1),'Color',COL,'Marker',on_r_marker,'MarkerSize',on_r_size,'LineWidth',lw,'LineStyle',on_r_style);
set(H(2),'Color',COL,'Marker',off_r_marker,'MarkerSize',off_r_size,'LineWidth',lw,'LineStyle',off_r_style);
X=[-0.25 +0.25];    
%H=area(X,[ywindow(2) ywindow(2)]-0.1,ywindow(1)+0.1);
%set(H,'FaceColor',[1 1 1],'FaceAlpha',1,'LineStyle','none');
H=plot(xx,ones(size(x)));
set(H,'LineWidth',3,'Color',0.5*[1 1 1]);
H=plot(xx,-ones(size(x)));
set(H,'LineWidth',3,'Color',0.5*[1 1 1]);
xmax=max(xx)+.25; %6.25
set(A,'Xlim',[-xmax xmax],'Ylim',ywindow);
set(A,'FontSize',fsize);

F=figure;
A=axes;
% ywindow=[min(min(pumpx_norm(3:4,:)))*1.2 max(max(pumpx_norm(3:4,:)))*1.2];
ywindow=[-6 11];
for i=1:numel(ind_even)
    X=[xx(ind_even(i))-0.25 xx(ind_even(i))+0.25];    
    H=area(X,[ywindow(2) ywindow(2)],ywindow(1)-0.5);
    set(H,'FaceColor',[0 0 0],'FaceAlpha',fa,'LineStyle','none');
    hold on;
end
X=[-0.25 +0.25];    
H=area(X,[ywindow(2) ywindow(2)],ywindow(1)-0.5);
set(H,'FaceColor',[0 0 0],'FaceAlpha',fa,'LineStyle','none');

%same phase
H=plot(xx(ind_even_neg),pumpx_norm(3:4,ind_even_neg)');
set(H(1),'Color',COL,'Marker',on_r_marker,'MarkerSize',on_r_size,'LineWidth',lw,'LineStyle',on_r_style);
set(H(2),'Color',COL,'Marker',off_r_marker,'MarkerSize',off_r_size,'LineWidth',lw,'LineStyle',off_r_style);
hold on;
H=plot(xx(ind_even_pos),pumpx_norm(3:4,ind_even_pos)');
set(H(1),'Color',COL,'Marker',on_r_marker,'MarkerSize',on_r_size,'LineWidth',lw,'LineStyle',on_r_style);
set(H(2),'Color',COL,'Marker',off_r_marker,'MarkerSize',off_r_size,'LineWidth',lw,'LineStyle',off_r_style);

%opposite phase
H=plot(xx(ind_odd_neg),pumpx_norm(3:4,ind_odd_neg)');
set(H(1),'Color',COL,'Marker',on_p_marker,'MarkerSize',on_p_size,'LineWidth',lw,'LineStyle',on_p_style);
set(H(2),'Color',COL,'Marker',off_p_marker,'MarkerSize',off_p_size,'LineWidth',lw,'LineStyle',off_p_style);
H=plot(xx(ind_odd_pos),pumpx_norm(3:4,ind_odd_pos)');
set(H(1),'Color',COL,'Marker',on_p_marker,'MarkerSize',on_p_size,'LineWidth',lw,'LineStyle',on_p_style);
set(H(2),'Color',COL,'Marker',off_p_marker,'MarkerSize',off_p_size,'LineWidth',lw,'LineStyle',off_p_style);
X=[-0.25 +0.25];    
%H=area(X,[ywindow(2) ywindow(2)]-0.1,ywindow(1)+0.1);
%set(H,'FaceColor',[1 1 1],'FaceAlpha',1,'LineStyle','none');
H=plot(xx,ones(size(x)));
set(H,'LineWidth',3,'Color',0.5*[1 1 1]);
H=plot(xx,-ones(size(x)));
set(H,'LineWidth',3,'Color',0.5*[1 1 1]);
set(A,'Xlim',[-xmax xmax],'Ylim',ywindow);
set(A,'FontSize',fsize);
end

function plot_pump_effects(struct)
COL=colormap('lines');
fsize=22;
pumps=struct.all_pump_stregth>0;
ret=struct.dur_rat(:,1)./(1+struct.dur_rat(:,2));
prot_ret=[struct.dur_rat(:,2).*ret ret]*1e3;
prot_ret=prot_ret+1*randn(size(prot_ret));
rise_fall=struct.rise_fall;
midpoints=struct.midpoints;
for i=1:2; 
    indp=find(pumps(:,i));
    indnp=find(~pumps(:,i));
    F=figure;
    A=axes;
    H=cdfplot(prot_ret(indnp,i));
    mnp=nanmean(prot_ret(indnp,i));
    set(H,'LineWidth',3,'LineStyle','--','Color',0*[1 1 1]);
    hold on;
    H=cdfplot(prot_ret(indp,i));
    mp=nanmean(prot_ret(indp,i));
    set(H,'LineWidth',3,'LineStyle','-','Color',0*[1 1 1]);
    set(A,'XLabel',[],'YLabel',[],'FontSize',fsize,'XLim',[40 160],'Title',[]);
    set(A,'XGrid','off','YGrid','off');
    display(['durtion: i=',num2str(i),'; mnp=',num2str(mnp),'; mp=',num2str(mp)]);
    %per database
%     idx1=indnp(indnp<=struct.idx(1));
%     H=cdfplot(prot_ret(idx1,i));
%     set(H,'LineWidth',2,'LineStyle',':','Color',COL(1,:));
%     idx2=indnp(indnp>struct.idx(1) & indnp<=struct.idx(2));
%     H=cdfplot(prot_ret(idx2,i));
%     set(H,'LineWidth',2,'LineStyle',':','Color',COL(2,:));
%     idx3=indnp(indnp>struct.idx(2) & indnp<=struct.idx(3));
%     H=cdfplot(prot_ret(idx3,i));
%     set(H,'LineWidth',2,'LineStyle',':','Color',COL(3,:));
%     F2=figure;
%     my_boxplot({prot_ret(idx1,i)',prot_ret(idx2,i)',prot_ret(idx3,i)'},{'hf','comp','ctl'},1e3,@mean,'double')
% 
%     figure(F);
%     idx=indp(indp<=struct.idx(1));
%     H=cdfplot(prot_ret(idx,i));
%     set(H,'LineWidth',2,'LineStyle','-','Color',COL(1,:));
%     idx=indp(indp>struct.idx(1) & indp<=struct.idx(2));
%     H=cdfplot(prot_ret(idx,i));
%     set(H,'LineWidth',2,'LineStyle','-','Color',COL(2,:));
%     idx=indp(indp>struct.idx(2) & indp<=struct.idx(3));
%     H=cdfplot(prot_ret(idx,i));
%     set(H,'LineWidth',2,'LineStyle','-','Color',COL(3,:));

%     figure;
%     axes;
%     cdfplot(rise_fall(indnp,i));
%     hold on;
%     cdfplot(rise_fall(indp,i));
    F=figure;
    A=axes;
    H=cdfplot(midpoints(indnp,i));
    mnp=nanmean(midpoints(indnp,i));
    set(H,'LineWidth',3,'LineStyle','--','Color',0*[1 1 1]);    
    hold on;
    H=cdfplot(midpoints(indp,i));
    mp=nanmean(midpoints(indp,i));
    set(H,'LineWidth',3,'LineStyle','-','Color',0*[1 1 1]);
    set(A,'XLabel',[],'YLabel',[],'FontSize',fsize,'XLim',[40 120],'Title',[]);    
    set(A,'XGrid','off','YGrid','off');
    display(['offset: i=',num2str(i),'; mnp=',num2str(mnp),'; mp=',num2str(mp)]);
end
end

function plot_pump_runs(struct)
fsize=22;
pumps=struct.all_pump_stregth>0;
segment_whisks=struct.segment_whisks;
M=15;
[residency,f,cond_prob]=get_pump_runs(pumps,segment_whisks,M);
N=500;
for n=1:N
    ctl=[pumps(randperm(size(pumps,1)),1) pumps(randperm(size(pumps,1)),2)];
    [residency_ctl,f_ctl,cond_prob_ctl]=get_pump_runs(ctl,segment_whisks,M);
    for i1=1:2
        for i2=1:2
            a(n,i1,i2)=f_ctl{i1,i2}.a;
            b(n,i1,i2)=f_ctl{i1,i2}.b;
        end
    end
end
x=[1:M];
for i1=1:2
    for i2=1:2
        [B,I]=sort(b(:,i1,i2));
        b_med(i1,i2)=B(N/2);
        a_med(i1,i2)=a(I(N/2),i1,i2);
        b_lq(i1,i2)=B(0.05*N);
        a_lq(i1,i2)=a(I(0.05*N),i1,i2);
        b_hq(i1,i2)=B(0.95*N);
        a_hq(i1,i2)=a(I(0.95*N),i1,i2);
        
        y_med(:,i1,i2)=a_med(i1,i2)*exp(b_med(i1,i2)*x);
        y_lq(:,i1,i2)=a_lq(i1,i2)*exp(b_lq(i1,i2)*x);
        y_hq(:,i1,i2)=a_hq(i1,i2)*exp(b_hq(i1,i2)*x);
    end
end

mrk={'o','square'};
edgs=[0 0 0];
fcs=[0 0 0;1 1 1];
szs=10;
for r=1:2
    F=figure;
    A=axes;
    for i1=1:2
        Q=sort([y_lq(:,i1,r)'; y_med(:,i1,r)'; y_hq(:,i1,r)']);
        my_plotWithConfQ(x,Q,0.5*[1 1 1]);
        hold on;
        H=plot(f{i1,r});
        set(H,'LineWidth',2,'Color','k');
        pval_res(i1,r)=sum(b(:,i1,r)>=f{i1,r}.b)/N;
    end
    H=plot(x,residency(:,:,r),'.');
    set(H(1),'Marker',mrk{2},'MarkerEdgeColor',edgs,'MarkerFaceColor',fcs(r,:),'MarkerSize',szs);
    set(H(2),'Marker',mrk{1},'MarkerEdgeColor',edgs,'MarkerFaceColor',fcs(r,:),'MarkerSize',szs);
    R=residency(:,:,r);
    set(A,'YScale','log','XLabel',[],'YLabel',[],'FontSize',fsize,'XLim',[1 M],'YLim',[min(R(R>0))/10 1]);
    set(A,'YTick',[1e-4 1e-3 1e-2 1e-1 1e0]);
    legend('off');
    %inset
%     F=figure;
%     A=axes;
%     X=[b(:,1,r);b(:,2,r)];
%     G=[zeros(N,1);ones(N,1)];
%     boxplot(X,G,'plotstyle','compact','symbol','','colorgroup',G,'median','target',...
%     'colors',0.5*[1 1 1],'Widths',.5,'FactorGap',[],'FactorSeparator',[],...
%     'orientation','vertical','LabelOrientation','inline');
%     hold on;
%     H=plot([f{1,r}.b f{2,r}.b],'.');
%     set(H,'Marker','.','Color',0*[1 1 1],'MarkerSize',32);
%     
end
pval_res
end

function [residency,f,cond_prob]=get_pump_runs(pumps,segment_whisks,M)
cond_prob=zeros(M,2,2);
counters=zeros(M,2,2);
residency=zeros(M,2,2);

for r=1:2
    k=1;
    for i=1:numel(segment_whisks)
        p=pumps(k:(k+segment_whisks(i)-1),:)';
        for m=1:M
            onev=[0 ones(1,m) 0];
            zerov=[1 zeros(1,m) 1];
            oneidx=strfind(p(r,:),onev);
            zeroidx=strfind(p(r,:),zerov);
            residency(m,1,r)=residency(m,1,r)+numel(zeroidx);
            residency(m,2,r)=residency(m,2,r)+numel(oneidx);
            for n=1:numel(zeroidx) %run of no pumps
                cond_prob(m,1,r)=cond_prob(m,1,r)+p(~(r-1)+1,zeroidx(n)+m);
                counters(m,1,r)=counters(m,1,r)+1;
            end
            for n=1:numel(oneidx) %run of no pumps
                cond_prob(m,2,r)=cond_prob(m,2,r)+p(~(r-1)+1,oneidx(n)+m);
                counters(m,2,r)=counters(m,2,r)+1;
            end
        end
        k=k+segment_whisks(i);
    end                 
    %normalize
    residency(:,1,r)=residency(:,1,r)/sum(residency(:,1,r));
    residency(:,2,r)=residency(:,2,r)/sum(residency(:,2,r));
    x=[1:M]';
    ind=find(residency(:,1,r)>0);
    f{1,r}=fit(x(ind),log(residency(ind,1,r)),'poly1');
    ind=find(residency(:,2,r)>0);
    f{2,r}=fit(x(ind),log(residency(ind,2,r)),'poly1');
end
cond_prob=cond_prob./counters;
f=conv_fits(f);
end

function fout=conv_fits(fin)
ft=fittype('a*exp(b*x)');
fout=fin;
for i=1:numel(fin)
    fout{i}=cfit(ft,exp(fin{i}.p2),fin{i}.p1);
end
end

function [corrs,pvals,all_c,all_p]=get_pump_smoothing(pumps,segment_whisks,M)
smooth_pump=[];
corrs=[];
pvals=[];
k=1;
for i=1:numel(segment_whisks)
    p=double(pumps(k:(k+segment_whisks(i)-1),:));
    sp=[];
    for m=1:M
        a=1;
        b=1/m*ones(1,m);
        prot=filter(b,a,p(:,1));
        ret=filter(b,a,p(:,2));
        sp(:,:,m)=[prot ret];
        [c(m),pv(m)]=corr(prot,ret);           
    end
    smooth_pump=[smooth_pump;sp];
    corrs=[corrs;c];
    pvals=[pvals;pv];
    k=k+segment_whisks(i);
end
for m=1:M
    [all_c(m),all_p(m)]=corr(smooth_pump(:,1,m),smooth_pump(:,2,m)); 
end
end
