%analyze behaving data
clear all;
% close all;
%% parameters
K=25;
glitch_th=15;

Fs=500;
d500=fdesign.lowpass('Fp,Fst,Ap,Ast',80,100,1,60,Fs);
Hd500 = design(d500);
% 
% Fs=1000;
% d1000=fdesign.lowpass('Fp,Fst,Ap,Ast',80,100,1,60,Fs);
% Hd1000 = design(d1000);


%% load data- take all left side C2 whiskers
% load('./COMP_DATA_STRUCTURE.mat');
% data=sComparative;
% colormatrix=colormap(lines);    
% 
% t=0;
% k=1;
% for i=1:length(data);
%     if(~strcmp(data(i).Strain,'Wis'))
%         continue;
%     end
%     ind_left=find(data(i).WhiskerSide==1);
%     for w=1:length(ind_left)
%         if(strcmp(data(i).WhiskerIdentity{ind_left(w)},'C2'))    
%             dat=data(i).Angle(:,ind_left(w));
%             dat=dat(~isnan(dat) & dat>0);
%             if(numel(dat)>500)
%                 wdata{k}=dat;
% 
%                 %remove glitches
%                 indg=find(abs(wdata{k}-smooth(wdata{k},K))>glitch_th);
%                 wdata{k}(indg)=(wdata{k}(indg-1)+wdata{k}(indg+1))/2;
% 
%                 %filter
%                 wdata{k}=filtfilt(Hd500.coeffs.Numerator,[1],wdata{k});
%                 session(k)=i;
%                 duration(k)=length(wdata{k})/500;            
%                 peak(k)=max(wdata{k});
%                 trough(k)=min(wdata{k});
%                 range(k)=peak(k)-trough(k);
% 
%                 H=plot(t+[1:length(wdata{k})]/500,wdata{k});
%                 set(H,'Color',colormatrix(mod(k,63)+1,:));
%                 text(t,100,num2str(i));
%                 t=t+duration(k);
%                 hold on;            
% 
%                 k=k+1;            
%             end
%         end
%     end
% end
% 
% save('comparative_c2_data','wdata','session','duration','peak','trough','range');
% 
%% create targe function
% 
% % level all trajectories
% for k=1:length(wdata)
%     wdata{k}=wdata{k}-trough(k)+mean(trough);
% end
% 
% for k=1:length(wdata)
%     for j=1:length(wdata)
%         if(j==k)
%             dist(k,j)=inf;
%         else
%             dist(k,j)=abs(wdata{k}(end)-wdata{j}(1));
%         end
%     end
% end
% 
% ind=find(range>70 | peak>115)
% dist(:,ind)=200;
% 
% for k=1:length(wdata)
%     order_mat(1,k)=k;
%     total_dist(k)=0;
%     for i=2:length(wdata)
%         vec=dist(order_mat(i-1,k),:);
%         vec(order_mat(1:i-1,k))=inf;
%         [M,j]=min(vec);
%         order_mat(i,k)=j;
%         total_dist(k)=total_dist(k)+M;
%     end
% end
%     
% [td,i]=nanmin(total_dist);
% order=order_mat(:,i);
%     
% t=0;
% figure;
% k=1;
% while(k<=length(order))
%     if(sum(order(k)==ind))        
%         order(k)=[];        
%         continue;
%     end
%     H=plot(t+[1:length(wdata{order(k)})]/500,wdata{order(k)});
%     set(H,'Color',colormatrix(mod(k,63)+1,:));
%     text(t,100,num2str(order(k)));
%     t=t+length(wdata{order(k)})/500;
%     hold on;
%     k=k+1;
% end
%     
% target=wdata{order(1)};%-min(wdata{order(1)});
% for k=2:length(order)
%     newd=wdata{order(k)};%-min(wdata{order(k)});
%     idx=find(abs(target-newd(1))<1);
%     target=[target(1:idx(end));newd];    
% end
% 
% target_slowed=interp(target,2);
% save('comparative_data','target','target_slowed');
% return;
%% compute slow vars
% filename='comparative_c2_data';
filename='headfixed_c2_data';
load(filename,'wdata','session','duration','peak','trough','range');
fs=500;
mode='spline';
maxlag=fs*5;
option='biased';
N_vector=zeros(1,2*maxlag+1);
all_amp=[];
all_offset=[];
all_whisk=[];
increasing=[];
decreasing=[];

for i=1:length(wdata)
    [amplitude,offset,phi,v,segmented_struct]=decompose_whisking(wdata{i},fs,mode);
    ind=find(~isnan(amplitude));
    
    %xcov
    [x_amp(i,:),lags]=xcov(amplitude(ind),maxlag,option);
    [x_offset(i,:),lags]=xcov(offset(ind),maxlag,option);
    [x_amp_offset(i,:),lags]=xcov(amplitude(ind),offset(ind),maxlag,option);    
    [x_doff(i,:),lags]=xcov(abs(diff(offset(ind))),maxlag,option);
    [x_amp_doffset(i,:),lags]=xcov(amplitude(ind(1:end-1)),abs(diff(offset(ind))),maxlag,option);
%     x_amp(i,:)=x_amp(i,:)*numel(ind);
    indout=find(abs(lags)>=numel(ind));
    indin=find(abs(lags)<numel(ind));
    x_amp(i,indout)=NaN;
    x_offset(i,indout)=NaN;
    x_amp_offset(i,indout)=NaN;
    x_doff(i,indout)=NaN;
    x_amp_doffset(i,indout)=NaN;
    N_vector(indin)=N_vector(indin)+numel(ind);
    
    %distributions
    all_amp=[all_amp amplitude(ind)];
    all_offset=[all_offset offset(ind)];
    
    %characteristic whisk
    Nsamp=100;
    for j=1:size(segmented_struct.whisks,1)
        ret=segmented_struct.whisks{j,2};
        prot=segmented_struct.whisks{j,1};        
        if(length(prot)<10 | length(ret)<10)
            continue;
        end
        ret_amp=ret(1)-ret(end);
        prot_amp=prot(end)-prot(1);
%         ret_amp=max(ret)-min(ret);
%         prot_amp=max(prot)-min(prot);
        if(prot_amp<5)
            continue;
        end
        scaled_ret=interp1(linspace(-1,0,length(ret)),ret,linspace(-1,0,Nsamp));
        scaled_prot=interp1(linspace(0,1,length(prot)),prot,linspace(0,1,Nsamp));
        whisk(j,:)=[scaled_prot scaled_ret];
        whisk(j,:)=(whisk(j,:)-min(prot))/(prot_amp);
        
        if(ret_amp<prot_amp)
            increasing=[increasing;whisk(j,:)];
        else
            decreasing=[decreasing;whisk(j,:)];  
        end
    end
    all_whisk=[all_whisk;whisk];
end
lags=lags/fs;
x_amp_m=nansum(x_amp)/size(x_amp,1);
sig_amp=sqrt(max(x_amp_m));

x_doff_m=nansum(x_doff)/size(x_doff,1);
sig_doff=sqrt(max(x_doff_m));

x_offset_m=nansum(x_offset)/size(x_offset,1);
sig_offset=sqrt(max(x_offset_m));

x_amp_offset_m=nansum(x_amp_offset)/size(x_amp_offset,1);
x_amp_doffset_m=nansum(x_amp_doffset)/size(x_amp_doffset,1);

figure;
plot(lags,x_amp_m/sig_amp^2);
hold on;
plot(lags,x_offset_m/sig_offset^2,'r');
plot(lags,x_amp_offset_m/(sig_amp*sig_offset),'k');

figure;
plot(lags,x_amp_m/sig_amp^2);
hold on;
plot(lags,x_doff_m/sig_doff^2,'r');
plot(lags,x_amp_doffset_m/(sig_amp*sig_doff),'k');

ctrs{1}=[0:50];
ctrs{2}=[20:120];
N=hist3([all_amp' all_offset'],ctrs);
figure, imagesc(ctrs{1},ctrs{2},N');


figure
sdec=std(decreasing)/sqrt(size(decreasing,1));
sinc=std(increasing)/sqrt(size(increasing,1));
my_plotWithConf(1:200,mean(decreasing),sdec,[1 0 0]);
hold on;
my_plotWithConf(1:200,mean(increasing),sinc,[0 0 1]);
d=mean(increasing)-mean(decreasing);
[m,i]=max(d(1:100));
x=increasing(:,i);
y=decreasing(:,i);
pvalue=boottest(x,y,1000,1e4);
text(i,0.75,['P=',num2str(pvalue)]);
% i=17;
% % target=wdata{i};
% [A,B,F,phi0,v]=decompose_whisking(target,500);
% t=[1:length(target)]/500;
% 
% figure;
% H=axes;
% plot(t,target);
% hold on;
% plot(t,v,':r');
% plot(t,B,'--k');
% plot(t,B+A,'k');
% plot(t,B-A,'-.k');
% set(H,'XLim',[t(1) t(end)]);
% legend('Experimental Data','Model','Bias','Max Protraction','Setpoint');
% ylabel('Protraction Angle (deg)');
% xlabel('Time (s)');
% 
% figure;
% H=subplot(3,1,1);
% plot(t,B);
% set(H,'XLim',[t(1) t(end)]);
% ylabel('Bias (deg)');
% H=subplot(3,1,2);
% plot(t,A);
% set(H,'XLim',[t(1) t(end)]);
% ylabel('Amplitude (deg)');
% subplot(3,1,3);
% plot(t,F);
% axis([t(1) t(end) 0 10]);
% ylabel('Inst. Frequency (Hz)');
% xlabel('Time (s)');
