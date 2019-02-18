close all;
%clear all;
%%
filename='comparative_c2_data';
% filename='control_c2_data';
% filename='headfixed_c2_data';
load(filename,'wdata','session','duration','peak','trough','range');
fs=500;
mode='spline';
% mode='zoh';
% mode='linear';
maxlag=fs*3;
option='biased';
N_vector=zeros(1,2*maxlag+1);
all_amp=[];
all_offset=[];
all_dur=[];
all_rat=[];
all_whisk=[];

amp=[];
off=[];
dur=[];
rise=[];
pr_ratio=[];
dur_1r=[];
dur_2r=[];

increasing=[];
model_increasing=[];
decreasing=[];
model_decreasing=[];
short=[];
model_short=[];
long=[];
model_long=[];
err=[];

%%
for i=1:length(wdata)
    [decomposed_data]=decompose_whisking(wdata{i},fs,mode);
    amplitude=decomposed_data.amplitude;
    offset=decomposed_data.offset;
    phi=decomposed_data.phi;
    v=decomposed_data.v;
    durations=decomposed_data.durations;   
    ratios=decomposed_data.ratios;   
    segmented_struct=decomposed_data.segmented_struct;
    amp=[amp decomposed_data.d_amp];
    off=[off decomposed_data.d_offset];
    pr_ratio=[pr_ratio segmented_struct.pr_ratio];
    dur=[dur segmented_struct.duration];
    rise=[rise segmented_struct.rise];
    [segmented_model]=segment_whisking(v(~isnan(v)),phi(~isnan(v)),pi/4,fs);
    ind=find(~isnan(amplitude));
    
    %model fitting error
    err=[err wdata{i}(ind)'-v(ind)];
    
    %first returns
    d1=[dur(1:end-1);dur(2:end)];
    d2=[dur(1:end-2);dur(3:end)];
    dur_1r=[dur_1r d1];
    dur_2r=[dur_2r d2];

    %xcov
    [x_amp(i,:),lags]=xcov(amplitude(ind),maxlag,option);
    [x_offset(i,:),lags]=xcov(offset(ind),maxlag,option);
    [x_amp_offset(i,:),lags]=xcov(amplitude(ind),offset(ind),maxlag,option);    
    [x_doff(i,:),lags]=xcov(abs(diff(offset(ind))),maxlag,option);
    [x_amp_doffset(i,:),lags]=xcov(amplitude(ind(1:end-1)),abs(diff(offset(ind))),maxlag,option);
    [x_duration(i,:),lags]=xcov(durations(ind),maxlag,option);
    [x_ratio(i,:),lags]=xcov(ratios(ind),maxlag,option);
%     x_amp(i,:)=x_amp(i,:)*numel(ind);
    indout=find(abs(lags)>=numel(ind));
    indin=find(abs(lags)<numel(ind));
    x_amp(i,indout)=NaN;
    x_offset(i,indout)=NaN;
    x_amp_offset(i,indout)=NaN;
    x_doff(i,indout)=NaN;
    x_amp_doffset(i,indout)=NaN;
    x_duration(i,indout)=NaN;
    x_ratio(i,indout)=NaN;
    x_amp(i,:)=x_amp(i,:)*duration(i);
    x_offset(i,:)=x_offset(i,:)*duration(i);
    x_amp_offset(i,:)=x_amp_offset(i,:)*duration(i);
    x_doff(i,:)=x_doff(i,:)*duration(i);
    x_amp_doffset(i,:)=x_amp_doffset(i,:)*duration(i);
    x_duration(i,:)=x_duration(i,:)*duration(i);
    x_ratio(i,:)=x_ratio(i,:)*duration(i);

%     N_vector(indin)=N_vector(indin)+numel(ind);
    
    %distributions
    all_amp=[all_amp amplitude(ind)];
    all_offset=[all_offset offset(ind)];
    all_dur=[all_dur durations(ind)];
    all_rat=[all_rat ratios(ind)];
    
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
        all_whisk=[all_whisk;whisk];

        if(length(ret)<length(prot)/1.6)
            ret=[ret;nan(round(length(prot)/1.6)-length(ret),1)];
            scaled_ret=interp1(linspace(-1,0,length(ret)),ret,linspace(-1,0,Nsamp));        
            w=[scaled_prot scaled_ret];
            w=(w-min(prot))/(prot_amp);
            short=[short ; w];
        else
            ret=ret(1:round(length(prot)/1.6));
            scaled_ret=interp1(linspace(-1,0,length(ret)),ret,linspace(-1,0,Nsamp));        
            w=[scaled_prot scaled_ret];
            w=(w-min(prot))/(prot_amp);
            long=[long;w];  
        end        

    end  
    
    for j=1:size(segmented_model.whisks,1)
        m_ret=segmented_model.whisks{j,2};
        m_prot=segmented_model.whisks{j,1};                
        if(length(m_prot)<10 | length(m_ret)<10)
            continue;
        end
        m_ret_amp=m_ret(1)-m_ret(end);
        m_prot_amp=m_prot(end)-m_prot(1);
        if(m_prot_amp<5)
            continue;
        end        
        scaled_m_ret=interp1(linspace(-1,0,length(m_ret)),m_ret,linspace(-1,0,Nsamp));
        scaled_m_prot=interp1(linspace(0,1,length(m_prot)),m_prot,linspace(0,1,Nsamp));
        m_whisk=[scaled_m_prot scaled_m_ret];
        m_whisk=(m_whisk-min(m_prot))/(m_prot_amp);
        if(m_ret_amp<m_prot_amp)
            model_increasing=[model_increasing;m_whisk];
        else
            model_decreasing=[model_decreasing;m_whisk];
        end
        
        pr_factor=1.6;

        if(length(m_ret)<length(m_prot)/pr_factor)
            m_ret=[m_ret nan(1,round(length(m_prot)/pr_factor)-length(m_ret))];
            scaled_m_ret=interp1(linspace(-1,0,length(m_ret)),m_ret,linspace(-1,0,Nsamp));        
            w=[scaled_m_prot scaled_m_ret];
            w=(w-min(m_prot))/(m_prot_amp);
            model_short=[model_short ; w];
        else
            m_ret=m_ret(1:round(length(m_prot)/pr_factor));
            scaled_m_ret=interp1(linspace(-1,0,length(m_ret)),m_ret,linspace(-1,0,Nsamp));        
            w=[scaled_m_prot scaled_m_ret];
            w=(w-min(m_prot))/(m_prot_amp);
            model_long=[model_long;w];  
        end        
        
    end
end

std(err)

lags=lags/fs;
x_amp_m=nansum(x_amp)/size(x_amp,1);
sig_amp=sqrt(max(x_amp_m));

x_doff_m=nansum(x_doff)/size(x_doff,1);
sig_doff=sqrt(max(x_doff_m));

x_offset_m=nansum(x_offset)/size(x_offset,1);
sig_offset=sqrt(max(x_offset_m));

x_duration_m=nansum(x_duration)/size(x_duration,1);
sig_duration=sqrt(max(x_duration_m));

x_ratio_m=nansum(x_ratio)/size(x_ratio,1);
sig_ratio=sqrt(max(x_ratio_m));

x_amp_offset_m=nansum(x_amp_offset)/size(x_amp_offset,1);
x_amp_doffset_m=nansum(x_amp_doffset)/size(x_amp_doffset,1);
%% temporal correlations

figure;
plot(lags,x_amp_m/sig_amp^2);
hold on;
plot(lags,x_offset_m/sig_offset^2,'r');
plot(lags,x_duration_m/sig_duration^2,'m');
plot(lags,x_ratio_m/sig_ratio^2,'g');
plot(lags,x_amp_offset_m/(sig_amp*sig_offset),'k');
legend({'amplitude','offset','duration','ratio','amp-offset'});

% figure;
% plot(lags,x_amp_m/sig_amp^2);
% hold on;
% plot(lags,x_doff_m/sig_doff^2,'r');
% plot(lags,x_amp_doffset_m/(sig_amp*sig_doff),'k');

%% dustributions

%bin=3;
clear ctrs;
ctrs{1}=linspace(0,50,50);
ctrs{2}=linspace(20,120,50);
Nx=hist3([amp' off'],ctrs);
Namp=hist(amp',ctrs{1});
Namp=Namp/sum(Namp);
Noff=hist(off',ctrs{2});
Noff=Noff/sum(Noff);
figure;
H=axes('position',[0.35 0.35 0.6 0.6]);
imagesc(ctrs{1},ctrs{2},Nx');
set(H,'YDir','normal');
set(H,'YTick',[],'YTickLabel',[]);
set(H,'XTick',[],'XTickLabel',[]);
H=axes('position',[0.1 0.35 0.2 0.6]);
barh(ctrs{2},Noff);
set(H,'XDir','reverse');
set(H,'YLim',[20,120]);
H=ylabel('Offset (deg)');
set(H,'FontSize',12);
H=axes('position',[0.35 0.1 0.6 0.2]);
bar(ctrs{1},Namp);
set(H,'XLim',[-1.5 50]);
H=xlabel('Amplitude (deg)');
set(H,'FontSize',12);

minamp=4;
ind=find(rise>minamp);
ctrs=linspace(3,18,50);
Nd=hist(1./dur(ind),ctrs);
Nd=Nd/sum(Nd);
figure;
subplot(2,1,1);
bar(ctrs,Nd);
xlabel('Frequency (Hz)');
ctrs=linspace(0,5,50);
Nr=hist(pr_ratio(ind),ctrs);
Nr=Nr/sum(Nr);
subplot(2,1,2);
bar(ctrs,Nr);
xlabel('Protraction-Retraction Ratio');

% ctrs{1}=[0:1:35];
% ctrs{2}
% N=hist3([all_amp' 1./all_dur'],ctrs);
% figure, imagesc(ctrs{1},ctrs{2},N');


%% cross-whisk interaction
T=linspace(-1,1,200);

sdec=std(decreasing)/sqrt(size(decreasing,1));
sinc=std(increasing)/sqrt(size(increasing,1));
mdec=mean(decreasing);
minc=mean(increasing);

smdec=std(model_decreasing)/sqrt(size(model_decreasing,1));
sminc=std(model_increasing)/sqrt(size(model_increasing,1));
mmdec=mean(model_decreasing);
mminc=mean(model_increasing);

figure
my_plotWithConf(T,mdec,sdec,[1 0 0]);
hold on;
my_plotWithConf(T,minc,sinc,[0 0 1]);
plot(T,mmdec,'r-');
plot(T,mminc,'b-');

d=minc-mdec;
[m,i]=max(d(1:100));
x=increasing(:,i);
y=decreasing(:,i);
pvalue=boottest(x,y,1e4);
tx=T(i);
text(tx-0.2,0.85,['Data P=',num2str(pvalue)]);

d=mminc-mmdec;
[m,i]=max(d(1:100));
x=model_increasing(:,i);
y=model_decreasing(:,i);
pvalue=boottest(x,y,1e4);

text(tx-0.2,0.75,['Model P=',num2str(pvalue)]);

figure
sshort=std(short)/sqrt(size(short,1));
slong=std(long)/sqrt(size(long,1));
mshort=mean(short);
mlong=mean(long);
T=linspace(-1,1,200);

my_plotWithConf(T,mshort,sshort,[1 0 0]);
hold on;
my_plotWithConf(T,mlong,slong,[0 0 1]);
plot(T,mean(model_short),'r-');
plot(T,mean(model_long),'b-');
% % d=minc-mdec;
% [m,i]=max(d(1:100));
% x=increasing(:,i);
% y=decreasing(:,i);
% pvalue=boottest(x,y,min(length(x),length(y)),1e4);
% text(T(i),0.75,['P=',num2str(pvalue)]);


%% save
% save([filename,'_',mode])