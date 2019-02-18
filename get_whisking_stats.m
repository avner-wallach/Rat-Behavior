function analyzed_sturct=get_whisking_stats(wdata,plotting)

if(nargin==1)
    plotting=0;
end

%%
fs=500;
mode='spline';
maxlag=fs;
option='coeff';
minamp=3; %cutoff between whisking and non-whisking amplitudes

N_vector=zeros(1,2*maxlag+1);
all_amp=[];
all_offset=[];
all_dur=[];
all_rat=[];
prot_pump=[];
ret_pump=[];

amp=[];
off=[];
dur=[nan];
rise=[];
pr_ratio=[nan];
dur_1r=[];
dur_2r=[];
k=1;
m=1;
%%
for i=1:length(wdata)
    i
    ddata=wdata{i};
    M=numel(ddata);
    while(numel(ddata)>maxlag)
        data=ddata;%(1:min(1e3,numel(ddata)));        
        duration=numel(data)/fs;
        [decomposed_data]=decompose_whisking(data,fs,mode);
        amplitude=decomposed_data.amplitude;
        ind=find(~isnan(amplitude) & amplitude>=minamp);
        if(~numel(ind))
            break;
        end
        d(m)=duration;
        a(m)=nanmedian(amplitude);
        offset=decomposed_data.offset;
        phi=decomposed_data.phi;
        v=decomposed_data.v;
        durations=decomposed_data.durations;   
        ratios=decomposed_data.ratios;   
        segmented_struct=decomposed_data.segmented_struct;
        amp=[amp decomposed_data.d_amp];
        off=decomposed_data.d_offset;
        pr_ratio=[pr_ratio segmented_struct.pr_ratio NaN];
        dur=[dur segmented_struct.duration NaN];
        rise=[rise segmented_struct.rise];
%         if(numel(ind)<maxlag)
%             continue;
%         end

        %distributions
        all_amp=[all_amp amplitude];
        all_offset=[all_offset offset];
    %     all_dur=[all_dur durations];
    %     all_rat=[all_rat ratios];

        %xcov
        [x_amp(m,:),lags]=xcov(amplitude(ind),maxlag,option);
%         x_amp(m,:)=x_amp(m,:)/d(m);
        [x_offset(m,:),lags]=xcov(offset(ind),maxlag,option);
        [x_amp_offset(m,:),lags]=xcov(amplitude(ind),offset(ind),maxlag,option);    
        [x_doff(m,:),lags]=xcov(abs(diff(offset(ind))),maxlag,option);
        [x_amp_doffset(m,:),lags]=xcov(amplitude(ind(1:end-1)),abs(diff(offset(ind))),maxlag,option);
        [x_duration(m,:),lags]=xcov(durations(ind),maxlag,option);
        [x_ratio(m,:),lags]=xcov(ratios(ind),maxlag,option);
        indout=find(abs(lags)>=numel(ind));
        indin=find(abs(lags)<numel(ind));
        x_amp(m,indout)=NaN;
        x_offset(m,indout)=NaN;
        x_amp_offset(m,indout)=NaN;
        x_doff(m,indout)=NaN;
        x_amp_doffset(m,indout)=NaN;
        x_duration(m,indout)=NaN;
        x_ratio(m,indout)=NaN;
%         x_amp(m,:)=x_amp(m,:)*duration;
        x_offset(m,:)=x_offset(m,:)*duration;
        x_amp_offset(m,:)=x_amp_offset(m,:)*duration;
        x_doff(m,:)=x_doff(m,:)*duration;
        x_amp_doffset(m,:)=x_amp_doffset(m,:)*duration;
%         x_duration(m,:)=x_duration(m,:)*duration;
        x_ratio(m,:)=x_ratio(m,:)*duration;
        
        mdur(m)=nanmean(durations(ind));
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
            if(prot_amp<5)
                continue;
            end
            scaled_ret=interp1(linspace(-1,0,length(ret)),ret,linspace(-1,0,Nsamp));
            scaled_prot=interp1(linspace(0,1,length(prot)),prot,linspace(0,1,Nsamp));

            whisk(k,:)=[scaled_prot scaled_ret];
            whisk(k,:)=(whisk(k,:)-min(prot))/(prot_amp);
            k=k+1;
            dprot=diff(segmented_struct.wphase{j,1});
            proty=sum(dprot(dprot<0));
            dret=diff(segmented_struct.wphase{j,2});
            rety=sum(dret(dret<0));        
            [y,proti]=min(dprot);
            [y,reti]=min(diff(segmented_struct.wphase{j,2}));
            prot_pump=[prot_pump;[proti proty]];
            ret_pump=[ret_pump;[reti rety]];
        end  
        
%         if(numel(ddata)>1e3+maxlag)
%             ddata=ddata(1e3+1:end);
%         else
            ddata=[];
%         end
        m=m+1;
    end
end

lags=lags/fs;

x_amp_m=nansum(x_amp)/size(x_amp,1);
sig_amp=sqrt(max(x_amp_m));


% for i=1:size(x_amp,1)
% ind=find(x_amp(i,500:end-1)>0 & x_amp(i,501:end)<=0);
% x_zc_amp(i)=ind(1)*2;
% x_int_amp(i)=sum(x_amp(i,500:(500+ind(1))));
% end
% 
% for i=1:size(x_offset,1)
% ind=find(x_offset(i,500:end-1)>0 & x_offset(i,501:end)<=0);
% x_zc_off(i)=ind(1)*2;
% end
% 
x_doff_m=nansum(x_doff)/size(x_doff,1);
sig_doff=sqrt(max(x_doff_m));

x_offset_m=nansum(x_offset)/size(x_offset,1);
sig_offset=sqrt(max(x_offset_m));

x_duration_m=nansum(x_duration)/size(x_duration,1);
sig_duration=sqrt(max(x_duration_m));

% for i=1:size(x_duration,1)
% ind=find(x_duration(i,500:end-1)>0 & x_duration(i,501:end)<=0);
% x_zc_dur(i)=ind(1)*2;
% end

x_ratio_m=nansum(x_ratio)/size(x_ratio,1);
sig_ratio=sqrt(max(x_ratio_m));

% for i=1:size(x_ratio,1)
% ind=find(x_ratio(i,500:end-1)>0 & x_ratio(i,501:end)<=0);
% x_zc_rat(i)=ind(1)*2;
% end

x_amp_offset_m=nansum(x_amp_offset)/size(x_amp_offset,1);
x_amp_doffset_m=nansum(x_amp_doffset)/size(x_amp_doffset,1);
%% temporal correlations
if(plotting)
    figure;
    plot(lags,x_amp_m/sig_amp^2);
    hold on;
    plot(lags,x_offset_m/sig_offset^2,'r');
    plot(lags,x_duration_m/sig_duration^2,'m');
    plot(lags,x_ratio_m/sig_ratio^2,'g');
    plot(lags,x_amp_offset_m/(sig_amp*sig_offset),'k');
    legend({'amplitude','offset','duration','ratio','amp-offset'});
end

analyzed_sturct.x_lags=lags;
analyzed_sturct.x_amp=x_amp;
analyzed_sturct.x_amp_m=x_amp_m;%/sig_amp^2;
analyzed_sturct.med_amp=a;
% analyzed_sturct.x_int_amp=x_int_amp;
% analyzed_sturct.x_zc_amp=x_zc_amp;
% analyzed_sturct.x_zc_off=x_zc_off;
% analyzed_sturct.x_zc_dur=x_zc_dur;
% analyzed_sturct.x_zc_rat=x_zc_rat;
analyzed_sturct.d=d;
analyzed_sturct.x_offset=x_offset_m/sig_offset^2;
analyzed_sturct.x_duration_m=x_duration_m/sig_duration^2;
analyzed_sturct.x_duration=x_duration;
analyzed_sturct.x_ratio=x_ratio_m;
analyzed_sturct.mdur=mdur;
%% dustributions

%bin=3;
clear ctrs;
ctrs{1}=linspace(0,50,50);
ctrs{2}=linspace(20,120,50);
Nx=hist3([all_amp' all_offset'],ctrs);
Nx=Nx/sum(Nx(:));
Namp=hist(all_amp',ctrs{1});
Namp=Namp/sum(Namp);
Noff=hist(all_offset',ctrs{2});
Noff=Noff/sum(Noff);

if(plotting)
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
end

minamp=4;
ind=find(amp(2:2:end)>minamp);
analyzed_sturct.amp_off=[all_amp' all_offset'];
analyzed_sturct.dur_rat=[dur(ind)' pr_ratio(ind)'];

if(plotting)
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
    % bar(ctrs,Nr);
    cdfplot(pr_ratio(ind));
    xlabel('Protraction-Retraction Ratio');
end
%% pumps
th=linspace(-0.5,0,200);
for i=1:200
pp(i)=sum(prot_pump(:,2)<=th(i))/numel(prot_pump(:,2));
pr(i)=sum(ret_pump(:,2)<=th(i))/numel(ret_pump(:,2));
end

analyzed_sturct.prot_pump=pp;
analyzed_sturct.ret_pump=pr;

if(plotting)
    figure, plot(th,pp);
    hold on;
    plot(th,pr,'r');
    xlabel('pump threshold')
    ylabel('precentage of whisks w/pumps')
    legend({'protraction','retraction'});

    indpump=find(prot_pump(:,2)<-0.05);
    pump_whisk=mean(whisk(indpump,:));
    nopump_whisk=mean(whisk(prot_pump(:,2)>=-0.05,:));
    figure, plot(pump_whisk);
    hold on;
    plot(nopump_whisk,'r');

    indpump=find(ret_pump(:,2)<-0.05);
    pump_whisk=mean(whisk(indpump,:));
    nopump_whisk=mean(whisk(ret_pump(:,2)>=-0.05,:));
    figure, plot(pump_whisk);
    hold on;
    plot(nopump_whisk,'r');

    [Y,ind]=sort(prot_pump(indpump,1));
    figure, plot(whisk(indpump(ind),:)');
end
end
