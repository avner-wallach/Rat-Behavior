function [decomposed_struct]=decompose_whisking_new(wsig,hsig,fs,mode)
% written by Avner Wallach
% decompose whisking signal using phase-amplitude-offset model
% variables:
%   wsig- whisking signal, vector of angle samples.
%   hsig- head motion matrix, [direction X Y]
%   fs- sampling frequency (default is 500) 
%   mode- interpolation mode of slow variables:
%       'linear'- piecwise linear interpolation (default)
%       'spline'- cubic spline interpolation (best)
%       'zoh'   - zero order hold, causal piecwise constant interpolation
%       'nearest'- non-causal piecwise constant interpolation
%
%%%%% IMPORTANT NOTE: partial whisks at the begining/end of wsig are
%%%%% discraded in the segmentation process. all the interpolated signals
%%%%% (amplitude, offset etc) are assigned the NaN value in these segments.
%
% decomposed_struct- a structure with the dcomposed variables in fields:
%   t- time vector of signal
%   amplitude- slow interpolated amplitude singal
%   offset- slow interpolated offset singal
%   phi- phase singal
%   durations- slow interpolated whisk duration singal
%   ratios- slow interpolated protraction/retraction duration ratio singal
%   v- estimation of wsig according to model
%           v(t)=amplitude(t)*cos(phi(t))+offset(t)
%   d_amp- amplitude values at peak and trough times
%   d_offset- offset values at peak and trough times
%   segmented_struct- output of segmentation algorithm


if nargin<4
    mode='pchip';    
end
if nargin<3
	fs=500;
end
if nargin<2
    hsig=[];
end
if(numel(hsig)==0)
    hsig=zeros(numel(wsig),3);
% else
%     hsig(:,1)=hsig(:,1)-nanmedian(hsig(:,1));
end
amp_hilbert=0; %0=use interpolation for amplitude; 1=use hilbert trans. abs value for ampltidue
%pre-filter
lp = [20];
lp = lp * 2 / fs; % convert Hz to radians/S
[N, Wn] = buttord( lp, lp .* [1.5], 3, 20); 
[B,A] = butter(N,Wn);
wsig = filtfilt(B,A,wsig); % zero-phase LPF

lp = [4];
lp = lp * 2 / fs; % convert Hz to radians/S
[N, Wn] = buttord( lp, lp .* [1.5], 3, 20); 
[B,A] = butter(N,Wn);
woff = filtfilt(B,A,wsig)'; % zero-phase LPF
wsig_hpf=wsig-woff';
Xh=hilbert(wsig_hpf);
phi=unwrap(angle(Xh));
wphi=mod(phi+pi,2*pi)-pi;
wamp_h=abs(Xh);
wamp= filtfilt(B,A,wamp_h)'; % zero-phase LPF

V=sqrt(hsig(:,2).^2+hsig(:,3).^2); %head velocity
% x=unique(phi);
% y=1:length(x);
% X_troughs=((ceil((phi(1)-pi)/(2*pi))*2+1)*pi):2*pi:max(x);
% ind_troughs=round(interp1(x,y,X_troughs));
% wfreq=nan(size(wsig));
% wtime=nan(size(wsig));
% for i=1:length(ind_troughs)-1
%     wfreq(ind_troughs(i):ind_troughs(i+1))=fs./(ind_troughs(i+1)-ind_troughs(i));
%     wtime(ind_troughs(i):ind_troughs(i+1))=0:(ind_troughs(i+1)-ind_troughs(i));
% end
% wtime=wtime/fs;

bp = [4 20];
bp = bp * 2 / fs; % convert Hz to radians/S
swind=pi/4;
t=[1:length(wsig)]/fs;

[N, Wn] = buttord( bp(2), bp(2).*1.5, 3, 20); 
[B,A] = butter(N,Wn);
X = filtfilt(B,A,wsig); % zero-phase filtering
[segmented_struct]=segment_whisking_new(X,phi,swind,fs);

%ind=(max(segmented_struct.ind_peak(1),segmented_struct.ind_trough(1)):...
%    min(segmented_struct.ind_peak(end),segmented_struct.ind_trough(end)));
ind=1:numel(t);
tops=nan(size(t));
bottoms=nan(size(t));
middles=nan(size(t));
durations=nan(size(t));
ratios=nan(size(t));
amplitude=nan(size(t));
offset=nan(size(t));
% if(~strcmp(mode,'zoh'))
tops(ind)=interp1(segmented_struct.ind_peak/fs,segmented_struct.peak,t(ind),mode);
bottoms(ind)=interp1(segmented_struct.ind_trough/fs,segmented_struct.trough,t(ind),mode);
middles(ind)=interp1(segmented_struct.ind_crossings/fs,segmented_struct.crossing,t(ind),mode);
% durations(ind)=interp1(segmented_struct.wmiddle/fs,segmented_struct.duration,t(ind),mode);
% ratios(ind)=interp1(segmented_struct.wmiddle/fs,segmented_struct.pr_ratio,t(ind),mode);
offset=(tops+bottoms)/2;
amplitude=(tops-bottoms)/2;
if( amp_hilbert==0)
    wamp=amplitude;
    woff=offset;
end
% else %zero order hold
%     tops(ind)=zoh(segmented_struct.ind_peak,segmented_struct.peak,ind);
%     bottoms(ind)=zoh(segmented_struct.ind_trough,segmented_struct.trough,ind);
% 
%     X=reshape([segmented_struct.wstart' segmented_struct.wmiddle']',[],1);
%     Y=reshape([segmented_struct.rise' segmented_struct.fall']',[],1);
%     amplitude(ind)=zoh(X,Y,ind)/2;
%     
%     off=nan(size(t));
%     off(segmented_struct.ind_trough)=segmented_struct.trough+amplitude(segmented_struct.ind_trough);
%     off(segmented_struct.ind_peak)=segmented_struct.peak-amplitude(segmented_struct.ind_peak);    
%     offset(ind)=zoh(X,off(X),ind);
%     
%     durations(ind)=zoh(segmented_struct.wstart,segmented_struct.duration,ind);
%     ratios(ind)=zoh(segmented_struct.wstart,segmented_struct.pr_ratio,ind);
% end

% ind=reshape([segmented_struct.wstart' segmented_struct.wmiddle']',[],1);
minamp=2; %minimal amplitude for active whisking

d_amp=wamp(segmented_struct.ind_crossings);
d_offset=woff(segmented_struct.ind_crossings);
d_hoffset=d_offset(:)+hsig(segmented_struct.ind_crossings,1);
d_hoffset_centered=d_offset(:)-nanmedian(hsig(:,1))+hsig(segmented_struct.ind_crossings,1);
d_V=V(segmented_struct.ind_crossings);
durations=segmented_struct.duration;
ratios=segmented_struct.pr_ratio;
rises=segmented_struct.rise;
falls=segmented_struct.fall;
durations(wamp(segmented_struct.wmiddle)<minamp)=NaN;
ratios(wamp(segmented_struct.wmiddle)<minamp)=NaN;
prots=durations.*ratios./(1+ratios);
rets=durations./(1+ratios);
% v=offset+amplitude.*cos(phi');
%% bout analysis
bout_durations=[];
bout_whisks=[];
inbout=(wamp>=minamp);
ind_bstart=find(diff(inbout)==1);
if(numel(ind_bstart)) %bout started in tracked segmented
    ind_bend=find(diff(inbout)==-1);
    ind_bend=ind_bend(ind_bend>ind_bstart(1));
    for i=1:numel(ind_bend)
        bout_durations(i)=(ind_bend(i)-ind_bstart(i))/fs;
        bout_whisks(i)=numel(find(segmented_struct.ind_crossings>ind_bstart(i) & ...
            segmented_struct.ind_crossings<ind_bend(i)))/2;
    end
end

%% correlations and dilutions
N=1e3;
minamp=2;
% [diluted_amp,X_amp,lag_amp]=dilute_series_new(d_amp,N,find(d_amp>=minamp),100);
% [diluted_off,X_off,lag_off]=dilute_series_new(d_offset,N,find(d_amp>=minamp),100);
% lag_amp_off=max(lag_amp,lag_off);
% ind_dilute=[ceil(lag_amp_off/2):lag_amp_off:numel(d_amp)];
% ind_dilute=intersect(ind_dilute,find(d_amp>=minamp));
% diluted_amp_off=[d_amp(ind_dilute)' d_offset(ind_dilute)'];
% [diluted_dur,X_dur,lag_dur]=dilute_series_new(durations,N,find(wamp(segmented_struct.wmiddle)>=minamp),50);
% [diluted_pr,X_pr,lag_pr]=dilute_series_new(ratios,N,find(wamp(segmented_struct.wmiddle)>=minamp),50);
% [diluted_prot,X_prot,lag_prot]=dilute_series_new(prots,N,find(wamp(segmented_struct.wmiddle)>=minamp),50);
% [diluted_ret,X_ret,lag_ret]=dilute_series_new(rets,N,find(wamp(segmented_struct.wmiddle)>=minamp),50);
[diluted_amp,X_amp,lag_amp]=dilute_series(d_amp,N,find(d_amp>=minamp));
[diluted_off,X_off,lag_off]=dilute_series(d_offset,N,find(d_amp>=minamp));
[diluted_hoff,X_hoff,lag_hoff]=dilute_series(d_hoffset,N,find(d_amp>=minamp & ~isnan(d_hoffset')));
[diluted_hoff_ctr,X_hoff_ctr,lag_hoff_ctr]=dilute_series(d_hoffset_centered,N,find(d_amp>=minamp& ~isnan(d_hoffset_centered')));
lag_amp_off=max(lag_amp,lag_off);
ind_dilute=[ceil(lag_amp_off/2):lag_amp_off:numel(d_amp)];
ind_dilute=intersect(ind_dilute,find(d_amp>=minamp));
diluted_amp_off=[d_amp(ind_dilute)' d_offset(ind_dilute)'];
[diluted_dur,X_dur,lag_dur]=dilute_series(durations,N,find(wamp(segmented_struct.wmiddle)>=minamp));
[diluted_pr,X_pr,lag_pr]=dilute_series(ratios,N,find(wamp(segmented_struct.wmiddle)>=minamp));
[diluted_prot,X_prot,lag_prot]=dilute_series(prots,N,find(wamp(segmented_struct.wmiddle)>=minamp));
[diluted_ret,X_ret,lag_ret]=dilute_series(rets,N,find(wamp(segmented_struct.wmiddle)>=minamp));
%% 
decomposed_struct.t=t;
decomposed_struct.wsig=wsig;
decomposed_struct.wamp=wamp;
decomposed_struct.woff=woff;
decomposed_struct.wphi=wphi;

decomposed_struct.bout_dirations=bout_durations;
decomposed_struct.bout_whisks=bout_whisks;

decomposed_struct.durations=durations;
decomposed_struct.ratios=ratios;
decomposed_struct.rises=rises;
decomposed_struct.falls=falls;
decomposed_struct.midpoint=segmented_struct.midpoint;
decomposed_struct.amplitudes=d_amp(:);
decomposed_struct.offsets=d_offset(:);
decomposed_struct.hoffsets=d_hoffset(:);
decomposed_struct.hoffsets_ctr=d_hoffset_centered(:);
decomposed_struct.V=d_V(:);
decomposed_struct.med_V=nanmedian(V);
decomposed_struct.diluted_amp=diluted_amp;
decomposed_struct.diluted_off=diluted_off;
decomposed_struct.diluted_hoff=diluted_hoff;
decomposed_struct.diluted_hoff_ctr=diluted_hoff_ctr;
decomposed_struct.diluted_amp_off=diluted_amp_off;
decomposed_struct.diluted_dur=diluted_dur;
decomposed_struct.diluted_pr=diluted_pr;
decomposed_struct.diluted_prot=diluted_prot;
decomposed_struct.diluted_ret=diluted_ret;
decomposed_struct.X_amp=X_amp;
decomposed_struct.X_off=X_off;
decomposed_struct.X_hoff=X_hoff;
decomposed_struct.X_hoff_ctr=X_hoff_ctr;
decomposed_struct.X_amp_off=[];
decomposed_struct.X_dur=X_dur;
decomposed_struct.X_pr=X_pr;
decomposed_struct.X_prot=X_prot;
decomposed_struct.X_ret=X_ret;
decomposed_struct.lag_amp=lag_amp;
decomposed_struct.lag_off=lag_off;
decomposed_struct.lag_hoff=lag_hoff;
decomposed_struct.lag_hoff_ctr=lag_hoff_ctr;
decomposed_struct.lag_amp_off=lag_amp_off;
decomposed_struct.lag_dur=lag_dur;
decomposed_struct.lag_pr=lag_pr;
decomposed_struct.lag_prot=lag_prot;
decomposed_struct.lag_ret=lag_ret;

decomposed_struct.pump_val=segmented_struct.pump_val;
decomposed_struct.pump_phase=segmented_struct.pump_phase;
decomposed_struct.pump_prom=segmented_struct.pump_prom;
decomposed_struct.pump_width=segmented_struct.pump_width;
decomposed_struct.ind_pumps=segmented_struct.ind_pumps;
decomposed_struct.maxspeed=segmented_struct.maxspeed;
decomposed_struct.rise_fall=[segmented_struct.rise' segmented_struct.fall'];
decomposed_struct.jerk_costs=segmented_struct.jerk_costs;

decomposed_struct.whisks=segmented_struct.whisks;
decomposed_struct.segmented_struct=segmented_struct;

end
