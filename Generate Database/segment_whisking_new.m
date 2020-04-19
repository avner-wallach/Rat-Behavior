function [segmented_struct]=segment_whisking_new(wsig,phi,swind,fs)
% written by Avner Wallach
% segment whisking signal using phase information
% variables:
%   wsig- whisking signal, vector of angle samples.
%   phi- phase of wsig, extracted using Hilbert transform
%   swind- window (in radians) around multiples of pi in which exrema are
%       looked for. 
%   fs- sampling frequency (default is 500) 
%
% segmented_struct- a structure with the segmented variables in fields:
%
% for the entire signal:
% ind_crossing are indices of mid-points in whisking signal
% crossing are angles at mid-points in whisking signal
% ind_trough are indices of troughs in whisking signal
% trough are angles at troughs (='setpoint')
% ind_peak are indices of peaks
% peak are angles at peaks
%
%%%%% IMPORTANT NOTE: partial whisks at the begining/end of wsig are
%%%%% discraded in the segmentation process. 
% for complete whisks only (discarding partial whisks):
% wstart are indices of whisk start (trough)
% wmiddle are indices of whisk middle (peak)
% whisks is an N-by-2 cell array, N being the number of complete whisks in signal.
%   first column are angle dynamics in protraction
%   second column are angle dynamics in retraction.
% duration are total whisk duration
% pr_ratio are protraction/retraction duration ratio (1 for equal duration)
% rise are protraction 'amplitude' (trough to subsequent peak)
% fall are retraction 'amplitude' (peak to subsequent trough)


if nargin<3
    swind=pi/4;
    fs=500;
end
if nargin<4
    fs=500;
end
jerk=[diff(wsig,3);nan(3,1)];

x=unique(phi);
y=1:length(x);
X_crossings=(ceil(phi(1)*2/pi)/2*pi+pi/2):pi:max(x);
X_troughs=((ceil((phi(1)-pi)/(2*pi))*2+1)*pi):2*pi:max(x);
X_peaks=(ceil(phi(1)/(2*pi))*2*pi):2*pi:max(x);
ind_crossings=round(interp1(x,y,X_crossings));
crossings=wsig(ind_crossings);

for i=1:length(X_troughs)
    ind=find(x>=(X_troughs(i)-swind/2) & x<(X_troughs(i)+swind/2));
    if(~numel(ind))
        ind_trough(i)=round(interp1(x,y,X_troughs(i)));
        trough(i)=wsig(ind_trough(i));
    else
        [trough(i) idx]=min(wsig(ind));
        ind_trough(i)=ind(idx);
    end
end

for i=1:length(X_peaks)
    ind=find(x>=(X_peaks(i)-swind/2) & x<(X_peaks(i)+swind/2));
    if(~numel(ind))
        ind_peak(i)=round(interp1(x,y,X_peaks(i)));
        peak(i)=wsig(ind_peak(i));
    else    
        [peak(i) idx]=max(wsig(ind));
        ind_peak(i)=ind(idx);
    end
end

wstart=[];
wmiddle=[];
duration=[];
pr_ratio=[];
rise=[];
fall=[];
midpoint=[];
% ind_pump=[];
whisks=cell(length(trough)-1,2);
jerk_costs=[];
iind=find(ind_peak>ind_trough(1));
for i=1:length(trough)-1
    wstart(i)=ind_trough(i);
    wmiddle(i)=ind_peak(iind(i));
    whisks{i,1}=wsig(ind_trough(i):ind_peak(iind(i))-1);
    whisks{i,2}=wsig(ind_peak(iind(i)):ind_trough(i+1)-1);
    wphase{i,1}=phi(ind_trough(i):ind_peak(iind(i))-1);
    wphase{i,2}=phi(ind_peak(iind(i)):ind_trough(i+1)-1);
    duration(i)=(ind_trough(i+1)-ind_trough(i))/fs;
    pr_ratio(i)=(ind_peak(iind(i))-ind_trough(i))/(ind_trough(i+1)-ind_peak(iind(i)));
    rise(i)=peak(iind(i))-trough(i);
    fall(i)=peak(iind(i))-trough(i+1);
    midpoint(i,:)=[(peak(iind(i))+trough(i))/2 (peak(iind(i))+trough(i+1))/2];
    j{1}=jerk(ind_trough(i):ind_peak(iind(i))-1)/rise(i);
    j{2}=jerk(ind_peak(iind(i)):ind_trough(i+1)-1)/fall(i);
    if(duration(i)<0 || pr_ratio(i)<0)
        display('error');
    end
    
    %pump analysis- protraction    
    dw=diff(whisks{i,1})*fs; %angular speed
    if(numel(dw)>10)
        maxspeed(i,1)=max(dw);
        [pks,locs,w,p]=findpeaks(-dw);
        if(numel(pks)) %trough exists=pump
           [pump_prom(i,1),I]=max(p);
           pump_val(i,1)=-pks(I);
           pump_phase(i,1)=locs(I)/numel(dw);
           pump_width(i,1)=w(I)/numel(dw);
           pump_ind(i,1)=ind_trough(i)+locs(I)-1;
%            ind_pump=[ind_pump ind_trough(i)+locs(I)-1];
        else
           pump_prom(i,1)=0;
           pump_val(i,1)=0;
           pump_phase(i,1)=0;
           pump_width(i,1)=0;
           pump_ind(i,1)=ind_trough(i)+(ind_peak(iind(i))-ind_trough(i))/2;
        end
    else
       pump_prom(i,1)=0;
       pump_val(i,1)=0;
       pump_phase(i,1)=0;
       pump_width(i,1)=0;
       pump_ind(i,1)=ind_trough(i)+(ind_peak(iind(i))-ind_trough(i))/2;
    end
    %pump analysis- retraction
    dw=-diff(whisks{i,2})*fs; %angular speed
    if(numel(dw)>10)
        maxspeed(i,2)=max(dw);
        [pks,locs,w,p]=findpeaks(-dw);
        if(numel(pks)) %trough exists=pump
           [pump_prom(i,2),I]=max(p);
           pump_val(i,2)=-pks(I);
           pump_phase(i,2)=locs(I)/numel(dw);
           pump_width(i,2)=w(I)/numel(dw);
           pump_ind(i,2)=ind_peak(iind(i))+locs(I)-1;
%            ind_pump=[ind_pump ind_peak(iind(i))+locs(I)-1];
        else
           pump_prom(i,2)=0;
           pump_val(i,2)=0;
           pump_phase(i,2)=0;
           pump_width(i,2)=0;
           pump_ind(i,2)=ind_peak(iind(i))+(ind_trough(i+1)-ind_peak(iind(i)))/2;
        end
    else
       pump_prom(i,2)=0;
       pump_val(i,2)=0;
       pump_phase(i,2)=0;
       pump_width(i,2)=0;
       pump_ind(i,2)=ind_peak(iind(i))+(ind_trough(i+1)-ind_peak(iind(i)))/2;
    end
    
    %jerk cost analysis
    for k=1:2
        jerk_costs(i,k)=nanmean(j{k}.^2);
    end
        
end

segmented_struct.ind_crossings=ind_crossings;
segmented_struct.crossing=crossings;
segmented_struct.ind_trough=ind_trough;
segmented_struct.trough=trough;
segmented_struct.ind_peak=ind_peak;
segmented_struct.peak=peak;
segmented_struct.wstart=wstart;
segmented_struct.wmiddle=wmiddle;
segmented_struct.whisks=whisks;
segmented_struct.wphase=wphase;
segmented_struct.duration=duration;
segmented_struct.pr_ratio=pr_ratio;
segmented_struct.rise=rise;
segmented_struct.fall=fall;
segmented_struct.midpoint=midpoint;
segmented_struct.pump_prom=pump_prom;
segmented_struct.pump_val=pump_val;
segmented_struct.pump_width=pump_width;
segmented_struct.pump_phase=pump_phase;
% segmented_struct.ind_pumps=ind_pump;
segmented_struct.pump_ind=round(pump_ind);
segmented_struct.maxspeed=maxspeed;
segmented_struct.jerk_costs=jerk_costs;
end
