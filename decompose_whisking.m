function [decomposed_struct]=decompose_whisking_new(wsig,fs,mode)
% written by Avner Wallach
% decompose whisking signal using phase-amplitude-offset model
% variables:
%   wsig- whisking signal, vector of angle samples.
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


if nargin<2
    fs=500;
    mode='linear';    
end

if nargin<3
    mode='linear';    
end

lp = [4];
lp = lp * 2 / fs; % convert Hz to radians/S
[N, Wn] = buttord( lp, lp .* [1.5], 3, 20); 
[B,A] = butter(N,Wn);
woff = filtfilt(B,A,wsig); % zero-phase LPF
wsig_hpf=wsig-woff;
Xh=hilbert(wsig_hpf);
phi=unwrap(angle(Xh));
wphi=mod(phi+pi,2*pi)-pi;
wamp=abs(Xh);

x=unique(phi);
y=1:length(x);
X_troughs=((ceil((phi(1)-pi)/(2*pi))*2+1)*pi):2*pi:max(x);
ind_troughs=round(interp1(x,y,X_troughs));
wfreq=nan(size(wsig));
wtime=nan(size(wsig));
for i=1:length(ind_troughs)-1
    wfreq(ind_troughs(i):ind_troughs(i+1))=fs./(ind_troughs(i+1)-ind_troughs(i));
    wtime(ind_troughs(i):ind_troughs(i+1))=0:(ind_troughs(i+1)-ind_troughs(i));
end
wtime=wtime/fs;

bp = [4 20];
bp = bp * 2 / fs; % convert Hz to radians/S
swind=pi/4;
t=[1:length(wsig)]/fs;

[N, Wn] = buttord( bp(2), bp(2).*1.5, 3, 20); 
[B,A] = butter(N,Wn);
X = filtfilt(B,A,wsig); % zero-phase filtering
[segmented_struct]=segment_whisking(X,phi,swind,fs);
ind=(max(segmented_struct.ind_peak(1),segmented_struct.ind_trough(1)):...
    min(segmented_struct.ind_peak(end),segmented_struct.ind_trough(end)));
tops=nan(size(t));
bottoms=nan(size(t));
middles=nan(size(t));
durations=nan(size(t));
ratios=nan(size(t));
amplitude=nan(size(t));
offset=nan(size(t));
if(~strcmp(mode,'zoh'))
    tops(ind)=interp1(segmented_struct.ind_peak/fs,segmented_struct.peak,t(ind),mode);
    bottoms(ind)=interp1(segmented_struct.ind_trough/fs,segmented_struct.trough,t(ind),mode);
    middles(ind)=interp1(segmented_struct.ind_crossing/fs,segmented_struct.crossing,t(ind),mode);
    durations(ind)=interp1(segmented_struct.wmiddle/fs,segmented_struct.duration,t(ind),mode);
    ratios(ind)=interp1(segmented_struct.wmiddle/fs,segmented_struct.pr_ratio,t(ind),mode);
    offset=(tops+bottoms)/2;
    amplitude=(tops-bottoms)/2;
else %zero order hold
    tops(ind)=zoh(segmented_struct.ind_peak,segmented_struct.peak,ind);
    bottoms(ind)=zoh(segmented_struct.ind_trough,segmented_struct.trough,ind);

    X=reshape([segmented_struct.wstart' segmented_struct.wmiddle']',[],1);
    Y=reshape([segmented_struct.rise' segmented_struct.fall']',[],1);
    amplitude(ind)=zoh(X,Y,ind)/2;
    
    off=nan(size(t));
    off(segmented_struct.ind_trough)=segmented_struct.trough+amplitude(segmented_struct.ind_trough);
    off(segmented_struct.ind_peak)=segmented_struct.peak-amplitude(segmented_struct.ind_peak);    
    offset(ind)=zoh(X,off(X),ind);
    
    durations(ind)=zoh(segmented_struct.wstart,segmented_struct.duration,ind);
    ratios(ind)=zoh(segmented_struct.wstart,segmented_struct.pr_ratio,ind);
end

ind=reshape([segmented_struct.wstart' segmented_struct.wmiddle']',[],1);

d_amp=amplitude(ind);
d_offset=offset(ind);

v=offset+amplitude.*cos(phi');

%% 
decomposed_struct.t=t;
decomposed_struct.amplitude=amplitude;
decomposed_struct.offset=offset;
decomposed_struct.phi=phi;
decomposed_struct.durations=durations;
decomposed_struct.ratios=ratios;
decomposed_struct.v=v;
decomposed_struct.d_amp=d_amp;
decomposed_struct.d_offset=d_offset;
decomposed_struct.segmented_struct=segmented_struct;
decomposed_struct.wsig_hpf=wsig_hpf;
end
