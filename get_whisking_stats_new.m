function [analyzed_sturct,pdata]=get_whisking_stats_new(wdata,hdata,plotting)

if(nargin==1)
    plotting=0;
end

%%
fs=500;
mode='spline';
maxlag=fs;
option='coeff';
minamp=2; %cutoff between whisking and non-whisking amplitudes

all_amp=[];
all_offset=[];
all_hoffset=[];
all_hoffset_ctr=[];
all_V=[];
all_dur=[];
all_rat=[];
all_rise=[];
all_fall=[];
all_d_amp=[];
all_d_offset=[];
all_d_dur=[];
all_d_rat=[];
all_wamp=[];

all_bout_dur=[];
all_bout_whisks=[];

all_pump_val=[];
all_pump_phase=[];
all_pump_prom=[];
all_pump_width=[];
all_pump_maxspeed=[];
all_rise_fall=[];
all_midpoints=[];
all_jerk_costs=[];

pdata=cell(0);

var={'amp','off','hoff','hoff_ctr','amp_off','dur','pr','prot','ret'};
field={'lag_','X_','diluted_'};
for k=1:numel(field)
    for j=1:numel(var)
        eval(['all_',field{k},var{j},'=[];']);
    end
end
k=1;
m=1;
%%
for i=1:length(wdata)
    i    
    ddata=wdata{i};    
    hhdata=hdata{i};
    M=numel(ddata);
    while(numel(ddata)>maxlag) % >1s
        duration=numel(ddata)/fs;
        [decomposed_data]=decompose_whisking_new(ddata,hhdata,fs,mode);
        d(m)=duration;
        a(m)=nanmedian(decomposed_data.amplitudes);
        v(m)=decomposed_data.med_V;
        
        %get measures per whisk
        all_amp=[all_amp ; decomposed_data.amplitudes];
        all_offset=[all_offset ; decomposed_data.offsets];
        all_hoffset=[all_hoffset ; decomposed_data.hoffsets];
        all_hoffset_ctr=[all_hoffset_ctr ; decomposed_data.hoffsets_ctr];
        all_V=[all_V; decomposed_data.V];
        all_dur=[all_dur ; decomposed_data.durations'];
        all_rat=[all_rat ; decomposed_data.ratios'];
        all_wamp=[all_wamp ; decomposed_data.wamp(:)];
        s(m)=numel(decomposed_data.amplitudes);
        w(m)=size(decomposed_data.maxspeed,1);
        
        %get bout measures
        all_bout_dur=[all_bout_dur;decomposed_data.bout_dirations'];
        all_bout_whisks=[all_bout_whisks;decomposed_data.bout_whisks'];

        %get pump measures
        all_pump_val=[all_pump_val;decomposed_data.pump_val];
        all_pump_phase=[all_pump_phase;decomposed_data.pump_phase];
        all_pump_prom=[all_pump_prom;decomposed_data.pump_prom];
        all_pump_width=[all_pump_width;decomposed_data.pump_width];
        all_pump_maxspeed=[all_pump_maxspeed;decomposed_data.maxspeed];
        all_rise_fall=[all_rise_fall;decomposed_data.rise_fall];
        all_midpoints=[all_midpoints;decomposed_data.midpoint];
        all_jerk_costs=[all_jerk_costs;decomposed_data.jerk_costs];
        %descrete correlations
        var={'amp','off','hoff','hoff_ctr','amp_off','dur','pr','prot','ret'};
        field={'lag_','X_','diluted_'};
        for k=1:numel(field)
            for j=1:numel(var)
                eval(['all_',field{k},var{j},'=[','all_',field{k},var{j},';decomposed_data.',field{k},var{j},'];']);
            end
        end
        sd(m)=numel(decomposed_data.diluted_amp);
        
        %pump 'raster'
        ind_pumps=decomposed_data.ind_pumps;
        SP=(decomposed_data.maxspeed-decomposed_data.pump_val)./decomposed_data.maxspeed;
        PH=[zeros(size(SP,1),1),ones(size(SP,1),1)]; %0=protraction, 1=retraction
        SP1=reshape(SP',[],1);
        PH1=reshape(PH',[],1);
        PV=reshape(decomposed_data.pump_val',[],1);
        MSP=reshape(decomposed_data.maxspeed',[],1);
        RF=reshape(decomposed_data.rise_fall',[],1);
        str_pumps=SP1(PV~=0);
        phs_pumps=PH1(PV~=0);
        msp_pumps=MSP(PV~=0);
        rf_pumps=RF(PV~=0);
        pdata{i}=[ind_pumps(:) str_pumps(:) phs_pumps(:) msp_pumps(:) rf_pumps(:)]; %time index, strength, phase, maxspeed,rise/fall
        
        ddata=[];
        m=m+1;
    end
end


% pump strength measure: 0=no pump; 0-1=delayed pump; >1=double pump
pumpstr=(all_pump_maxspeed-all_pump_val)./all_pump_maxspeed;
pumpstr(all_pump_val==0)=0;
% pumpstr(all_rise_fal<minamp)=NaN;
% all_jerk_costs(all_rise_fal<minamp)=NaN;
% lags=lags/fs;

%% export
analyzed_sturct.segment_dur=d;
analyzed_sturct.segment_amp=a;
analyzed_sturct.segment_length=s;
analyzed_sturct.segment_whisks=w;
analyzed_sturct.segmented_diluted_length=sd;
analyzed_sturct.segment_velocity=v;
%diluted data
analyzed_sturct.d_amp=all_d_amp;
analyzed_sturct.d_off=all_d_offset;
analyzed_sturct.d_dur=all_d_dur;
analyzed_sturct.d_rat=all_d_rat;

%cont. data
analyzed_sturct.amp_off=[all_amp all_offset];
analyzed_sturct.amp_hoff=[all_amp all_hoffset];
analyzed_sturct.amp_hoff_ctr=[all_amp all_hoffset_ctr];
analyzed_sturct.d_V=all_V;
analyzed_sturct.dur_rat=[all_dur all_rat];
analyzed_sturct.rise_fall=all_rise_fall;
analyzed_sturct.midpoints=all_midpoints;
analyzed_sturct.wamp=all_wamp;

%bout
analyzed_sturct.all_bout_dur=all_bout_dur;
analyzed_sturct.all_bout_whisks=all_bout_whisks;

%pumps
analyzed_sturct.all_pump_val=all_pump_val;
analyzed_sturct.all_pump_phase=all_pump_phase;
analyzed_sturct.all_pump_prom=all_pump_prom;
analyzed_sturct.all_pump_width=all_pump_width;
analyzed_sturct.all_pump_maxspeed=all_pump_maxspeed;
analyzed_sturct.all_pump_stregth=pumpstr;
analyzed_sturct.all_jerk_costs=all_jerk_costs;
%xcov
var={'amp','off','hoff','hoff_ctr','amp_off','dur','pr','prot','ret'};
field={'lag_','X_','diluted_'};
for k=1:numel(field)
    for j=1:numel(var)
        eval(['analyzed_sturct.all_',field{k},var{j},'=all_',field{k},var{j},';']);
    end
end


end
