% generate structures and datasets

%% extract whisker and head data
extract_C2L_head_data('headfixed','L');
extract_C2L_head_data('comparative','R');
extract_C2L_head_data('comparative','L');
extract_C2L_head_data('control','R');
extract_C2L_head_data('control','L');

%% get whisking stats
load('../Data/headfixed_c2h_data.mat');
[hf_struct,pdata]=get_whisking_stats_new(wdata,hdata,10,0);

load('../Data/comparative_c2h_r_data.mat');
[comp_struct,pdata]=get_whisking_stats_new(wdata,hdata,headsize,0);
save('../Data/comparative_c2h_r_data.mat','pdata','-append')
load('../Data/comparative_c2h_data.mat');
[comp_struct,pdata]=get_whisking_stats_new(wdata,hdata,headsize,0);
save('../Data/comparative_c2h_data.mat','pdata','-append')


load('../Data/control_c2h_r_data.mat');
[ctl_struct,pdata]=get_whisking_stats_new(wdata,hdata,headsize,0);
save('../Data/control_c2h_r_data.mat','pdata','-append')
load('../Data/control_c2h_data.mat');
[ctl_struct,pdata]=get_whisking_stats_new(wdata,hdata,headsize,0);
save('../Data/control_c2h_data.mat','pdata','-append')
%% save to file
save('../Data/all_analyzed_structs_head.mat','hf_struct','comp_struct','ctl_struct');