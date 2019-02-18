% generate structures and datasets

%% extract whisker and head data
% extract_C2L_head_data('comparative','R');
% extract_C2L_head_data('comparative','L');
% extract_C2L_head_data('control','R');
% extract_C2L_head_data('control','L');

%% get whisking stats
load('headfixed_c2h_data.mat');
[hf_struct,pdata]=get_whisking_stats_new(wdata,hdata,0);

load('comparative_c2h_r_data.mat');
[comp_struct,pdata]=get_whisking_stats_new(wdata,hdata,0);
save('comparative_c2h_r_data.mat','pdata','-append')
load('comparative_c2h_data.mat');
[comp_struct,pdata]=get_whisking_stats_new(wdata,hdata,0);
save('comparative_c2h_data.mat','pdata','-append')


load('control_c2h_r_data.mat');
[ctl_struct,pdata]=get_whisking_stats_new(wdata,hdata,0);
save('control_c2h_r_data.mat','pdata','-append')
load('control_c2h_data.mat');
[ctl_struct,pdata]=get_whisking_stats_new(wdata,hdata,0);
save('control_c2h_data.mat','pdata','-append')
%% save to file
save('all_analyzed_structs_head.mat','hf_struct','comp_struct','ctl_struct');