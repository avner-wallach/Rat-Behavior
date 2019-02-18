function create_db()
load('all_analyzed_structs_head.mat');
% load('headfixed_c2_data');
% hdata=cell(size(wdata));
% hf_struct=get_whisking_stats_new(wdata,hdata,0)
load('comparative_c2h_data');
comp_struct=get_whisking_stats_new(wdata,hdata,0);
load('control_c2h_data');
ctl_struct=get_whisking_stats_new(wdata,hdata,0);
save('all_analyzed_structs_head.mat','hf_struct','comp_struct','ctl_struct');
end

