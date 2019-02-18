date=[];
file=[];
m=1;
for i=1:numel(su_struct_total)
    if(~isfield(su_struct_total{i},'generated'))
        continue;
    end
    if(~isfield(su_struct_total{i}.generated,'wsig'))
        continue;
    end
    if(su_struct_total{i}.date==date & su_struct_total{i}.gen_datafile==file)
        continue;
    end
    date=su_struct_total{i}.date;
    file=su_struct_total{i}.gen_datafile;    
    wdata{m}=su_struct_total{i}.generated.wsig+60;
    m=m+1;
end

