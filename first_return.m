function [c,p]=first_return(data)
data=data(~isnan(data));
% figure;
% plot(data(1:end-1),data(2:end),'.');
% hold on;
% plot([min(data) max(data)],[min(data) max(data)],':k');
% axis([min(data) max(data) min(data) max(data)]);
[c,p]=corr(data(1:end-1)',data(2:end)')
end
