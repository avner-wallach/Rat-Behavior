function my_plotWithConfQ(x,Q,colorvec)
y=Q(2,:)';
hiconf=Q(3,:)';
loconf = Q(1,:)';

% check inputs
if length(x)~=length(y)
    error('length(x)~=length(y)');
end

if length(x)~=length(hiconf)
    error('length(x)~=length(hiconf)');
end

x=reshape(x,[],1);

hold('all');

ind=find(isnan(y)==0);
area1 = area(x(ind),[loconf(ind) -loconf(ind)+hiconf(ind)],'LineStyle','none');
set(area1(1),'FaceColor','none');
% set(area1(2),'FaceColor',1-(1-colorvec)/2);
set(area1(2),'FaceColor',colorvec);
set(area1(2),'FaceAlpha',0.5);

% Create plot
plot(x(ind),y(ind),'LineWidth',2,'Color',colorvec);

xlabel('x');
ylabel('y');

% axis([min(x) max(x) min(y-max(loconf)) max(y+max(hiconf))])