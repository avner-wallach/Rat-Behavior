function my_plotWithConf(x,y,hiconf,colorvec)

% check inputs
if length(x)~=length(y)
    error('length(x)~=length(y)');
end

if length(x)~=length(hiconf)
    error('length(x)~=length(hiconf)');
end

x=reshape(x,[],1);
y=reshape(y,[],1);
hiconf=reshape(hiconf,[],1);
loconf = hiconf;

hold('all');

ind=find(isnan(y)==0);
area1 = area(x(ind),[y(ind)-loconf(ind) loconf(ind)+hiconf(ind)],'LineStyle','none');
set(area1(1),'FaceColor','none');
set(area1(2),'FaceColor',1-(1-colorvec)/2);

% Create plot
plot(x(ind),y(ind),'LineWidth',2,'Color',colorvec);

xlabel('x');
ylabel('y');

% axis([min(x) max(x) min(y-max(loconf)) max(y+max(hiconf))])