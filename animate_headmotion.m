function animate_headmotion(LE,RE,N,hdat)
alpha=hdat(:,1)*pi/180;
X=hdat(:,2);
Y=hdat(:,3);
x=[RE(:,1) LE(:,1) N(:,1)];
y=[RE(:,2) LE(:,2) N(:,2)];
R=25;
U=R*cos(alpha);
V=R*sin(alpha);
F=figure;
A=axes;
loops=size(N,1);
%         F(loops) = struct('cdata',[],'colormap',[]);
for k=1:loops
    fill(x(k,:),y(k,:),[0.5 0.5 0.5]);
    hold on;
    H=plot(x(k,3),y(k,3),'k.');
    set(H,'MarkerSize',30);
    H=plot(x(k,1),y(k,1),'r.'); % red==right eye
    set(H,'MarkerSize',40);
    H=plot(x(k,2),y(k,2),'b.'); % blue==left eye
    set(H,'MarkerSize',40);
    H=quiver(X(k),Y(k),U(k),V(k));
    hold off;
    axis([min(x(:)) max(x(:)) min(y(:)) max(y(:))]);
%     set(A,'YDir','reverse');
    drawnow;
%             F(i)=getframe;
end
end
