function [hdata]=extract_head_data(database)
% this script extracts the head direction and location from all free moving
% videos

%% parameters
colormatrix=colormap(lines);    
if(ismac)
    comparative_path='/Users/avner_wallach/Documents/Whisking Behavioral data/WIS';
    control_path='/Users/avner_wallach/Documents/Whisking Behavioral data/WIS- C2-C1-D1 whiskers';
    sep='/';
else    
    comparative_path='C:\Users\awallach\Documents\comparative\WIS\';
    control_path='C:\Users\awallach\Documents\comparative\WIS- C2-C1-D1 whiskers';
    sep='\';
end
plotmotion=1;
if(plotmotion)
    F=figure;
    A1=subplot(1,2,1);
    A2=subplot(1,2,2);
    
end
%% load data
if(strcmp(database,'comparative'))
    k=1;
    D=dir(comparative_path);
    for d=3:numel(D)
        d
        DD=dir([comparative_path,sep,D(d).name]);
        for dd=3:numel(DD) %go over each file
            dd
            load([comparative_path,sep,D(d).name,sep,DD(dd).name]);
            data(k).RightEye=g_tMovieInfo.RightEye;
            data(k).LeftEye=g_tMovieInfo.LeftEye;
            data(k).Nose=g_tMovieInfo.Nose;
%             
%             data(k)=g_tMovieInfo;
            k=k+1
        end
    end
elseif(strcmp(database,'control'))
    k=1;
    D=dir(control_path);
    for d=3:numel(D)
        d
        DD=dir([control_path,sep,D(d).name]);
        for dd=3:numel(DD) %go over each file
            dd
            load([control_path,sep,D(d).name,sep,DD(dd).name]);
            if(~isfield(g_tMovieInfo,'Curvature'))
                g_tMovieInfo.Curvature=[];
            end
            if(size(g_tMovieInfo.Angle,1)>500) %at least one second
                data(k)=g_tMovieInfo;
                k=k+1
            end
        end
    end    
end   
%% get all eyes and nose data from tracked data
framerate=500;
for i=1:length(data);
    RE=data(i).RightEye;
    LE=data(i).LeftEye;
    N=data(i).Nose;
    x=[RE(:,1) LE(:,1) N(:,1)];
    y=[RE(:,2) LE(:,2) N(:,2)];
    X=(RE(:,1)+LE(:,1))/2;
    Y=(RE(:,2)+LE(:,2))/2;
    alpha=atan2((N(:,2)-Y),(N(:,1)-X)); %head angle
    Vx=diff(X)*framerate;
    Vy=diff(Y)*framerate;
    beta=(atan2(Vy,Vx)); %head translational motion
    da=diff(alpha)*framerate; %head rotation
    da=[da(1);da(:)];
    gamma=beta-alpha(2:end); %head translation egocentric direction
    
    
    hdata{i}=[rad2deg(alpha) X Y];
    R=25;
    U=R*cos(alpha);
    V=R*sin(alpha);
    W=R*cos(beta);
    Z=R*sin(beta);
    K=R*cos(gamma);
    L=R*sin(gamma);
    
    
    if(plotmotion)
        if i==15
        loops=size(N,1)-1;
%         F(loops) = struct('cdata',[],'colormap',[]);
        for k=1:loops
            axes(A1);
            fill(x(k,:),y(k,:),[0.5 0.5 0.5]);
            hold on;
            H=plot(x(k,3),y(k,3),'k.');
            set(H,'MarkerSize',30);
            H=plot(x(k,1),y(k,1),'r.'); % red==right eye
            set(H,'MarkerSize',40);
            H=plot(x(k,2),y(k,2),'b.'); % blue==left eye
            set(H,'MarkerSize',40);
            H=quiver(X(k),Y(k),U(k),V(k));
            H=quiver(X(k),Y(k),W(k),Z(k));
            hold off;
            axis([min(x(:)) max(x(:)) min(y(:)) max(y(:))]);
%             set(A,'YDir','reverse');
            axis('ij');
            axes(A2);
            H=polarplot([gamma(k) gamma(k)],[0 1]);  
            hold on;
            H=polarplot([pi/2 pi/2],[0 da(k)]); 
            hold off;
            A2=gca;
            set(A2,'ThetaDir','clockwise','ThetaZeroLocation','top');
            drawnow;
%             F(i)=getframe;
        end
    end

        
end
end    
    

    