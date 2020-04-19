function extract_C2L_head_data(database,side)
% this function extracts only analyzable whisking sections of 
% whisker C2 on one side from a database. it also removes glitches and downsample to
% 500fps if needed.
% also the head direction and location from all free moving videos
% function saves extracted data into file

%% parameters
K=10;

fs=500;
lp = [200];
lp = lp * 2 / 1e3; % convert Hz to radians/S
[N, Wn] = buttord( lp, lp .* [1.5], 3, 20); 
[B,A] = butter(N,Wn);

headfixed_path='../Data/sAvner_HeadFixed.mat';
comparative_path='../Data/Free_1row';
control_path='../Data/Free_trimmed';

%% load data
if(strcmp(database,'headfixed'))
    load(headfixed_path);
    data=sAvner_HeadFixed;
elseif(strcmp(database,'comparative'))
    k=1;
    D=dir(comparative_path);
    for d=3:numel(D)
        d
        DD=dir([comparative_path,'/',D(d).name]);
        for dd=3:numel(DD) %go over each file
            dd
            load([comparative_path,'/',D(d).name,'/',DD(dd).name]);
            if(numel(g_tMovieInfo.Angle))
                data(k)=g_tMovieInfo;
                k=k+1
            end
        end
    end
elseif(strcmp(database,'control'))
    k=1;
    D=dir(control_path);
    for d=3:numel(D)
        d
        DD=dir([control_path,'/',D(d).name]);
        for dd=3:numel(DD) %go over each file
            dd
            load([control_path,'/',D(d).name,'/',DD(dd).name]);
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

%% take all left side C2 whiskers
t=0;
k=1;
for i=1:length(data);
    ind_left=find(data(i).WhiskerSide==(1+(side=='L')));
    for w=1:length(ind_left)
        if(strcmp(data(i).WhiskerIdentity{ind_left(w)},'C2') | strcmp(data(i).WhiskerIdentity{ind_left(w)},[side,'C3']))    
            dat=data(i).Angle(:,ind_left(w));
            indok=find(~isnan(dat) & dat>0);
            dat=dat(indok);
            if(numel(dat)>500)
               
                %remove glitches
                fdat=medfilt1(dat,K);
                dat=[dat(1);fdat(2:end-1);dat(end)];
                
                if(data(i).FramesPerSecond==1e3)
                    dat= filtfilt(B,A,dat); % zero-phase LPF
                    wdata{k}=decimate(dat,2); %down-sample to 500 hz
                else
                    wdata{k}=dat;
                end

                % head motion
                if(strcmp(database,'control') | strcmp(database,'comparative'))
                    RE=data(i).RightEye;
                    LE=data(i).LeftEye;
                    N=data(i).Nose;
                    indok=find(~isnan(RE(:,1)) & ~isnan(LE(:,1)) & ~isnan(N(:,1)));
                    RE=RE(indok,:);   LE=LE(indok,:);   N=N(indok,:);
                    x=[RE(:,1) LE(:,1) N(:,1)];
                    y=[RE(:,2) LE(:,2) N(:,2)];
                    X=(RE(:,1)+LE(:,1))/2;
                    Y=(RE(:,2)+LE(:,2))/2;
                    alpha=atan2((N(:,2)-Y),(N(:,1)-X));
                    Vx=diff(X)*fs;
                    Vy=diff(Y)*fs;
                    beta=(atan2(Vy,Vx));
                    da=diff(alpha)*fs;
                    da=[da(1);da(:)];
                    gamma=beta-alpha(2:end);


                    hdata{k}=[rad2deg(alpha) X Y];
                    H=sqrt((X-N(:,1)).^2+(Y-N(:,2)).^2);
                else
                    hdata{k}=[];
                    H=0;
                end

                %stats
                session(k)=i;
                duration(k)=length(wdata{k})/500;            
                peak(k)=max(wdata{k});
                trough(k)=min(wdata{k});
                range(k)=peak(k)-trough(k);
                headsize(k)=nanmean(H);

%                 H=plot(t+[1:length(wdata{k})]/500,wdata{k});
%                 set(H,'Color',colormatrix(mod(k,63)+1,:));
%                 hold on;    
%                 if(numel(hdata{k}))
%                     H=plot(t+[1:length(hdata{k}(:,1))]/500,hdata{k}(:,1));
%                     set(H,'Color',colormatrix(mod(k,63)+1,:));
%                     V=sqrt(diff(hdata{k}(:,2)).^2+diff(hdata{k}(:,3)).^2)*500;
%                     H=plot(t+[1:length(V)]/500,V);
%                     set(H,'Color',colormatrix(mod(k,63)+1,:));
%                 end

                text(t,100,num2str(i));                
                t=t+duration(k);
                k=k+1;            
            end
        end
    end
end
if(side=='L')
    suffix='_c2h_data';
else
    suffix='_c2h_r_data';
end
save(['../Data/',database,suffix],'wdata','hdata','session','duration','peak','trough','range','headsize');
end
