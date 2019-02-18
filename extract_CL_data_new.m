function all_angle=extract_CL_data_new(database)
% this script extracts only analyzable whisking sections of left side
% whisker C2 from a database. it also removes glitches and downsample to
% 500fps if needed.

%% parameters
K=10;
% glitch_th=15;

fs=500;
lp = [200];
lp = lp * 2 / 1e3; % convert Hz to radians/S
[N, Wn] = buttord( lp, lp .* [1.5], 3, 20); 
[B,A] = butter(N,Wn);

colormatrix=colormap(lines);    

headfixed_path='C:\Users\awallach\Dropbox\postdoc\work\new behaving analysis\sAvner_HeadFixed.mat';
comparative_path='C:\Users\awallach\Documents\comparative\WIS\';
control_path='C:\Users\awallach\Documents\comparative\WIS- C2-C1-D1 whiskers';
%% load data
if(strcmp(database,'headfixed'))
    load(headfixed_path);
    data=sAvner_HeadFixed;
elseif(strcmp(database,'comparative'))
    k=1;
    D=dir(comparative_path);
    for d=3:numel(D)
        d
        DD=dir([comparative_path,'\',D(d).name]);
        for dd=3:numel(DD) %go over each file
            dd
            load([comparative_path,'\',D(d).name,'\',DD(dd).name]);
            data(k)=g_tMovieInfo;
            k=k+1
        end
    end
elseif(strcmp(database,'control'))
    k=1;
    D=dir(control_path);
    for d=3:numel(D)
        d
        DD=dir([control_path,'\',D(d).name]);
        for dd=3:numel(DD) %go over each file
            dd
            load([control_path,'\',D(d).name,'\',DD(dd).name]);
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

%% take all left side whiskers
t=0;
k=1;
all_angle{1}=[];
all_angle{2}=[];
all_angle{3}=[];
for i=1:length(data);
%     if(~strcmp(data(i).Strain,'Wis'))
%         continue;
%     end
    ind_left=find(data(i).WhiskerSide==2);
    for w=1:length(ind_left)
        % in control: LC1=D1;LC2=C1;LC3=C2
        if(strcmp(data(i).WhiskerIdentity{ind_left(w)},'C2')...
                | strcmp(data(i).WhiskerIdentity{ind_left(w)},'LC3'))    
            indw=2;
        elseif(strcmp(data(i).WhiskerIdentity{ind_left(w)},'C1')...
                | strcmp(data(i).WhiskerIdentity{ind_left(w)},'LC2'))
            indw=1;
        elseif(strcmp(data(i).WhiskerIdentity{ind_left(w)},'C3'))
            indw=3;
        else
            continue;
        end
        dat=data(i).Angle(:,ind_left(w));
        dat=dat(~isnan(dat) & dat>0);
        if(numel(dat)>500)

            %remove glitches
            fdat=medfilt1(dat,K);
            dat=[dat(1);fdat(2:end-1);dat(end)];


            if(data(i).FramesPerSecond==1e3)
                dat= filtfilt(B,A,dat); % zero-phase LPF
                all_angle{indw}=[all_angle{indw};decimate(dat,2)]; %down-sample to 500 hz
            else
                all_angle{indw}=[all_angle{indw};dat];
            end

%             H=plot(t+[1:length(wdata{k})]/500,wdata{k});
%             set(H,'Color',colormatrix(mod(k,63)+1,:));
%             text(t,100,num2str(i));
%             t=t+duration(k);
%             hold on;            
% 
%             k=k+1;            
        end
    end
    end
end

