function extract_C2L_data(database)
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

%% take all left side C2 whiskers
t=0;
k=1;
for i=1:length(data);
%     if(~strcmp(data(i).Strain,'Wis'))
%         continue;
%     end
    ind_left=find(data(i).WhiskerSide==2);
    for w=1:length(ind_left)
        if(strcmp(data(i).WhiskerIdentity{ind_left(w)},'C2') | strcmp(data(i).WhiskerIdentity{ind_left(w)},'LC3'))    
            dat=data(i).Angle(:,ind_left(w));
            dat=dat(~isnan(dat) & dat>0);
            if(numel(dat)>500)
               
                %remove glitches
%                 indg=find(abs(wdata{k}-smooth(wdata{k},K))>glitch_th);
%                 wdata{k}(indg)=(wdata{k}(indg-1)+wdata{k}(indg+1))/2;
                fdat=medfilt1(dat,K);
                dat=[dat(1);fdat(2:end-1);dat(end)];
                

                if(data(i).FramesPerSecond==1e3)
                    dat= filtfilt(B,A,dat); % zero-phase LPF
                    wdata{k}=decimate(dat,2); %down-sample to 500 hz
                else
                    wdata{k}=dat;
                end

                %filter
%                 wdata{k}=filtfilt(Hd500.coeffs.Numerator,[1],wdata{k});
                session(k)=i;
                duration(k)=length(wdata{k})/500;            
                peak(k)=max(wdata{k});
                trough(k)=min(wdata{k});
                range(k)=peak(k)-trough(k);

                H=plot(t+[1:length(wdata{k})]/500,wdata{k});
                set(H,'Color',colormatrix(mod(k,63)+1,:));
                text(t,100,num2str(i));
                t=t+duration(k);
                hold on;            

                k=k+1;            
            end
        end
    end
end

save([database,'_c2_data'],'wdata','session','duration','peak','trough','range');
end
