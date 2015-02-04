%Updated 7/6/09 to subtract offsets from I and Q data!  We are very close
%to noise floor, so it is important to do this.
%Updated 7/6/09 to include calculation of mua changes and to work for both
%ISS and Joel's Instruments

close all
clear all

%% Enter patient-specific identifiers here
patID='SVC169';
patdate='060214';
probeoffmark=[999 1]; % Separate pairs of marks where probe was not placed correctly, 999 1 limits to measurement range
baselinemarks=[5 7];
%If you want something added to your file name
%Extension you want added to filename, if nothing, just do ext='';
ext='';

usedISS=1;%0=Joel's Instrument, 1=ISS
usednewsettings=1;%Dennis made new settings file for us to turn off PMTs while DCS laser is on.  
useddet=2;%If using ISS instrument, DetA=1, DetB=2
threelambda=1;%Sometimes we have trouble with 685 nm laser, if this is the case, just use 785 and 830 for DPF calcs
syncmark=1;

% approximate dpf factors for the above wavelengths, AO Vol32(418)-Fig 8a
% 41weeks gestation, dont cite for lower than 720nm
% Healthy neonate DPF values in Duncan paper, Phys Med Biol, Vol 40, 1995,
% p. 295-304.  Her values are ~20% higher
DPF=[4.35 4.2 4.4];

if usedISS
    %ISS Analysis
    sourcesused=[7 8 6]; % which sources, out of 16
    sourcesnotused=[1 2 4 3 5];
    lambdas=[826 688 786];

    %%LOAD DATA
    fdir=['../' patID '/' patID '_' patdate '/'];
    % Must sort the rest
    files=dir([ fdir '*.txt']);
    for numfile=1:size(files,1)
        names(numfile,:)=files(numfile).name(12:end-4);
    end
    names=sortrows(names);

    data=[];
    for i=1:size(names,1)
        %%Find time each file was saved "HH:MM:SS"
        timetmp{i}=[ num2str(names(i,2:3)) ':' num2str(names(i,4:5)) ':' num2str(names(i,6:7)) ];
        
        fname=['../' patID '/' patID '_' patdate '/Data_' patdate names(i,:) '.txt'];
        fid = fopen([ fname ], 'r');
        if usednewsettings
            [tmpdata, count]=fscanf(fid,'%c %s',1010);
            clear tmpdata count
            [tmpdata, count]=fscanf(fid,'%g %g',[79 inf]);
        else
            [tmpdata, count]=fscanf(fid,'%c %s',1683);
            clear tmpdata count
            [tmpdata, count]=fscanf(fid,'%g %g',[151 inf]);
        end
        fclose(fid);
        tmpdata=tmpdata.';
        datalength(i)=size(tmpdata,1)/12;
        data=cat(1,data,tmpdata);

    end

    % The columns of data:
    % Column 1: Time (in seconds) *NOT RELIABLE! DON'T USE FOR TIME AXIS!*
    % Column 2: Data group
    % Column 3: Step (has to do with external GPIB devices)
    % Column 4: Mark
    % Column 5: Flag (turn into binary for message about data)
    % Columns 6-21: AC data for Detector A, Sources 1-16
    % Columns 22-37: DC data for Detector A, Sources 1-16
    % Columns 38-53: Phase data for Detector A, Sources 1-16
    % Columns 54-69: AC data for Detector B, Sources 1-16
    % Columns 70-85: DC data for Detector B, Sources 1-16
    % Columns 86-101: Phase data for Detector B, Sources 1-16
    % Columns 102: Analog output data 1
    % Column 103: Digital output data

    if usednewsettings
        updatetime=1/6.3776;
    else
        updatetime=1/6.25; % Update time, set to 6.25 Hz.
    end
    Marks=find(data(:,4));
    Marks=floor(Marks/12);
    
    ISStime4frame=datenum(timetmp,'HH:MM:SS')-floor(datenum(timetmp,'HH:MM:SS'));%In arbitrary units--1a.u.=24hrs, counting from 1/1/2000.
    ISStime4frame(find(datalength==0))=[];
    datalength(find(datalength==0))=[];
    %Subtract off whole days since year 2000 and you are left with a fraction
    %of a day that has elapsed since the start of this data
    
    %Zero time axis so that t=0 is first frame
    ISStime4frame=(ISStime4frame-ISStime4frame(1))*1.44e3;%Convert time in a.u. to minutes
    
    %Now, we only have time points for when each file was saved.  However,
    %each file has approximately 5 frames.  So we will interpolate the time
    %points for the frames within each file.
    
    %%Get frame number that each ISS time corresponds to:
    for i=1:length(ISStime4frame)
        frame(i+1)=sum(datalength(1:i));
    end
    %Since first time point corresponds to the first file, we need to add a fake file in the
    %beginning in order for interpolation to catch frames 1-5
    frame(1)=0;
    ISStime4frame(2:end+1)=ISStime4frame(1:end);
    ISStime4frame(1)=ISStime4frame(2)-(ISStime4frame(3)-ISStime4frame(2));
    
    %%Interpolate times for the rest of the frames
    %First find index of jumps in time (we dont want to interpolate in the
    %jumps)
    ind(1)=1;
    ind(2:length(find(diff(ISStime4frame)>3))+1)=find(diff(ISStime4frame)>3);
    ind(end+1)=length(ISStime4frame);
    for j=1:size(ind,2)-1
        %Define chunk of continuous data without jump in time
        chunk{j}=ind(j):1:ind(j+1);
        if j>1
            ISStime4frame(ind(j))=ISStime4frame(ind(j)+1)-(ISStime4frame(ind(j)+2)-ISStime4frame(ind(j)+1));
        end
        %Interpolate time data within this chunk
        ISStime(min(frame(chunk{j}))+1:1:max(frame(chunk{j}))+1)=interp1(frame(chunk{j}),ISStime4frame(chunk{j}),min(frame(chunk{j})):1:max(frame(chunk{j})));
    end
    ISStimetmp=ISStime(2:end);%Ignore first point, just used it to interpolate 
    clear ind ISStime
    ISStime=ISStimetmp;%Now we have timepoints for each datapoint!!!
    
    %%Set t=0 to the sync mark
    ISStime=ISStime-ISStime(Marks(syncmark));
    
    %%Plot times just to make sure they are accurate
    figure,plot(1:length(ISStime),ISStime,'r.','MarkerSize',30)
    hold on,plot(frame,ISStime4frame,'b.','MarkerSize',25)
    axis tight
    xlabel('Frame')
    ylabel('Time (min)')
    set(findall(gcf,'-property','FontSize'),'FontSize',25)
    legend('Time Used','Raw Frame Time',2)
    set(gcf,'PaperPositionMode','Auto')
    maxwindows(gcf)
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/ISSTimes_' patID '_' patdate '.fig'],'fig')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/ISSTimes_' patID '_' patdate '.jpg'],'jpg')
    
    %% Plot noncalibrated AC, DC and phase data for Detectors A & B
    for j=1:length(sourcesused)
        if useddet==1
            ind=5;
        elseif useddet==2
            ind=29;
        end
        ACtmp(:,j)=data(:,ind+sourcesused(j));
        DCtmp(:,j)=data(:,ind+8+sourcesused(j));
        phasetmp(:,j)=data(:,ind+8*2+sourcesused(j));
        %Convert phase to radians
        phasetmp(:,j)=phasetmp(:,j)*pi/180;
        phasetmp(:,j)=unwrap(phasetmp(:,j));
        %Convert back to degrees
        phasetmp(:,j)=phasetmp(:,j)*180/pi;
    end
    
    %%Find data for unused channels
    for j=1:length(sourcesnotused)
        if useddet==1
            ind=5;
        elseif useddet==2
            ind=29;
        end
        ACtmp_offsource(:,j)=data(:,ind+sourcesnotused(j));
        DCtmp_offsource(:,j)=data(:,ind+8+sourcesnotused(j));
        
    end
    
    %Find the mean and standard deviation of each frame of data
    for i=1:length(data(:,2))/12
        AC(i,:)=mean(ACtmp((i-1)*12+1:12*i,:),1);
        DC(i,:)=mean(DCtmp((i-1)*12+1:12*i,:),1);
        
        ACstd(i,:)=std(ACtmp((i-1)*12+1:12*i,:),0,1);
        DCstd(i,:)=std(DCtmp((i-1)*12+1:12*i,:),0,1);
        
        AC_offsource(i,:)=mean(ACtmp_offsource((i-1)*12+1:12*i,:),1);
        DC_offsource(i,:)=mean(DCtmp_offsource((i-1)*12+1:12*i,:),1);
        
        ACstd_offsource(i,:)=std(ACtmp_offsource((i-1)*12+1:12*i,:),0,1);
        DCstd_offsource(i,:)=std(DCtmp_offsource((i-1)*12+1:12*i,:),0,1);
        
    end
    
    figure,subplot(2,1,1)
    plot(AC)
    %xlabel('Frame','FontSize',20)
    ylabel('Amplitude','FontSize',20)
    axis tight
    grid on
    set(gcf,'PaperPositionMode','Auto')
    set(gca,'FontSize',20)
    tmplim=get(gca,'YLim');
    for kkkk=1:length(Marks)
        if isempty(find(baselinemarks==kkkk))
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','b','LineWidth',1)
        else
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','k','LineWidth',2)
        end
        text(Marks(kkkk),tmplim(2)-5,num2str(kkkk),'Color',[0 0 0],'FontSize',16)
    end
    subplot(2,1,2)
    plot(AC_offsource)
    xlabel('Frame','FontSize',20)
    ylabel('Amp (unused sources)','FontSize',20)
    axis tight
    ylim([tmplim])
    grid on
    set(gcf,'PaperPositionMode','Auto')
    set(gca,'FontSize',20)
    for kkkk=1:length(Marks)
        if isempty(find(baselinemarks==kkkk))
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','b','LineWidth',1)
        else
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','k','LineWidth',2)
        end
        
    end
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Amplitude_' patID '_' patdate ext '.fig'],'fig')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Amplitude_' patID '_' patdate ext '.jpg'],'jpg')
    
    figure,subplot(2,1,1)
    plot(DC)
    ylabel('DC Amplitude','FontSize',20)
    axis tight
    grid on
    set(gcf,'PaperPositionMode','Auto')
    set(gca,'FontSize',20)
    tmplim=get(gca,'YLim');
    for kkkk=1:length(Marks)
        if isempty(find(baselinemarks==kkkk))
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','b','LineWidth',1)
        else
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','k','LineWidth',2)
        end
        text(Marks(kkkk),tmplim(2)-5,num2str(kkkk),'Color',[0 0 0],'FontSize',16)
    end
    subplot(2,1,2)
    plot(DC_offsource)
    xlabel('Frame','FontSize',20)
    ylabel('DC Amp (unused sources)','FontSize',20)
    axis tight
    ylim([tmplim])
    grid on
    set(gcf,'PaperPositionMode','Auto')
    set(gca,'FontSize',20)
    for kkkk=1:length(Marks)
        if isempty(find(baselinemarks==kkkk))
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','b','LineWidth',1)
        else
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','k','LineWidth',2)
        end
        
    end
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DCAmplitude_' patID '_' patdate ext '.fig'],'fig')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DCAmplitude_' patID '_' patdate ext '.jpg'],'jpg')
    
    %load the probe map, sdlist etc.
    plott=0;
    flatbabyprobe;%Loads s-d seps in cm
    
    %Every time probe moves, start new baselines
    for k=1:size(baselinemarks,1)
        baselineoxy=Marks(baselinemarks(k,1)):Marks(baselinemarks(k,2)); %baseline frames
        if k==1
            if threelambda
                [hbseries,hbo2series,muaseries]=threewavelengthdpf(squeeze(AC(:,1)),squeeze(AC(:,2)),squeeze(AC(:,3)),...
                    mean(squeeze(AC(baselineoxy,1))),mean(squeeze(AC(baselineoxy,2))),mean(squeeze(AC(baselineoxy,3))),...
                    lambdas,sdlist(1,1),DPF);
            else
                [hbseries,hbo2series,muaseries]=twowavelengthdpf(squeeze(AC(:,1)),squeeze(AC(:,2)),...
                    mean(squeeze(AC(baselineoxy,1))),mean(squeeze(AC(baselineoxy,2))),...
                    lambdas(1:2),sdlist(1,1),DPF(1:2));
            end
        else
            if threelambda
                [hbseries(baselineoxy(1):length(AC)),hbo2series(baselineoxy(1):length(AC)),muaseries(:,baselineoxy(1):length(AC))]=threewavelengthdpf(squeeze(AC(baselineoxy(1):end,1)),squeeze(AC(baselineoxy(1):end,2)),squeeze(AC(baselineoxy(1):end,3)),...
                    mean(squeeze(AC(baselineoxy,1))),mean(squeeze(AC(baselineoxy,2))),mean(squeeze(AC(baselineoxy,3))),...
                    lambdas,sdlist(1,1),DPF);
            else
                [hbseries(baselineoxy(1):length(AC)),hbo2series(baselineoxy(1):length(AC)),muaseries(:,baselineoxy(1):length(AC))]=twowavelengthdpf(squeeze(AC(baselineoxy(1):end,1)),squeeze(AC(baselineoxy(1):end,2)),...
                    mean(squeeze(AC(baselineoxy,1))),mean(squeeze(AC(baselineoxy,2))),...
                    lambdas(1:2),sdlist(1,1),DPF(1:2));
            end

        end
        clear baselineoxy
    end
    %If probe was off head, set values to NaN
    for m=1:size(probeoffmark,1)
        if probeoffmark(m,1)==999
            hbseries(1:Marks(probeoffmark(m,2)))=NaN;
            hbo2series(1:Marks(probeoffmark(m,2)))=NaN;
            muaseries(:,1:Marks(probeoffmark(m,2)))=NaN;
        
        elseif probeoffmark(m,2)==999
            hbseries(Marks(probeoffmark(m,1)):end)=NaN;
            hbo2series(Marks(probeoffmark(m,1)):end)=NaN;
            muaseries(:,Marks(probeoffmark(m,1)):end)=NaN;
            
        else
            hbseries(Marks(probeoffmark(m,1)):Marks(probeoffmark(m,2)))=NaN;
            hbo2series(Marks(probeoffmark(m,1)):Marks(probeoffmark(m,2)))=NaN;
            muaseries(:,Marks(probeoffmark(m,1)):Marks(probeoffmark(m,2)))=NaN;
        end
    end
    
   
    %If the change in mua goes below -0.09 or above 0.09, this doesnt seem
    %realistic, so ignore these data points
    ind=find(muaseries>0.09);
    muaseries(ind)=NaN;
    clear ind 
    ind=find(muaseries<-0.09);
    muaseries(ind)=NaN;
    
    figure,plot((hbo2series+hbseries).*1000,'.-g'), hold on
    plot(hbseries.*1000,'.-r'), hold on
    plot(hbo2series.*1000,'.-k'), hold on
    legend({['THC'],['Hb'],['HbO_{2}']})
    tmplim=get(gca,'YLim');
    for kkkk=1:length(Marks)
        if isempty(find(baselinemarks==kkkk))
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','b','LineWidth',1)
        else
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','k','LineWidth',2)
        end
    end
    grid on
    if threelambda
        title('3 \lambda DPF')
    else
        title('2 \lambda DPF')
    end
    maxwindows(gcf)
    %Plot mua changes
    figure,plot(squeeze(muaseries(1,:)),'b.-','MarkerSize',25,'LineWidth',3)
    hold on,plot(squeeze(muaseries(2,:)),'r.-','MarkerSize',25,'LineWidth',3)
    if threelambda
        hold on,plot(squeeze(muaseries(3,:)),'.-','Color',[0 0.5 0],'MarkerSize',25,'LineWidth',3)
        legend({[ num2str(lambdas(1)) 'nm'],[ num2str(lambdas(2)) 'nm'],[ num2str(lambdas(3)) 'nm']})
    else
        legend({[ num2str(lambdas(1)) 'nm'],[ num2str(lambdas(2)) 'nm']})
    end
    axis tight
    xlabel('Frame','FontSize',35)
    ylabel('\Delta\mu_{a} (cm^{-1})','FontSize',35)
    set(gca,'FontSize',35)
    set(gcf,'PaperPositionMode','Auto')
    grid on
    ylim([-0.1 0.1])
    tmplim=get(gca,'YLim');
    for kkkk=1:length(Marks)
        if isempty(find(baselinemarks==kkkk))
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','b','LineWidth',1)
        else
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','k','LineWidth',2)
        end
    end
    maxwindows(gcf)
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/mua_' patID '_' patdate ext '.fig'],'fig')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/mua_' patID '_' patdate ext '.jpg'],'jpg')
    tt=['save ' patID '_' patdate ext '_dpfout.mat ISStime AC ACstd Marks hbo2series hbseries muaseries useddet DPF usedISS threelambda probeoffmark baselinemarks'];
    eval(tt);

else
    %Analysis for Joel's Instrument
    numsources=1+1;%1 O2 source + 1 for offset measurements
    numlambda=3;

    lambdas=[830 786 686];
    usedlambda=[1 2 3];
    lambdas=lambdas(usedlambda);
    numdets=4; %Even though we are only using 1 of the 4 dets

    measnum=0;
    numfiles=0;
    morefiles=1; %1:more files to read, 0:quit

    while (morefiles)

        fnametmp=[fdir fname sprintf('%01d',numfiles) '.dat'];

        if exist(fnametmp)==2 %file exists

            tmp=readdpdw081704(fnametmp);
            fg(:,measnum+1:(measnum+size(tmp,2)))=tmp;
            measnum=measnum+size(tmp,2);
            numfiles=numfiles+1;
        elseif exist(fname)==0 %file does not exist
            morefiles=0;
        end

    end

    %get only the amplitudes for now
    fga=fg(3:5:22,:);
    %get rid of half-frame stops
    fga=fga(:,1:(size(fga,2)-rem(size(fga,2),numlambda.*(numsources))));
    %reshape it as dets wavelength sources frames
    fg2a=reshape(fga,[numdets numlambda (numsources) prod(size(fga))./(numdets.*(numsources).*numlambda)]);

    %get phases
    fgp=fg(4:5:23,:);
    %get rid of half-frame stops
    fgp=fgp(:,1:(size(fgp,2)-rem(size(fgp,2),numlambda.*(numsources))));
    %reshape it as dets wavelength sources frames
    fg2p=reshape(fgp,[numdets numlambda (numsources) prod(size(fgp))./(numdets.*(numsources).*numlambda)]);

    %get I
    fgI=fg(5:5:22,:);
    %get rid of half-frame stops
    fgI=fgI(:,1:(size(fgI,2)-rem(size(fgI,2),numlambda.*(numsources))));
    %reshape it as dets wavelength sources frames
    fg2I=reshape(fgI,[numdets numlambda (numsources) prod(size(fgI))./(numdets.*(numsources).*numlambda)]);

    %get Q
    fgQ=fg(6:5:22,:);
    %get rid of half-frame stops
    fgQ=fgQ(:,1:(size(fgQ,2)-rem(size(fgQ,2),numlambda.*(numsources))));
    %reshape it as dets wavelength sources frames
    fg2Q=reshape(fgQ,[numdets numlambda (numsources) prod(size(fgQ))./(numdets.*(numsources).*numlambda)]);

    sources=fg(1,:);
    %get rid of half-frame stops
    sources=sources(:,1:(size(sources,2)-rem(size(sources,2),numlambda.*(numsources))));
    %reshape it as dets wavelength sources frames
    sources2=reshape(sources,[1 numlambda (numsources) prod(size(sources))./(1.*(numsources).*numlambda)]);

    %get the Marks
    Marks=ceil((find(fg(end,:)>0))./numlambda./(numsources));

    %%%%%%%%%%%%%

    %Subtract offsets from data (source 2=offsets)
    for i=1:length(fg2I)
        fg2I(:,:,1,i)=fg2I(:,:,1,i)-fg2I(:,:,2,i);
        fg2Q(:,:,1,i)=fg2Q(:,:,1,i)-fg2Q(:,:,2,i);
    end
    %Calculate amp and phase from offset subtracted data
    amp=sqrt(fg2I.^2+fg2Q.^2);
    pha=atan(fg2Q./fg2I);

    %load the probe map, sdlist etc.
    plott=0;
    flatbabyprobe;%Loads s-d seps in cm
    for p=1:size(probeoffmark,2)
        amp(:,:,:,Marks(probeoffmark(p,1):probeoffmark(p,2)))=NaN;
    end
    %Plot amplitudes for reference
    figure,plot(squeeze(amp(useddet,1,1,:))*1000,'.-r','MarkerSize',30,'LineWidth',3)
    hold on,plot(squeeze(amp(useddet,2,1,:))*1000,'.-b','MarkerSize',30,'LineWidth',3)
    hold on,plot(squeeze(amp(useddet,3,1,:))*1000,'.-g','MarkerSize',30,'LineWidth',3)
    legend({['830nm'],['785nm'],['685nm']},'FontSize',24)
    ylabel('Amplitude (mV)','FontSize',24)
    xlabel('Frame','FontSize',24)
    set(gca,'FontSize',24)
    ylim([0 20])
    maxwindows(gcf)
    tmplim=get(gca,'YLim');
    for kkkk=1:length(Marks)
        if isempty(find(baselinemarks==kkkk))
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','b','LineWidth',1)
        else
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','k','LineWidth',2)
        end
        text(Marks(kkkk),tmplim(2)-5,num2str(kkkk),'Color',[0 0 0],'FontSize',16)
    end
    set(gcf,'PaperPositionMode','Auto')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Amp_' patID '_' patdate ext '.fig'],'fig')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Amp_' patID '_' patdate ext '.jpg'],'jpg')
    fig1=figure;

    numdets=4;
    ct=1;
    for kk=1:numdets
        for kkk=1:numsources-1
            %Reset baseline every time probe moves
            for k=1:size(baselinemarks,1)
                baselineoxy=Marks(baselinemarks(k,1)):Marks(baselinemarks(k,2)); %baseline frames
                if k==1
                    if threelambda
                        [hbseriestmp,hbo2seriestmp,muaseriestmp]=threewavelengthdpf(squeeze(amp(kk,1,kkk,:)),squeeze(amp(kk,2,kkk,:)),squeeze(amp(kk,3,kkk,:)),...
                            mean(squeeze(amp(kk,1,kkk,baselineoxy))),mean(squeeze(amp(kk,2,kkk,baselineoxy))),mean(squeeze(amp(kk,3,kkk,baselineoxy))),...
                            lambdas,sdlist(kkk,kk),DPF);
                    else
                        [hbseriestmp,hbo2seriestmp,muaseriestmp]=twowavelengthdpf(squeeze(amp(kk,1,kkk,:)),squeeze(amp(kk,2,kkk,:)),...
                            mean(squeeze(amp(kk,1,kkk,baselineoxy))),mean(squeeze(amp(kk,2,kkk,baselineoxy))),...
                            lambdas(1:2),sdlist(kkk,kk),DPF(1:2));
                    end
                else
                    if threelambda
                        [hbseriestmp(baselineoxy(1):length(amp)),hbo2seriestmp(baselineoxy(1):length(amp)),muaseriestmp(:,baselineoxy(1):length(amp))]=...
                            threewavelengthdpf(squeeze(amp(kk,1,kkk,baselineoxy(1):end)),squeeze(amp(kk,2,kkk,baselineoxy(1):end)),...
                            squeeze(amp(kk,3,kkk,baselineoxy(1):end)),mean(squeeze(amp(kk,1,kkk,baselineoxy))),mean(squeeze(amp(kk,2,kkk,baselineoxy))),...
                            mean(squeeze(amp(kk,3,kkk,baselineoxy))),lambdas,sdlist(kkk,kk),DPF);
                    else
                        [hbseriestmp(baselineoxy(1):length(amp)),hbo2seriestmp(baselineoxy(1):length(amp)),muaseriestmp(:,baselineoxy(1):length(amp))]=...
                            twowavelengthdpf(squeeze(amp(kk,1,kkk,baselineoxy(1):end)),squeeze(amp(kk,2,kkk,baselineoxy(1):end)),...
                            mean(squeeze(amp(kk,1,kkk,baselineoxy))),mean(squeeze(amp(kk,2,kkk,baselineoxy))),...
                            lambdas(1:2),sdlist(kkk,kk),DPF(1:2));
                    end
                end
                clear baselineoxy
            end

            figure(fig1),subplot(numdets,numsources-1,ct),plot((hbo2seriestmp+hbseriestmp).*1000,'.g'), hold on
            figure(fig1),subplot(numdets,numsources-1,ct),plot(hbseriestmp.*1000,'.r'), hold on
            figure(fig1),subplot(numdets,numsources-1,ct),plot(hbo2seriestmp.*1000,'.k'), hold on
            legend({['THC'],['Hb'],['HbO_{2}']})
            hbo2series(kk,kkk,:)=hbo2seriestmp;
            hbseries(kk,kkk,:)=hbseriestmp;
            muaseries(kk,kkk,:,:)=muaseriestmp;
            ylim([-30 30])
            tmplim=get(gca,'YLim');
            for kkkk=1:length(Marks)
                if isempty(find(baselinemarks==kkkk))
                    line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','b','LineWidth',1)
                else
                    line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','k','LineWidth',2)
                end
                
            end
            maxwindows(gcf)

            ct=ct+1;
        end
    end

    for m=1:size(probeoffmark,1)
        if probeoffmark(m,1)==999
            hbseries(:,:,1:Marks(probeoffmark(m,2)))=NaN;
            hbo2series(:,:,1:Marks(probeoffmark(m,2)))=NaN;
            muaseries(:,:,:,1:Marks(probeoffmark(m,2)))=NaN;
            
        elseif probeoffmark(m,2)==999
            hbseries(:,:,Marks(probeoffmark(m,1)):end)=NaN;
            hbo2series(:,:,Marks(probeoffmark(m,1)):end)=NaN;
            muaseries(:,:,:,Marks(probeoffmark(m,1)):end)=NaN;
            
        else
            hbseries(:,:,Marks(probeoffmark(m,1)):Marks(probeoffmark(m,2)))=NaN;
            hbo2series(:,:,Marks(probeoffmark(m,1)):Marks(probeoffmark(m,2)))=NaN;
            muaseries(:,:,:,Marks(probeoffmark(m,1)):Marks(probeoffmark(m,2)))=NaN;

        end
    end
        %If the change in mua goes below -0.09 or above 0.09, this doesnt seem
    %realistic, so ignore these data points
    ind=find(muaseries>0.09);
    muaseries(ind)=NaN;
    clear ind 
    ind=find(muaseries<-0.09);
    muaseries(ind)=NaN;
    clear ind
    %Plot mua changes
    figure,plot(squeeze(muaseries(useddet,1,1,:)),'b.-','MarkerSize',25,'LineWidth',3)
    hold on,plot(squeeze(muaseries(useddet,1,2,:)),'r.-','MarkerSize',25,'LineWidth',3)
    if threelambda==1
        hold on,plot(squeeze(muaseries(useddet,1,3,:)),'.-','Color',[0 0.5 0],'MarkerSize',25,'LineWidth',3)
        legend({[ num2str(lambdas(1)) 'nm'],[ num2str(lambdas(2)) 'nm'],[ num2str(lambdas(3)) 'nm']})
    else
        legend({[ num2str(lambdas(1)) 'nm'],[ num2str(lambdas(2)) 'nm']})
    end
    axis tight
    xlabel('Frame','FontSize',35)
    ylabel('\Delta\mu_{a} (cm^{-1})','FontSize',35)
    set(gca,'FontSize',35)
    set(gcf,'PaperPositionMode','Auto')
    grid on
    ylim([-0.1 0.1])
    tmplim=get(gca,'YLim');
    for kkkk=1:length(Marks)
        if isempty(find(baselinemarks==kkkk))
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','b','LineWidth',1)
        else
            line([Marks(kkkk) Marks(kkkk)],[tmplim(1) tmplim(2)],'Color','k','LineWidth',2)
        end
    end
    maxwindows(gcf)
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/mua_' patID '_' patdate ext '.fig'],'fig')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/mua_' patID '_' patdate ext '.jpg'],'jpg')
    tt=['save ' patID '_' patdate ext '_dpfout.mat amp Marks hbo2series hbseries muaseries useddet DPF threelambda probeoffmark baselinemarks usedISS'];
    eval(tt);

end



