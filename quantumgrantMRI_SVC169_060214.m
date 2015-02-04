%%Time-stamp: "2007-05-15 15:14:13 matlabuser"

%%% EDITED TO MAKE A PLOT FOR QUANTUM GRANT
%% AFTER PLOTTING THE GRAPH , IT LETS ME TAKE POINTS FOR A BAR PLOT

close all
clear all
    
subjectid='SVC169';
patdate='060214';
ext='';

%Copy checkfitframes from Checkfit.m

markstoshow=[7 8 11 12 13 15 19 20]+1;
markstolabel={'pCASL start','pCASL end','Diffusion start','Diffusion end','CO2 start','CO2 seen','pCASL start','pCASL end'};

%If using the Moberg CNS vitals monitor:
usedCNS=1; % Change to 1 if using vitals data
usedCompuRecord=0;
%Set parameters for excluding erroneous vitals data
hrmin=50;
bpmin=15;
bpmax=130;
spo2min=20;
rapmin=1;
%Using invivo or GE vitals monitor?
invivo=1;%0=GE, 1=invivo
savematfile=1;
timesyncmark_dcs=[1]+1;%Note the mark in optics data where time synced

tt=['load ' subjectid '_' patdate ext '_1_flow_output_fitavg.mat'];
eval(tt);
%Remove marks which exceed timeaxis_flow
Marksflow = setdiff(Marksflow,Marksflow(find(Marksflow>length(timeaxis_flow))));

excelfiletoread=['../' subjectid '/' subjectid 'notes/' subjectid '_' patdate ext 'VitalSigns_EBedits.xls'];

%oxygen data good or not?
oxydatagood=1;
flowdatagood=1;

%plot oxygen or not?
plotoxygen=1;
plotflow=1;
subplotsize=2;
secondplot=plotflow+plotoxygen;

if plotoxygen
    tt=['load ' subjectid '_' patdate ext '_dpfout.mat'];
    eval(tt);
     %if first mark is missing
    tmp=Marks;
    Marks(2:length(tmp)+1)=tmp;
    Marks(1)=1;
    if usedISS==0
        hbseries=squeeze(hbseries(useddet,1,:));
        hbo2series=squeeze(hbo2series(useddet,1,:));
    end
else
    Marks=Marksflow;
end

%Make first frame be mark 1 for flow
tmp=Marksflow;
Marksflow(2:length(tmp)+1)=tmp;
Marksflow(1)=1;


%%for saving figures
savefigures=1; %=1 saves

%filter for gaussian
myfiltersize=1;

labelshift=0; %in seconds
labelshift=labelshift./60;
labelshiftyflow=15;
labelshiftyoxy=3;
legendlocation=2;
rotatelabel=1;%If want labels rotated to save space

%TO CUT DATA
cutdata=0; %to cut the data a some mark
if cutdata==0
    if plotoxygen
        Marks(length(Marks)+1)=length(hbo2series);
        cutpointo2=length(Marks);
    end
    if fitavg
        Marksflow(length(Marksflow)+1)=length(Dbfitavg);
    else
        Marksflow(length(Marksflow)+1)=length(Dbfit);
    end
    cutpointflow=length(Marksflow);
end

%plot ranges
oxylimrange=[-10 20];
flowlimrange=[50 300];
hrlim=[80 200];
bplim=[0 130];
SpO2lim=[50 100];


%calculate time axis for flow and oxygenation
%if plotoxygen
    timeaxis=ISStime;
    Markstime=timeaxis(Marks); %time corrected
    
%end
if plotflow
    Markstime_flow=timeaxis_flow(Marksflow); %time corrected
    baselinerangeflow=Marksflow(baselinemarks+1);
    
end

if ~oxydatagood
    hbseries(:)=nan;
    hbo2series(:)=nan;
end

fig1=figure;
%To plot oxygen data
if plotoxygen
    %Set Hbsmooth=Hb
    hbseries_smooth=hbseries(1:Marks(cutpointo2));
    hbo2series_smooth=hbo2series(1:Marks(cutpointo2));
    %Filter data, however, if NaN, filtfilt doesnt work, so break into chunks
    %where there are no NaNs
    nans=find(isnan(hbseries(1:Marks(cutpointo2))));
    if isempty(nans)
        hbseries_smooth=strokefilter(hbseries,myfiltersize);
        hbo2series_smooth=strokefilter(hbo2series,myfiltersize);
    elseif length(nans)==1
        %Smooth the inital part of the data until reaching first NaN
        if length(1:nans-1)>3*length(-2*(ceil(myfiltersize./2)/2):1:(ceil(myfiltersize./2)/2)*2)
            hbseries_smooth(1:nans-1)=strokefilter(hbseries(1:nans-1),myfiltersize);
            hbo2series_smooth(1:nans-1)=strokefilter(hbo2series(1:nans-1),myfiltersize);
        else
            hbseries_smooth(1:nans-1)=hbseries(1:nans-1);
            hbo2series_smooth(1:nans-1)=hbo2series(1:nans-1);
        end
        %Smooth the final part of the data after last NaN
        if length(nans+1:length(hbseries))>3*length(-2*(ceil(myfiltersize./2)/2):1:(ceil(myfiltersize./2)/2)*2)
            hbseries_smooth(nans+1:end)=strokefilter(hbseries(nans+1:end),myfiltersize);
            hbo2series_smooth(nans+1:end)=strokefilter(hbo2series(nans+1:end),myfiltersize);
        else
            hbseries_smooth(nans+1:length(hbseries))=hbseries(nans+1:end);
            hbo2series_smooth(nans+1:length(hbo2series))=hbo2series(nans+1:end);
        end
    else
        for i=1:length(nans)-1
            diff(i)=nans(i+1)-nans(i);
        end
        minchunk=nans(find(diff>5));
        maxchunk=nans(find(diff>5)+1);
        %Smooth the inital part of the data until reaching first NaN
        if length(1:min(nans)-1)>3*length(-2*(ceil(myfiltersize./2)/2):1:(ceil(myfiltersize./2)/2)*2)
            hbseries_smooth(1:min(nans)-1)=strokefilter(hbseries(1:min(nans)-1),myfiltersize);
            hbo2series_smooth(1:min(nans)-1)=strokefilter(hbo2series(1:min(nans)-1),myfiltersize);
        else
            hbseries_smooth(1:min(nans)-1)=hbseries(1:min(nans)-1);
            hbo2series_smooth(1:min(nans)-1)=hbo2series(1:min(nans)-1);
        end
        %Smooth the final part of the data after last NaN
        if length(max(nans)+1:length(hbseries))>3*length(-2*(ceil(myfiltersize./2)/2):1:(ceil(myfiltersize./2)/2)*2)
            hbseries_smooth(max(nans)+1:end)=strokefilter(hbseries(max(nans)+1:end),myfiltersize);
            hbo2series_smooth(max(nans)+1:end)=strokefilter(hbo2series(max(nans)+1:end),myfiltersize);
        else
            hbseries_smooth(max(nans)+1:length(hbseries))=hbseries(max(nans)+1:end);
            hbo2series_smooth(max(nans)+1:length(hbseries))=hbo2series(max(nans)+1:end);
            
        end
        %Smooth chunks, but only if they are at least 3 times greater than the
        %filter size, otherwise filtfilt doesnt work.
        for l=1:length(minchunk)
            chunk=minchunk(l)+1:maxchunk(l)-1;
            if length(chunk)>3*length(-2*(ceil(myfiltersize./2)/2):1:(ceil(myfiltersize./2)/2)*2)
                hbseries_smooth(chunk)=strokefilter(hbseries(chunk),myfiltersize);
                hbo2series_smooth(chunk)=strokefilter(hbo2series(chunk),myfiltersize);
            else
                hbseries_smooth(chunk)=hbseries(chunk);
                hbo2series_smooth(chunk)=hbo2series(chunk);
            end
            clear chunk
        end
        clear nans diff minchunk maxchunk
    end
    
    figure(fig1),subplot('Position',[0.13    0.5    0.775    0.32])
    plot(timeaxis,hbo2series_smooth*1000,'x-r','LineWidth',2)
    hold on, plot(timeaxis,hbseries_smooth.*1000,'o-b','LineWidth',2)
    t=title(['Patient ID=' subjectid '\_' patdate ],'FontSize',24);
    pos=get(t,'Position');
    pos(2)=pos(2)+5;
    set(t,'Position',pos)
    ylim(oxylimrange)
    xlim([min(timeaxis_flow) max(timeaxis_flow)])
    set(gca,'FontSize',24)
    grid on
    set(gca,'XTick',[])
    legend({'HbO2','Hb'},legendlocation,'FontSize',16)
    ylabel('\Delta\muM','FontSize',24)
    %xlabel('Min','FontSize',24)
    tmplim3=get(gca,'YLim');

    for kkkk=1:length(markstoshow)
        line([Markstime_flow(markstoshow(kkkk)) Markstime_flow(markstoshow(kkkk))],[tmplim3(1) tmplim3(2)],'Color',[0 0 0])
        if rotatelabel==0
            ht=text(Markstime_flow(markstoshow(kkkk))+labelshift, tmplim3(2)-mod(labelshiftyoxy-kkkk*3,tmplim3(2)-tmplim3(1)),markstolabel(kkkk),'Color',[0 0 0],'FontSize',16);
        else
            ht=text(Markstime_flow(markstoshow(kkkk))+labelshift, tmplim3(2),markstolabel(kkkk),'Color',[0 0 0],'FontSize',16);
            set(ht,'Rotation',45)
        end
    end
    for kkkk=1:size(baselinemarks,1)*2
        line([Markstime_flow(baselinemarks(kkkk)+1) Markstime_flow(baselinemarks(kkkk)+1)],[tmplim3(1) tmplim3(2)],'Color',[0 0 0],'LineWidth',3)
    end
    for kkkk=1:size(baselinemarks,1)
        hold on,plot([Markstime_flow(baselinemarks(kkkk,1)+1) Markstime_flow(baselinemarks(kkkk,2)+1)],[tmplim3(1) tmplim3(1)],'-k','LineWidth',4)
        hold on,plot([Markstime_flow(baselinemarks(kkkk,1)+1) Markstime_flow(baselinemarks(kkkk,2)+1)],[tmplim3(2) tmplim3(2)],'-k','LineWidth',4)
    end

    

end

if ~flowdatagood
Dbfit(:,usedflowdets)=nan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

if numframestoavg==1
    if ~fitavg
        Dbfitavg=nanmean(Dbfit(:,usedflowdets),2);
    end
    %Filter data, however, if NaN, filtfilt doesnt work, so break into chunks
    %where there are no NaNs
    %Set DbFitsmooth=Dbfit
    Dbfitavg_smooth=Dbfitavg(1:Marksflow(cutpointflow));
    %timeaxis_smooth=timeaxis;
    nans=find(isnan(Dbfitavg(1:Marksflow(cutpointflow))));
    if isempty(nans)
        Dbfitavg_smooth=strokefilter(Dbfitavg,myfiltersize);
    elseif length(nans)==1
        %Smooth the inital part of the data until reaching first NaN
        if length(1:nans-1)>3*length(-2*(ceil(myfiltersize./2)/2):1:(ceil(myfiltersize./2)/2)*2)
            Dbfitavg_smooth(1:nans-1)=strokefilter(Dbfitavg(1:nans-1),myfiltersize);
        else
            Dbfitavg_smooth(1:nans-1)=Dbfitavg(1:nans-1);
        end
        %Smooth the final part of the data after last NaN
        if length(nans+1:length(Dbfitavg_smooth))>3*length(-2*(ceil(myfiltersize./2)/2):1:(ceil(myfiltersize./2)/2)*2)
            Dbfitavg_smooth(nans+1:end)=strokefilter(Dbfitavg(nans+1:end),myfiltersize);
        else
            Dbfitavg_smooth(nans+1:length(Dbfitavg))=Dbfitavg(nans+1:end);
        end
    else
        for i=1:length(nans)-1
            diff(i)=nans(i+1)-nans(i);
        end
        minchunk=nans(find(diff>5));
        maxchunk=nans(find(diff>5)+1);
        %Smooth the inital part of the data until reaching first NaN
        if length(1:min(nans)-1)>3*length(-2*(ceil(myfiltersize./2)/2):1:(ceil(myfiltersize./2)/2)*2)
            Dbfitavg_smooth(1:min(nans)-1)=strokefilter(Dbfitavg(1:min(nans)-1),myfiltersize);
        else
            Dbfitavg_smooth(1:min(nans)-1)=Dbfitavg(1:min(nans)-1);
        end
        %Smooth the final part of the data after last NaN
        if length(max(nans)+1:length(Dbfitavg))>3*length(-2*(ceil(myfiltersize./2)/2):1:(ceil(myfiltersize./2)/2)*2)
            Dbfitavg_smooth(max(nans)+1:end)=strokefilter(Dbfitavg(max(nans)+1:end),myfiltersize);
        else
            Dbfitavg_smooth(max(nans)+1:length(Dbfitavg))=Dbfitavg(max(nans)+1:end);
        end
        %Smooth chunks, but only if they are at least 3 times greater than the
        %filter size, otherwise filtfilt doesnt work.
        for l=1:length(minchunk)
            chunk=minchunk(l)+1:maxchunk(l)-1;
            if length(chunk)>3*length(-2*(ceil(myfiltersize./2)/2):1:(ceil(myfiltersize./2)/2)*2)
                Dbfitavg_smooth(chunk)=strokefilter(Dbfitavg(chunk),myfiltersize);
            else
                Dbfitavg_smooth(chunk)=Dbfitavg(chunk);
            end
            clear chunk
        end    
    end
else 
    %If I averaged multiple frames in my data analysis, every other frame
    %is NaN.  Still want to be able to smooth the data though.
    Dbfitavg_smoothtmp=Dbfitavg(1:numframestoavg:Marksflow(cutpointflow));
    %timeaxis_smooth=timeaxis(1:numframestoavg:Marksflow(cutpointflow));
    nans=find(isnan(Dbfitavg_smoothtmp));
    Dbfitavg_smooth(nans)=NaN;
    if isempty(nans)
        Dbfitavg_smooth=strokefilter(Dbfitavg_smoothtmp,myfiltersize);
    elseif length(nans)==1
        %Smooth the inital part of the data until reaching first NaN
        if length(1:nans-1)>3*length(-2*(ceil(myfiltersize./2)/2):1:(ceil(myfiltersize./2)/2)*2)
            Dbfitavg_smooth(1:nans-1)=strokefilter(Dbfitavg_smoothtmp(1:nans-1),myfiltersize);
        else
            Dbfitavg_smooth(1:nans-1)=Dbfitavg_smoothtmp(1:nans-1);
        end
        %Smooth the final part of the data after last NaN
        if length(nans+1:length(Dbfitavg_smoothtmp))>3*length(-2*(ceil(myfiltersize./2)/2):1:(ceil(myfiltersize./2)/2)*2)
            Dbfitavg_smooth(nans+1:length(Dbfitavg_smoothtmp))=strokefilter(Dbfitavg_smoothtmp(nans+1:end),myfiltersize);
        else
            Dbfitavg_smooth(nans+1:length(Dbfitavg_smoothtmp))=Dbfitavg_smoothtmp(nans+1:end);
        end
    else
        Dbfitavg_smooth(nans)=NaN;
        for i=1:length(nans)-1
            diff(i)=nans(i+1)-nans(i);
        end
        %Find frame that is the start of each chunk of non-NaN data sets
        minchunk=nans(find(diff>5));
        %Find frame that is the end of each chunk of non-NaN data sets
        maxchunk=nans(find(diff>5)+1);
        %Smooth the inital part of the data until reaching first NaN
        if length(1:min(nans)-1)>3*length(-2*(ceil(myfiltersize./2)/2):1:(ceil(myfiltersize./2)/2)*2)
            Dbfitavg_smooth(1:min(nans)-1)=strokefilter(Dbfitavg_smoothtmp(1:min(nans)-1),myfiltersize);
        else
            Dbfitavg_smooth(1:min(nans)-1)=Dbfitavg_smoothtmp(1:min(nans)-1);
        end
        %Smooth the final part of the data after last NaN
        if length(max(nans)+1:length(Dbfitavg_smoothtmp))>3*length(-2*(ceil(myfiltersize./2)/2):1:(ceil(myfiltersize./2)/2)*2)
            Dbfitavg_smooth(max(nans)+1:length(Dbfitavg_smoothtmp))=strokefilter(Dbfitavg_smoothtmp(max(nans)+1:end),myfiltersize);
        else
            Dbfitavg_smooth(max(nans)+1:length(Dbfitavg_smoothtmp))=Dbfitavg_smoothtmp(max(nans)+1:end);
        end
        %Smooth chunks, but only if they are at least 3 times greater than the
        %filter size, otherwise filtfilt doesnt work.
        for l=1:length(minchunk)
            chunk=minchunk(l)+1:maxchunk(l)-1;
            if length(chunk)>3*length(-2*(ceil(myfiltersize./2)/2):1:(ceil(myfiltersize./2)/2)*2)
                Dbfitavg_smooth(chunk)=strokefilter(Dbfitavg_smoothtmp(chunk),myfiltersize);
            else
                Dbfitavg_smooth(chunk)=Dbfitavg_smoothtmp(chunk);
            end
            clear chunk
        end    
    end
    var=ones(size(Dbfitavg))*NaN;
    var(1:numframestoavg:length(Dbfitavg_smooth)*numframestoavg)=Dbfitavg_smooth;
    Dbfitavg_smooth=var;
end

for i=1:size(baselinerangeflow,1)
    baselineflow=baselinerangeflow(i,1):baselinerangeflow(i,2); %baseline frames
    baselines=squeeze(nanmean(Dbfitavg(baselineflow)));
    baselines_smooth=squeeze(nanmean(Dbfitavg_smooth(baselineflow)));
    if i==1
        flowavg=Dbfitavg./baselines;
        flowavg_smooth=Dbfitavg_smooth./baselines_smooth;
    else
        flowavg(baselineflow(1):end)=Dbfitavg(baselineflow(1):end)./baselines;
        flowavg_smooth(baselineflow(1):end)=Dbfitavg_smooth(baselineflow(1):end)./baselines_smooth;
    end
    clear baselineflow
end

if ~flowdatagood
    flowavg_smooth(:)=nan;
end

if plotflow
    figure(fig1)
    subplot(subplotsize,1,secondplot)
    % Shift time axis to match smoothed frame data    
    indexdiff = abs(length(flowavg_smooth) - length(timeaxis_flow));
    timeaxis_flow = timeaxis_flow(1:end-indexdiff);
    h1=plot(timeaxis_flow,flowavg_smooth.*100,'.-','Color',[0 0.5 0],'LineWidth',2,'MarkerSize',25);
    axis tight
    ylim(flowlimrange)
    xlim([min(timeaxis_flow) max(timeaxis_flow)])
    set(gca,'FontSize',24)
    ylabel('rCBF(%)')
    xlabel('Min')
    tmplim3=get(gca,'YLim');
    tmplimx=get(gca,'XLim');
    for kkkk=1:length(markstoshow)
        line([Markstime_flow(markstoshow(kkkk)) Markstime_flow(markstoshow(kkkk))],[tmplim3(1) tmplim3(2)],'Color',[0 0 0])
        if rotatelabel==0
            ht=text(Markstime_flow(markstoshow(kkkk))+labelshift,tmplim3(2)-mod(labelshiftyflow-kkkk*15,tmplim3(2)-tmplim3(1)),markstolabel(kkkk),'Color',[0 0 0],'FontSize',16);
        end
    end
    for kkkk=1:size(baselinemarks,1)*2
        line([Markstime_flow(baselinemarks(kkkk)+1) Markstime_flow(baselinemarks(kkkk)+1)],[tmplim3(1) tmplim3(2)],'Color',[0 0 0],'LineWidth',3)
    end
    for kkkk=1:size(baselinemarks,1)
        hold on,plot([Markstime_flow(baselinemarks(kkkk,1)+1) Markstime_flow(baselinemarks(kkkk,2)+1)],[tmplim3(1) tmplim3(1)],'-k','LineWidth',4)
        hold on,plot([Markstime_flow(baselinemarks(kkkk,1)+1) Markstime_flow(baselinemarks(kkkk,2)+1)],[tmplim3(2) tmplim3(2)],'-k','LineWidth',4)
    
    end
    grid on

end

set(fig1,'PaperPositionMode','Auto')
maxwindows(fig1);

if savefigures
    figure(fig1)
    saveas(gcf,['../' subjectid '/' subjectid 'notes/savedfigs/' patdate '_' subjectid ext '_flowoxy.fig'],'fig')
    saveas(gcf,['../' subjectid '/' subjectid 'notes/savedfigs/' patdate '_' subjectid ext '_flowoxy.eps'],'epsc2')
    saveas(gcf,['../' subjectid '/' subjectid 'notes/savedfigs/' patdate '_' subjectid ext '_flowoxy.png'],'png')

end

if plotoxygen
    tt=['save quantum' subjectid '_' patdate ext ' hbseries_smooth hbo2series_smooth flowavg flowavg_smooth Dbfitavg_smooth baselinerangeflow timeaxis usedCNS timeaxis_flow Marks Marksflow markstoshow markstolabel myfiltersize Markstime_flow Markstime'];
    eval(tt);
else
    tt=['save quantum' subjectid '_' patdate ext ' flowavg flowavg_smooth Dbfitavg_smooth baselinerangeflow timeaxis Marks Marksflow markstoshow  timeaxis_flow markstolabel myfiltersize usedCNS Markstime_flow Markstime'];
    eval(tt);
end

if usedCNS
    CNSvitals(subjectid,patdate,ext,savefigures,timeaxis_flow,Markstime_flow,timesyncmark_dcs,savematfile,rotatelabel,invivo,hrmin,bpmin,bpmax,spo2min,rapmin);
elseif usedCompuRecord
    [vitals vitalheader]=xlsread(['../' subjectid '/' subjectid 'notes/' subjectid '_' patdate ext 'vitals.xls']);
   
    noteminutes= vitals(:,1).*1440-vitals(1,1).*1440;
    Marksvitals=vitals(:,8);
    indsync=find(Marksvitals==timesyncmark_dcs-1);
    
    noteminutes=noteminutes-noteminutes(indsync)+timeaxis(Marksflow(timesyncmark_dcs)); %Mark 1 is alignment point + 1 =2

    
    plotlist=1:size(vitals,1);

    %HR
    figure,
    plot(noteminutes(1:length(plotlist)),vitals(plotlist,2),'.k-','MarkerSize',40,'LineWidth',3)
    set(gca,'FontSize',24)
    xlim([min(timeaxis) max(timeaxis)])
    ylim(hrlim)
    tmplim3=get(gca,'YLim');
    for kkkk=1:length(markstoshow)
        line([Markstime(markstoshow(kkkk)) Markstime(markstoshow(kkkk))],[tmplim3(1) tmplim3(2)],'Color',[0 0 0])
    end
    for kkkk=1:size(baselinemarks,1)*2
        line([Markstime(baselinemarks(kkkk)+1) Markstime(baselinemarks(kkkk)+1)],[tmplim3(1) tmplim3(2)],'Color',[0 0 0],'LineWidth',3)
    end
    labelshifty=-(tmplim3(2)-tmplim3(1)).*0.03;
    for kkkk=1:length(markstolabel)
        if rotatelabel==0
            ht=text(Markstime(markstoshow(kkkk))+labelshift, tmplim3(2)-mod(labelshifty-kkkk*5,tmplim3(2)-tmplim3(1)),markstolabel(kkkk),'Color',[0 0 0],'FontSize',16);
        else
            ht=text(Markstime(markstoshow(kkkk))+labelshift, tmplim3(2),markstolabel(kkkk),'Color',[0 0 0],'FontSize',16);
            set(ht,'Rotation',45)
        end
    end
    grid on
    h1=legend('HR');
    set(h1,'FontSize',16)
    ylabel('Heart Rate (bpm)')
    xlabel('min')
    set(gcf,'PaperPositionMode','Auto')
    maxwindows(gcf);
    if savefigures

        saveas(gcf,['../' subjectid '/' subjectid 'notes/savedfigs/' patdate '_' subjectid ext '_HR.fig'],'fig')
        saveas(gcf,['../' subjectid '/' subjectid 'notes/savedfigs/' patdate '_' subjectid ext '_HR.jpg'],'jpg')

    end
    %SpO2
    figure,
    plot(noteminutes(1:length(plotlist)),vitals(plotlist,3),'.k-','MarkerSize',40,'LineWidth',3)
    set(gca,'FontSize',24)
    xlim([min(timeaxis) max(timeaxis)])
    ylim(SpO2lim)
    tmplim3=get(gca,'YLim');
    for kkkk=1:length(markstoshow)
        line([Markstime(markstoshow(kkkk)) Markstime(markstoshow(kkkk))],[tmplim3(1) tmplim3(2)],'Color',[0 0 0])
    end
    for kkkk=1:size(baselinemarks,1)*2
        line([Markstime(baselinemarks(kkkk)+1) Markstime(baselinemarks(kkkk)+1)],[tmplim3(1) tmplim3(2)],'Color',[0 0 0],'LineWidth',3)
    end
    labelshifty=-(tmplim3(2)-tmplim3(1)).*0.03;
    for kkkk=1:length(markstolabel)
        if rotatelabel==0
            ht=text(Markstime(markstoshow(kkkk))+labelshift, tmplim3(2)-mod(labelshifty-2*(kkkk-1),tmplim3(2)-tmplim3(1)),markstolabel(kkkk),'Color',[0 0 0],'FontSize',16);
        else
            ht=text(Markstime(markstoshow(kkkk))+labelshift, tmplim3(2),markstolabel(kkkk),'Color',[0 0 0],'FontSize',16);
            set(ht,'Rotation',45)
        end
    end

    grid on
    h1=legend('O_2 Sat');
    set(h1,'FontSize',16)
    ylabel('O_2 Sat(%)')
    xlabel('min')
    set(gcf,'PaperPositionMode','Auto')
    maxwindows(gcf);
    if savefigures

        saveas(gcf,['../' subjectid '/' subjectid 'notes/savedfigs/' patdate '_' subjectid ext '_SpO2.fig'],'fig')
        saveas(gcf,['../' subjectid '/' subjectid 'notes/savedfigs/' patdate '_' subjectid ext '_SpO2.jpg'],'jpg')

    end
    %BP
    figure,
    plot(noteminutes(1:length(plotlist)),vitals(plotlist,4),'.b-','MarkerSize',40,'LineWidth',3)
    hold on
    plot(noteminutes(1:length(plotlist)),vitals(plotlist,5),'.-','Color',[0 0.5 0],'MarkerSize',40,'LineWidth',3)
    hold on
    plot(noteminutes(1:length(plotlist)),vitals(plotlist,6),'r.-','MarkerSize',40,'LineWidth',3)
    hold on
    set(gca,'FontSize',24)
    ylim(bplim)
    tmplim3=get(gca,'YLim');
    xlim([min(timeaxis) max(timeaxis)])
    for kkkk=1:length(markstoshow)
        line([Markstime(markstoshow(kkkk)) Markstime(markstoshow(kkkk))],[tmplim3(1) tmplim3(2)],'Color',[0 0 0])
    end
    for kkkk=1:size(baselinemarks,1)*2
        line([Markstime(baselinemarks(kkkk)+1) Markstime(baselinemarks(kkkk)+1)],[tmplim3(1) tmplim3(2)],'Color',[0 0 0],'LineWidth',3)
    end
    labelshifty=-(tmplim3(2)-tmplim3(1)).*0.03;
    for kkkk=1:length(markstolabel)
        if rotatelabel==0
            ht=text(Markstime(markstoshow(kkkk))+labelshift, tmplim3(2)-mod(labelshifty-3*(kkkk-1),tmplim3(2)-tmplim3(1)),markstolabel(kkkk),'Color',[0 0 0],'FontSize',16);
        else
            ht=text(Markstime(markstoshow(kkkk))+labelshift, tmplim3(2),markstolabel(kkkk),'Color',[0 0 0],'FontSize',16);
            set(ht,'Rotation',45)
        end
    end
    grid on
    set(gca,'FontSize',24)
    h1=legend('Systolic','Diastolic','MAP');
    set(h1,'FontSize',16)
    ylabel('mmHg')
    xlabel('min')
    set(gcf,'PaperPositionMode','Auto')
    maxwindows(gcf);
    if savefigures

        saveas(gcf,['../' subjectid '/' subjectid 'notes/savedfigs/' patdate '_' subjectid ext '_BP.fig'],'fig')
        saveas(gcf,['../' subjectid '/' subjectid 'notes/savedfigs/' patdate '_' subjectid ext '_BP.jpg'],'jpg')

    end
    %InCO2
    figure,plot(noteminutes(1:length(plotlist)),vitals(plotlist,7),'-','LineWidth',3,'Color','b','Marker','.','MarkerSize',40);
    set(gca,'FontSize',24)
    tmplim3=get(gca,'YLim');
    for kkkk=1:length(markstoshow)
        line([Markstime(markstoshow(kkkk)) Markstime(markstoshow(kkkk))],[tmplim3(1) tmplim3(2)],'Color',[0 0 0])
    end
    for kkkk=1:size(baselinemarks,1)*2
        line([Markstime(baselinemarks(kkkk)+1) Markstime(baselinemarks(kkkk)+1)],[tmplim3(1) tmplim3(2)],'Color',[0 0 0],'LineWidth',3)
    end
    labelshifty=-(tmplim3(2)-tmplim3(1)).*0.03;
    for kkkk=1:length(markstolabel)
        if rotatelabel==0
            ht=text(Markstime(markstoshow(kkkk))+labelshift, tmplim3(2)-mod(labelshifty-0.1*(kkkk-1),tmplim3(2)-tmplim3(1)),markstolabel(kkkk),'Color',[0 0 0],'FontSize',16);
        else
            ht=text(Markstime(markstoshow(kkkk))+labelshift, tmplim3(2),markstolabel(kkkk),'Color',[0 0 0],'FontSize',16);
            set(ht,'Rotation',45)
        end
    end
    
    grid on
    ylabel('InCO_{2} (mmHg)','FontSize',24)
    xlabel('Min.')
    xlim([min(timeaxis) max(timeaxis)])
    set(gcf,'PaperPositionMode','Auto')
    maxwindows(gcf);
    if savefigures

        saveas(gcf,['../' subjectid '/' subjectid 'notes/savedfigs/' patdate '_' subjectid ext '_inco2.fig'],'fig')
        saveas(gcf,['../' subjectid '/' subjectid 'notes/savedfigs/' patdate '_' subjectid ext '_inco2.jpg'],'jpg')

    end
end

