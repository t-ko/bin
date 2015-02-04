%Updated 7/6/09 to subtract offsets from I and Q data!  We are very close
%to noise floor, so it is important to do this.
%Updated 7/6/09 to include calculation of mua changes and to work for both
%ISS and Joel's Instruments

%CORRECTED ON 10/12/12 SO NOW USING N=1.406 FOR NA OF PHANTOM!!!!

%close all
clear all

patID='SVC169';
patdate='060214';
%If you want something added to your file name
%Extension you want added to filename, if nothing, just do ext='';
ext='';

% Enter indices of respective measurements HERE
% NOTE: SVC116_060614 Bad calibration absorption/scattering values
% post-measurement
groups={[1 3 5];[2 4 6];[7 8 9 10];[11 12 13 14]};
grouplabels={'Calibration block';'Check block 1';'Left forehead';'Right forehead'};

%% Automatically generate legacy {labels} variable (updated 08-08-2013)
n = 0;
for i = 1:length(groups)
	n = n+length(groups{i});
end

for i = 1:n
	for j = 1:length(groups)
		if ~isempty(find(groups{j}==i))
			labels{i} = grouplabels{j};
		end
	end
end
%% Measurement number of the reference we will use for calibration
ref=1;
framesperfile=100;

plotfigs=0;
savefigs=1;

%PC
fdir=['..\' patID '\' patID '_' patdate '\absolute\']

%LINUX
%fdir=['../' patID '/' patID '_' patdate '/absolute/'];

usedISS=1;%0=Joel's Instrument, 1=ISS
usednewsettings=0;%Dennis made new settings file for us to turn off PMTs while DCS laser is on.
useddet=2;%If using ISS instrument, DetA=1, DetB=2
appliedcalibration=1;

%Adult Absolute Probe
r=[2.0 2.51 3.01 3.45; 1.94 2.45 2.96 3.41];
%Custom Probe
%r=[1.44 1.94 2.41 2.90; 1.51 1.99 2.49 2.98];
%Neonate Probe
%r=[1.51 2.01 2.5 3.01; 1.45 1.95 2.44 2.95];
%Custom Probe 2
%r=[1.38 1.87 2.38 2.86; 1.46 1.93 2.46 2.93];


sourceindex=[1 2 3 4; 5 6 7 8];



n0=1.4;%index of refraction for water
c=2.99792458e10; %speed of light in vacuum,cm/s
v=c/n0; %speed of light in medium
w=2*pi*110e6;
R=-1.440./n0^2+0.710/n0+0.668+0.0636.*n0;
ze=2./3.*(1+R)./(1-R); %check this number, no need for 1/musp as Cecil takes ca\re of it in the files

lambdas=[690 830];
numlambda=length(lambdas);
%ISS Analysis
sourcesused=[4 8 13 14 7 11 12 16]; % which sources, out of 16

%LOAD DATA
% Must sort the rest
files=dir([ fdir 'Data_' patdate '*.txt'])
for numfile=1:size(files,1)
    names(numfile,:)=files(numfile).name(12:end-4);
end
names=sortrows(names);
% print out filenames and associated analysis labelto ensure accuracy
filecheck = {};
for i=1:length(labels)
    filecheck{i} = [names(i,:) ' - ' labels{i}];
end
fprintf('%20s\n',filecheck{:})

data=[];
for i=1:size(names,1)
    fname=[ fdir 'Data_' patdate names(i,:) '.txt'];
    fid = fopen([ fname ], 'r');
    if usednewsettings
        [tmpdata, count]=fscanf(fid,'%c %s',1600);
        clear tmpdata count
        [tmpdata, count]=fscanf(fid,'%g %g',[151 inf]);
    else
        if appliedcalibration
            [tmpdata, count]=fscanf(fid,'%c %s',657);
            clear tmpdata count
            [tmpdata, count]=fscanf(fid,'%c %s',96);
            phacalibcoeffs=str2num(tmpdata);
            clear tmpdata count
            [tmpdata, count]=fscanf(fid,'%c %s',2);
            clear tmpdata count
            [tmpdata, count]=fscanf(fid,'%c %s',64);
            ACDCcalibcoeffs=str2num(tmpdata);
            clear tmpdata count
            [tmpdata, count]=fscanf(fid,'%c %s',390);
            clear tmpdata count
            [tmpdata, count]=fscanf(fid,'%g %g',[103 inf]);
        else
            [tmpdata, count]=fscanf(fid,'%c %s',1209);
            clear tmpdata count
            [tmpdata, count]=fscanf(fid,'%g %g',[103 inf]);
        end
    end
    fclose(fid);
    tmpdata=tmpdata.';
    data=cat(1,data,tmpdata);
    l(i)=size(tmpdata,1);
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
    updatetime=1/10;
else
    updatetime=1/6.25; % Update time, set to 6.25 Hz.
end

% Plot noncalibrated AC, DC and phase data for Detectors A & B
for j=1:length(sourcesused)
    if useddet==1
        ind=5;
    elseif useddet==2
        ind=53;
    end
    if appliedcalibration
        ACtmp(:,j)=data(:,ind+sourcesused(j))/ACDCcalibcoeffs(sourcesused(j));
        DCtmp(:,j)=data(:,ind+16+sourcesused(j))/ACDCcalibcoeffs(16+sourcesused(j));
        phasetmp(:,j)=(data(:,ind+16*2+sourcesused(j))+phacalibcoeffs(32+sourcesused(j)))*pi/180;
    else
        ACtmp(:,j)=data(:,ind+sourcesused(j));
        DCtmp(:,j)=data(:,ind+16+sourcesused(j));
        phasetmp(:,j)=data(:,ind+16*2+sourcesused(j))*pi/180;
    end    
end

phaseref=nanmean(phasetmp(:,5));
for j=1:length(sourcesused)
    phase_new(:,j)=phasetmp(:,j)-phaseref;
end

load colors.mat

%Find the mean and standard deviation of each frame of data
for i=1:length(l)
    if i==1
        AC(i,:)=nanmean(ACtmp(1:sum(l(1:i)),:),1);
        DC(i,:)=nanmean(DCtmp(1:sum(l(1:i)),:),1);
        phase(i,:)=nanmean(phase_new(1:sum(l(1:i)),:),1);
        for j=1:numlambda
            phase(i,sourceindex(j,:))=phase(i,sourceindex(j,:))-phase(i,sourceindex(j,1));%Shift phase axis so that phi(r=0.8)=0
            tmp=find(phase(i,sourceindex(j,:))<=-pi/2);
            if ~isempty(tmp)
                for k=1:length(tmp)
                    phase(i,sourceindex(j,tmp(k)))=phase(i,sourceindex(j,tmp(k)))+2*pi;%Unwrap phase if need be
                end
            end
            clear tmp
        end
    else
        AC(i,:)=nanmean(ACtmp(sum(l(1:i-1))+1:sum(l(1:i)),:),1);
        DC(i,:)=nanmean(DCtmp(sum(l(1:i-1))+1:sum(l(1:i)),:),1);
        phase(i,:)=nanmean(phase_new(sum(l(1:i-1))+1:sum(l(1:i)),:),1);
        ACstd(i,:)=nanstd(ACtmp(sum(l(1:i-1))+1:sum(l(1:i)),:),1);
        DCstd(i,:)=nanstd(DCtmp(sum(l(1:i-1))+1:sum(l(1:i)),:),1);
        phasestd(i,:)=nanstd(phase_new(sum(l(1:i-1))+1:sum(l(1:i)),:),1);
        for j=1:numlambda
            phase(i,sourceindex(j,:))=phase(i,sourceindex(j,:))-phase(i,sourceindex(j,1));%Shift phase axis so that phi(r=0.8)=0
            tmp=find(phase(i,sourceindex(j,:))<=-pi/2);
            if ~isempty(tmp)
                for k=1:length(tmp)
                    phase(i,sourceindex(j,tmp(k)))=phase(i,sourceindex(j,tmp(k)))+2*pi;%Unwrap phase if need be
                end
            end
            clear tmp
        end
    end
end

figure,subplot(1,3,1)
for i=1:length(sourcesused)
    hold on,plot(phasetmp(:,i)*180/pi,'.-','LineWidth',3,'Color',colors(i,:))
end
xlabel('Frame')
ylabel('Raw Phase (deg.)')
axis tight
set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(gcf,'PaperPositionMode','Auto')
subplot(1,3,2)
for i=1:length(sourcesused)
    hold on,plot(ACtmp(:,i),'.-','LineWidth',3,'Color',colors(i,:))
end
axis tight
%ylim([-20 20])
%legend('Source 4','Source 8','Source 13','Source 14','Source 7','Source 11','Source 12','Source 16')
xlabel('Frame')
ylabel('AC Amplitude')
subplot(1,3,3)
for i=1:length(sourcesused)
    hold on,plot(DCtmp(:,i),'.-','LineWidth',3,'Color',colors(i,:))
end
axis tight
%ylim([-20 20])
legend('1.5cm-690','2.0cm-690','2.5cm-690','3.0cm-690','1.5cm-830','2.0cm-830','2.5cm-830','3.0cm-830')
xlabel('Frame')
ylabel('DC Amplitude')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(gcf,'PaperPositionMode','Auto')
maxwindows(gcf)
saveas(gcf,['../' patID '/' patID 'notes/savedfigs/ACDCamp_phase.jpg'])
saveas(gcf,['../' patID '/' patID 'notes/savedfigs/ACDCamp_phase.fig'])
    

figure,
for i=1:length(sourcesused)
    hold on,errorbar(AC(:,i),ACstd(:,i),'LineWidth',3,'Color',colors(i,:))
end
%xlabel('Frame','FontSize',20)
ylabel('Amplitude','FontSize',20)
axis tight
grid on
legend([ num2str(sourcesused.') ],'Location','EastOutside')
set(gcf,'PaperPositionMode','Auto')
set(gca,'FontSize',20)
if savefigs
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Amplitudemeans_' patID '_' patdate ext '.fig'],'fig')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Amplitudemeans_' patID '_' patdate ext '.jpg'],'jpg')
end
figure,
for i=1:length(sourcesused)
    hold on,errorbar(DC(:,i),DCstd(:,i),'LineWidth',3,'Color',colors(i,:))
end
%xlabel('Frame','FontSize',20)
ylabel('DC','FontSize',20)
axis tight
grid on
legend([ num2str(sourcesused.') ],'Location','EastOutside')
set(gcf,'PaperPositionMode','Auto')
set(gca,'FontSize',20)
if savefigs
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DCmeans_' patID '_' patdate ext '.fig'],'fig')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DCmeans_' patID '_' patdate ext '.jpg'],'jpg')
end
figure,
for i=1:length(sourcesused)
    hold on,errorbar(phase(:,i)*180/pi,phasestd(:,i)*180/pi,'LineWidth',3,'Color',colors(i,:))
end
xlabel('Frame','FontSize',25)
ylabel('Phase (deg)','FontSize',25)
set(gca,'FontSize',25)
axis tight
legend([ num2str(sourcesused.') ],'Location','EastOutside')
grid on
set(gcf,'PaperPositionMode','Auto')
if savefigs
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Phasemeans_' patID '_' patdate ext '.fig'],'fig')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Phasemeans_' patID '_' patdate ext '.jpg'],'jpg')
end


%Calibrate measurements based on "known" optical properties
mua=[0.137 0.133];
musp=[5.19 4.32];
D=v/3./musp;
amp=AC;
pha=phase;

for i=1:numlambda
    %Calculate slopes of amp and phase based on these assumptions
    kr(i)=sqrt(v.*mua(i)./(2.*D(i)).*(sqrt(1+(w./v./mua(i)).^2)-1));%Slope of phase vs. r
    ki(i)=sqrt(v.*mua(i)./(2.*D(i)).*(sqrt(1+(w./v./mua(i)).^2)+1));%Slope of log(Amp*r2) vs. r
    %Calculate amp and phase corrections for each separation
    for rr=1:length(groups{1})
        ref=groups{1}(rr);
        %From intercepts, find expected amplitude and phase for each r
        expLNampr2tmp(rr,i,:)=-ki(i).*(r(i,:)-r(i,1))+log(AC(ref,sourceindex(i,1)).*r(i,1)^2);%Expected LN(Ar^2)
        expamp(rr,i,:)=squeeze(exp(expLNampr2tmp(rr,i,:))).'./r(i,:).^2;%Expected A
        expphatmp(rr,i,:)=kr(i).*(r(i,:)-r(i,1))+phase(ref,sourceindex(i,1));%Expected phi
        expphatmp(rr,i,:)=expphatmp(rr,i,:)-expphatmp(rr,i,1);%Shift phase axis so that phi(r=1.5)=0
        
        ampcorrectiontmp(rr,i,:)=squeeze(expamp(rr,i,:)).'./AC(ref,sourceindex(i,:));
        phacorrectiontmp(rr,i,:)=squeeze(expphatmp(rr,i,:)).'-phase(ref,sourceindex(i,:));
    end
    expLNampr2(i,:)=nanmean(squeeze(expLNampr2tmp(:,i,:)),1);
    exppha(i,:)=nanmean(squeeze(expphatmp(:,i,:)),1);

    ampcorrection(i,:)=nanmean(squeeze(ampcorrectiontmp(:,i,:)),1);
    phacorrection(i,:)=nanmean(squeeze(phacorrectiontmp(:,i,:)),1);
end

%Plot corrected data from calibration source
figure,subplot(2,1,1)
plot(r(1,:),log(amp(ref,sourceindex(1,:)).*r(1,:).^2),'b.','MarkerSize',30,'LineWidth',3)
hold on,plot(r(2,:),log(amp(ref,sourceindex(2,:)).*r(2,:).^2),'.','Color',[0 0.5 0],'MarkerSize',30,'LineWidth',3)
hold on,plot(r(1,:),log(amp(ref,sourceindex(1,:)).*ampcorrection(1,:).*r(1,:).^2),'bo','MarkerSize',10,'LineWidth',3)
hold on,plot(r(2,:),log(amp(ref,sourceindex(2,:)).*ampcorrection(2,:).*r(2,:).^2),'o','MarkerSize',10,'Color',[0 0.5 0],'LineWidth',3)
hold on,plot(r(1,:),expLNampr2(1,:),'b--','LineWidth',3)
hold on,plot(r(2,:),expLNampr2(2,:),'--','Color',[0 0.5 0],'LineWidth',3)
title([ labels(ref) ])
ylabel('LN(A*r^{2})')
legend('688nm','830nm')
set(gca,'XTickLabel',[])
axis tight
xlim([1.3 3.2])
xlim1=get(gca,'XLim');
ylim1=get(gca,'YLim');
text(xlim1(1)+0.25,(ylim1(2)-ylim1(1))*-1/8,['Expect \mu_{a}=' num2str(mua(1),'%6.2f') 'cm^{-1}, \mu_{sp}=' num2str(musp(1),'%6.2f') 'cm^{-1}' ],'Color','b')
text(xlim1(1)+0.25,(ylim1(2)-ylim1(1))*-1/4,['Expect \mu_{a}=' num2str(mua(2),'%6.2f') 'cm^{-1}, \mu_{sp}=' num2str(musp(2),'%6.2f') 'cm^{-1}' ],'Color',[0 0.5 0])
grid on
subplot(2,1,2)
plot(r(1,:),pha(ref,sourceindex(1,:)),'b.','MarkerSize',30,'LineWidth',3)
hold on,plot(r(2,:),pha(ref,sourceindex(2,:)),'.','Color',[0 0.5 0],'MarkerSize',30,'LineWidth',3)
hold on,plot(r(1,:),pha(ref,sourceindex(1,:))+phacorrection(1,:),'bo','MarkerSize',10,'LineWidth',3)
hold on,plot(r(2,:),pha(ref,sourceindex(2,:))+phacorrection(2,:),'o','MarkerSize',10,'Color',[0 0.5 0],'LineWidth',3)
hold on,plot(r(1,:),exppha(1,:),'b--','LineWidth',3)
hold on,plot(r(2,:),exppha(2,:),'--','Color',[0 0.5 0],'LineWidth',3)
xlabel('Separation (cm)')
ylabel('Phase (rad)')
axis tight
xlim([1.3 3.2])
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',40)
maxwindows(gcf)
set(gcf,'PaperPositionMode','Auto')

for m=1:size(AC,1)
  
    n=1.4;
    v=c/n;
    for i=1:numlambda
        [R,p]=corrcoef(r(i,:),log(r(i,:).^2.*amp(m,sourceindex(i,:)).*ampcorrection(i,:)));
        R2amp(m,i)=R(2)^2;
        pamp(m,i)=p(2);
        p_amp(i,:)=polyfit(r(i,:),log(r(i,:).^2.*amp(m,sourceindex(i,:)).*ampcorrection(i,:)),1);
        fitamp(i,:)=polyval(p_amp(i,:),r(i,:));
        [b_amp(m,i,:),b_amp_int(m,i,:,:)] = regress((log(r(i,:).^2.*amp(m,sourceindex(i,:)).*ampcorrection(i,:))).',[r(i,:).' ones(size(r(i,:).'))]);
        %{
            if i==2
                p_pha(i,:)=polyfit(r(i,1:3),pha(m,sourceindex(i,1:3))+phacorrection(i,1:3),1);
                [R,p]=corrcoef(r(i,1:3),pha(m,sourceindex(i,1:3))+phacorrection(i,1:3));
                R2pha(m,i)=R(2)^2;
                ppha(m,i)=p(2);
            else
        %}
        p_pha(i,:)=polyfit(r(i,:),pha(m,sourceindex(i,:))+phacorrection(i,:),1);
        [R,p]=corrcoef(r(i,:),pha(m,sourceindex(i,:))+phacorrection(i,:));
        R2pha(m,i)=R(2)^2;
        ppha(m,i)=p(2);
        %            end
        fitpha(i,:)=polyval(p_pha(i,:),r(i,:));

        fitmua(m,i)=abs(w/(2*v)*(p_amp(i,1)/p_pha(i,1)-p_pha(i,1)/p_amp(i,1)));
        fitmusp(m,i)=abs(2*v/3/w*p_amp(i,1)*p_pha(i,1));
        if (R2amp(m,i) <0.97) || (R2pha(m,i)<0.97)
            fitmua(m,i)=NaN;
            fitmusp(m,i)=NaN;
        end
    end
    
    if plotfigs
        figure,subplot(2,1,1)
        plot(r(1,:),log(amp(m,sourceindex(1,:)).*ampcorrection(1,:).*r(1,:).^2),'b.','MarkerSize',30,'LineWidth',3)
        hold on,plot(r(2,:),log(amp(m,sourceindex(2,:)).*ampcorrection(2,:).*r(2,:).^2),'.','Color',[0 0.5 0],'MarkerSize',30,'LineWidth',3)
        hold on,plot(r(1,:),fitamp(1,:),'b--','LineWidth',3)
        hold on,plot(r(2,:),fitamp(2,:),'--','Color',[0 0.5 0],'LineWidth',3)
        title([ labels(m) ])
        ylabel('LN(A*r^{2})')
        legend('688nm','830nm')
        axis tight
        set(gca,'XTickLabel',[])
        xlim1=get(gca,'XLim');
        ylim1=get(gca,'YLim');
        text(xlim1(1)+0.25,(ylim1(2)-ylim1(1))*-1/8,['Fit \mu_{a}=' num2str(fitmua(m,1),'%6.2f') 'cm^{-1}, \mu_{sp}=' num2str(fitmusp(m,1),'%6.2f') 'cm^{-1}' ],'Color','b')
        text(xlim1(1)+0.25,(ylim1(2)-ylim1(1))*-1/4,['Fit \mu_{a}=' num2str(fitmua(m,2),'%6.2f') 'cm^{-1}, \mu_{sp}=' num2str(fitmusp(m,2),'%6.2f') 'cm^{-1}' ],'Color',[0 0.5 0])
        grid on
        subplot(2,1,2)
        plot(r(1,:),pha(m,sourceindex(1,:))+phacorrection(1,:),'b.','MarkerSize',30,'LineWidth',3)
        hold on,plot(r(2,:),pha(m,sourceindex(2,:))+phacorrection(2,:),'.','Color',[0 0.5 0],'MarkerSize',30,'LineWidth',3)
        hold on,plot(r(1,:),fitpha(1,:),'b--','LineWidth',3)
        hold on,plot(r(2,:),fitpha(2,:),'--','Color',[0 0.5 0],'LineWidth',3)
        xlabel('Separation (cm)')
        ylabel('Phase (rad)')
        axis tight
        grid on
        set(findall(gcf,'-property','FontSize'),'FontSize',40)
        maxwindows(gcf)
        set(gcf,'PaperPositionMode','Auto')
        if savefigs
            saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Absolute_' patID '_' patdate ext '_' char(labels(m)) '.fig'],'fig')
            saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Absolute_' patID '_' patdate ext '_' char(labels(m)) '.jpg'],'jpg')
        end
    end
end

%Plot mua
figure,
for i=1:size(groups,1)
    subplot(2,2,i)
    plot(fitmua(groups{i},1),'b.','MarkerSize',30)
    hold on,plot(fitmua(groups{i},2),'.','Color',[0 0.5 0],'MarkerSize',30)
    title([ grouplabels(i) ])
    ylabel('\mu_{a} (cm-1)')
    ylim([0.005 0.3])
    xlim([0.75 size(groups{i},2)+0.25])
    grid on
    %legend('688nm','830nm')
    if i>2
        xlabel('Trial No.')
    end
end
set(findall(gcf,'-property','FontSize'),'FontSize',20)
maxwindows(gcf)
set(gcf,'PaperPositionMode','Auto')
if savefigs
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Absolute_' patID '_' patdate ext '_muaALL_usingRef_' char(labels(ref)) '.fig'],'fig')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Absolute_' patID '_' patdate ext '_muaALL_usingRef_' char(labels(ref)) '.jpg'],'jpg')
end
%Plot musp
figure,
for i=1:size(groups,1)
    subplot(2,2,i)
    plot(fitmusp(groups{i},1),'b.','MarkerSize',30)
    hold on,plot(fitmusp(groups{i},2),'.','Color',[0 0.5 0],'MarkerSize',30)
    title([ grouplabels(i) ])
    ylim([0 15])
    xlim([0.75 size(groups{i},2)+0.25])
    grid on
    ylabel('\mu_{sp} (cm-1)')
    %legend('688nm','830nm')
    if i>2
        xlabel('Trial No.')
    end
end
set(findall(gcf,'-property','FontSize'),'FontSize',20)
maxwindows(gcf)
set(gcf,'PaperPositionMode','Auto')
if savefigs
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Absolute_' patID '_' patdate ext '_muspALL_usingRef_' char(labels(ref)) '.fig'],'fig')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Absolute_' patID '_' patdate ext '_muspALL_usingRef_' char(labels(ref)) '.jpg'],'jpg')
end
%}
tt=['save ' patID '_' patdate ext '_absolute_dpfout.mat AC ACstd DC DCstd phase phasestd data useddet fitmua fitmusp usedISS groups grouplabels labels'];
eval(tt);
