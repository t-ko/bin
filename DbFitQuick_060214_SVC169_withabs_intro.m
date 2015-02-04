%Dbfit code template- EB 6/26/09
%see labnote book on why I chose some of these parameters
%Right now only works for 1 s-d separation

clear all
close all

patID='SVC169';
patdate='060214';
%Extension you want added to filename, if nothing, just use ext='';
exten={'_1'};


usedflowdets=[2 8];%must be consecutive, cannot be ex. 1,3,4 or 2,4.  must change code if this is the case 
marksforgamma=[5 7];%Set a range of marks in which the curve decay is ~ the same order of magnitude.  For hypercapnia study, would want to choose only baseline
collectedo2=1;%Was oxygenation data collected in this study?  If so, will use changes in mua in our Dbfit
syncmark=1;
%Only for use with ISS, bc the ISS and DCS data are not lined up and will affect the choice of mua for fitting Db

n0=1.4;%index of refraction for tissue
c=2.99792458e10; %speed of light in vacuum,cm/s
vo=c/n0; %speed of light in medium
lambda=7.85e-5;
R=-1.440./n0^2+0.710/n0+0.668+0.0636.*n0;
ze=2./3.*(1+R)./(1-R); %check this number, no need for 1/musp as Cecil takes ca\re of it in the files
muao=0.1;
muspo=7;
Dbr0=1e-8; %background brownian diffusion coeff
D=vo./(3.0.*muspo); %background diffusion coeff
k0=2*pi*n0/lambda; %this is the k0 for flow!
thickness = 1; %thickness of extra-cerebral tissue (cm)

%%%% FITTING SETTINGS:
fitbeta=1; %=0 do not fit beta,==1 fit beta, DO NOT RECOMMEND FITTING BETA UNLESS INTENSITY > 20kHz
betafitmin=0.4;
betafitmax=0.6;
op=optimset('fminsearch');
options=optimset(op,'MaxIter',500,'MaxFunEvals',500,'TolFun',1.000e-7,'TolX',1.000e-9,'Display','Final');
startcorr=1;% if you are fitting beta, want to use a higher number here (e.g. 9), if not fitting beta, use smallest number possible
datalength=90; %chao's way of cutting data for decay.
avgnum=10;%How many points to average in each curve before fitting
cutoff=0;%where to cut correlation curve off when fitting
setbeta=0.45; %used when calculating sigma
x0=[1e-8 setbeta];%Initial guess for Db and beta
cutoffmeanerror_perc=2; %mean error suggested cutoff, upper threshold; whole integer percentage removed

plott=0; %=1 plots the probe map

%%%%PLOT SETTINGS
intenymax=70;%ylim for intensity plot

%Get probe information
flatbabyprobe;
Ns=size(fsources,1);
Nd=size(fdetectors,1);
r=fsdlist(1,1);

%Load flow data
minfiles=0;
measnum=0;
morefiles=1; %1:more files to read, 0:quit
snum=1;

%Step 1, plot amplitudes for each file extension and assess probe off marks
%and baseline marks
measnumtmp=0;
measnumtmp1=0;
for e=1:length(exten)
    fdir1=['../' patID '/' patID '_' patdate '/'];
    fname1=[ patID '_' patdate  char(exten(e)) '_'];

    %Load flow data
    minfiles=0;
    measnum=0;
    morefiles=1; %1:more files to read, 0:quit
    snum=1;

    %Load all data
    while (morefiles)

        if snum>Ns
            snum=1;
            measnum=measnum+1;
        end

        numfiles=measnum.*Ns+snum+minfiles-1;
        fname=[fdir1 fname1 'flow_' sprintf('%01d',numfiles) '.dat'];

        if exist(fname)==2 %file exists
            %Find time of frame:
            fid = fopen([ fname ], 'r');
            [tmpdata, count]=fscanf(fid,'%c %s',3);
            clear tmpdata
            [tmpdata, count]=fscanf(fid,'%c %s',2);
            time{measnum+1}=tmpdata;
            fclose(fid);
            
            %Load correlator data
            data=load(fname);
            
            %Record marks
            if data(end,1)>0
                Marksflow(data(end,1))=measnum+1;
            end

            if size(data,2)==9
                intensitydata(measnum+1,:)=data(1,2:9);
                corrs(measnum+1,:,:)=data(3:end,2:9);
            else
                intensitydata(measnum+1,:)=data(1,2:5);
                corrs(measnum+1,:,:)=data(3:end,2:5);
            end
            %Find frames with light leakage
            for d=1:size(corrs,3)
                if mean(corrs(measnum+1,datalength:datalength+20,d))>1.01
                    corrs(measnum+1,:,d)=NaN;
                end
            end
            taus=data(3:end,1);
            snum=snum+1;
        elseif exist(fname)==0 %file does not exist
            morefiles=0;
        end
    end
    marksperext(e)=length(Marksflow);
    filelength(e)=size(intensitydata,1)+measnumtmp;
    %Make running amplitude and Marks variables
    intensitydatatmp(measnumtmp+1:filelength(e),:)=intensitydata;
    corrstmp(measnumtmp+1:filelength(e),:,:)=corrs;
    timetmp(measnumtmp+1:filelength(e))=time;
    Marksflowtmp(measnumtmp1+1:length(Marksflow)+measnumtmp1)=Marksflow+measnumtmp+1;
    measnumtmp=filelength(e);
    measnumtmp1=length(Marksflow)+measnumtmp1;
    %Get integration time (sec)
    t=data(1,1)/1000;

    clear data Marksflow corrs intensitydata time
end
Marksflow=Marksflowtmp;
intensitydata=intensitydatatmp;
corrs=corrstmp;

%Calculate time points corresponding to each flow curve
timeaxis_flow=datenum(timetmp,'HH:MM:SS')-floor(datenum(timetmp,'HH:MM:SS'));%In arbitrary units--1a.u.=24hrs, counting from 1/1/2000.

%Data in DCS files is not in military time, so add 12 hours to all data if
%need be
if timeaxis_flow(1)<0.5
    timeaxis_flow=timeaxis_flow+0.5;
end

%Unwrap time vector
for i=1:length(timeaxis_flow)-1
    if abs(timeaxis_flow(i+1)-timeaxis_flow(i))>0.4
        timeaxis_flow(i+1)=timeaxis_flow(i+1)+0.5;
    end
end

%Zero time axis so that t=0 is first flow curve
timeaxis_flow=(timeaxis_flow-timeaxis_flow(1))*1.44e3;%Convert time to minutes

%Correct for time since off bypass
timeaxis_flow=timeaxis_flow-timeaxis_flow(Marksflow(syncmark));

figure,plot(timeaxis_flow)

%Remove marks which exceed timeaxis_flow
Marksflow = setdiff(Marksflow,Marksflow(find(Marksflow>length(timeaxis_flow))));

%Import absolute mua and musp, if measured (TK 2014-8-2)
if collectedo2
    if exist([ patID '_' patdate '_dpfout_withabs.mat'],'file')
        load([ patID '_' patdate '_dpfout_withabs.mat']);
         
        %Sync ISS and flow 

        if threelambda
            mua=interp1(ISStime,muaseries(3,:),timeaxis_flow);
        else
            mua=interp1(ISStime,muaseries(2,:),timeaxis_flow);
        end
        %muspo=mean(muspmean);%Mean across three wavelenghtsave musp for 830 and 690nm, just take mean to get musp at 790 for now
        muspo=muspmean(3);%select scattering at wavelength closest to 785 (TK 2014-08)
    else
        %Initialize mua vector
        mua=ones(size(corrs,1),1)*muao;
        %Sync ISS and flow 
        diff=Marks(syncmark)-Marksflow(syncmark);
        %Marksflow(syncmark)=# of flow frames recorded since start of data
        if threelambda
            if diff<0
                %Frames 1:abs(diff)==muao, rest of frames will be
                %muao+deltamua
                mua(abs(diff)+1:abs(diff)+length(muaseries))=muao+squeeze(muaseries(3,:));
            elseif diff>=0
                mua=muao+squeeze(muaseries(3,diff+1:end));
            end
        else
            if diff<0
                %Frames 1:abs(diff)==muao, rest of frames will be
                %muao+deltamua
                mua(abs(diff)+1:abs(diff)+length(muaseries))=muao+squeeze(muaseries(2,:));
            elseif diff>=0
                mua=muao+squeeze(muaseries(2,diff+1:end));
            end
                        
        end
    end
end

%Determine bin width for each tau
T=zeros(size(taus));
for indt=1:length(T)-1
    T(indt)=taus(indt+1)-taus(indt);
end
%If we didnt collect oxygen data, need to keep muao constant:
if ~collectedo2
    mua=muao+zeros(1,size(intensitydata,1));
end
%Plot intensities to determine whether to use the avg of dets or each
%individually
if length(usedflowdets)==1
    figure,plot(intensitydata(:,usedflowdets),'.-','MarkerSize',20,'LineWidth',3)
else
    figure,plot(intensitydata(:,usedflowdets(1):usedflowdets(2)),'.-','MarkerSize',20,'LineWidth',3)
end
xlabel('Frame','FontSize',25)
ylabel('Intensity (kHz)','FontSize',25)
axis tight
ylim([0 intenymax])
grid on
set(gca,'FontSize',20)
set(gcf,'PaperPositionMode','Auto')
saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Intensity_' patID '_' patdate '.fig'],'fig')
saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Intensity_' patID '_' patdate '.jpg'],'jpg')

if length(usedflowdets)==1
    fitavg=0;
else
    %Ask for user input as to whether to fit average or individual curves
    fitavg=input('Fit average of corr. curves? 0=no, 1=yes  ');
end

fitmultframes=input('Fit average of multiple frames?  0=no, 1=yes (this option is for very noisy data)  ');
if fitmultframes==1
    numframestoavg=input('How many frames would you like to average?   ');
else
    numframestoavg=1;
end

%Define a cutoff intensity to get rid of bad frames
cutoffintensity=input('Cutoff intensity value (kHz)? (All frames with intensity above this value will be discarded)  ');

%Plot all correlation curves to get a sense of what beta should be
if length(usedflowdets)==1
    d=usedflowdets;
    %Find frames where I is bigger than cutoff
    int=find(intensitydata(:,d)>cutoffintensity);
    %Set frames with I above cutoff = NaN
    corrs(int,:,d)=NaN;

    figure,semilogx(taus,corrs(1,:,d),'r-','LineWidth',1)
    for m=2:size(corrs,1)
        hold on,semilogx(taus,corrs(m,:,d),'r-','LineWidth',1)
    end
    hold on,semilogx(taus,nanmean(corrs(:,:,d),1),'k-','LineWidth',5)
    xlabel('\tau','FontSize',25)
    ylabel('g2','FontSize',25)
    title(['Detector ' num2str(d) ],'FontSize',25)
    ylim([0.9 1.6])
    xlim([0 1e-2])
    set(gca,'FontSize',20)
    h=line([taus(datalength) taus(datalength)],[0.9 1.6]);
    set(h,'Color',[0 0 0]);
    set(gcf,'PaperPositionMode','Auto')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/corrs_det' num2str(d) '_' patID '_' patdate '.fig'],'fig')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/corrs_det' num2str(d) '_' patID '_' patdate '.jpg'],'jpg')
    clear int
else
    for d=usedflowdets(1):1:usedflowdets(2)
        %Find frames where I is bigger than cutoff
        int=find(intensitydata(:,d)>cutoffintensity);
        %Set frames with I above cutoff = NaN
        corrs(int,:,d)=NaN;

        figure,semilogx(taus,corrs(1,:,d),'r-','LineWidth',1)
        for m=2:size(corrs,1)
            hold on,semilogx(taus,corrs(m,:,d),'r-','LineWidth',1)
        end
        hold on,semilogx(taus,nanmean(corrs(:,:,d),1),'k-','LineWidth',5)
        xlabel('\tau','FontSize',25)
        ylabel('g2','FontSize',25)
        title(['Detector ' num2str(d) ],'FontSize',25)
        ylim([0.9 1.6])
        xlim([0 1e-2])
        set(gca,'FontSize',20)
        h=line([taus(datalength) taus(datalength)],[0.9 1.6]);
        set(h,'Color',[0 0 0]);
        set(gcf,'PaperPositionMode','Auto')
        saveas(gcf,['../' patID '/' patID 'notes/savedfigs/corrs_det' num2str(d) '_' patID '_' patdate '.fig'],'fig')
        saveas(gcf,['../' patID '/' patID 'notes/savedfigs/corrs_det' num2str(d) '_' patID '_' patdate '.jpg'],'jpg')
        clear int
    end
end


%Calculate gamma for sigma calculation
%Will use the mean of all correlation curves taken over the range set above.  In
%this way, the curve will be smooth, and it will provide us with
%approximately the right order of magnitude for gamma
if length(usedflowdets==1)
    corrsavgtmp=nanmean(corrs(:,:,usedflowdets),3);%First take avg over all dets
else
    corrsavgtmp=nanmean(corrs(:,:,usedflowdets(1):usedflowdets(2)),3);%First take avg over all dets
end
diff=abs(nanmean(corrsavgtmp(Marksflow(marksforgamma(1)):Marksflow(marksforgamma(2)),:),1)-(1+setbeta*1/exp(1)));
ind=find(diff==min(diff));
gamma=1/taus(min(ind));%use min(ind) in case ind is not a 1x1 vector.

%Now fit each frame of data
% for measnum=0:numframestoavg:size(intensitydata,1)-numframestoavg
% 
%     if fitavg==1
%         %Since we take data from multiple SM flow fibers at the same
%         %location, fits are more accurate if we take the mean of all these
%         %curves
%         
%         intensityavg(measnum+1)=nanmean(nanmean(intensitydata(measnum+1:measnum+numframestoavg,usedflowdets(1):usedflowdets(2)),2),1);
%         corrsavg(measnum+1,:)=nanmean(nanmean(corrs(measnum+1:measnum+numframestoavg,:,usedflowdets(1):usedflowdets(2)),3),1);
%         
%         if isempty(find(isnan(corrsavg(measnum+1,:)))) & ~isnan(mua(measnum+1))
%             foo=min(find(slidingavg(corrsavg(measnum+1,:),avgnum)<=cutoff));
%             if isempty(foo) || foo<startcorr
%                 tmpf=datalength;
%             else
%                 tmpf=min(find(slidingavg(corrsavg(measnum+1,:),avgnum)<=cutoff));
%             end
%             corrstmpavg=squeeze(corrsavg(measnum+1,startcorr:tmpf));
%             corrstmpavg=slidingavg(corrstmpavg,avgnum);
%             taustmp=taus(startcorr:tmpf);
%             %Calculate noise from Chao's noise model
%             sigma=1./intensityavg(measnum+1).*sqrt(1./T./t).*sqrt(1+mean(setbeta).*exp(-gamma*taus));
% 
% %             if fitbeta==0
% %                 %Dont fit beta, bound Db fit to reasonable values
% %                 Betasaveavg(measnum+1)=mean(setbeta);
% %                 [Dbfitavg(measnum+1),fvalavg(measnum+1),exitflagavg(measnum+1)]=fminsearchbnd(@xg2fitx_new1_withsigma,x0(1),[0],[1e-2],options,...
% %                     Betasaveavg(measnum+1),r,taustmp,muspo,mua(measnum+1),k0,ze,corrstmpavg,length(taustmp),sigma);
% %             elseif fitbeta==2
% %                 %Dont fit beta, bound Db fit to reasonable values
% %                 Betasaveavg(measnum+1)=nanmean([1.5*corrsavg(measnum+1,1) corrsavg(measnum+1,2) 0.5*corrsavg(measnum+1,3)])-1;
% %                 [Dbfitavg(measnum+1),fvalavg(measnum+1),exitflagavg(measnum+1)]=fminsearchbnd(@xg2fitx_new1_withsigma,x0(1),[0],[1e-2],options,...
% %                     Betasaveavg(measnum+1),r,taustmp,muspo,mua(measnum+1),k0,ze,corrstmpavg,length(taustmp),sigma);
% %             else
% %                 %Fit both Beta and BFI at the same time, bound Beta fit to
% %                 %reasonable values
% %                 [betaDbfitavg(:,measnum+1),fvalavg(measnum+1),exitflagavg(measnum+1)]=fminsearchbnd(@xg2fitx_betaandDB_new1_withsigma,x0,...
% %                     [0 betafitmin],[1e-1 betafitmax],options,r,taustmp,muspo,mua(measnum+1),k0,ze,corrstmpavg,length(taustmp),sigma);
% %                 Betasaveavg(measnum+1)=squeeze(betaDbfitavg(2,measnum+1));
% %                 Dbfitavg(measnum+1)=squeeze(betaDbfitavg(1,measnum+1));
% %             end
% 
%             % Modified to NOT fit beta, set to average of first 5 pts of
%             % all frames
%             
%             g2smooth = smooth(corrsavg(measnum+1,:), avgnum);
%             Betasaveavg(measnum+1)=mean(g2smooth(1:5))-1;
%             
%             %Semi-inf fit of g2
%             [Dbfitavg(measnum+1),fvalavg(measnum+1),exitflagavg(measnum+1)]=fminsearchbnd(@xg2fitx_new1_withsigma,x0(1),[0],[1e-2],options,...
%                      Betasaveavg(measnum+1),r,taustmp,muspo,mua(measnum+1),k0,ze,corrstmpavg,length(taustmp),sigma);
%             
%             %2-layer fit of g1
%             zd=1./muspo;
%             g1avg(:,measnum+1)=sqrt(abs((corrsavg(measnum+1,:)-1)./Betasaveavg(measnum+1)));
%             g1tmp = g1avg(1:length(taustmp),measnum+1);
% % 
% %             try
% %                [Dbfit_2layer(measnum+1,i,:),fval_2layer(measnum+1,i),exitflag_2layer(measnum+1,i)]=fminsearchbnd(@(x)xg1fit2layer(x,n0,zd,[mua(measnum+1) mua(measnum+1)],[muspo muspo],lambda*10^7,thickness,r,taustmp,g1tmp),x0,[0 0],[],options);
% %                %xg1fit2layer_betaandDB_withsigma(x0,n0,R,mua,musp,lambda,ell,rho,taustmp,g2,sigma)
% %             catch E
% %                 Dbfit_2layer(measnum+1,i,:)=NaN;
% %                 fval_2layer(measnum+1,i)=NaN;
% %                 exitflag_2layer(measnum+1,i)=NaN;
% %             end
%                         
%             Curvefitavg(:,measnum+1)=g1fitx(Dbfitavg(measnum+1),r,taus,muspo,mua(measnum+1),k0,ze);
%             Curvefitg2avg(:,measnum+1)=g2fitx([Dbfitavg(measnum+1) Betasaveavg(measnum+1)],r,taus,muspo,mua(measnum+1),k0,ze);
%             %Calculate error in fit
%             indtmp(measnum+1)=min(find(abs(squeeze(Curvefitavg(:,measnum+1))-0.3)==min(abs(squeeze(Curvefitavg(:,measnum+1))-0.3))));%Use min in case size(ind)>1
%             errorfitavg(:,measnum+1)=(g1avg(:,measnum+1)-squeeze(Curvefitavg(:,measnum+1)))./squeeze(Curvefitavg(:,measnum+1))*100;
%             meanerroravg(measnum+1)=mean(errorfitavg(1:indtmp(measnum+1),measnum+1));
%             stderroravg(measnum+1)=std(errorfitavg(1:indtmp(measnum+1),measnum+1));
%         else
%             indtmp(measnum+1)=NaN;
%             errorfitavg(:,measnum+1)=ones(1,length(taus)).*NaN;
%             meanerroravg(measnum+1)=NaN;
%             stderroravg(measnum+1)=NaN;
%             Betasaveavg(measnum+1)=NaN;
%             Dbfitavg(measnum+1)=NaN;
%             fvalavg(measnum+1)=NaN;
%             exitflagavg(measnum+1)=NaN;
%             g1avg(:,measnum+1)=ones(1,length(taus)).*NaN;
%             Curvefitavg(:,measnum+1)=ones(1,length(taus)).*NaN;
%             Curvefitg2avg(:,measnum+1)=ones(1,length(taus)).*NaN;
%         end
%         clear corrstmpavg taustmp
% 
%     else
%         %Fit data from each detector
%         for d=1:length(usedflowdets)
%             i=usedflowdets(d);
%             if isempty(find(isnan(corrs(measnum+1,:,i)))) & ~isnan(mua(measnum+1))
%                 corrsmean=nanmean(corrs(measnum+1:measnum+numframestoavg,:,i),1);
%                 foo=min(find(slidingavg(corrsmean,avgnum)<=cutoff));
%                 if isempty(foo) || foo<startcorr
%                     tmpf(i)=datalength;
%                 else
%                     tmpf(i)=min(find(slidingavg(corrsmean,avgnum)<=cutoff));
%                 end
%                 corrstmp=squeeze(corrsmean(startcorr:tmpf(i)));
%                 corrstmp=slidingavg(corrstmp,avgnum);
%                 taustmp=taus(startcorr:tmpf(i));
%                 %Calculate noise from Chao's noise model
%                 sigma(:,i)=1./intensitydata(measnum+1,i).*sqrt(1./T./t).*sqrt(1+mean(setbeta).*exp(-gamma*taus));
% 
%                 if fitbeta==0
%                     %Dont fit beta, bound Db fit to reasonable values
%                     Betasave(measnum+1,i)=mean(setbeta);
%                     [Dbfit(measnum+1,i),fval(measnum+1,i),exitflag(measnum+1,i)]=fminsearchbnd(@xg2fitx_new1_withsigma,x0(1),[0],[1e-2],options,...
%                         Betasave(measnum+1,i),r,taustmp,muspo,mua(measnum+1),k0,ze,corrstmp,length(taustmp),sigma(:,i));
%                 else
%                     %Fit both Beta and BFI at the same time, bound Beta fit to
%                     %reasonable values
%                     [betaDbfit(:,measnum+1,i),fval(measnum+1,i),exitflag(measnum+1,i)]=fminsearchbnd(@xg2fitx_betaandDB_new1_withsigma,x0,...
%                         [0 betafitmin],[1e-1 betafitmax],options,r,taustmp,muspo,mua(measnum+1),k0,ze,corrstmp,length(corrstmp),sigma(:,i));
%                     Betasave(measnum+1,i)=squeeze(betaDbfit(2,measnum+1,i));
%                     Dbfit(measnum+1,i)=squeeze(betaDbfit(1,measnum+1,i));
%                 end
%                 g1(:,i,measnum+1)=sqrt(abs((corrs(measnum+1:measnum+numframestoavg,:,i)-1)./Betasave(measnum+1,i)));
%                 Curvefit(:,measnum+1,i)=g1fitx(Dbfit(measnum+1,i),r,taus,muspo,mua(measnum+1),k0,ze);
%                 Curvefitg2(:,measnum+1,i)=g2fitx([Dbfit(measnum+1,i) Betasave(measnum+1,i)],r,taus,muspo,mua(measnum+1),k0,ze);
%                 %Calculate error in fit
%                 indtmp(measnum+1,i)=min(find(abs(squeeze(Curvefit(:,measnum+1,i))-0.3)==min(abs(squeeze(Curvefit(:,measnum+1,i))-0.3))));%Use min in case size(ind)>1
%                 errorfit(:,measnum+1,i)=(g1(:,i,measnum+1)-squeeze(Curvefit(:,measnum+1,i)))./squeeze(Curvefit(:,measnum+1,i))*100;
%                 meanerror(measnum+1,i)=mean(errorfit(1:indtmp(measnum+1,i),measnum+1,i));
%                 stderror(measnum+1,i)=std(errorfit(1:indtmp(measnum+1,i),measnum+1,i));
%             else
%                 indtmp(measnum+1,i)=NaN;
%                 errorfit(:,measnum+1,i)=ones(1,length(taus)).*NaN;
%                 meanerror(measnum+1,i)=NaN;
%                 stderror(measnum+1,i)=NaN;
%                 g1(:,i,measnum+1)=ones(1,length(taus)).*NaN;
%                 Curvefit(:,measnum+1,i)=ones(1,length(taus)).*NaN;
%                 Curvefitg2(:,measnum+1,i)=ones(1,length(taus)).*NaN;
%                 Betasave(measnum+1,i)=NaN;
%                 Dbfit(measnum+1,i)=NaN;
%             end
%             clear corrstmp taustmp
%         end
% 
%     end
% end
% 
% if fitavg==0
% 
%     for i=1:length(usedflowdets)
%         %Make sure fit converged
%         ind=find(exitflag(:,usedflowdets(i))==0);
%         Dbfit(ind,usedflowdets(i))=NaN;
%         %Make sure fval is small (0.35 is an arbitrary cutoff)
%         ind1=find(fval(:,usedflowdets(i))>0.35);
%         Dbfit(ind1,usedflowdets(i))=NaN;
%         frames=1:1:length(meanerror);
%         if exist('probeoffmark')
%             for m=1:size(probeoffmark,1)
%                 if probeoffmark(m,1)==999
%                     Dbfit(1:Marksflow(probeoffmark(m,2)),:)=NaN;
%                 elseif probeoffmark(m,2)==999
%                     Dbfit(Marksflow(probeoffmark(m,1)):end,:)=NaN;
%                 else
%                     Dbfit(Marksflow(probeoffmark(m,1)):Marksflow(probeoffmark(m,2)),:)=NaN;
%                 end
%             end
%         end
%                 
%         satisfied=0;
%         while satisfied==0
%             %Plot errors in fits
%             figure,
%             [AX,h1,h2]=plotyy(frames,Dbfit(:,usedflowdets(i)),frames,meanerror(:,usedflowdets(i)));
%             set(get(AX(1),'Ylabel'),'String','BFI','FontSize',25)
%             set(get(AX(2),'Ylabel'),'String','Mean Error (%)','FontSize',25)
%             set(AX(2),'YTick',-50:10:50,'YLim',[-50 50],'FontSize',20)
%             grid on
%             set(AX(1),'FontSize',20)
%             set(h1,'LineWidth',3)
%             set(h2,'LineWidth',3)
%             xlabel('Frame','FontSize',25)
%             title(['Detector ' num2str(usedflowdets(i)) ],'FontSize',30)
%             set(gca,'FontSize',20)
%             set(gcf,'PaperPositionMode','Auto')
%             maxwindows(gcf)
%             saveas(gcf,['../' patID '/' patID 'notes/savedfigs/errorinfit_indivD' num2str(usedflowdets(i)) '_' patID '_' patdate '.fig'],'fig')
%             saveas(gcf,['../' patID '/' patID 'notes/savedfigs/errorinfit_indivD' num2str(usedflowdets(i)) '_' patID '_' patdate '.jpg'],'jpg')
%             %If notice a correlation between erroneous Db values and mean error,
%             %define cutoff for mean error
%             cutoffmeanerror(i)=input('Cutoff mean error value? (All frames with average error above this value will be discarded)  ');
%             goners=find(abs(meanerror(:,i))>cutoffmeanerror(i));
%             Dbfit(goners,i)=NaN;
%             satisfied=input('Happy with this cutoff value? 0=no, 1=yes  ');
%             clear goners
%         end
%     end
% 
%     %Final plot of Dbfit
%     figure,plot(Dbfit,'.-','MarkerSize',20,'LineWidth',3)
%     xlabel('Frame','FontSize',25)
%     ylabel('BFI','FontSize',25)
%     ylim([0 5e-7])
%     axis tight
%     tmplim=get(gca,'YLim');
%     for kkkk=1:length(Marksflow)
%         h=line([Marksflow(kkkk) Marksflow(kkkk)],[tmplim(1) tmplim(2)]);
%         set(h,'Color',[0 0 0]);
%     end
%     set(gca,'FontSize',20)
%     set(gcf,'PaperPositionMode','Auto')
%     maxwindows(gcf)
%     if fitbeta==0
%         saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DbFit_' patID '_' patdate '_fixbeta.fig'],'fig')
%         saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DbFit_' patID '_' patdate '_fixbeta.jpg'],'jpg')
%         
%     else
%         saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DbFit_' patID '_' patdate '_fitbeta.fig'],'fig')
%         saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DbFit_' patID '_' patdate '_fitbeta.jpg'],'jpg')
%         
%     end
%     muao=mua;
%     ff=['save ' patID '_' patdate '_MRI_flow_output_fitindiv_abs_quick.mat muao muspo taus Dbfit Curvefit Curvefitg2 corrs g1 intensitydata Marksflow Betasave numframestoavg fval exitflag usedflowdets fitbeta fitavg startcorr datalength avgnum cutoff'];
%     eval(ff);
% 
% else
%     %Make sure fit converged
%     ind=find(exitflagavg==0);
%     Dbfitavg(ind)=NaN;
%     %Make sure fval is small
%     ind1=find(fvalavg>0.35);
%     Dbfitavg(ind1)=NaN;
%     %FOR THIS PATIENT ONLY!
%     %Dbfitavg(1:Marksflow(1))=NaN;
%     %Dbfitavg(Marksflow(probeoffmark(2,1)):end)=NaN;
%     %{
%     if exist('probeoffmark')
%         for m=1:size(probeoffmark,1)
%             if probeoffmark(m,1)==999
%                 Dbfitavg(1:Marksflow(probeoffmark(m,2)))=NaN;
%             elseif probeoffmark(m,2)==999
%                 Dbfitavg(Marksflow(probeoffmark(m,1)):end)=NaN;
%             else
%                 Dbfitavg(Marksflow(probeoffmark(m,1)):Marksflow(probeoffmark(m,2)))=NaN;
%             end
%         end
%     end
%     %}
% 
%     %Plot errors in fits
%     frames=1:1:length(meanerroravg);
%     figure,
%     [AX,h1,h2]=plotyy(frames,Dbfitavg,frames,meanerroravg);%,'b.-','MarkerSize',20,'LineWidth',3)
%     set(get(AX(1),'Ylabel'),'String','BFI','FontSize',25)
%     set(get(AX(2),'Ylabel'),'String','Mean Error (%)','FontSize',25)
%     set(AX(2),'YTick',-50:10:50,'YLim',[-50 50],'FontSize',20)
%     grid on
%     set(AX(1),'FontSize',20)
%     set(h1,'LineWidth',3, 'LineStyle','.', 'MarkerSize', 20)
%     set(h2,'LineWidth',3)
%     xlabel('Frame','FontSize',25)
%     title('Average of all Detectors','FontSize',30)
%     set(gca,'FontSize',20)
%     set(gcf,'PaperPositionMode','Auto')
%     maxwindows(gcf)
%     axis tight
%     
%     satisfied=0;
%     suggest=0;
%     while satisfied==0
%         %If notice a correlation between erroneous Db values and mean error,
%         %define cutoff for mean error
%         if ~suggest
%             errorquartiles = quantile(abs(nonzeros(meanerroravg)), 99);
%             if cutoffmeanerror_perc
%                 cutoffmeanerror = errorquartiles(100-round(cutoffmeanerror_perc));
%             else
%                 cutoffmeanerror = max(abs(nonzeros(meanerroravg)));
%             end
%             ['Suggested cutoff mean error value: ' num2str(cutoffmeanerror)]
%             suggest = 1;
%         else
%             cutoffmeanerror=input('Cutoff mean error value? (All frames with average error above this value will be discarded)  ');            
%         end
%         goners=find(abs(meanerroravg)>cutoffmeanerror);
%         Dbfitavgtemp  = Dbfitavg;
%         Dbfitavgtemp(goners)=NaN;
%         [AX,h1,h2]=plotyy(frames,Dbfitavgtemp,frames,meanerroravg);
%         set(AX(2),'YTick',-50:10:50,'YLim',[-50 50],'FontSize',20)
%         grid on
%         set(AX(1),'FontSize',20)
%         set(h1,'LineWidth',3, 'LineStyle','.', 'MarkerSize', 20)
%         set(h2,'LineWidth',3)
%         xlabel('Frame','FontSize',25)
%         title('Average of all Detectors','FontSize',30)
%         set(gca,'FontSize',20)
%         set(gcf,'PaperPositionMode','Auto')
%         maxwindows(gcf)
%         axis tight
%         satisfied=input('Happy with this cutoff value? 0=no, 1=yes  ');
%         if satisfied
%             Dbfitavg  = Dbfitavgtemp;
%             saveas(gcf,['../' patID '/' patID 'notes/savedfigs/errorinfit_avg_' patID '_' patdate '.fig'],'fig')
%             saveas(gcf,['../' patID '/' patID 'notes/savedfigs/errorinfit_avg_' patID '_' patdate '.jpg'],'jpg')
%         end
%     end
% 
%     %Final plot of Dbfit
%     figure,plot(Dbfitavg,'.-','MarkerSize',20,'LineWidth',3)
%     xlabel('Frame','FontSize',25)
%     ylabel('BFI','FontSize',25)
%     axis tight
%     ylim([0 1e-7])
%     tmplim=get(gca,'YLim');
%     for kkkk=1:length(Marksflow)
%         h=line([Marksflow(kkkk) Marksflow(kkkk)],[tmplim(1) tmplim(2)]);
%         set(h,'Color',[0 0 0]);
%     end
%     set(gca,'FontSize',20)
%     maxwindows(gcf)
%     set(gcf,'PaperPositionMode','Auto')
%     if fitbeta==0
%         
%         saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DbFitavg_' patID '_' patdate '_fixbeta.fig'],'fig')
%         saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DbFitavg_' patID '_' patdate '_fixbeta.jpg'],'jpg')
%         
%     else
%         saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DbFitavg_' patID '_' patdate '_fitbeta.fig'],'fig')
%         saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DbFitavg_' patID '_' patdate '_fitbeta.jpg'],'jpg')
%     end
%     muao=mua;
%     ff=['save ' patID '_' patdate '_1_flow_output_fitavg_abs_quick.mat timeaxis_flow muao muspo taus Dbfitavg Curvefitavg Curvefitg2avg corrsavg g1avg intensityavg Marksflow numframestoavg fvalavg exitflagavg Betasaveavg usedflowdets fitbeta fitavg startcorr datalength avgnum cutoff'];
%     eval(ff);
% 
% end
% 
% if fitbeta==1
%     if fitavg==0
%         figure,plot(Betasave,'.','MarkerSize',10)
%     elseif fitavg==1
%         figure,plot(Betasaveavg,'.','MarkerSize',10)
%     end
%     xlabel('Frame','FontSize',25)
%     ylabel('\beta','FontSize',25)
%     axis tight
%     ylim([0 1])
%     set(gca,'FontSize',20)
%     set(gcf,'PaperPositionMode','Auto')
%     if fitbeta==1 & fitavg==0
%         saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Beta_' patID '_' patdate '_MRI.fig'],'fig')
%         saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Beta_' patID '_' patdate '_MRI.jpg'],'jpg')
%         
%     elseif fitbeta==1 & fitavg==1
%         saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Betaavg_' patID '_' patdate '_MRI.fig'],'fig')
%         saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Betaavg_' patID '_' patdate '_MRI.jpg'],'jpg')
%         
%     end
% end
% 
