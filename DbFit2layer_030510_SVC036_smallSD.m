%Dbfit code template- EB 6/26/09
%see labnote book on why I chose some of these parameters
%Right now only works for 1 s-d separation
%Flow source 8 SMALL S-D!!!
clear all
close all

patID='SVC169';
patdate='060214';
%Extension you want added to filename, if nothing, just use ext='';
ext='';

fdir1=['../' patID '/' patID '_' patdate '/'];
fname1=[ patID '_' patdate ext '_1_'];
probeoffmark=[999 1];
baselinemarks=[1 7];
thickness=0.67;

usedflowdets=[2 8];%must be consecutive, cannot be ex. 1,3,4 or 2,4.  must change code if this is the case
marksforgamma=[5 7];%Set a range of marks in which the curve decay is ~ the same order of magnitude.  For hypercapnia study, would want to choose only baseline
collectedo2=1;%Was oxygenation data collected in this study?  If so, will use changes in mua in our Dbfit
syncmark=1;%Only for use with ISS, bc the ISS and DCS data are not lined up and will affect the choice of mua for fitting Db
usesigma=0;

n0=1.4;%index of refraction for tissue
c=2.99792458e10; %speed of light in vacuum,cm/s
vo=c/n0; %speed of light in medium
lambda=7.85e-5;
R=-1.440./n0^2+0.710/n0+0.668+0.0636.*n0;
ze=2./3.*(1+R)./(1-R); %check this number, no need for 1/musp as Cecil takes ca\re of it in the files
muao=0.1;
muspo=10;
Dbr0=1e-8; %background brownian diffusion coeff
D=vo./(3.0.*muspo); %background diffusion coeff
k0=2*pi*n0/lambda; %this is the k0 for flow!
DPF=[4.35 4.2 4.4];
zd=1./muspo;



%%%% FITTING SETTINGS:
fitbeta=2; %=0 do not fit beta,==1 fit beta,==2 use first 3 pts
betafitmin=0.3;
betafitmax=0.6;
op=optimset('fminsearch');
options=optimset(op,'MaxIter',30000,'MaxFunEvals',30000,'TolFun',1.000e-16,'TolX',1.000e-16,'Display','Final');
startcorr=1;% if you are fitting beta, want to use a higher number here (e.g. 9), if not fitting beta, use smallest number possible
datalength=100; %chao's way of cutting data for decay.
avgnum=1;%How many points to average in each curve before fitting
cutoff=1;%where to cut correlation curve off when fitting
setbeta=0.45; %used when calculating sigma
x0=[1e-8 1e-7];%Initial guess for Db and beta

plott=0; %=1 plots the probe map

%%%%PLOT SETTINGS
intenymax=1100;%ylim for intensity plot

%Get probe information
flatbabyprobe;
Ns=size(fsources,1);
Nd=size(fdetectors,1);
sdsep=0.8;
sourceo=[0 0 zd].';
detectoro=[sdsep 0 0].';

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
            if mean(corrs(measnum+1,datalength:datalength+20,d))>1.05
               corrs(measnum+1,:,d)=NaN;
            end
        end
        taus=data(3:end,1);
        snum=snum+1;
    elseif exist(fname)==0 %file does not exist
        morefiles=0;
    end
end
%Get integration time (sec)
t=data(1,1)/1000;
Markstmp=Marksflow;
Marksflow(8)=Marksflow(7)+13;
Marksflow(9:end+1)=Markstmp(8:end);

%Determine bin width for each tau
T=zeros(size(taus));
for indt=1:length(T)-1
    T(indt)=taus(indt+1)-taus(indt);
end


%Define a cutoff intensity to get rid of bad frames
fitavg=0;
numframestoavg=1;
cutoffintensity=1000;

Ibaseline=nanmean(intensitydata(Marksflow(baselinemarks(1)):Marksflow(baselinemarks(2)),usedflowdets),1);
muaseries=-log(intensitydata(:,usedflowdets)./Ibaseline)./(sdsep*DPF(3));
 %If the change in mua goes below -0.09 or above 0.09, this doesnt seem
%realistic, so ignore these data points
ind=find(muaseries>0.09);
muaseries(ind)=NaN;
clear ind 
ind=find(muaseries<-0.09);
muaseries(ind)=NaN;
clear ind

%We have intensity changes from DCS laser, so we can calculate mua changes
%at this wavelength
muao=muao+muaseries;
    
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
for measnum=0:numframestoavg:size(intensitydata,1)-numframestoavg

    i=usedflowdets;
    if isempty(find(isnan(corrs(measnum+1,:,i)))) & ~isnan(muao(measnum+1))
        corrsmean=nanmean(corrs(measnum+1,:,i),1);
        tmpf(i)=datalength;
        corrstmp=squeeze(corrsmean(startcorr:tmpf(i)));
        corrstmp=slidingavg(corrstmp,avgnum);
        taustmp=taus(startcorr:tmpf(i));
        %Fit for Db now
        if fitbeta==2
            %Dont fit beta, bound Db fit to reasonable values
            Betasave(measnum+1,i)=mean([1.5*corrstmp(1) corrstmp(2) 0.5*corrstmp(3)])-1;
            %Convert to g1 using g2 and beta
            g1(:,i,measnum+1)=sqrt(abs((corrs(measnum+1,:,i)-1)./Betasave(measnum+1,i)));
            g1tmp=squeeze(g1(1:length(taustmp),i,measnum+1)).';
            %2-layer fit
            %[Dbfit_2layer(measnum+1,i,:),fval(measnum+1,i),exitflag(measnum+1,i)]=fminsearchbnd(@xg1fit2layer,x0,[0 0],[],options,...
                %n0,zd,muao(measnum+1),muspo,vo,lambda,thickness,sdsep,taustmp,g1tmp);
            [Dbfit_2layer(measnum+1,i,:),fval(measnum+1,i),exitflag(measnum+1,i)]=fminsearchbnd(@(x)xg1fit2layer(x,n0,zd,[muao(measnum+1) muao(measnum+1)],[muspo muspo],lambda*10^7,thickness,sdsep,taustmp,g1tmp),x0,[0 0],[],options);
            %X2=xg1fit2layer(Db0,n0,R,mua,musp,lambda,ell,rho,taustmp,g1)
            %Semi-inf fit
            [Dbfit_semiinf(measnum+1,i),fval(measnum+1,i),exitflag(measnum+1,i)]=fminsearchbnd(@xg2fitx_new1,x0(1),[0],[1e-2],options,...
                Betasave(measnum+1,i),sdsep,taustmp,muspo,muao(measnum+1),k0,ze,corrstmp,length(taustmp));
        end
        %Using fit values for Dbtop and Dbbottom, generate g1 curve from
        %2layer model
        for tindex2=1:length(taus)
            detectedsignal(tindex2)=twolayerdynamicjorgegeneric(1, muao(measnum+1), muspo, muao(measnum+1), muspo,thickness,taus(tindex2),Dbfit_2layer(measnum+1,i,1),Dbfit_2layer(measnum+1,i,2),n0,0,vo,lambda,sourceo,detectoro);
        end
        g1fit_2layer(:,measnum+1,i)=detectedsignal./detectedsignal(1);
        g2fit_2layer(:,measnum+1,i)=(detectedsignal./detectedsignal(1)).^2.*Betasave(measnum+1,i)+1;
        g1fit_semiinf(:,measnum+1,i)=g1fitx(Dbfit_semiinf(measnum+1,i),sdsep,taus,muspo,muao(measnum+1),k0,ze);
        g2fit_semiinf(:,measnum+1,i)=g2fitx([Dbfit_semiinf(measnum+1,i) Betasave(measnum+1,i)],sdsep,taus,muspo,muao(measnum+1),k0,ze);
        
        %{
        figure,subplot('Position',[0.125    0.4    0.35    0.42])
        maxwindows(gcf)
        semilogx(taustmp,g1tmp,'b-','LineWidth',3)
        hold on,semilogx(taustmp,detectedsignal./detectedsignal(1),'r--','LineWidth',3)
        hold on,semilogx(taus,g1fit_semiinf(:,measnum+1,i),'k--','LineWidth',3)
        ylabel('g_{1}(\tau)')
        xlabel('\tau')
        xlim([0 1e-2])
        ylim([-0.1 1.1])
        legend('Data','2layer fit','Semi-inf fit')
        set(findall(gcf,'-property','FontSize'),'FontSize',30)
        subplot('Position',[0.6    0.4    0.35    0.42])
        semilogx(taustmp,corrstmp,'b-','LineWidth',3)
        hold on,semilogx(taustmp,(detectedsignal./detectedsignal(1)).^2.*Betasave(measnum+1,i)+1,'r--','LineWidth',3)
        hold on,semilogx(taus,g2fit_semiinf(:,measnum+1,i),'k--','LineWidth',3)
        ylabel('g_{2}(\tau)')
        xlabel('\tau')
        xlim([0 1e-2])
        ylim([0.9 1.6])
        legend('Data','2layer fit','Semi-inf fit')
        gtext({['\alpha D_{B}(top)=' num2str(Dbfit(measnum+1,i,1),'%6.2e') ' cm^{2}/s' ],['\alpha D_{B}(bottom)=' num2str(Dbfit(measnum+1,i,2),'%6.2e') ' cm^{2}/s' ]},'FontSize',30)
        set(findall(gcf,'-property','FontSize'),'FontSize',30)
        set(gcf,'PaperPositionMode','Auto')
        saveas(gcf,['../' patID '/' patID 'notes/savedfigs/fit2layerandsemiinf_Det' num2str(usedflowdets) '_' patID '_' patdate ext '_smallSD_' num2str(measnum+1) '.fig'],'fig')
        saveas(gcf,['../' patID '/' patID 'notes/savedfigs/fit2layerandsemiinf_Det' num2str(usedflowdets) '_' patID '_' patdate ext '_smallSD_' num2str(measnum+1) '.jpg'],'jpg')
        %}
        
        %Calculate error in fit
        indtmp_semiinf(measnum+1,i)=min(find(abs(squeeze(g1fit_semiinf(:,measnum+1,i))-0.3)==min(abs(squeeze(g1fit_semiinf(:,measnum+1,i))-0.3))));%Use min in case size(ind)>1
        errorfit_semiinf(:,measnum+1,i)=(g1(:,i,measnum+1)-squeeze(g1fit_semiinf(:,measnum+1,i)))./squeeze(g1fit_semiinf(:,measnum+1,i))*100;
        meanerror_semiinf(measnum+1,i)=mean(errorfit_semiinf(1:indtmp_semiinf(measnum+1,i),measnum+1,i));
        stderror_semiinf(measnum+1,i)=std(errorfit_semiinf(1:indtmp_semiinf(measnum+1,i),measnum+1,i));
        indtmp_2layer(measnum+1,i)=min(find(abs(squeeze(g1fit_2layer(:,measnum+1,i))-0.3)==min(abs(squeeze(g1fit_2layer(:,measnum+1,i))-0.3))));%Use min in case size(ind)>1
        errorfit_2layer(:,measnum+1,i)=(g1(:,i,measnum+1)-squeeze(g1fit_2layer(:,measnum+1,i)))./squeeze(g1fit_2layer(:,measnum+1,i))*100;
        meanerror_2layer(measnum+1,i)=mean(errorfit_2layer(1:indtmp_2layer(measnum+1,i),measnum+1,i));
        stderror_2layer(measnum+1,i)=std(errorfit_2layer(1:indtmp_2layer(measnum+1,i),measnum+1,i));
    else
        indtmp_semiinf(measnum+1,i)=NaN;
        errorfit_semiinf(:,measnum+1,i)=ones(1,length(taus)).*NaN;
        meanerror_semiinf(measnum+1,i)=NaN;
        stderror_semiinf(measnum+1,i)=NaN;
        indtmp_2layer(measnum+1,i)=NaN;
        errorfit_2layer(:,measnum+1,i)=ones(1,length(taus)).*NaN;
        meanerror_2layer(measnum+1,i)=NaN;
        stderror_2layer(measnum+1,i)=NaN;
        
        g1(:,i,measnum+1)=ones(1,length(taus)).*NaN;
        g1fit_semiinf(:,measnum+1,i)=ones(1,length(taus)).*NaN;
        g2fit_semiinf(:,measnum+1,i)=ones(1,length(taus)).*NaN;
        g1fit_2layer(:,measnum+1,i)=ones(1,length(taus)).*NaN;
        g2fit_2layer(:,measnum+1,i)=ones(1,length(taus)).*NaN;
        Betasave(measnum+1,i)=NaN;
        Dbfit_semiinf(measnum+1,i)=NaN;
        Dbfit_2layer(measnum+1,i,1:2)=NaN;
    end
    clear corrstmp taustmp g1tmp
end

for i=1:length(usedflowdets)
    %Make sure fit converged
    ind=find(exitflag(:,usedflowdets(i))==0);
    Dbfit_2layer(ind,usedflowdets(i),:)=NaN;
    Dbfit_semiinf(ind,usedflowdets(i))=NaN;
    %Make sure fval is small (0.35 is an arbitrary cutoff)
    %ind1=find(fval(:,usedflowdets(i))>0.35);
    %Dbfit(ind1,usedflowdets(i))=NaN;
    frames=1:1:length(meanerror_2layer);
    if exist('probeoffmark')
        for m=1:size(probeoffmark,1)
            if probeoffmark(m,1)==999
                Dbfit_2layer(1:Marksflow(probeoffmark(m,2)),:,:)=NaN;
                Dbfit_semiinf(1:Marksflow(probeoffmark(m,2)),:)=NaN;
            elseif probeoffmark(m,2)==999
                Dbfit_2layer(Marksflow(probeoffmark(m,1)):end,:,:)=NaN;
                Dbfit_semiinf(Marksflow(probeoffmark(m,1)):end,:)=NaN;
            else
                Dbfit_2layer(Marksflow(probeoffmark(m,1)):Marksflow(probeoffmark(m,2)),:,:)=NaN;
                Dbfit_semiinf(Marksflow(probeoffmark(m,1)):Marksflow(probeoffmark(m,2)),:)=NaN;
            end
        end
    end
    %Error in semiinf fit
    satisfied=0;
    while satisfied==0
        %Plot errors in fits
        figure,
        [AX,h1,h2]=plotyy(frames,Dbfit_semiinf(:,usedflowdets(i)),frames,meanerror_semiinf(:,usedflowdets(i)));
        set(get(AX(1),'Ylabel'),'String','BFI','FontSize',25)
        set(get(AX(2),'Ylabel'),'String','Mean Error (%)','FontSize',25)
        set(AX(2),'YTick',-50:10:50,'YLim',[-50 50],'FontSize',20)
        grid on
        set(AX(1),'FontSize',20)
        set(h1,'LineWidth',3)
        set(h2,'LineWidth',3)
        xlabel('Frame','FontSize',25)
        title(['Detector ' num2str(usedflowdets(i)) ],'FontSize',30)
        set(gca,'FontSize',20)
        set(gcf,'PaperPositionMode','Auto')
        maxwindows(gcf)
        saveas(gcf,['../' patID '/' patID 'notes/savedfigs/errorinsemiinffit_indivD' num2str(usedflowdets(i)) '_' patID '_' patdate ext '_smallSD.fig'],'fig')
        saveas(gcf,['../' patID '/' patID 'notes/savedfigs/errorinsemiinffit_indivD' num2str(usedflowdets(i)) '_' patID '_' patdate ext '_smallSD.jpg'],'jpg')
        %If notice a correlation between erroneous Db values and mean error,
        %define cutoff for mean error
        cutoffmeanerror(i)=input('Cutoff mean error value? (All frames with average error above this value will be discarded)  ');
        goners=find(abs(meanerror_semiinf(:,i))>cutoffmeanerror(i));
        Dbfit_semiinf(goners,i)=NaN;
        satisfied=input('Happy with this cutoff value? 0=no, 1=yes  ');
        clear goners
    end
    %Error in 2layer fit
    satisfied=0;
    while satisfied==0
        %Plot errors in fits
        figure,
        [AX,h1,h2]=plotyy(frames,Dbfit_2layer(:,usedflowdets(i),1),frames,meanerror_2layer(:,usedflowdets(i)));
        hold on,plot(frames,Dbfit_2layer(:,usedflowdets(i),2),'r');
        set(get(AX(1),'Ylabel'),'String','BFI','FontSize',25)
        set(get(AX(2),'Ylabel'),'String','Mean Error (%)','FontSize',25)
        set(AX(2),'YTick',-50:10:50,'YLim',[-50 50],'FontSize',20)
        grid on
        set(AX(1),'FontSize',20)
        set(h1,'LineWidth',3)
        set(h2,'LineWidth',3)
        xlabel('Frame','FontSize',25)
        title(['Detector ' num2str(usedflowdets(i)) ],'FontSize',30)
        set(gca,'FontSize',20)
        set(gcf,'PaperPositionMode','Auto')
        maxwindows(gcf)
        saveas(gcf,['../' patID '/' patID 'notes/savedfigs/errorin2layerfit_indivD' num2str(usedflowdets(i)) '_' patID '_' patdate ext '_smallSD.fig'],'fig')
        saveas(gcf,['../' patID '/' patID 'notes/savedfigs/errorin2layerfit_indivD' num2str(usedflowdets(i)) '_' patID '_' patdate ext '_smallSD.jpg'],'jpg')
        %If notice a correlation between erroneous Db values and mean error,
        %define cutoff for mean error
        cutoffmeanerror(i)=input('Cutoff mean error value? (All frames with average error above this value will be discarded)  ');
        goners=find(abs(meanerror_2layer(:,i))>cutoffmeanerror(i));
        Dbfit_2layer(goners,i,:)=NaN;
        satisfied=input('Happy with this cutoff value? 0=no, 1=yes  ');
        clear goners
    end
end

%Final plot of Dbfit
figure,plot(Dbfit_2layer(:,usedflowdets,1),'b.-','MarkerSize',20,'LineWidth',3)
hold on,plot(Dbfit_2layer(:,usedflowdets,2),'r.-','MarkerSize',20,'LineWidth',3)
hold on,plot(Dbfit_semiinf(:,usedflowdets),'.-','Color',[0 0.5 0],'MarkerSize',20,'LineWidth',3)
xlabel('Frame','FontSize',25)
ylabel('BFI','FontSize',25)
ylim([0 5e-7])
axis tight
legend('Top','Bottom','Semiinf')
tmplim=get(gca,'YLim');
for kkkk=1:length(Marksflow)
    h=line([Marksflow(kkkk) Marksflow(kkkk)],[tmplim(1) tmplim(2)]);
    set(h,'Color',[0 0 0]);
end
set(gca,'FontSize',20)
set(gcf,'PaperPositionMode','Auto')
maxwindows(gcf)
saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DbFit2layer_' patID '_' patdate ext '_smallSD.fig'],'fig')
saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DbFit2layer_' patID '_' patdate ext '_smallSD.jpg'],'jpg')


ff=['save ' fname1 'smallSD_flow_output_fitindiv_semiinf2layer.mat baselinemarks muao muspo taus Dbfit_2layer Dbfit_semiinf g1fit_semiinf g2fit_semiinf g1fit_2layer g2fit_2layer corrs g1 intensitydata Marksflow Betasave numframestoavg fval exitflag usedflowdets fitbeta fitavg startcorr datalength avgnum cutoff'];
eval(ff);

