close all
clear all

patID='SVC169';
patdate='060214';
ext='';

% Remove all data when diffusion sequence occurs.  Enter start/stop marks
diffusionmarks = [11 12];

fname1=[ patID '_' patdate ext '_1_'];
load([ fname1 'flow_output_fitavg.mat']);
if fitbeta==0
    h=open(['../' patID '/' patID 'notes/savedfigs/DbFitavg_' patID '_' patdate ext '_fixbeta.fig']);
else
    h=open(['../' patID '/' patID 'notes/savedfigs/DbFitavg_' patID '_' patdate ext '_fitbeta.fig']);
end

checkingfits=1;
while checkingfits==1
    h1=get(h,'Children'); %get the axes, 3 = flow, 4= oxy
    ['Click on Line to Analyze and press a key']
    pause;

    ['Click on Point to see Correlation Curve and Fit.']
    
    [x,y,r]=selectpoints(gco,1,0);


    figure,semilogx(taus,corrsavg(x,:),'b-','LineWidth',3)
    hold on,semilogx(taus,Curvefitg2avg(:,x),'k--','LineWidth',3)
    ylim([0.9 1.6])
    xlim([0 1e-2])
    text(5e-4,1.3,['I=' num2str(intensityavg(x),'%6.2f') ],'FontSize',16 )
    text(5e-4,1.4,['\beta=' num2str(Betasaveavg(x),'%6.2f') ],'FontSize',16 )
    xlabel('\tau','FontSize',20)
    ylabel('g2','FontSize',20)
    set(gca,'FontSize',20)
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    set(gcf,'PaperPositionMode','Auto')
    %saveas(gcf,['../' patID '/' patID 'notes/savedfigs/fit_taupower_' patID '_' patdate ext '_' num2str(x) '.fig'],'fig')
    %saveas(gcf,['../' patID '/' patID 'notes/savedfigs/fit_taupower_' patID '_' patdate ext '_' num2str(x) '.jpg'],'jpg')

    likefit=input('Satisfied with this fit? 0=no, 1=yes  ');

    if likefit==0
        Dbfitavg(x)=NaN;
    end
    checkingfits=input('Want to look at more fits? 0=no, 1=yes  ');
end
close all

% Remove all data when diffusion sequence occurs.  Enter start/stop marks
Dbfitavg(Marksflow(diffusionmarks(1)):Marksflow(diffusionmarks(2)))=NaN;


%Final plot of Dbfit
figure,plot(Dbfitavg,'.-','MarkerSize',20,'LineWidth',3)
xlabel('Frame','FontSize',25)
ylabel('BFI','FontSize',25)
axis tight
%ylim([0 1e-7])
tmplim=get(gca,'YLim');
for kkkk=1:length(Marksflow)
    h=line([Marksflow(kkkk) Marksflow(kkkk)],[tmplim(1) tmplim(2)]);
    set(h,'Color',[0 0 0]);
end
set(gca,'FontSize',20)
maxwindows(gcf)
set(gcf,'PaperPositionMode','Auto')
if fitbeta==0
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DbFitavg_' patID '_' patdate ext '_fixbeta.fig'],'fig')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DbFitavg_' patID '_' patdate ext '_fixbeta.jpg'],'jpg')
else
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DbFitavg_' patID '_' patdate ext '_fitbeta.fig'],'fig')
    saveas(gcf,['../' patID '/' patID 'notes/savedfigs/DbFitavg_' patID '_' patdate ext '_fitbeta.jpg'],'jpg')
end

ff=['save ' fname1 'flow_output_fitavg.mat timeaxis_flow muao muspo taus Dbfitavg Curvefitavg Curvefitg2avg corrsavg g1avg intensityavg Marksflow numframestoavg fvalavg exitflagavg Betasaveavg usedflowdets fitbeta fitavg startcorr datalength avgnum cutoff'];
eval(ff);

