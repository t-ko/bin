%Assume 75% water
close all
clear all

patID='SVC169';
patdate='060214';
ext='';

waterconc=0.75;%Assume 75% water

left=[];
right=[];

fname=[patID '_' patdate ext '_absolute_dpfout.mat'];
load([ fname ])

for j=1:length(grouplabels)
    tmp_left=findstr(char(grouplabels{j}),'eft');
    if ~isempty(tmp_left)
        left=groups{j};
    end
    tmp_right=findstr(char(grouplabels{j}),'ight');
    if ~isempty(tmp_right)
        right=groups{j};
    end
    tmp_calib=findstr(char(grouplabels{j}),'alibration');
    if ~isempty(tmp_calib)
        calib=groups{j};
    end
    tmp_check=findstr(char(grouplabels{j}),'heck');
    if ~isempty(tmp_check)
        check=groups{j};
    end
end

muawater=[0.0047 0.0288 0.0212];%mua h2o 688, 826nm, 786nm
%Adjust mua for water contribution
if ~isempty(left)
    for j=1:length(left)
        leftmua(j,:)=fitmua(left(j),:)-muawater(1:2)*waterconc;
    end
    leftmusp=fitmusp(left,:);
else
    leftmua=NaN(length(right),3);
    leftmusp=NaN(length(right),3);
    hbtmp_left=NaN(length(right),2);
end

if ~isempty(right)
    for j=1:length(right)
        rightmua(j,:)=fitmua(right(j),:)-muawater(1:2)*waterconc;
    end
    rightmusp=fitmusp(right,:);
else
    rightmua=NaN(length(left),3);
    rightmusp=NaN(length(left),3);
    hbtmp_right=NaN(length(left),2);
end
calibmua=nanmean(fitmua(calib,:),1);
checkmua=nanmean(fitmua(check,:),1);
calibmusp=nanmean(fitmusp(calib,:),1);
checkmusp=nanmean(fitmusp(check,:),1);

ext1=[276 2051.96; 974 693.04]*2.303/1e6;%W. B. Gratzer, converted to uM; from Prahl's website
ext786=[740	957.36]*2.303/1e6;%W. B. Gratzer, converted to uM; from Prahl's website [HbO2(786) Hb(786)]

for i=1:size(rightmua)
    
    hbtmp_right(i,:)=inv(ext1)*rightmua(i,1:2).';
   
    Hb_right(i)=hbtmp_right(i,2);
    HbO2_right(i)=hbtmp_right(i,1);
    
    THC_right(i)=sum(hbtmp_right(i,:));
    StO2_right(i)=hbtmp_right(i,1)./THC_right(i)*100;
    
    %Calculate mua at 785 from extinction coef's and HbO2/Hb concentrations
    rightmua(i,3)=[HbO2_right(i) Hb_right(i)]*ext786.'+muawater(3)*waterconc;
    %Calculate musp at 785 by averaging 685 and 830 (crude, I know...)
    rightmusp(i,3)=(rightmusp(i,1)+rightmusp(i,2))/2;
    DPF_right(i,:)=3.*rightmusp(i,:).*2.5./2./(2.5*sqrt(3*rightmua(i,:).*rightmusp(i,:))+1);

end
for i=1:size(leftmua)
    
    hbtmp_left(i,:)=inv(ext1)*leftmua(i,1:2).';
    
    Hb_left(i)=hbtmp_left(i,2);
    HbO2_left(i)=hbtmp_left(i,1);
    
    THC_left(i)=sum(hbtmp_left(i,:));
    StO2_left(i)=hbtmp_left(i,1)./THC_left(i)*100;
    
    %Calculate mua at 785 from extinction coef's and HbO2/Hb concentrations
    leftmua(i,3)=[HbO2_left(i) Hb_left(i)]*ext786.'+muawater(3)*waterconc;
    %Calculate musp at 785 by averaging 685 and 830 (crude, I know...)
    leftmusp(i,3)=(leftmusp(i,1)+leftmusp(i,2))/2;
    DPF_left(i,:)=3.*leftmusp(i,:).*2.5./2./(2.5*sqrt(3*leftmua(i,:).*leftmusp(i,:))+1);
    
end

DPF_left=nanmean(DPF_left,1);
DPF_right=nanmean(DPF_right,1);


figure,subplot(1,2,1)
plot(StO2_left,'.','MarkerSize',40,'LineWidth',3)
hold on,plot(StO2_right,'o','MarkerSize',15,'LineWidth',3)
legend('Left','Right')
ylabel('StO2(%)')
xlim([0.5 length(StO2_right)+0.5])
set(gca,'XTick',1:length(THC_right))
%set(gca,'XTick',t)
xlabel('Repetition No.')
grid on
ylim([35 90])
subplot(1,2,2)
plot(THC_left,'.','MarkerSize',40,'Color',[1 0.6 0],'LineWidth',3)
hold on,plot(Hb_left,'b.','MarkerSize',40,'LineWidth',3)
hold on,plot(HbO2_left,'r.','MarkerSize',40,'LineWidth',3)
hold on,plot(THC_right,'o','Color',[1 0.6 0],'MarkerSize',15,'LineWidth',3)
hold on,plot(Hb_right,'bo','MarkerSize',15,'LineWidth',3)
hold on,plot(HbO2_right,'ro','MarkerSize',15,'LineWidth',3)
xlim([0.5 length(THC_right)+0.5])
ylim([0 150])
set(gca,'XTick',1:length(THC_right))
ylabel('\mu M ')
legend('THC-Left','Hb-Left','HbO2-Left','THC-Right','Hb-Right','HbO2-Right')
xlabel('Repetition No.')
grid on
maxwindows(gcf)
set(findall(gcf,'-property','FontSize'),'FontSize',30)
set(gcf,'PaperPositionMode','Auto')
saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Absoluteoxygenation_' patID '_' patdate ext '.fig'],'fig')
saveas(gcf,['../' patID '/' patID 'notes/savedfigs/Absoluteoxygenation_' patID '_' patdate ext '.jpg'],'jpg')

tt=['save ' patID '_' patdate ext '_baselines.mat leftmua rightmua leftmusp rightmusp HbO2_left HbO2_right Hb_left Hb_right THC_left THC_right StO2_left StO2_right DPF_left DPF_right'];
eval(tt);



ttt=cat(2,nanmean(StO2_right),nanstd(StO2_right),nanmean(StO2_left),nanstd(StO2_left),nanmean(THC_right),nanstd(THC_right),nanmean(THC_left),nanstd(THC_left));
cat(1,ttt,nan(size(ttt)))

muaforxl=nanmean(cat(1,leftmua(:,1:2),rightmua(:,1:2)),1);
muspforxl=nanmean(cat(1,leftmusp(:,1:2),rightmusp(:,1:2)),1);
ttt=cat(2,muaforxl,muspforxl);
cat(1,ttt,nan(size(ttt)))

allHbO2=cat(2,HbO2_right,HbO2_left);
allHb=cat(2,Hb_right,Hb_left);

ttt=cat(2,nanmean(allHb),nanstd(allHb),nanmean(allHbO2),nanstd(allHbO2));
cat(1,ttt,nan(size(ttt)))
