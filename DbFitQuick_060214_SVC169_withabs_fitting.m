
%Now fit each frame of data
for measnum=0:numframestoavg:size(intensitydata,1)-numframestoavg

    if fitavg==1
        %Since we take data from multiple SM flow fibers at the same
        %location, fits are more accurate if we take the mean of all these
        %curves
        
        intensityavg(measnum+1)=nanmean(nanmean(intensitydata(measnum+1:measnum+numframestoavg,usedflowdets(1):usedflowdets(2)),2),1);
        corrsavg(measnum+1,:)=nanmean(nanmean(corrs(measnum+1:measnum+numframestoavg,:,usedflowdets(1):usedflowdets(2)),3),1);
        
        if isempty(find(isnan(corrsavg(measnum+1,:)))) & ~isnan(mua(measnum+1))
            foo=min(find(slidingavg(corrsavg(measnum+1,:),avgnum)<=cutoff));
            if isempty(foo) || foo<startcorr
                tmpf=datalength;
            else
                tmpf=min(find(slidingavg(corrsavg(measnum+1,:),avgnum)<=cutoff));
            end
            corrstmpavg=squeeze(corrsavg(measnum+1,startcorr:tmpf));
            corrstmpavg=slidingavg(corrstmpavg,avgnum);
            taustmp=taus(startcorr:tmpf);
            %Calculate noise from Chao's noise model
            sigma=1./intensityavg(measnum+1).*sqrt(1./T./t).*sqrt(1+mean(setbeta).*exp(-gamma*taus));

%             if fitbeta==0
%                 %Dont fit beta, bound Db fit to reasonable values
%                 Betasaveavg(measnum+1)=mean(setbeta);
%                 [Dbfitavg(measnum+1),fvalavg(measnum+1),exitflagavg(measnum+1)]=fminsearchbnd(@xg2fitx_new1_withsigma,x0(1),[0],[1e-2],options,...
%                     Betasaveavg(measnum+1),r,taustmp,muspo,mua(measnum+1),k0,ze,corrstmpavg,length(taustmp),sigma);
%             elseif fitbeta==2
%                 %Dont fit beta, bound Db fit to reasonable values
%                 Betasaveavg(measnum+1)=nanmean([1.5*corrsavg(measnum+1,1) corrsavg(measnum+1,2) 0.5*corrsavg(measnum+1,3)])-1;
%                 [Dbfitavg(measnum+1),fvalavg(measnum+1),exitflagavg(measnum+1)]=fminsearchbnd(@xg2fitx_new1_withsigma,x0(1),[0],[1e-2],options,...
%                     Betasaveavg(measnum+1),r,taustmp,muspo,mua(measnum+1),k0,ze,corrstmpavg,length(taustmp),sigma);
%             else
%                 %Fit both Beta and BFI at the same time, bound Beta fit to
%                 %reasonable values
%                 [betaDbfitavg(:,measnum+1),fvalavg(measnum+1),exitflagavg(measnum+1)]=fminsearchbnd(@xg2fitx_betaandDB_new1_withsigma,x0,...
%                     [0 betafitmin],[1e-1 betafitmax],options,r,taustmp,muspo,mua(measnum+1),k0,ze,corrstmpavg,length(taustmp),sigma);
%                 Betasaveavg(measnum+1)=squeeze(betaDbfitavg(2,measnum+1));
%                 Dbfitavg(measnum+1)=squeeze(betaDbfitavg(1,measnum+1));
%             end

            % Modified to NOT fit beta, set to average of first 5 pts of
            % all frames
            
            g2smooth = smooth(corrsavg(measnum+1,:), avgnum);
            Betasaveavg(measnum+1)=mean(g2smooth(1:5))-1;
            
            %Semi-inf fit of g2
            [Dbfitavg(measnum+1),fvalavg(measnum+1),exitflagavg(measnum+1)]=fminsearchbnd(@xg2fitx_new1_withsigma,x0(1),[0],[1e-2],options,...
                     Betasaveavg(measnum+1),r,taustmp,muspo,mua(measnum+1),k0,ze,corrstmpavg,length(taustmp),sigma);
            
            %2-layer fit of g1
            zd=1./muspo;
            g1avg(:,measnum+1)=sqrt(abs((corrsavg(measnum+1,:)-1)./Betasaveavg(measnum+1)));
            g1tmp = g1avg(1:length(taustmp),measnum+1);
% 
%             try
%                [Dbfit_2layer(measnum+1,i,:),fval_2layer(measnum+1,i),exitflag_2layer(measnum+1,i)]=fminsearchbnd(@(x)xg1fit2layer(x,n0,zd,[mua(measnum+1) mua(measnum+1)],[muspo muspo],lambda*10^7,thickness,r,taustmp,g1tmp),x0,[0 0],[],options);
%                %xg1fit2layer_betaandDB_withsigma(x0,n0,R,mua,musp,lambda,ell,rho,taustmp,g2,sigma)
%             catch E
%                 Dbfit_2layer(measnum+1,i,:)=NaN;
%                 fval_2layer(measnum+1,i)=NaN;
%                 exitflag_2layer(measnum+1,i)=NaN;
%             end
                        
            Curvefitavg(:,measnum+1)=g1fitx(Dbfitavg(measnum+1),r,taus,muspo,mua(measnum+1),k0,ze);
            Curvefitg2avg(:,measnum+1)=g2fitx([Dbfitavg(measnum+1) Betasaveavg(measnum+1)],r,taus,muspo,mua(measnum+1),k0,ze);
            %Calculate error in fit
            indtmp(measnum+1)=min(find(abs(squeeze(Curvefitavg(:,measnum+1))-0.3)==min(abs(squeeze(Curvefitavg(:,measnum+1))-0.3))));%Use min in case size(ind)>1
            errorfitavg(:,measnum+1)=(g1avg(:,measnum+1)-squeeze(Curvefitavg(:,measnum+1)))./squeeze(Curvefitavg(:,measnum+1))*100;
            meanerroravg(measnum+1)=mean(errorfitavg(1:indtmp(measnum+1),measnum+1));
            stderroravg(measnum+1)=std(errorfitavg(1:indtmp(measnum+1),measnum+1));
        else
            indtmp(measnum+1)=NaN;
            errorfitavg(:,measnum+1)=ones(1,length(taus)).*NaN;
            meanerroravg(measnum+1)=NaN;
            stderroravg(measnum+1)=NaN;
            Betasaveavg(measnum+1)=NaN;
            Dbfitavg(measnum+1)=NaN;
            fvalavg(measnum+1)=NaN;
            exitflagavg(measnum+1)=NaN;
            g1avg(:,measnum+1)=ones(1,length(taus)).*NaN;
            Curvefitavg(:,measnum+1)=ones(1,length(taus)).*NaN;
            Curvefitg2avg(:,measnum+1)=ones(1,length(taus)).*NaN;
        end
        clear corrstmpavg taustmp

    else
        %Fit data from each detector
        for d=1:length(usedflowdets)
            i=usedflowdets(d);
            if isempty(find(isnan(corrs(measnum+1,:,i)))) & ~isnan(mua(measnum+1))
                corrsmean=nanmean(corrs(measnum+1:measnum+numframestoavg,:,i),1);
                foo=min(find(slidingavg(corrsmean,avgnum)<=cutoff));
                if isempty(foo) || foo<startcorr
                    tmpf(i)=datalength;
                else
                    tmpf(i)=min(find(slidingavg(corrsmean,avgnum)<=cutoff));
                end
                corrstmp=squeeze(corrsmean(startcorr:tmpf(i)));
                corrstmp=slidingavg(corrstmp,avgnum);
                taustmp=taus(startcorr:tmpf(i));
                %Calculate noise from Chao's noise model
                sigma(:,i)=1./intensitydata(measnum+1,i).*sqrt(1./T./t).*sqrt(1+mean(setbeta).*exp(-gamma*taus));

                if fitbeta==0
                    %Dont fit beta, bound Db fit to reasonable values
                    Betasave(measnum+1,i)=mean(setbeta);
                    [Dbfit(measnum+1,i),fval(measnum+1,i),exitflag(measnum+1,i)]=fminsearchbnd(@xg2fitx_new1_withsigma,x0(1),[0],[1e-2],options,...
                        Betasave(measnum+1,i),r,taustmp,muspo,mua(measnum+1),k0,ze,corrstmp,length(taustmp),sigma(:,i));
                else
                    %Fit both Beta and BFI at the same time, bound Beta fit to
                    %reasonable values
                    [betaDbfit(:,measnum+1,i),fval(measnum+1,i),exitflag(measnum+1,i)]=fminsearchbnd(@xg2fitx_betaandDB_new1_withsigma,x0,...
                        [0 betafitmin],[1e-1 betafitmax],options,r,taustmp,muspo,mua(measnum+1),k0,ze,corrstmp,length(corrstmp),sigma(:,i));
                    Betasave(measnum+1,i)=squeeze(betaDbfit(2,measnum+1,i));
                    Dbfit(measnum+1,i)=squeeze(betaDbfit(1,measnum+1,i));
                end
                g1(:,i,measnum+1)=sqrt(abs((corrs(measnum+1:measnum+numframestoavg,:,i)-1)./Betasave(measnum+1,i)));
                Curvefit(:,measnum+1,i)=g1fitx(Dbfit(measnum+1,i),r,taus,muspo,mua(measnum+1),k0,ze);
                Curvefitg2(:,measnum+1,i)=g2fitx([Dbfit(measnum+1,i) Betasave(measnum+1,i)],r,taus,muspo,mua(measnum+1),k0,ze);
                %Calculate error in fit
                indtmp(measnum+1,i)=min(find(abs(squeeze(Curvefit(:,measnum+1,i))-0.3)==min(abs(squeeze(Curvefit(:,measnum+1,i))-0.3))));%Use min in case size(ind)>1
                errorfit(:,measnum+1,i)=(g1(:,i,measnum+1)-squeeze(Curvefit(:,measnum+1,i)))./squeeze(Curvefit(:,measnum+1,i))*100;
                meanerror(measnum+1,i)=mean(errorfit(1:indtmp(measnum+1,i),measnum+1,i));
                stderror(measnum+1,i)=std(errorfit(1:indtmp(measnum+1,i),measnum+1,i));
            else
                indtmp(measnum+1,i)=NaN;
                errorfit(:,measnum+1,i)=ones(1,length(taus)).*NaN;
                meanerror(measnum+1,i)=NaN;
                stderror(measnum+1,i)=NaN;
                g1(:,i,measnum+1)=ones(1,length(taus)).*NaN;
                Curvefit(:,measnum+1,i)=ones(1,length(taus)).*NaN;
                Curvefitg2(:,measnum+1,i)=ones(1,length(taus)).*NaN;
                Betasave(measnum+1,i)=NaN;
                Dbfit(measnum+1,i)=NaN;
            end
            clear corrstmp taustmp
        end

    end
end