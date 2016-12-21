close all;
clear all;
clc;
dur=42;
stepsize=6;
imdim=75;
numpoints=20;
nper=1;
csize=4;
tic

pointsperpoint=10;

totalpoints=numpoints*pointsperpoint;
percentage = 0;

for i=1:numpoints
   
    for j=1:pointsperpoint
        
        ratiopoisson=(i)/numpoints;

        [total,parangs,centers] = generatefakedata(dur,imdim,stepsize,ratiopoisson,imdim/4.5,nper);
        center = centers*2/sqrt(2*log(2));

        [finalKLIP(:,:,i),KLIP_res,KLIP_parangs]=KLIP(dur,imdim,total,parangs);
        finalICA(:,:,i)=ICA(KLIP_res,KLIP_parangs,imdim);
        finalSSA(:,:,i)=SSAstuff(dur,imdim,total,parangs);
        finalCSP(:,:,i)=cspfunction(dur,imdim,total,parangs);
        clc
        disp(horzcat(['The Analysis is roughly ' num2str(percentage) '% complete.']))

        img_scaled = imresize(finalKLIP(:,:,i),2/sqrt(2*log(2)),'bicubic');
        [klipA(:,:,i),D] = snrMap2(img_scaled,2,center);        
        centKLIP = gaussCent2(finalKLIP(:,:,i), 0.5, [centers(1)-csize centers(1)+csize], [centers(2)-csize centers(2)+csize]);
        errKLIP(j)=norm(centKLIP-centers);
        clc
        disp(horzcat(['The Analysis is roughly ' num2str(percentage) '% complete.']))

        img_scaled = imresize(finalICA(:,:,i),2/sqrt(2*log(2)),'bicubic');
        [icaA(:,:,i),D] = snrMap2(img_scaled,2,center);
        centICA = gaussCent2(finalICA(:,:,i), 0.5, [centers(1)-csize centers(1)+csize], [centers(2)-csize centers(2)+csize]);
        errICA(j)=norm(centICA-centers);
        clc
        disp(horzcat(['The Analysis is roughly ' num2str(percentage) '% complete.']))
        
        img_scaled = imresize(finalSSA(:,:,i),2/sqrt(2*log(2)),'bicubic');
        [ssaA(:,:,i),D] = snrMap2(img_scaled,2,center);
        centSSA = gaussCent2(finalSSA(:,:,i), 0.5, [centers(1)-csize centers(1)+csize], [centers(2)-csize centers(2)+csize]);
        errSSA(j)=norm(centSSA-centers);
        clc
        disp(horzcat(['The Analysis is roughly ' num2str(percentage) '% complete.']))
        
        img_scaled = imresize(finalCSP(:,:,i),2/sqrt(2*log(2)),'bicubic');
        [cspA(:,:,i),D] = snrMap2(img_scaled,2,center);
        centCSP = gaussCent2(finalCSP(:,:,i), 0.5, [centers(1)-csize centers(1)+csize], [centers(2)-csize centers(2)+csize]);
        errCSP(j)=norm(centCSP-centers);
        clc
        disp(horzcat(['The Analysis is roughly ' num2str(percentage) '% complete.']))
        
        klipSNR(j)=max(max(klipA(:,:,i)));
        ssaSNR(j)=max(max(ssaA(:,:,i)));
        icaSNR(j)=max(max(icaA(:,:,i)));
        cspSNR(j)=max(max(cspA(:,:,i)));
        clc
        disp(horzcat(['The Analysis is roughly ' num2str(percentage) '% complete.']))
        
        donepoints=((i-1)*pointsperpoint+j);
        percentage= 100*donepoints/totalpoints;
        disp(horzcat(['The Analysis is ' num2str(percentage) '% complete.']))
        
    end
    
    klipsnr(i)=sum(klipSNR)/pointsperpoint;
    ssasnr(i)=sum(ssaSNR)/pointsperpoint;
    icasnr(i)=sum(icaSNR)/pointsperpoint;
    cspsnr(i)=sum(cspSNR)/pointsperpoint;
    
    errorKLIP(i)=sum(errKLIP)/pointsperpoint;
    errorSSA(i)=sum(errSSA)/pointsperpoint;
    errorICA(i)=sum(errICA)/pointsperpoint;
    errorCSP(i)=sum(errCSP)/pointsperpoint;
    
end

percents=100/numpoints:100/numpoints:100;
figure(1)
plot(percents,klipsnr,'r',percents,ssasnr,'k',percents,icasnr,'b',percents,cspsnr,'g')
xlim([0 100])
xlabel('Percent Poisson-dominated Noise')
ylabel('Average SNR')
legend('KLIP','SSA','ICA','CSP')

figure(2)
plot(percents,errorKLIP,'r',percents,errorSSA,'k',percents,errorICA,'b',percents,errorCSP,'g')
xlim([0 100])
xlabel('Percent Poisson-dominated Noise')
ylabel('Average Astrometric Error (pixels)')
legend('KLIP','SSA','ICA','CSP')

toc