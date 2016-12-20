close all;
clear all;
clc;
dur=42;
stepsize=6;
imdim=81;
numpoints=2;
nper=1;
csize=4;
tic

pointsperpoint=2;

for i=1:numpoints
   
    for j=1:pointsperpoint
        
        ratiopoisson=(i)/numpoints;

        [total,parangs,centers] = generatefakedata(dur,imdim,stepsize,ratiopoisson,imdim/4.5,nper);

        [finalKLIP(:,:,i),KLIP_res,KLIP_parangs]=KLIP(dur,imdim,total,parangs);
        finalICA(:,:,i)=ICA(KLIP_res,KLIP_parangs,imdim);
        finalSSA(:,:,i)=SSAstuff(dur,imdim,total,parangs);
        finalCSP(:,:,i)=cspfunction(dur,imdim,total,parangs);
        clc
        disp('The analyses are done for one point, once')
        toc

        img_scaled = imresize(finalKLIP(:,:,i),2/sqrt(2*log(2)),'bicubic');
        [klipA(:,:,i),D] = snrMap2(img_scaled,2);
        disp('first snr measurement complete')
        toc
        centKLIP = gaussCent2(finalKLIP(:,:,i), 0.5, [centers(1)-csize centers(1)+csize], [centers(2)-csize centers(2)+csize]);
        disp('first centroiding done')
        toc
        errKLIP(j)=norm(centKLIP-centers);

        img_scaled = imresize(finalICA(:,:,i),2/sqrt(2*log(2)),'bicubic');
        [icaA(:,:,i),D] = snrMap2(img_scaled,2);
        centICA = gaussCent2(finalICA(:,:,i), 0.5, [centers(1)-csize centers(1)+csize], [centers(2)-csize centers(2)+csize]);
        errICA(j)=norm(centICA-centers);

        img_scaled = imresize(finalSSA(:,:,i),2/sqrt(2*log(2)),'bicubic');
        [ssaA(:,:,i),D] = snrMap2(img_scaled,2);
        centSSA = gaussCent2(finalSSA(:,:,i), 0.5, [centers(1)-csize centers(1)+csize], [centers(2)-csize centers(2)+csize]);
        errSSA(j)=norm(centSSA-centers);

        img_scaled = imresize(finalCSP(:,:,i),2/sqrt(2*log(2)),'bicubic');
        [cspA(:,:,i),D] = snrMap2(img_scaled,2);
        centCSP = gaussCent2(finalCSP(:,:,i), 0.5, [centers(1)-csize centers(1)+csize], [centers(2)-csize centers(2)+csize]);
        errCSP(j)=norm(centCSP-centers);

        klipSNR(j)=max(max(klipA(:,:,i)));
        ssaSNR(j)=max(max(ssaA(:,:,i)));
        icaSNR(j)=max(max(icaA(:,:,i)));
        cspSNR(j)=max(max(cspA(:,:,i))); 
    end
    
    klipSNR(i)=sum(klipSNR)/pointsperpoint;
    ssaSNR(i)=sum(ssaSNR)/pointsperpoint;
    icaSNR(i)=sum(icaSNR)/pointsperpoint;
    cspSNR(i)=sum(cspSNR)/pointsperpoint;
    
    errorKLIP(i)=sum(errKLIP)/pointsperpoint;
    errorSSA(i)=sum(errSSA)/pointsperpoint;
    errorICA(i)=sum(errICA)/pointsperpoint;
    errorCSP(i)=sum(errCSP)/pointsperpoint;
    
end

percents=100/numpoints:100/numpoints:100;
figure(1)
plot(percents,klipSNR,'r',percents,ssaSNR,'k',percents,icaSNR,'b',percents,cspSNR,'g')
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