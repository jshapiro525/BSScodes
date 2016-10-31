close all;
clear all;
clc;
dur=60;
step=2;
imdim=121;
numpoints=10;

for i=1:numpoints
    ratiopoisson=(i)/numpoints;
    [total,rottotal,Sum,parangs] = makefakedatathistime(dur,imdim,step,ratiopoisson);
    [finalKLIP(:,:,i),KLip_res,KLIP_parangs]=KLIPstuff(dur,imdim,total,parangs);
    %finalICA(:,:,i)=ICAstuff(KLIP_res,KLIP_parangs);
    finalSSA(:,:,i)=SSAstuff(dur,imdim,step,total,parangs);
    img_scaled = imresize(finalKLIP(:,:,i),2,'bicubic');
    [klipA(:,:,i),D] = snrMapUpd(img_scaled,4);
%     img_scaled = imresize(finalICA(:,:,i),2,'bicubic');
%     [icaA(:,:,i),D] = snrMapUpd(img_scaled,4);
    img_scaled = imresize(finalSSA(:,:,i),2,'bicubic');
    [ssaA(:,:,i),D] = snrMapUpd(img_scaled,4);
    klipSNR(i)=max(max(klipA));
    ssaSNR(i)=max(max(ssaA));
%     icaSNR(i)=max(max(icaA));    
    
end
percents=100/numpoints:numpoints:100;
plot(klipSNR,percents,'r',ssaSNR,percents,'k')%,icaSNR,percents,'b')
xlabel('Percent Poisson-dominated Noise')
Ylabel('SNR')
legend('KLIP','SSA')%,'ICA')