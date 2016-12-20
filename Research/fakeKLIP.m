function [images,centers] = fakeKLIP(ntot,nper)

dur=120;
imdim=81;
stepsize=1;

IWA=imdim/14+5;
OWA=imdim/2.5;
range=OWA-IWA;

reqimages=ceil(ntot/nper);
radii=IWA:range/(ntot-1):OWA;

if rem(ntot,nper)~=0
    radii=[radii zeros(1,nper-rem(ntot,nper))];
end

radii=reshape(radii,reqimages,nper);

disp(horzcat(['There will be ' num2str(reqimages) ' total images.']))

for k=1:reqimages    
    [total, parangs, centers(:,:,k)]=generatefakedata(dur,imdim,stepsize,.5,radii(k,:),nper);
    disp(horzcat(['Finished generating dataset ' num2str(k)]))
    images(:,:,k)=KLIP(dur,imdim,total,parangs);
    disp(horzcat(['Finished  analysing dataset ' num2str(k)]))
end

% for i = 1:reqimages
%     figure(i)
%     image(images(:,:,i),'CDataMapping','scaled')
%     colorbar
% end

end