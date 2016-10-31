function [images,centers] = fakeKLIP(n)

dur=120;
imdim=121;
step=1;

IWA=imdim/14+5;
OWA=imdim/2.5;
range=OWA-IWA;

reqimages=ceil(n/6);
radii=[IWA:range/(n-1):OWA zeros(1,6-rem(n,6))];
radii=reshape(radii,reqimages,6);

disp(horzcat(['There will be ' num2str(reqimages) ' total images.']))

for k=1:reqimages    
    [total, parangs, centers(:,:,k)]=generatefakedata(dur,imdim,step,.5,radii(k,:));
    disp(horzcat(['Finished generating dataset ' num2str(k)]))
    images(:,:,k)=KLIP(dur,imdim,total,parangs);
    disp(horzcat(['Finished  analysing dataset ' num2str(k)]))
end

for i = 1:reqimages
    figure(i)
    image(images(:,:,i),'CDataMapping','scaled')
    colorbar
end

end