close all;
clear all;
clc;

load('gpiPCAresult.mat')
%load('myPCAresults.mat')

channels=3;
delta=dur/channels;

for i = 1:channels
    for j=1:delta
        pos=delta*(i-1)+j;
        intermediatederotate(:,:,j)=imrotate(Signal(:,:,pos),parangs(pos-j+1)-parangs(pos),'bicubic','crop');
    end
    baseimages(:,:,i)=sum(intermediatederotate,3);
    baseimages(:,:,i)=baseimages(:,:,i)-mean(mean(baseimages(:,:,i)));
    baseimages(:,:,i)=imrotate(baseimages(:,:,i),-1*parangs(delta*(i-1)+1),'bicubic','crop');
    X(i,:)=reshape(baseimages(:,:,i),1,imdim^2);
end

M=1/(imdim^2-1)*X*X';
determin=det(M);

mixnorm=M^(-1/2)*X;
covmixnorm=1/(imdim^2-1)*mixnorm*mixnorm.';

seperated=fastica(mixnorm);


for i=1:channels
  signal(:,:,i)=reshape(seperated(i,:),imdim,imdim);
  figure(i)
  image(signal(:,:,i),'CDataMapping','scaled')
  colorbar
  thistitle=horzcat('Independent Channel ',num2str(i));
  title(thistitle)
  subsig(:,:,i)=signal(:,:,i)-imrotate(signal(:,:,i),180,'bicubic');
end

TOTAL=sum(subsig,3);

for i=1:channels
  figure(i+channels)
  image(baseimages(:,:,i),'CDataMapping','scaled')
  colorbar
  thistitle=horzcat('Base Image (combined PCA result) ',num2str(i));
  title(thistitle)
end

figure(2*channels+1)
image(TOTAL,'CDataMapping','Scaled')
colorbar
title('Combined self-subtracted Channels')
% fitswrite(total,'derotatedica.fits')

figure(2*channels+2)
image(Sum,'CDataMapping','Scaled')
colorbar
title('KLIP Results')

figure(2*channels+3)
image(Sum-imrotate(Sum,180,'bicubic','crop'),'CDataMapping','Scaled')
colorbar
title('KLIP Results Self-Subtracted')

fitswrite(signal,'manimnotsure.fits')
