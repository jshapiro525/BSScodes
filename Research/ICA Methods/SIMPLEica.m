close all;
clear all;
clc;

dur=6;
imdim=121;
step=24;

[stat,planetsignal,nonstatnoise,total,rottotal,Sum,speckle,parangs] = makefakedata(dur,imdim,step,0);
%[total,parangs]=getobsdata(dur,imdim,step);

for i = 1:dur
    total(:,:,i)=total(:,:,i)-mean(mean(total(:,:,i)));
    X(i,:)=reshape(total(:,:,i),1,imdim^2);
end

M=1/(imdim^2-1)*X*X';
determin=det(M);

mixnorm=M^(-1/2)*X;
covmixnorm=1/(imdim^2-1)*mixnorm*mixnorm.';

seperated=fastica(mixnorm);


for i=1:dur
  signal(:,:,i)=reshape(seperated(i,:),imdim,imdim);
  figure(i)
  image(signal(:,:,i),'CDataMapping','scaled')
  colorbar
  subsig(:,:,i)=signal(:,:,i)-imrotate(signal(:,:,i),180,'bicubic');
end

SUM=sum(subsig,3);

figure(dur+1)
image(SUM,'CDataMapping','Scaled')
colorbar

figure(dur+2)
image(total(:,:,1),'CDataMapping','Scaled')
colorbar
