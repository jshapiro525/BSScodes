close all;
clear all;
clc;

dur=40;
imdim=121;
step=3;
p=imdim^2;

Sum=zeros(imdim);

[stat,planetsignal,nonstatnoise,total,rottotal,sum,speckle,parangs] = makefakedatatrue(dur,imdim,step,0);
%[total,parangs]=getobsdata(dur,imdim,step);

channels=120;
images='All'; %All or stk
Prepro='not'; %PCA or not
datatype='obs0'; %Mine or obs0

%[total,rottotal,parangs,Sum] = getHR2563(1,imdim);
% [parangs,total] = getdata(images,Prepro,datatype,channels);
 
K=1;

for k=1:dur    
    t=reshape(total(:,:,k),1,p)';
    tbar=t-mean(t);
    r=[];
    rbar=[];
    for i = 1:k-5
        r(:,i)=reshape(total(:,:,i),1,p)';
    end
    for i = k+5:dur
        r=[r reshape(total(:,:,i),1,p)'];
    end
    
    thiset=length(r(1,:));
    
    for i=1:thiset
        rbar(:,i)=r(:,i)-mean(r(:,i));
    end
    
    
    
    Rbar=rbar';
    S=1/(imdim^2-1)*Rbar*Rbar';
    
    [Phi,lambda] = eig(S);
    
    Z=sqrt(inv(lambda))/sqrt(p-1)*Phi'*Rbar;
    
    ZK=[eye(K) zeros(K,thiset-K)]*Z;
    
    sig(:,k)=(eye(p)-ZK'*ZK)*tbar;
    
    psf=(ZK'*ZK)*tbar;
    
    PSF(:,:,k)=reshape(psf',imdim,imdim);
    
    Signal(:,:,k)=reshape(sig(:,k)',imdim,imdim);
    
    rotsig(:,:,k)=imrotate(Signal(:,:,k),-1*parangs(k),'bicubic','crop');
    
    Sum=Sum+rotsig(:,:,k);
    
    k
   
end
 
figure(1)
image(total(:,:,1),'CDataMapping','Scaled')
colorbar
title('Sample Total Image')

figure(2)
image(Signal(:,:,1),'CDataMapping','Scaled')
colorbar
title('Single Signal')

figure(3)
image(Sum,'CDataMapping','Scaled')
fitswrite(Sum,'PCAKLIP.fits')
colorbar
title('Summed Signal')

figure(4)
image(PSF(:,:,1),'CDataMapping','Scaled')
colorbar
title('Single PSF')


fitswrite(Signal,'thistest2.fits')

