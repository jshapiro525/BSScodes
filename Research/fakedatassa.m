close all;
clear all;
clc;

imdim=121;
p=imdim^2;

channels=20;
images='stk'; %All or stk
Prepro='not'; %PCA or not
datatype='Mine'; %Mine or obs0
%[parangs,images] = getdata(images,Prepro,datatype,channels);
dur=120;
step=2;

channels=dur;
%[stat,planetsignal,images,rottotal,SUM,parangs] = makefakedata2(dur,imdim,step,0);
[stat,planetsignal,nonstatnoise,images,rottotal,SUM,speckle,parangs] = makefakedatatrue(dur,imdim,2,0);
%[images,rottotal,parangs,SUM] = getHR2563(1,imdim);
%dur=34;

figure(1)
image(images(:,:,1),'CDataMapping','Scaled')
colorbar
title('Original Image Example')

for i=1:channels
    %images(:,:,i)=imrotate(images(:,:,i),-1*parangs(i),'bicubic','crop');
    X(i,:)=reshape(images(:,:,i),1,imdim^2);
end

[Ps, Pn, As, An, ssa_results] = ssa(X, channels-1);%,'equal_epochs',floor(p/11));

S_hat=inv([As An])*X;

for i=1:channels
    frames(:,:,i)=reshape(S_hat(i,:),imdim,imdim);
end

for i=channels-3:channels
    figure()
    image(frames(:,:,i),'CDataMapping','Scaled')
    colorbar
    title('Reconstructed images')
end


b=S_hat(end,:);
b=b./norm(b);
psf=b'*b;

for i=1:channels
    blah=(eye(p)-psf)*reshape(images(:,:,i),imdim^2,1);
    sig(:,:,i)=reshape(blah,imdim,imdim);
    rotsig(:,:,i)=imrotate(sig(:,:,i),-1*parangs(i),'bicubic','crop');
end

figure()
image(sum(rotsig,3),'CDataMapping','Scaled')
colorbar
title('After SSA')

for k=1:channels
    i1=centers(1,k);
    j1=centers(2,k);
    i2=centers(3,k);
    j2=centers(4,k);
    closepic1(:,:,k)=sig(i1-6:i1+6,j1-6:j1+6,k);
    closepic2(:,:,k)=sig(i2-6:i2+6,j2-6:j2+6,k);
end

for i=k:channels
    X1(k,:)=reshape(closepic1(:,:,k),1,13^2);
    X2(k,:)=reshape(closepic2(:,:,k),1,13^2);
end

[Ps1, Pn1, As1, An1, ssa_results1] = ssa(X1, channels-1);
[Ps2, Pn2, As2, An2, ssa_results2] = ssa(X2, channels-1);

S_hat1=inv([As1 An1])*X1;
S_hat2=inv([As2 An2])*X2;

for i=1:channels
    frames1(:,:,i)=reshape(S_hat1(i,:),13,13);
    frames2(:,:,i)=reshape(S_hat2(i,:),13,13);
end

for i=channels-3:channels
    figure()
    image(frames2(:,:,i),'CDataMapping','Scaled')
    colorbar
    title('Reconstructed Focused Images')
end
