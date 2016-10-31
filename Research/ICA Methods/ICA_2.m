clear all;
close all;
clc;

imdim=281;
p=imdim^2;


channels=4;
images='stk'; %All or stk
Prepro='not'; %PCA or not
datatype='obs0'; %Mine or obs0

%[parangs,data] = getdata(images,Prepro,datatype,channels);
%[data,rottotal,parangs,Sum] = getHR2563(0,imdim);
[stat,planetsignal,data,rottotal,Sum,parangs] = makefakedata2(10,imdim,10,0);

for i=1:length(data(1,1,:))
    data(:,:,i)=data(:,:,i)-mean(mean(data(:,:,i)));
    data(:,:,i)=imrotate(data(:,:,i),-1*parangs(i),'bicubic','crop');
    X(i,:)=reshape(data(:,:,i),1,p);
end

M=1/(p-1)*X*X';
determin=det(M);
mixnorm=M^(-1/2)*X;
covmixnorm=1/(p-1)*mixnorm*mixnorm.';
seperated=fastica(mixnorm);

for i=1:length(data(1,1,:))
  signal(:,:,i)=reshape(seperated(i,:),imdim,imdim);
end

if images=='stk'
    for i=1:length(data(1,1,:))
        figure(i)
        image(signal(:,:,i),'CDataMapping','scaled')
        colorbar
        thistitle=horzcat('Independent Channel ',num2str(i));
        title(thistitle)
    end
end

fitswrite(signal,horzcat(['HR ICA_2.fits']))