clear all;
close all;
clc;

imdim=121;
p=imdim^2;


channels=4;
images='stk'; %All or stk
Prepro='PCA'; %PCA or not
datatype='Mine'; %Mine or obs0

%[data,rottotal,parangs,Sum] = getHR2563(0,imdim);
%[stat,planetsignal,data,rottotal,Sum,parangs] = makefakedata2(4,255,15,0);

[parangs,data] = getdata(images,Prepro,datatype,channels);

for i=1:length(data(1,1,:))
    data(:,:,i)=data(:,:,i)-mean(mean(data(:,:,i)));
    %data(:,:,i)=imrotate(data(:,:,i),-1*parangs(i),'bicubic','crop');
    X(i,:)=reshape(data(:,:,i),1,p);
end

M=1/(p-1)*X*X';
determin=det(M);
mixnorm=M^(-1/2)*X;
covmixnorm=1/(p-1)*mixnorm*mixnorm.';
seperated=fastica(mixnorm);

for i=1:length(data(1,1,:))
  signal(:,:,i)=reshape(seperated(i,:),imdim,imdim);
  signal(:,:,i)=signal(:,:,i)-imrotate(signal(:,:,i),180,'bicubic','crop');
end

currentmax=zeros(imdim);
for i=1:2
    if j==2
        addsig(:,:,1)=-1*signal(:,:,1);
    else addsig(:,:,1)=signal(:,:,1);
    end
    
    for j=1:2
       if j==2
           addsig(:,:,2)=-1*signal(:,:,2);
       else addsig(:,:,2)=signal(:,:,2);
       end
       
       for k=1:2
           if k==2
               addsig(:,:,3)=-1*signal(:,:,3);
           else addsig(:,:,3)=signal(:,:,3);
           end
           
           for m=1:2
               if m==2
                   addsig(:,:,4)=-1*signal(:,:,4);
               else addsig(:,:,4)=signal(:,:,4);
               end
               
               pos=m+2*(k-1)+4*(j-1)+8*(i-1);
               
               Addedtotal(:,:,pos)=sum(addsig,3);
               
               if max(max(Addedtotal(:,:,pos)))>max(max(currentmax))
                   currentmax=Addedtotal(:,:,pos);
               end
           end
       end
    end
end
Addedtotal(:,:,17)=currentmax;

for i=1:length(data(1,1,:))
    figure(i)
    image(signal(:,:,i),'CDataMapping','scaled')
    colorbar
    thistitle=horzcat('Independent Channel ',num2str(i));
    title(thistitle)
end

fitswrite(Addedtotal,horzcat(['ICA_5' images Prepro datatype '.fits']))