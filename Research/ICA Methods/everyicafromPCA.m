close all;
clear all;
clc;

%load('gpiPCAresult.mat')
load('myPCAresults.mat')
dur=120;
imdim=121;

channels=6;
delta=dur/channels;


for i = 1:channels
    for j=1:delta
        pos=delta*(i-1)+j;
        intermediatederotate(:,:,j)=imrotate(Signal(:,:,pos),(parangs(pos-j+1)-parangs(pos)),'bicubic','crop');
    end
    baseimages(:,:,i)=sum(intermediatederotate,3);
    baseimages(:,:,i)=baseimages(:,:,i)-mean(mean(baseimages(:,:,i)));
    X(i,:)=reshape(baseimages(:,:,i),1,imdim^2);
    newparangs(i)=parangs(delta*(i-1)+1);
end

M=1/(imdim^2-1)*X*X';
determin=det(M);

mixnorm=M^(-1/2)*X;
covmixnorm=1/(imdim^2-1)*mixnorm*mixnorm.';

seperated=fastica(mixnorm);


for i=1:channels
  signal(:,:,i)=reshape(seperated(i,:),imdim,imdim);
end

currentmax=zeros(2);
vals=zeros(1,dur);

progress=0;



startime=clock;

for a=1:2
    if mod(a,2)
        addsig(:,:,1)=signal(:,:,1)*-1;
    else addsig(:,:,1)=signal(:,:,1);
    end
        
    for b=1:4
        if mod(b,2)
            addsig(:,:,2)=signal(:,:,2)*-1;
        else addsig(:,:,2)=signal(:,:,2);
        end
            
        for c=1:6
            if mod(c,2)
                addsig(:,:,3)=signal(:,:,3)*-1;
            else addsig(:,:,3)=signal(:,:,3);
            end
           
            for d=1:8
                if mod(d,2)
                    addsig(:,:,4)=signal(:,:,4)*-1;
                else addsig(:,:,4)=signal(:,:,4);
                end
                
                for e=1:10
                    if mod(e,2)
                        addsig(:,:,5)=signal(:,:,5)*-1;
                    else addsig(:,:,5)=signal(:,:,5);
                    end
            
                    for f=1:12
                        if mod(f,2)
                            addsig(:,:,6)=signal(:,:,6)*-1;
                        else addsig(:,:,6)=signal(:,:,6);
                        end
                        par=zeros(1,channels);
                        pars=newparangs;
                        
                        par(6)=pars(ceil(f/2));
                        pars(ceil(f/2))=[];
                        
                        par(5)=pars(ceil(e/2));
                        pars(ceil(e/2))=[];
                        
                        par(4)=pars(ceil(d/2));
                        pars(ceil(d/2))=[]; 
                        
                        par(3)=pars(ceil(c/2));
                        pars(ceil(c/2))=[];
                        
                        par(2)=pars(ceil(b/2));
                        pars(ceil(b/2))=[];
                        
                        para=pars;
                                    
                        for z=1:channels
                            Sig(:,:,z)=imrotate(addsig(:,:,z),par(z)*-1,'bicubic','crop');
                        end
                        
                        pos=(12*10*8*6*4)*(a-1)+(12*10*8*6)*(b-1)+(12*10*8)*(c-1)+(12*10)*(d-1)+12*(e-1)+f;
                        
                        Addedtotal(:,:,pos)=sum(Sig,3);
                        
                        if max(max(Addedtotal(:,:,pos)))>max(max(currentmax))
                            currentmax=Addedtotal(:,:,pos);
                            vals=[a b c d];% e f];
                        end
                        
                        currentpos=((12)*(e-1)+(12*10)*(d-1)+(12*10*8)*(c-1)+(12*10*8*6)*(b-1)+(12*10*8*6*4)*(a-1));
                        final=2*4*6*8*10*12;
                        
                        if currentpos/final>progress+.001                     
                            disp(currentpos/final);
                            progress=currentpos/final;
                        end
                    end
                end
            end
        end
    end
end


now=clock;

duration=startime-now;
    
figure(1)
image(currentmax,'CDataMapping','Scaled')
colorbar
title('EveryICA Combined')

figure(2)
image(Signal(:,:,1),'CDataMapping','Scaled')
colorbar
title('Original Observation')

% figure(3)
% image(Sum,'CDataMapping','Scaled')
% colorbar
% title('KLIP Results')

for i=1:channels
  figure(i+3)
  image(signal(:,:,i),'CDataMapping','scaled')
  colorbar
  thistitle=horzcat('Independent Channel ',num2str(i));
  title(thistitle)
end

for i=1:channels
  figure(i+3+channels)
  image(baseimages(:,:,i),'CDataMapping','scaled')
  colorbar
  thistitle=horzcat('Base Image (combined PCA result) ',num2str(i));
  title(thistitle)
end