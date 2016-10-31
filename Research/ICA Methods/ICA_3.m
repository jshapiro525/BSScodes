clear all;
close all;
clc;

imdim=121;
p=imdim^2;


channels=4;
images='stk'; %All or stk
Prepro='PCA'; %PCA or not
datatype='Mine'; %Mine or obs0

[parangs,data] = getdata(images,Prepro,datatype,channels);

for i=1:length(data(1,1,:))
    data(:,:,i)=data(:,:,i)-mean(mean(data(:,:,i)));
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

currentmax=zeros(2);
vals=zeros(1,length(data(1,1,:)));

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
                
%                 for e=1:10
%                     if mod(e,2)
%                         addsig(:,:,5)=signal(:,:,5)*-1;
%                     else addsig(:,:,5)=signal(:,:,5);
%                     end
%             
%                     for f=1:12
%                         if mod(f,2)
%                             addsig(:,:,6)=signal(:,:,6)*-1;
%                         else addsig(:,:,6)=signal(:,:,6);
%                         end
                        par=zeros(1,channels);
                        pars=parangs;
                        
%                         par(6)=pars(ceil(f/2));
%                         pars(ceil(f/2))=[];
%                         
%                         par(5)=pars(ceil(e/2));
%                         pars(ceil(e/2))=[];
                        
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
                        
                        pos=(8*6*4)*(a-1)+(8*6)*(b-1)+8*(c-1)+d;
                        
                        Addedtotal(:,:,pos)=sum(Sig,3);
                        
                        if max(max(Addedtotal(:,:,pos)))>max(max(currentmax))
                            currentmax=Addedtotal(:,:,pos);
                            vals=[a b c d];% e f];
                        end
                        
%                         currentpos=((12)*(e-1)+(12*10)*(d-1)+(12*10*8)*(c-1)+(12*10*8*6)*(b-1)+(12*10*8*6*4)*(a-1));
%                         final=2*4*6*8*10*12;
%                         
%                         if currentpos/final>progress+.001                     
%                             disp(currentpos/final);
%                             progress=currentpos/final;
%                         end
%                     end
%                 end
            end
        end
    end
end

addsig(:,:,end+1)=currentmax;
now=clock;

duration=startime-now;


if images=='stk'
    for i=1:length(data(1,1,:))
        figure(i)
        image(signal(:,:,i),'CDataMapping','scaled')
        colorbar
        thistitle=horzcat('Independent Channel ',num2str(i));
        title(thistitle)
    end
end


fitswrite(Addedtotal,horzcat(['ICA_3' images Prepro datatype '.fits']))
  