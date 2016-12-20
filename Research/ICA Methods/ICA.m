function [finalICA,signal,Addedtotal]=ICA(data,Parangs,imdim)
p=imdim^2;

for i=1:length(data(1,1,:))
    data(:,:,i)=data(:,:,i)-mean(mean(data(:,:,i)));
    X(i,:)=reshape(data(:,:,i),1,p);
end

M=1/(p-1)*X*X';
mixnorm=M^(-1/2)*X;
seperated=fastica(mixnorm);

for i=1:length(data(1,1,:))
  signal(:,:,i)=reshape(seperated(i,:),imdim,imdim);
end

currentmax=zeros(imdim);
blah=0;
addsig(:,:,1)=signal(:,:,1);
       
for b=1:2
%     if mod(b,2)
%         addsig(:,:,2)=signal(:,:,2)*-1;
%     else addsig(:,:,2)=signal(:,:,2);
    if max(max(signal(:,:,2)))<abs(min(min(signal(:,:,2))))
        addsig(:,:,2)=signal(:,:,2)*-1;
    else addsig(:,:,2)=signal(:,:,2);
    end
    
    for c=1:3
%         if mod(c,2)
%             addsig(:,:,3)=signal(:,:,3)*-1;
%         else addsig(:,:,3)=signal(:,:,3);
        if max(max(signal(:,:,3)))<abs(min(min(signal(:,:,3))))
            addsig(:,:,3)=signal(:,:,3)*-1;
        else addsig(:,:,3)=signal(:,:,3);
        end
        
        par=zeros(1,3);
        pars=Parangs;
        
%         par(3)=pars(ceil(c/2));
%         pars(ceil(c/2))=[];
%         
%         par(2)=pars(ceil(b/2));
%         pars(ceil(b/2))=[];
%             
%         par(1)=pars;
        par(3)=pars(c);
        pars(c)=[];
        
        par(2)=pars(b);
        pars(b)=[];
            
        par(1)=pars;
                                    
        for z=1:3
            Sig(:,:,z)=imrotate(addsig(:,:,z),par(z),'bicubic','crop');
        end
        
        pos=3*(b-1)+c;
            
        Addedtotal(:,:,pos)=sum(Sig,3);
        
        if max(max(Addedtotal(:,:,pos)))>max(max(currentmax))
            currentmax=Addedtotal(:,:,pos);
        end
        blah=blah+1;
                       
    end
end

finalICA=currentmax;

% for i=1:length(data(1,1,:))
%     figure(i)
%     image(signal(:,:,i),'CDataMapping','scaled')
%     colorbar
%     thistitle=horzcat('Independent Channel ',num2str(i));
%     title(thistitle)
% end

end