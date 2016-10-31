function [Sum]=KLIP(dur,imdim,total,parangs)

p=imdim^2;
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
    k
end

Sum=sum(rotsig,3);

end