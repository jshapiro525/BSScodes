function [Sum,KLIP_res,KLIP_parangs,lambda]=KLIP(dur,imdim,total,parangs)

p=imdim^2;
K=2;

parfor k=1:dur    
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
    rotsig(:,:,k)=imrotate(Signal(:,:,k),parangs(k),'bicubic','crop');
end

Sum=sum(rotsig,3);

delta=floor(dur/3);

    for i=1:3
        for j=1:delta
            pos=delta*(i-1)+j;
            intermediatederotate(:,:,j)=imrotate(Signal(:,:,pos),parangs(pos)-parangs(pos-j+1),'bicubic','crop');
        end
        tempparangs(i)=parangs(delta*(i-1)+1);
        tempdata(:,:,i)=sum(intermediatederotate,3);
    end
    KLIP_parangs=tempparangs;
    KLIP_res=tempdata;
    finalKLIP=Sum;

end