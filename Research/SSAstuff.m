function [finalSSA,b,stat]=SSAstuff(dur,imdim,total,parangs)

p=imdim^2;
for i=1:dur
    X(i,:)=reshape(total(:,:,i),1,imdim^2);
end

nonstat=1;

[Ps, Pn, As, An, ssa_results] = ssa(X, dur-nonstat);%,'equal_epochs',floor(p/11));

S_hat=inv([As An])*X;

for i=1:dur
    frames(:,:,i)=reshape(S_hat(i,:),imdim,imdim);
end

b=S_hat(end-nonstat+1:end,:);
b=b./norm(b);
psf=b'*b;

stat = S_hat(end-nonstat,:);

for i=1:dur
    blah=(eye(p)-psf)*reshape(total(:,:,i),imdim^2,1);
    sig(:,:,i)=reshape(blah,imdim,imdim);
    rotsig(:,:,i)=imrotate(sig(:,:,i),parangs(i),'bicubic','crop');
end

finalSSA=sum(rotsig,3);

end