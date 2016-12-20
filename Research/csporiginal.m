% close all;
% clear all;
% clc;

imdim = 81;
p = imdim^2;
dur = 37;

stepsize = 120/dur;
ratiopoisson = 0.2;
nper = 2;

[total,parangs,Sum] = getBetaPic(1,imdim);
r=[];
for i = 1:dur
    r(:,i)=reshape(total(:,:,i),p,1);
    r(:,i) = r(:,i)-mean(r(:,i));
end

r=r';

Xh = r(:,1:round(p/2));
Xf = r(:,round(p/2)+1:end);

Rh = Xh*Xh'*1/trace(Xh*Xh');
Rf = Xf*Xf'*1/trace(Xf*Xf');

R = Rh+Rf;

[U0,Sig] = eig(R);
[Sig,inds] = sort(diag(Sig),1,'descend');
Sig = diag(Sig,0);
U0 = U0(:, inds);

invSig = Sig;

P = Sig^(-.5)*U0';

Sh = P*Rh*P';
Sf = P*Rf*P';

[Uh,eigh] = eig(Sh);
[eigh,indsh] = sort(diag(eigh),1,'descend');
eigh = diag(eigh,0);
Uh = Uh(:, indsh);

[Uf,eigf] = eig(Sf);
[eigf,indsf] = sort(diag(eigf),1,'ascend');
eigf = diag(eigf,0);
Uf = Uf(:, indsf);

figure()
image(Uh-Uf,'CDataMapping','scaled')
colorbar
        
blah = [diag(eigh) diag(eigf) diag(eigh)+diag(eigf)]

W = Uh'*P;

sig = -(W'*W)*r; 

for i = 1:dur
    Signal(:,:,i) = reshape(sig(i,:),imdim,imdim);
%     figure(i)
%     image(Signal(:,:,i),'CDataMapping','scaled')
%     colorbar
    rotated = imrotate(Signal(:,:,i),parangs(i),'bicubic','crop');
end

Summed=sum(rotated,3);

figure()
image(Summed,'CDataMapping','scaled')
colorbar

