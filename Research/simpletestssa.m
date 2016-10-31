close all;
clear all;
clc;

dur = 10;

imdim=13;
p=imdim^2;

mover=[.05 .05;.05 .05];

images(:,:,1)=zeros(imdim);

for i=1:imdim
    images(:,i,1)=10/i;
end
images(:,:,1)=images(:,:,1)/(max(max(images(:,:,1))))*.5;

for i=1:dur
    images(:,:,i)=images(:,:,1);
    if i==1
        figure(i)
        image(images(:,:,i),'CDataMapping','Scaled')
        colorbar
        title(horzcat(['Pre Signal Image ' num2str(i)]))
    end
end

for i=2:dur
    images(i-1:i,4:5,i)=images(i-1:i,4:5,i)+mover;
end

 for i=1:dur
    images(:,:,i)=imnoise(images(:,:,i),'gaussian',0,.0005);
    figure(i+1)
    image(images(:,:,i),'CDataMapping','Scaled')
    colorbar
    title(horzcat(['After Signal Image ' num2str(i)]))
 end

for i=1:dur
    X(i,:)=reshape(images(:,:,i),1,imdim^2);
end


[Ps, Pn, As, An, ssa_results] = ssa(X, dur-1,'equal_epochs',15);

S_hat=inv([As An])*X;

for i=1:dur
    frames(:,:,i)=reshape(S_hat(i,:),imdim,imdim);
end

for i=1:dur
    figure(dur+1+i)
    image(frames(:,:,i),'CDataMapping','Scaled')
    colorbar
    title('Reconstructed images')
end

psf=S_hat(10,:)/norm(S_hat(10,:));

for i=2:dur
    newimage=(eye(p)-psf'*psf)*X(2,:)';
end
blah=reshape(newimage,imdim,imdim);
figure()
image(blah,'CDataMapping','Scaled')
colorbar
