function [total,rottotal,sum,parangs] = makefakedatathistime(dur,imdim,step,ratiopoisson)

maxsigstr=.08;

mu = [imdim/2 imdim/2];
Sigma = [5*imdim 0; 0 5*imdim];
x1 = 1:1:imdim; x2 = 1:1:imdim;
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu,Sigma);
F = reshape(F,length(x2),length(x1));

base(:,:,1)=F;
for i=1:imdim
    for j=1:imdim
        theta=atan((ceil(imdim/2)-i)/(ceil(imdim/2)-j));
        if sqrt((i-ceil(imdim/2))^2+(j-ceil(imdim/2))^2)<imdim/14
            base(i,j,1)=0;
        end
    end
end
base(:,:,1)=base(:,:,1)/max(max(base(:,:,1)));

for k=2:dur
    base(:,:,k)=base(:,:,1);
end

mu=[0 0];
Sigma = [.015 0; 0 .015];
x1 = -.4:.1:.4; x2 = -.4:.1:.4;
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu,Sigma);
G = mvnpdf([X1(:) X2(:)],mu,Sigma./2);
F = reshape(F,length(x2),length(x1));
G = reshape(G,length(x2),length(x1));
F = F./max(max(F))*maxsigstr;
G = G./max(max(G))*maxsigstr/2;

move=360/(24*60)*pi/180*step;
theta1=2*pi/3;
theta2=10*pi/2;
theta3=pi/3;
sum=zeros(imdim);

for k=1:dur
    
    i1(k)=ceil(imdim/2)+floor(imdim/10*cos(theta1+(k-1)*move));
    j1(k)=ceil(imdim/2)+floor(imdim/10*sin(theta1+(k-1)*move));
   
    base(i1(k)-4:i1(k)+4,j1(k)-4:j1(k)+4,k)=base(i1(k)-4:i1(k)+4,j1(k)-4:j1(k)+4,k)+F;
    
    pois(:,:,k)=imnoise(base(:,:,k),'poisson');
    multipl(:,:,k)=base(:,:,k)+(pois(:,:,k)-base(:,:,k))*100*ratiopoisson;
    speck(:,:,k)=addspeckle(multipl(:,:,k),imdim);
    total(:,:,k)=multipl(:,:,k)+(speck(:,:,k)-base(:,:,k))*(1-ratiopoisson);
    
    parang=move*180/pi*(k-1);
    parangs(k)=parang;
    rottotal(:,:,k) = imrotate(total(:,:,k),parang*-1,'nearest','crop');
    sum=sum+rottotal(:,:,k);
end

end