function [stat,planetsignal,nonstatnoise,total,rottotal,sum,speckle,parangs,centers] = makefakedatatruetestspots(dur,imdim,step,showfigs)

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
        %if abs(cos(theta)) < sqrt(3)/2;
            base(i,j,1)=0;
        %end
        if sqrt((i-ceil(imdim/2))^2+(j-ceil(imdim/2))^2)<imdim/14
            base(i,j,1)=0;
            if abs(sin(theta)) < sqrt(3)/2;
                base(i,j,1)=0;
            end
        end
    end
end
stat=base(:,:,1);

for k=2:dur
    base(:,:,k)=base(:,:,1);
end

mu=[0 0];
Sigma = [.015 0; 0 .015];
x1 = -.1:.1:.1; x2 = -.1:.1:.1;
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu,Sigma);
F = reshape(F,length(x2),length(x1));
F = F./max(max(F))*maxsigstr;


move=360/(24*60)*pi/180*step;
theta1=2*pi/3;
theta2=10*pi/2;
theta3=pi/3;
sum=zeros(imdim);

for k=1:dur
    
    i1(k)=ceil(imdim/2)+floor(imdim/10*cos(theta1+(k-1)*move));
    j1(k)=ceil(imdim/2)+floor(imdim/10*sin(theta1+(k-1)*move));
    i2(k)=ceil(imdim/2)+floor(imdim/4*cos(theta2+(k-1)*move));
    j2(k)=ceil(imdim/2)+floor(imdim/4*sin(theta2+(k-1)*move));
    i3(k)=ceil(imdim/2)-floor(imdim/4.5*cos(theta1+(k-1)*move));
    j3(k)=ceil(imdim/2)-floor(imdim/4.5*sin(theta1+(k-1)*move));
    i4(k)=ceil(imdim/2)+floor(imdim/3*cos(theta3+(k-1)*move));
    j4(k)=ceil(imdim/2)+floor(imdim/3*sin(theta3+(k-1)*move));
    i5(k)=ceil(imdim/2)-floor(imdim/6*cos(theta2+(k-1)*move));
    j5(k)=ceil(imdim/2)-floor(imdim/6*sin(theta2+(k-1)*move));
    i6(k)=ceil(imdim/2)-floor(imdim/8*cos(theta3+(k-1)*move));
    j6(k)=ceil(imdim/2)-floor(imdim/8*sin(theta3+(k-1)*move));
    
    base(i1(k)-1:i1(k)+1,j1(k)-1:j1(k)+1,k)=base(i1(k)-1:i1(k)+1,j1(k)-1:j1(k)+1,k)+F;
    base(i2(k)-1:i2(k)+1,j2(k)-1:j2(k)+1,k)=base(i2(k)-1:i2(k)+1,j2(k)-1:j2(k)+1,k)+F;
    base(i3(k)-1:i3(k)+1,j3(k)-1:j3(k)+1,k)=base(i3(k)-1:i3(k)+1,j3(k)-1:j3(k)+1,k)+F;
    base(i4(k)-1:i4(k)+1,j4(k)-1:j4(k)+1,k)=base(i4(k)-1:i4(k)+1,j4(k)-1:j4(k)+1,k)+F;
    base(i5(k)-1:i5(k)+1,j5(k)-1:j5(k)+1,k)=base(i5(k)-1:i5(k)+1,j5(k)-1:j5(k)+1,k)+F;
    base(i6(k)-1:i6(k)+1,j6(k)-1:j6(k)+1,k)=base(i6(k)-1:i6(k)+1,j6(k)-1:j6(k)+1,k)+F;

    
    planetsignal(:,:,k)=base(:,:,k)-stat;
    
    pois(:,:,k)=imnoise(base(:,:,k),'poisson');
    multipl(:,:,k)=imnoise(pois(:,:,k),'speckle',.01);
    total(:,:,k)=addspeckle(multipl(:,:,k),imdim);
    
    nonstatnoise(:,:,k)=total(:,:,k)-base(:,:,k);
    speckle(:,:,k)=total(:,:,k)-multipl(:,:,k);
    
    parang=move*180/pi*(k-1);
    parangs(k)=parang;
    rottotal(:,:,k) = imrotate(total(:,:,k),parang*-1,'nearest','crop');
    sum=sum+rottotal(:,:,k);
end
centers=[i1(1) j1(1);
         i2(1) j2(1);  
         i3(1) j3(1);  
         i4(1) j4(1);  
         i5(1) j5(1);  
         i6(1) j6(1)];  
centers=centers-0.5;
if showfigs
    close all;
    
    figure(1)
    imshow(mat2gray(stat))
    title('Stationary Noise Source')
    
    figure(2)
    imshow(mat2gray(base(:,:,1)))
    title('Stationary Noise and Signal')
    
    figure(3)
    imshow(mat2gray(planetsignal(:,:,1)))
    title('Just Planet Signal')
        
    figure(4)
    imshow(mat2gray(total(:,:,1)-base(:,:,1)))
    title('Just Speckle and Poisson Noise')
    
    for k=1:dur
        figure(k+4)
        imshow(mat2gray(total(:,:,k)));
        title(strcat('Total Image',num2str(k)))
        if showfigs~=2
            figure(6)
            imshow(mat2gray(sum))
            title('All Images Summed')
            break
        end
    end
    
    if showfigs==2        
        figure(dur+5)
        imshow(mat2gray(sum))
        title('All Images Summed')
    end
end

end