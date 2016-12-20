function [total, parangs, centers]=generatefakedata(dur,imdim,stepsize,rpois,radii,nper);

maxsigstr=.15;
radii=radii(radii~=0);
nplanets=length(radii);

mu = [imdim/2 imdim/2];
Sigma = [5*imdim 0; 0 5*imdim];
x1 = 1:1:imdim; x2 = 1:1:imdim;
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu,Sigma);
F = reshape(F,length(x2),length(x1));

base=F;

for i=1:imdim
    for j=1:imdim
        theta=atan((ceil(imdim/2)-i)/(ceil(imdim/2)-j));
        if sqrt((i-ceil(imdim/2))^2+(j-ceil(imdim/2))^2)<imdim/14
            base(i,j)=0;
        end
    end
end

base=base/max(max(base));
move=(1/4)*pi/180*stepsize;

thetas=[0 pi/3 2*pi/3];

for k=1:dur
    for m=1:nplanets
        if m<4
            i(m)=ceil(imdim/2)+radii(m)*cos(thetas(m)+(k-1)*move);
            j(m)=ceil(imdim/2)+radii(m)*sin(thetas(m)+(k-1)*move);
        else
            i(m)=ceil(imdim/2)-radii(m)*cos(thetas(m-3)+(k-1)*move);
            j(m)=ceil(imdim/2)-radii(m)*sin(thetas(m-3)+(k-1)*move);
        end
    end
    
    [injected]=injectplanets(base,i,j,imdim,maxsigstr);
    
    pois=imnoise(injected,'poisson');
    multipl=imnoise(pois,'speckle',.02*rpois);
    total(:,:,k)=addspeckle(multipl,imdim,1-rpois);
    
    parangs(k)=move*180/pi*(k-1);
    
    if k==1
        centers=[i' j'];
    end
end

if nplanets~=nper
    centers=[centers;zeros(nper-nplanets,2)];
end



end