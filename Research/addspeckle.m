function [newimage] = addspeckle(tot,imdim)


OWA=50;
[E,u]=ft(tot,1,1,1,2*OWA,2*OWA,1,imdim,imdim);
PSF=E.*conj(E);

for i=1:imdim
    for j=1:imdim
        radius=sqrt((i-round(imdim/2))^2+(j-round(imdim/2))^2);
        if radius<imdim/14 
            PSF(i,j)=0;
%         elseif i==imdim/2+.5
%             PSF(i,j)=PSF(i+1,j);
%         elseif j==imdim/2+.5
%             PSF(i,j)=PSF(i,j+1);
        end
    end
end

PSF=PSF*200000/6;

newimage=PSF+tot;
% figure(1)
% image(PSF,'CDataMapping','Scaled')
% colorbar
% title('psf')
% 
% figure(2)
% image(log10(PSF),'CDataMapping','Scaled')
% title('logpsf')
% colorbar
% 
% figure(3)
% image(newimage,'CDataMapping','Scaled')
% title('newimage')
% colorbar

end
