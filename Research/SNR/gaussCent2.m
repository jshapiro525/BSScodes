function centroid = gaussCent2(img, sigma, row_bounds, col_bounds)
[rows,cols] = size(img);
img(1:row_bounds(1),:) = 0;
img(row_bounds(2):rows,:) = 0;
img(:,1:col_bounds(1)) = 0;
img(:,col_bounds(2):cols) = 0;
%figure
%image(img,'CDataMapping','Scaled')
f = fittype('a1*exp(-(x-x0)^2/(2*sigmax^2)-(y-y0)^2/(2*sigmay^2))','independent',{'x','y'},'dependent','z');
[c,r] = meshgrid(1:rows,1:cols);
r = r(:); c = c(:); img_alt = img(:);
[M,I] = max(img_alt);
[I1,I2] = ind2sub(size(img),I);
%I1 = round(mean(row_bounds));
%I2 = round(mean(col_bounds));
[sf,gof] = fit([r,c],img_alt,f,'StartPoint',[img(I1,I2),sigma,sigma,I1,I2],'Exclude',img_alt ==0);
%plot(sf,[r,c],img_alt)
coefficients = coeffvalues(sf);
centroid = [coefficients(4),coefficients(5)];