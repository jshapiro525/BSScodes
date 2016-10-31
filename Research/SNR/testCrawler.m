img_scaled = imresize(Sum,2,'bicubic');
[A,D] = snrMapUpd(img_scaled,4);
figure(1)
image(Sum,'CDataMapping','Scaled')
title('Original Image')
colorbar
figure(2)
image(A,'CDataMapping','Scaled')
colorbar
title('SNR Map')