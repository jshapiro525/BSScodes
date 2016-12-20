 close all;
clear all;
clc;

imdim = 101;
dur=37;

[DATA,parangs,Sum] = getBetaPic(1,imdim);
figure(1)
image(Sum,'CDataMapping','scaled')
colorbar
title('Beta Pic b')

 [KLIPSum,KLIP_res,KLIP_parangs,lambda]=KLIP(length(DATA(1,1,:)),imdim,DATA,parangs);

figure(2)
image(KLIPSum,'CDataMapping','scaled')
colorbar
title('KLIP Result')

figure(3)
plot(diag(lambda))


[finalICA,intermeds,bunch]=ICA(KLIP_res,KLIP_parangs,imdim);

for i=1:length(bunch(1,1,:))
    figure(2+i)
    image(bunch(:,:,i),'CDataMapping','scaled')
    colorbar
    title(horzcat(['Derotation ' num2str(i)]))
end

figure(i+3)
image(finalICA,'CDataMapping','scaled')
colorbar
title('Temporary ICA Result')
%     
% [finalSSA,b,stat]=SSAstuff(dur,imdim,DATA,parangs);
% figure(i+4)
% image(finalSSA,'CDataMapping','scaled')
% colorbar
% title('SSA Result')
% 
% for j=1:length(b(:,1))
%     figure(j)
%     image(reshape(b(j,:),imdim,imdim),'CDataMapping','scaled')
%     colorbar
%     title('Stationary Parts')
% end
% 
% figure()
% plot(b(j,:),'.')
% xlabel('Index')
% ylabel('Value')
% 
% figure()
% image(reshape(stat,imdim,imdim),'CDataMapping','scaled')
% 
% figure()
% plot(stat,'.')
% xlabel('Index')
% ylabel('Value')
% 
% 
% 
% 
% 
