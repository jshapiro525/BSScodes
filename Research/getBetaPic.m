function [DATA,parangs,Sum] = getBetaPic(trunc,imdim)

parangs=[];
DATA=[];
rotDATA=[];

for i = 291:328
    
    if i~=304
       
        name='S20151106S0';
        exten='_spdc_distorcorr.fits';
        filenumb=i;
        if filenumb<10
            numb=strcat('00',num2str(filenumb));
        elseif filenumb<100
            numb=strcat('0',num2str(filenumb));
        else
            numb=num2str(filenumb);
        end
        filename=strcat(name,numb,exten);
        f{i}=filename;
        
        thisparang = fitsheader(filename,'PAR_ANG');
        parangs =[parangs thisparang];

        data=fitsread(f{i},'image');
        if trunc      
            dim=max(size(data));
            c=floor(dim/2)+1;
            w=floor(imdim/2);
            data=data(c-w:c+w,c-w:c+w,:);
        end
        thissum=sum(data(:,:,24:25),3);
        DATA=cat(3,DATA,thissum);
        rotDATA=cat(3,rotDATA,imrotate(thissum,thisparang,'bicubic','crop'));
    end
    
end

for i=1:imdim
    for j=1:imdim
        for k=1:34
            if isnan(DATA(i,j,k))
                DATA(i,j,k)=0;
            end
        end
    end
end

Sum=sum(rotDATA,3);
end