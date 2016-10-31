function [ images dim ] = imdiv( fullimages, mindim, maxdim )

p=length(fullimages(:,1,1))^2;

for i=1:maxdim-mindim
    if floor(p/(maxdim+1-i)^2)==p/(maxdim+1-i)^2
        pmin=(maxdim+1-i)^2;
        break
    end
end

if not(exist('pmin'))
    disp('no easy image divide, sorry, but it just sucks')
    return
end

dim=sqrt(pmin);

for i = 1:p/pmin
    for k=1:length(fullimages(1,1,:))
        blah=reshape(fullimages(:,:,k),1,p);
        images(i).section(:,:,k)=reshape(blah((i-1)*pmin+1:pmin*i),dim,dim)
    end
end


end

