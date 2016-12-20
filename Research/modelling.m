dur=12;
imdim=81;

[total, parangs, centers]=generatefakedata(dur,imdim,10,.5,imdim/3,1);

for i=1:dur
    A(:,i)=reshape(total(:,:,i),imdim^2,1);
end

cvx_begin
    variable P(imdim^2,dur)
    variable B(imdim^2)
    OBJ = 0;
    for i=1:dur
        tempvect=3-(A(:,i)-B-P(:,i))*log(3);
        thissum=0;
        for k=1:imdim^2
            thissum=thissum+tempvect(k);
        end
        OBJ=OBJ+thissum;
    end
    minimize(OBJ)
cvx_end

for i=1:dur
    planets(:,:,i)=reshape(P(:,i),imdim,imdim);
end
    