function [ EM ] = emat( P,i1,i2,i3 )
    EM(:,1) = P(:,i1);
    EM(:,2) = P(:,i2);
    EM(:,3) = P(:,i3);
end

