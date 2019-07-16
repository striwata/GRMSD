function [T] = permatrix(i,k,n)
    if k==1
        T = eye(n);
    else
        j = ceil(i/factorial(k-1));
        T = permatrix(j*factorial(k-1)-i+1,k-1,n);
        T(:,k) = T(:,j);
        T(:,j) = zeros(n,1);
        T(k,j) = 1;
    end
end

