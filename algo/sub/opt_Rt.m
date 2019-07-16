function[opt_val]=opt_Rt(AA,BB,data_num,permit_mirror)

A=repmat(AA,1);
B=repmat(BB,1);

Ag=sum(A,2)/data_num;
Bg=sum(B,2)/data_num;
A=A-Ag;
B=B-Bg;

[U,~,V]=svd(B*A');
if permit_mirror
    R=V*U';
else
    hoge=eye(3);
    hoge(3,3)=det(U)*det(V);
    R=V*hoge*U';
end

opt_val=norm(A-R*B)/(data_num^0.5);
end




