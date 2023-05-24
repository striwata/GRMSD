function[RMSD]=ICP(A,B,label_A,label_B,permit_mirror)
a=size(A);
b=size(B);
data_num_A=a(2);
data_num_B=b(2);
%data_num_A<=data_num_B‚ð‘z’è

penalty=zeros(data_num_B,data_num_A);
for i=1:data_num_B
    for j=1:data_num_A
        if label_B(i)~=label_A(j)
            penalty(i,j)=-Inf;
        end
    end
end

Motion20=s20mentai();
umeyama=eye(3);
umeyama(3,3)=-1;
if data_num_A==data_num_B
    tempF = 0;
    best_T=eye(data_num_A);
    
    Ag=sum(A,2)/data_num_A;
    Bg=sum(B,2)/data_num_A;
    A=A-Ag;
    B=B-Bg;


    for i=1:120
        R0 = Motion20(:,:,i);
        if(det(R0)<0) && ~permit_mirror
            continue
        end
        temp = eye(data_num_A);
        R = R0;
        for iteration=1:10000
            T = max_tr_h(A'*R*B+penalty);
            [U,S,V] = svd(B * T * A');
            if permit_mirror || det(U)*det(V)>0
                R = V * U';
            else
                R=V*umeyama*U';
            end
            if(temp == T)
                break;
            end
            temp = T;
        end
        maxF = trace(S);
        if maxF>tempF
            best_T=temp;
        end
        tempF = max(tempF,maxF);
    end
    RMSD=opt_Rt(A,B*best_T,data_num_A,permit_mirror);
else
    tempF = inf;
    A=A-sum(A,2)/data_num_A;
    for i0=1:data_num_A
        for j0=1:data_num_B
            BB=B-B(:,j0)+A(:,i0);
            for i=1:120
                R0 = Motion20(:,:,i);
                if(det(R0)<0) && ~permit_mirror 
                    continue
                end
                
                R = R0;
                for iteration=1:10000
                    BBB=R*BB;
                    cost=sqrt(ones(data_num_B,3) * (A .* A) - 2 * BBB' * A + (BBB .* BBB)' * ones(3,data_num_A));
                    T= max_tr_h3(-1*cost+penalty);
                    if iteration==1
                        temp=T;
                    end
                    BBB=B*T;
                    BBB=BBB-sum(BBB,2)/data_num_A;
                    [U,~,V] = svd(BBB * A');
                    if permit_mirror || det(U)*det(V)>0
                        R = V * U';
                    else
                        R=V*umeyama*U';
                    end  

                    if(temp == T)
                        break;
                    end
                    temp = T;
                end
                val=opt_Rt(A,B*T,data_num_A,permit_mirror);
                tempF = min(tempF,val);
            end
        end
    end
    RMSD=tempF;
end
