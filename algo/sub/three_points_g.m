function[RMSD]=three_points_g(data_num,A,B,label_A,label_B,permit_mirror,reduction,iter_num)

Ag=sum(A,2)/data_num;
Bg=sum(B,2)/data_num;
A=A-Ag;
B=B-Bg;

A_norm_list=zeros(1,data_num);
B_norm_list=zeros(1,data_num);

for i=1:data_num
    A_norm_list(i)=norm(A(:,i));
    B_norm_list(i)=norm(B(:,i));
    fprintf('A_norm and B_norm');
end

penalty=zeros(data_num);
for i=1:data_num
    for j=1:data_num
        if label_A(i)~=label_B(j)
            penalty(i,j)=-Inf;
        end
    end
end

if reduction
    thre=ceil(data_num^0.5);
else
    thre=data_num;
end

min_val=Inf;

for try_num=1:iter_num
    
    j0=randi([1,data_num]);
    while true
        j1=randi([1,data_num]);
        if j0~=j1
            break
        end
    end

    hoge0=abs(A_norm_list-B_norm_list(j0));
    hoge1=abs(A_norm_list-B_norm_list(j1));
    [~,sort0]=sort(hoge0);
    [~,sort1]=sort(hoge1);

    b0=B(:,j0);
    b1=B(:,j1)-(B(:,j1)'*B(:,j0))*B(:,j0)/(B_norm_list(j0)^2);
    b2=cross(B(:,j0),B(:,j1));

    if norm(b2)<1e-10
      continue
    end

    RB=[b0'/norm(b0);b1'/norm(b1);b2'/norm(b2)];
    BB=RB*B;

    Rz=[1,0,0;0,1,0;0,0,-1];
    BBminus=Rz*B;


    i0_count=0;
    for i0=sort0
        if i0_count>thre
            break
        end
        if label_A(i0)~=label_B(j0)
            continue
        else
            i0_count=i0_count+1;
            i1_count=0;
        end
        for i1=sort1
            if i1_count>thre
                break
            end
            if label_A(i1)~=label_B(j1) ||i0==i1
                continue
            else
                i1_count=i1_count+1;
            end
            
            a0=A(:,i0);
            a1=A(:,i1)-(A(:,i1)'*A(:,i0))*A(:,i0)/(A_norm_list(i0)^2);
            a2=cross(A(:,i0),A(:,i1));

            if norm(a2)<1e-10
              continue
            endif

            RA=[a0'/norm(a0);a1'/norm(a1);a2'/norm(a2)];
            AA=RA*A;
            
            %pruning
            cost=ones(data_num,3) * (AA .* AA) - 2 * BB' * AA + (BB .* BB)' * ones(3,data_num);
            cost=cost';
            [many_M,many_I]=min(cost-penalty);
            many_val=sqrt(mean(many_M));
            if many_val>2*min_val
                continue;
            end
            if length(unique(many_I))==data_num
                 rmsd=opt_Rt(A(:,many_I),B,data_num,permit_mirror); 
            else
                T= max_tr_h(-1*cost+penalty);
                rmsd=opt_Rt(A,B*T,data_num,permit_mirror); 
            end
            %pruning
            
            if rmsd<min_val
                min_val=rmsd;
            end
            if permit_mirror
                %pruning
                cost=ones(data_num,3) * (AA .* AA) - 2 * BBminus' * AA + (BBminus .* BBminus)' * ones(3,data_num);
                cost=cost';
                [many_M,many_I]=min(cost-penalty);
                many_val=sqrt(mean(many_M));
                if many_val>2*min_val
                    continue;
                end
                if length(unique(many_I))==data_num
                    rmsd=opt_Rt(A,B(:,many_I),data_num,permit_mirror);
                else
                    T= max_tr_h(-1*cost+penalty);
                    rmsd=opt_Rt(A,B*T,data_num,permit_mirror); 
                end
                %pruning
                if rmsd<min_val
                    min_val=rmsd;
                end
            end
        end
    end
end
RMSD=min_val;
end
