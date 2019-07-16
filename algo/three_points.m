function[RMSD]=three_points(A,B,label_A,label_B,permit_mirror,reduction,iter_num)
data_num_A=size(A,2);
data_num_B=size(B,2);

if data_num_B<data_num_A
    C=B;
    label_C=label_B;
    data_num_C=data_num_B;
    
    B=A;
    label_B=label_A;
    data_num_B=data_num_A;
    A=C;
    label_A=label_C;
    data_num_A=data_num_C;
end

%data_num_A<=data_num_B‚ð‘z’è
if data_num_A==0 ||data_num_B==0
    RMSD=inf;
    return
end
for i=1:data_num_A
    for j=i+1:data_num_A
        if norm(A(:,i)-A(:,j))<1e-10
            RMSD=inf;
            return 
        end
    end
end
for i=1:data_num_B
    for j=i+1:data_num_B
        if norm(B(:,i)-B(:,j))<1e-10
            RMSD=inf;
            return 
        end
    end
end

if data_num_A==data_num_B
    RMSD=three_points_g(data_num_A,A,B,label_A,label_B,permit_mirror,reduction,iter_num);
    
else
    penalty=zeros(data_num_B,data_num_A);
    for i=1:data_num_B
        for j=1:data_num_A
            if label_B(i)~=label_A(j)
                penalty(i,j)=-Inf;
            end
        end
    end
    
    if reduction
        thre=ceil(data_num_B^0.5);
    else
        thre=data_num_B;
    end
    
    min_val=Inf;
    for try_num=1:iter_num
        i0=randi([1,data_num_A]);
        while true
            i1=randi([1,data_num_A]);
            if i1~=i0
                break
            end
        end
        while true
            i2=randi([1,data_num_A]);
            if i2~=i0 && i2~=i1
                break
            end
        end
        
        AA=A-A(:,i0);
        
        AA_norm_list=zeros(1,data_num_A);
        for i=1:data_num_A
            AA_norm_list(i)=norm(AA(:,i));
        end
        
        a1=AA(:,i1);
        %revision
        a2=AA(:,i2)-(AA(:,i2)'*AA(:,i1))*AA(:,i1)/AA_norm_list(i1)^2;
        %revision
        a3=cross(a1,a2);

        if norm(a3)<1e-10
            continue
        end
        RA=[a1'/norm(a1);a2'/norm(a2);a3'/norm(a3)];
        
        AA=RA*AA;
        Rz=[1,0,0;0,1,0;0,0,-1];
        AAminus=Rz*AA;
        
        for j0=1:data_num_B
            if label_A(i0)~=label_B(j0)
                continue
            end
            
            BB=B-B(:,j0);
            
            BB_norm_list=zeros(1,data_num_B);
            for i=1:data_num_B
                BB_norm_list(i)=norm(BB(:,i));
            end
            
            hoge0=abs(BB_norm_list-AA_norm_list(i1));
            hoge1=abs(BB_norm_list-AA_norm_list(i2));
            [~,sort0]=sort(hoge0);
            [~,sort1]=sort(hoge1);
            
            j1_count=0;
            
            for j1=sort0
                if j1_count>thre
                    break
                end
                if label_B(j1)~=label_A(i1) || j1==j0
                    continue
                else
                    j1_count=j1_count+1;
                    j2_count=0;
                end
                for j2=sort1
                    if j2_count> thre
                        break
                    end
                    if label_B(j2)~=label_A(i2) || j2==j1 || j2==j0
                        continue
                    else
                        j2_count=j2_count+1;
                    end
                    b1=BB(:,j1);
                    %revision
                    b2=BB(:,j2)-(BB(:,j2)'*BB(:,j1))*BB(:,j1)/BB_norm_list(j1)^2;
                    %revision
                    b3=cross(b1,b2);

                    if norm(b3)<1e-10
                       continue
                    end
                    RB=[b1'/norm(b1);b2'/norm(b2);b3'/norm(b3)];
                    
                    BBB=RB*BB;
                    
                    cost=ones(data_num_B,3) * (AA .* AA) - 2 * BBB' * AA + (BBB .* BBB)' * ones(3,data_num_A);
                    %pruning
                    [many_M,many_I]=min(cost-penalty);
                    many_val=sqrt(mean(many_M));
                    if many_val>2*min_val
                        continue;
                    end
                    if length(unique(many_I))==data_num_A
                        rmsd=opt_Rt(A,B(:,many_I),data_num_A,permit_mirror);
                    else
                        T= max_tr_h3(-1*cost+penalty);
                        rmsd=opt_Rt(A,B*T,data_num_A,permit_mirror); 
                    end
                    %pruning
                    if rmsd<min_val
                        min_val=rmsd;
                    end
                    if permit_mirror
                        cost=ones(data_num_B,3) * (AAminus .* AAminus) - 2 * BBB' * AAminus + (BBB .* BBB)' * ones(3,data_num_A);
                        %pruning
                        [many_M,many_I]=min(cost-penalty);
                        many_val=sqrt(mean(many_M));
                        if many_val>2*min_val
                            continue;
                        end
                        if length(unique(many_I))==data_num_A
                            rmsd=opt_Rt(A,B(:,many_I),data_num_A,permit_mirror);
                        else
                            T= max_tr_h3(-1*cost+penalty);
                            rmsd=opt_Rt(A,B*T,data_num_A,permit_mirror); 
                        end
                        %pruning
                        if rmsd<min_val
                            min_val=rmsd;
                        end
                    end
                end
            end
        end
    end
    
    RMSD=min_val;
end
end

        
        
        
