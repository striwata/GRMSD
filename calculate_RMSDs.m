function[]=calculate_RMSDs(options)
% path の追加
% add paths
addpath('algo');
addpath('algo/sub');

%% options の確認
%% get parameters from options
try
    method=options.method;
catch
    method='three_points';
end
if ~strcmp(method,'ICP') && ~strcmp(method,'three_points')
    method='three_points';
end
try
    use_label=options.use_label;
catch
    use_label=true;
end
try
    permit_mirror=options.permit_mirror;
catch
    permit_mirror=false;
end
try
    ignore_atom=options.ignore_atom;
catch
    ignore_atom=1;
%%    ignore_atom=[];  % consider H
end
try 
    clus_mode=options.clus_mode;
catch
    clus_mode=false;
end
try
    reduction=options.reduction;
catch
    reduction=true;
end
try
    iter_num=options.iter_num;
catch
    iter_num=4;
end
if ~reduction
    iter_num=1;
end
%% queryとtargetフォルダの中身を検索
%% get data from query and target folders
qInfo=dir('query'); %name属性にファイル名を格納 (store file names to name-attribute)
query_num=length(qInfo)-2; %ファイルの数 (the number of files)
tInfo=dir('target'); %name属性にファイル名を格納[ (store file names to name-attribute)
target_num=length(tInfo)-2; %ファイルの数 (the number of files)

if(query_num==0 || target_num==0)
    fprintf('queryフォルダ　または　targetフォルダにファイルがありません\n')
    fprintf('no files in query and/or target folder(s)\n')
    return
end

%% main処理
%% main process
for q=1:query_num
    for t=1:target_num
        query=qInfo(q+2).name;
        target=tInfo(t+2).name;
        if ~strcmp('csv',query(length(query)-2:length(query))) || ~strcmp('csv',target(length(target)-2:length(target)))
            continue
        end
        
        
        if clus_mode
            if ~strcmp(query,target)
                continue
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % データの読み込み
        % loading data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        data1=importdata(strcat('query/',query));
        data_num_1=data1(1,1);
        mol_num_1=(size(data1,1)-1)/data_num_1;

        data2=importdata(strcat('target/',target));
        data_num_2=data2(1,1);
        mol_num_2=(size(data2,1)-1)/data_num_2;

        data_struct_1(1:mol_num_1)=struct('pos',zeros(3,data_num_1),'label',zeros(1,data_num_1));
        data_struct_2(1:mol_num_2)=struct('pos',zeros(3,data_num_2),'label',zeros(1,data_num_2));

        for i=1:mol_num_1
            data_struct_1(i).pos=data1(data_num_1*(i-1)+2:data_num_1*i+1,1:3)';
            
            final_label=true(1,data_num_1);
            data_struct_1(i).label=data1(data_num_1*(i-1)+2:data_num_1*i+1,4)';
            for ind=ignore_atom
                final_label=final_label&(data_struct_1(i).label~=ind);
            end
            data_struct_1(i).pos=data_struct_1(i).pos(:,final_label);
            data_struct_1(i).label=data_struct_1(i).label(final_label);
            if ~use_label
                data_struct_1(i).label=ones(1,data_num_1);
            end
            
        end

        for i=1:mol_num_2
            data_struct_2(i).pos=data2(data_num_2*(i-1)+2:data_num_2*i+1,1:3)';
            final_label=true(1,data_num_2);
            data_struct_2(i).label=data2(data_num_2*(i-1)+2:data_num_2*i+1,4)';
            for ind=ignore_atom
                final_label=final_label&(data_struct_2(i).label~=ind);
            end
            data_struct_2(i).pos=data_struct_2(i).pos(:,final_label);
            data_struct_2(i).label=data_struct_2(i).label(final_label);
            if ~use_label
                data_struct_2(i).label=ones(1,data_num_2);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % RMSD計算の実行
        % run RMSD calculations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(strcat(query(1:length(query)-4),'//',target(1:length(target)-4),'\n'))
        result=zeros(mol_num_1,mol_num_2);
        tot=mol_num_1*mol_num_2; 
        now_num=0;
        progress=0.1; 
        for i=1:mol_num_1
            for j=1:mol_num_2
                now_num=now_num+1;
                if(now_num/tot>progress)
                    fprintf('*');
                    progress=progress+0.1;
                end
                if strcmp(method,'three_points')
                    result(i,j)=three_points(data_struct_1(i).pos,data_struct_2(j).pos,data_struct_1(i).label,data_struct_2(j).label,permit_mirror,reduction,iter_num);
                else
                    result(i,j)=three_points(data_struct_1(i).pos,data_struct_2(j).pos,data_struct_1(i).label,data_struct_2(j).label,permit_mirror);
                end
            end
        end
        fprintf('\n')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 計算結果の保存
        % store the results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        csvwrite(strcat(query(1:length(query)-4),'_',target(1:length(target)-4),'_result.csv'),result)
    end
end
