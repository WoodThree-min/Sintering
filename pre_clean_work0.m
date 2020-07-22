function [pre_clean_data_rev,pre_clean_data_rev_indx,indx,error_data,plat_like_data_indx,max_tem_300_indx]=pre_clean_work0(data)
%% 变量说明
% 输出变量----------------------------------------------------------------
% pre_clean_data_rev：需要修补的数据，nrows * 27，nrows=sum(pre_clean_data_rev_indx)
% pre_clean_data_rev_indx： 需要修补的数据的索引，n个元素的向量,n=size(data,1)
% indx:find_ideal_sinter脚本所得数据情况标记，n个元素的向量,n=size(data,1)
% error_data：find_ideal_sinter脚本所得数据需要修补的位置，cell类型，第一列是位置，第二列是对应温度最高点风箱
% plat_like_data_indx：多个风箱平缓的数据的索引 logical类型，n个元素的向量
% max_tem_300：最高温度在300度以下的数据的索引 logical类型，n个元素的向量
% 输入变量----------------------------------------------------------------
% data：源数据， nrows * 27 nrows=size(data,1)
% zb0：对应温度最高点风箱数据 ，n个元素的向量,n=size(data,1)
% V：台车速度
%% 对数据进行诊断，确定完美数据和待修数据

pre_data=data(size(data,1),:);
[indx,error_data]=find_ideal_sinter_convex(pre_data);


% 返回待修数据和实际索引
pre_clean=((indx<5)&indx~=0);
pre_clean_data=data(pre_clean,:);


%% 排除多个风箱平缓的数据
plat_threshold=40;
zb0 = find(pre_data == max(pre_data),1);

mean_bias = (abs(pre_data(zb0-1)-pre_data(zb0-3))+ ...
    abs(pre_data(zb0-2)-pre_data(zb0-4)))/2;
if mean_bias < plat_threshold %& ((indx<4)&indx~=0)
    plat_like_data_indx = 1;
else
    plat_like_data_indx = 0;
end

plat_like_data_indx=logical(plat_like_data_indx);
plat_like_data=data(plat_like_data_indx,:);
% 移除pre_clean数据中的plat_like_data
pre_clean_data_rev_indx=double(pre_clean)-double(plat_like_data_indx)';
   if pre_clean_data_rev_indx <0
       pre_clean_data_rev_indx=0;
   end

pre_clean_data_rev_indx=logical(pre_clean_data_rev_indx);
pre_clean_data_rev=data(pre_clean_data_rev_indx,:);

%% 排除最高温度在300度以下的
% 温度阈值
tem_threshold=300;
if pre_data(zb0) < tem_threshold %& ((indx<4)&indx~=0)
    max_tem_300_indx=1;
else
    max_tem_300_indx=0;
end

pre_clean_data_rev_indx=double(pre_clean_data_rev_indx)-double(max_tem_300_indx);
max_tem_300_indx=logical(max_tem_300_indx);
   if pre_clean_data_rev_indx<0
       pre_clean_data_rev_indx=0;
   end
pre_clean_data_rev_indx=logical(pre_clean_data_rev_indx);
pre_clean_data_rev=data(pre_clean_data_rev_indx,:);

end