function [indx,error_data]=find_ideal_sinter(data)
%% 变量说明
% 输出变量----------------------------------------------------------------
% indx:find_ideal_sinter脚本所得数据情况标记，n个元素的向量,n=size(data,1)
% error_data：find_ideal_sinter脚本所得数据需要修补的位置，cell类型，第一列是位置，第二列是对应温度最高点风箱
% 输入变量----------------------------------------------------------------
% data：源数据， nrows * 27 nrows=size(data,1)

[zb0,indx]=deal(zeros(1,size(data,1)));
error_data=cell(size(data,1),2);
datadf=diff(data,1,2);
datadf_2=diff(datadf,1,2);
for i=1:size(data,1)
    % 温度曲线应该是上凸的、温升点前的温度应该是<100°C的
    % 寻找最高点
    zb0(i) = find(data(i,:) == max(data(i,:)),1);
    % 前17个点的温度低于100°C
    temp_1_17_data=data(i,1:17);
    temp=find(temp_1_17_data>100);
%     if length(temp)==2 & (diff(temp)~=1)
%         con_1=1;
%     elseif length(temp)==1
%          con_1=1;
%     elseif length(temp)==0
%         con_1=1;
%     else
%         con_1=0;
%     end
    con_1=max(temp_1_17_data)<100;
    % 先看单增单减
    con_2=min(datadf(i,zb0(i)-4:zb0(i)-1))>=0;
    con_3=max(datadf(i,zb0(i):26))<=0;
    % 再看二次减
    con_4=max(datadf_2(i,zb0(i)-4:zb0(i)-2))<0;
    con_5=max(datadf_2(i,zb0(i):25))<0;
    % 逐个条件判断
    if isempty(con_5)
        con_5=0;
    end
    if con_1 & con_2 & con_3 & con_4 & con_5
        indx(i)=5;
        error_data{i,1}=[];
    end
    % 不满足上凸约束
    if con_1 & con_2 & con_3 & (con_4==0 | con_5==0)
        indx(i)=4;
        error_data{i,1}=[];
    end
    % 左右两侧均不满足单增单减约束
    if  con_1 & (con_2==0) & (con_3==0)
        indx(i)=3;
        error_data{i,1}=[find(datadf(i,zb0(i):26)>0)+zb0(i)-1 find(datadf(i,zb0(i)-4:zb0(i)-1)<0)+zb0(i)-5];
    end
    % 只右侧不满足单减约束
    if  con_1 & con_2 & (con_3==0)
        indx(i)=2;
        error_data{i,1}=find(datadf(i,zb0(i):26)>0)+zb0(i)-1;
    end
    % 只左侧不满足单增约束
    if  con_1 & (con_2==0) & con_3
        %原始数据中做标记
        indx(i)=1;
        %找到原始数据中不符合约束的点
        error_data{i,1}=find(datadf(i,zb0(i)-4:zb0(i)-1)<0)+(zb0(i)-4);
    end
    error_data{i,2}=zb0(i);
end
end