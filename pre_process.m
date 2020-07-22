function [cleaned_data2,cleaned_data]=pre_process(data,pos)
% 大致思路：
% 首先分侧进行处理
% 1.
% 2.最小二乘的上凸修补，如果不能修补，则整段数据都进行排除，后续清洗工作也不予进行修补
% 对两侧的处理索引取交集索引，以交集索引对应数据进行平均，再进行一次上述1、2操作

%% 修补-1000,+1000数据
% 温度为负的或大于LIMIT的跑一遍记录索引
LIMIT = 600;
abnormal_tem_indx=zeros(size(data,1),2);
for i=1:size(data,1)
    if max(data(i,:))>LIMIT
        abnormal_tem_indx(i,1)=1;
    end
    if  min(data(i,:))<=0
        abnormal_tem_indx(i,2)=1;
    end
end
abnormal_tem_indx=logical(abnormal_tem_indx);

i=size(data,1);
if abnormal_tem_indx==0
    
    % 温度超过LIMIT的或温度为负数的，以其他数据记录的均值进行代替
    for j=1:size(data,2)
        if data(i,j)>LIMIT
            data(i,j)=mean(data(abnormal_tem_indx(:,1)==0,j));
        end
        if data(i,j)<0
            data(i,j)=mean(data(abnormal_tem_indx(:,2)==0,j));
        end
    end
    
else
    
    if max(data(i,:))>LIMIT
        if data(i,1)>LIMIT
            data(i,1)=data(i,2);
        elseif data(i,size(data,2))>LIMIT
            data(i,size(data,2))=300+10*randn(1);
        end
        
        if max(data(i,2:size(data,2)-1))>LIMIT
            for j=2:size(data,2)-1
                if data(i,j)>LIMIT
                    data(i,j)=mean([data(i,j-1),data(i,j+1)])-10;
                end
            end
        end
    elseif  min(data(i,:))<=0
        if data(i,1)<=0
            data(i,1)=120;
        elseif data(i,size(data,2))<=0
            data(i,size(data,2))=300+10*randn(1);
        end
        if min(data(i,2:size(data,2)-1))<=0
            for j=2:size(data,2)-1
                if data(i,j)<=0
                    data(i,j)=mean([data(i,j-1),data(i,j+1)])-10;
                end
            end
        end
        
    end
end

%% 修补前17温度大于100的风箱的温度记录

% 前17温度大于100的数据,如果不超过3个，而且不是在最后风箱处，那么以前17均值进行补救
LIMIT2 = 100;
temp_1_17_data=data(i,1:17);
temp=find(temp_1_17_data>LIMIT2);
temp_1_17_data(temp)=[];
if max(temp)<17
    
    for j=1:17
        if data(i,j)>LIMIT2
            if j>=1
                data(i,j)=mean(temp_1_17_data);
            end
        end
    end
end

%% 修补数据的V型缺失和单调性
% 一补V型缺失：data_cleansed

data=data_cleansed0(data,1000);

% 补凹凸性：data_cleansed_convex
cleaned_data=data_cleansed_convex(data,pos);

% 二补单调性：clean_the_data0

[pre_clean2_data,pre_clean2_data_indx,indx_s,error_data_s]=pre_clean_work0(cleaned_data);
cleaned_data=clean_the_data0(cleaned_data,error_data_s,pre_clean2_data_indx,pos,indx_s);

[pre_clean_data,pre_clean_data_rev_indx,indx,error_data]=pre_clean_work0(cleaned_data);
cleaned_data2=clean_the_data0(cleaned_data,error_data,pre_clean_data_rev_indx,pos,indx);
%
% [pre_clean2_data,pre_clean_data_rev_indx,indx,error_data2]=pre_clean_work0(cleaned_data);
%
% [cleaned_data2]=clean_the_data0(cleaned_data,error_data2,pre_clean_data_rev_indx,pos,indx);

end

