function [cleaned_data2,cleaned_data]=pre_process(data,pos)
% ����˼·��
% ���ȷֲ���д���
% 1.
% 2.��С���˵���͹�޲�����������޲������������ݶ������ų���������ϴ����Ҳ��������޲�
% ������Ĵ�������ȡ�����������Խ���������Ӧ���ݽ���ƽ�����ٽ���һ������1��2����

%% �޲�-1000,+1000����
% �¶�Ϊ���Ļ����LIMIT����һ���¼����
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
    
    % �¶ȳ���LIMIT�Ļ��¶�Ϊ�����ģ����������ݼ�¼�ľ�ֵ���д���
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

%% �޲�ǰ17�¶ȴ���100�ķ�����¶ȼ�¼

% ǰ17�¶ȴ���100������,���������3�������Ҳ����������䴦����ô��ǰ17��ֵ���в���
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

%% �޲����ݵ�V��ȱʧ�͵�����
% һ��V��ȱʧ��data_cleansed

data=data_cleansed0(data,1000);

% ����͹�ԣ�data_cleansed_convex
cleaned_data=data_cleansed_convex(data,pos);

% ���������ԣ�clean_the_data0

[pre_clean2_data,pre_clean2_data_indx,indx_s,error_data_s]=pre_clean_work0(cleaned_data);
cleaned_data=clean_the_data0(cleaned_data,error_data_s,pre_clean2_data_indx,pos,indx_s);

[pre_clean_data,pre_clean_data_rev_indx,indx,error_data]=pre_clean_work0(cleaned_data);
cleaned_data2=clean_the_data0(cleaned_data,error_data,pre_clean_data_rev_indx,pos,indx);
%
% [pre_clean2_data,pre_clean_data_rev_indx,indx,error_data2]=pre_clean_work0(cleaned_data);
%
% [cleaned_data2]=clean_the_data0(cleaned_data,error_data2,pre_clean_data_rev_indx,pos,indx);

end

