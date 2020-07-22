function [indx,error_data]=find_ideal_sinter(data)
%% ����˵��
% �������----------------------------------------------------------------
% indx:find_ideal_sinter�ű��������������ǣ�n��Ԫ�ص�����,n=size(data,1)
% error_data��find_ideal_sinter�ű�����������Ҫ�޲���λ�ã�cell���ͣ���һ����λ�ã��ڶ����Ƕ�Ӧ�¶���ߵ����
% �������----------------------------------------------------------------
% data��Դ���ݣ� nrows * 27 nrows=size(data,1)

[zb0,indx]=deal(zeros(1,size(data,1)));
error_data=cell(size(data,1),2);
datadf=diff(data,1,2);
datadf_2=diff(datadf,1,2);
for i=1:size(data,1)
    % �¶�����Ӧ������͹�ġ�������ǰ���¶�Ӧ����<100��C��
    % Ѱ����ߵ�
    zb0(i) = find(data(i,:) == max(data(i,:)),1);
    % ǰ17������¶ȵ���100��C
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
    % �ȿ���������
    con_2=min(datadf(i,zb0(i)-4:zb0(i)-1))>=0;
    con_3=max(datadf(i,zb0(i):26))<=0;
    % �ٿ����μ�
    con_4=max(datadf_2(i,zb0(i)-4:zb0(i)-2))<0;
    con_5=max(datadf_2(i,zb0(i):25))<0;
    % ��������ж�
    if isempty(con_5)
        con_5=0;
    end
    if con_1 & con_2 & con_3 & con_4 & con_5
        indx(i)=5;
        error_data{i,1}=[];
    end
    % ��������͹Լ��
    if con_1 & con_2 & con_3 & (con_4==0 | con_5==0)
        indx(i)=4;
        error_data{i,1}=[];
    end
    % ��������������㵥������Լ��
    if  con_1 & (con_2==0) & (con_3==0)
        indx(i)=3;
        error_data{i,1}=[find(datadf(i,zb0(i):26)>0)+zb0(i)-1 find(datadf(i,zb0(i)-4:zb0(i)-1)<0)+zb0(i)-5];
    end
    % ֻ�Ҳ಻���㵥��Լ��
    if  con_1 & con_2 & (con_3==0)
        indx(i)=2;
        error_data{i,1}=find(datadf(i,zb0(i):26)>0)+zb0(i)-1;
    end
    % ֻ��಻���㵥��Լ��
    if  con_1 & (con_2==0) & con_3
        %ԭʼ�����������
        indx(i)=1;
        %�ҵ�ԭʼ�����в�����Լ���ĵ�
        error_data{i,1}=find(datadf(i,zb0(i)-4:zb0(i)-1)<0)+(zb0(i)-4);
    end
    error_data{i,2}=zb0(i);
end
end