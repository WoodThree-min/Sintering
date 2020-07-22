function [pre_clean_data_rev,pre_clean_data_rev_indx,indx,error_data,plat_like_data_indx,max_tem_300_indx]=pre_clean_work0(data)
%% ����˵��
% �������----------------------------------------------------------------
% pre_clean_data_rev����Ҫ�޲������ݣ�nrows * 27��nrows=sum(pre_clean_data_rev_indx)
% pre_clean_data_rev_indx�� ��Ҫ�޲������ݵ�������n��Ԫ�ص�����,n=size(data,1)
% indx:find_ideal_sinter�ű��������������ǣ�n��Ԫ�ص�����,n=size(data,1)
% error_data��find_ideal_sinter�ű�����������Ҫ�޲���λ�ã�cell���ͣ���һ����λ�ã��ڶ����Ƕ�Ӧ�¶���ߵ����
% plat_like_data_indx���������ƽ�������ݵ����� logical���ͣ�n��Ԫ�ص�����
% max_tem_300������¶���300�����µ����ݵ����� logical���ͣ�n��Ԫ�ص�����
% �������----------------------------------------------------------------
% data��Դ���ݣ� nrows * 27 nrows=size(data,1)
% zb0����Ӧ�¶���ߵ�������� ��n��Ԫ�ص�����,n=size(data,1)
% V��̨���ٶ�
%% �����ݽ�����ϣ�ȷ���������ݺʹ�������

pre_data=data(size(data,1),:);
[indx,error_data]=find_ideal_sinter_convex(pre_data);


% ���ش������ݺ�ʵ������
pre_clean=((indx<5)&indx~=0);
pre_clean_data=data(pre_clean,:);


%% �ų��������ƽ��������
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
% �Ƴ�pre_clean�����е�plat_like_data
pre_clean_data_rev_indx=double(pre_clean)-double(plat_like_data_indx)';
   if pre_clean_data_rev_indx <0
       pre_clean_data_rev_indx=0;
   end

pre_clean_data_rev_indx=logical(pre_clean_data_rev_indx);
pre_clean_data_rev=data(pre_clean_data_rev_indx,:);

%% �ų�����¶���300�����µ�
% �¶���ֵ
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