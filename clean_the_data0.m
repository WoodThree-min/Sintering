function [cleaned_data]=clean_the_data0(data,error_data,pre_clean_data_rev_indx,pos,indx)
% ��ϴ˼·��������Ҫ�޲�λ�õ�ǰ��������ݣ�����Լ����С���ˣ��޲���Ҫ�޲���λ��
%% ����˵��
% �������----------------------------------------------------------------
% cleaned_data����ϴ������ݣ�nrows * 27��nrows=sum(pre_clean_data_rev_indx)
% [a1 b1 c1]��[a2 b2 c2] H-LSE������ò���
% �������----------------------------------------------------------------
% data��Դ���ݣ� nrows * 27 nrows=size(data,1)
% error_data��find_ideal_sinter�ű�����������Ҫ�޲���λ�ã�cell���ͣ���һ����λ�ã��ڶ����Ƕ�Ӧ�¶���ߵ����
% pre_clean_data_rev_indx�� ��Ҫ�޲������ݵ�������n��Ԫ�ص�����,n=size(data,1)
% zb0����Ӧ�¶���ߵ�������� ��n��Ԫ�ص�����,n=size(data,1)
% pos:�¶ȵ��Ӧ��ʵλ��
% indx:find_ideal_sinter�ű��������������ǣ�n��Ԫ�ص�����,n=size(data,1)

data=data(size(data,1),:);
zb0 = find(data == max(data),1);

cleaned_data=data;
w=0.75;th_up=0.1;
if pre_clean_data_rev_indx==1
    % ֻ��಻���㵥��Լ��
    if indx==1
        % ȷ���޲�����
        raw_indx=-length((error_data{1,1}))+zb0-4:zb0-1;
        if length(error_data{1,1})==1
            rev_indx=raw_indx;
            error_data_etmp=error_data{1,1};
            temp_con_indx=find(raw_indx==error_data_etmp);
            rev_indx([temp_con_indx-1 temp_con_indx])=[];
            pos1=pos(rev_indx);
            YY1 = data(rev_indx);
            [beta_hat1] = panel_ols(YY1, pos1, w, th_up);
            [a1, b1, c1] = deal(beta_hat1(3), beta_hat1(2), beta_hat1(1));
            p=[a1, b1, c1];
            % ��������ߵ�ֵ��Ϊ��error���ݵ��޲�
            pos_fix=pos(error_data_etmp);
            temp_fix_data=polyval(p,pos_fix);temp_fix_data2=polyval(p,(pos_fix-3));
            temp_fix_data_bias=temp_fix_data-cleaned_data(error_data_etmp);
            temp_fix_data_bias2=temp_fix_data2-cleaned_data(error_data_etmp-1);
            if temp_fix_data_bias2>temp_fix_data_bias
                if temp_fix_data2>0 & temp_fix_data2<550
                cleaned_data(error_data_etmp-1)=temp_fix_data2;
                end
            else
                if temp_fix_data>0 & temp_fix_data<550
                cleaned_data(error_data_etmp)=temp_fix_data;
                end
            end
        else
            for j=1:length(error_data{1,1})
                rev_indx=raw_indx;
                error_data_etmp=error_data{1,1};
                temp_con_indx=raw_indx==error_data_etmp(j);
                rev_indx(temp_con_indx)=[];
            end
            pos1=pos(rev_indx);
            YY1 = data(rev_indx);
            [beta_hat1] = panel_ols(YY1, pos1, w, th_up);
            [a1, b1, c1] = deal(beta_hat1(3), beta_hat1(2), beta_hat1(1));
            p=[a1, b1, c1];
            % ��������ߵ�ֵ��Ϊ��error���ݵ��޲�
            pos_fix=pos(error_data_etmp);
            if polyval(p,pos_fix)>0 & polyval(p,pos_fix)<550
            cleaned_data(error_data_etmp)=polyval(p,pos_fix);
            end
        end
    end
    % ֻ�Ҳ಻���㵥��Լ��
    if indx==2
        % ȷ���޲�����
        raw_indx=zb0+1:zb0+length(error_data{1,1});
        if (zb0+length(error_data{1,1}))>27
            raw_indx=zb0+1:zb0+length(error_data{1,1})-1;
        end
        if length(raw_indx)<3
            raw_indx=[zb0 raw_indx];
        end
        for j=1:length(error_data{1,1})
            rev_indx=raw_indx;
            error_data_etmp=error_data{1,1};
            temp_con_indx=raw_indx==error_data_etmp(j);
            rev_indx(temp_con_indx)=[];
        end
        if length(rev_indx)<3 & (zb0==26)
            rev_indx=[25 zb0 rev_indx];
        end
        pos1=pos(rev_indx);
        YY1 = data(rev_indx);
        [beta_hat1] = panel_ols(YY1, pos1, w, th_up);
        [a2, b2, c2] = deal(beta_hat1(3), beta_hat1(2), beta_hat1(1));
        b=[a2, b2, c2];
        % ��������ߵ�ֵ��Ϊ��error���ݵ��޲�
        pos_fix=pos(error_data_etmp);
        if polyval(b,pos_fix)>0 & polyval(b,pos_fix)<550
        cleaned_data(error_data_etmp)=polyval(b,pos_fix);
        end
    end
    % ˫�಻���㵥��Լ��
    if indx==3
        %�ȷֱ��ȡ���������
        error_data_indx_double=error_data{1,1};
        error_data_indx_p=error_data_indx_double(error_data_indx_double<zb0);
        error_data_indx_b=error_data_indx_double(error_data_indx_double>zb0);
        % ȷ��B���޲�����
        raw_indx=zb0+1:zb0+length(error_data_indx_b);
        if (zb0+length(error_data_indx_b))>27
            raw_indx=zb0+1:zb0+length(error_data_indx_b)-1;
        end
        if length(raw_indx)<3
            raw_indx=[zb0 raw_indx];
        end
        for j=1:length(error_data_indx_b)
            rev_indx=raw_indx;
            error_data_etmp=error_data_indx_b;
            temp_con_indx=(raw_indx==error_data_etmp(j));
            rev_indx(temp_con_indx)=[];
        end
        if length(rev_indx)<3
            rev_indx=[zb0 rev_indx];
        end
        pos1=pos(rev_indx);
        YY1 = data(rev_indx);
        [beta_hat1] = panel_ols(YY1, pos1, w, th_up);
        [a2, b2, c2] = deal(beta_hat1(3), beta_hat1(2), beta_hat1(1));
        p=[a2, b2, c2];
        % ��������ߵ�ֵ��Ϊ��error���ݵ��޲�
        pos_fix=pos(error_data_etmp);
        if polyval(p,pos_fix)>0 & polyval(p,pos_fix)<550
            cleaned_data(error_data_etmp)=polyval(p,pos_fix);
        end
        % ȷ��P���޲�����
        if length(error_data_indx_p)==1
            raw_indx=-length(error_data_indx_p)+zb0-4:zb0-1;
            rev_indx=raw_indx;
            error_data_etmp=error_data_indx_p;
            temp_con_indx=find(raw_indx==error_data_etmp);
            rev_indx([temp_con_indx-1 temp_con_indx])=[];
            pos_fix=pos(error_data_etmp);
            temp_fix_data=polyval(p,pos_fix);temp_fix_data2=polyval(p,(pos_fix-3));
            temp_fix_data_bias=temp_fix_data-cleaned_data(error_data_etmp);
            temp_fix_data_bias2=temp_fix_data2-cleaned_data(error_data_etmp-1);
            if temp_fix_data_bias2>temp_fix_data_bias
                if polyval(p,(pos_fix-3))>0 & polyval(p,(pos_fix-3))<550
                    cleaned_data(error_data_etmp-1)=temp_fix_data2;
                end
            else
                if polyval(p,pos_fix)>0 & polyval(p,pos_fix) <550
                    
                    cleaned_data(error_data_etmp-1)=temp_fix_data;
                end
            end
            % ��������ߵ�ֵ��Ϊ��error���ݵ��޲�
        else
            raw_indx=-length(error_data_indx_p)+zb0-3:zb0-1;
            for j=1:length(error_data_indx_p)
                rev_indx=raw_indx;
                error_data_etmp=error_data_indx_p;
                temp_con_indx=raw_indx==error_data_etmp(j);
                rev_indx(temp_con_indx)=[];
            end
            pos1=pos(rev_indx);
            YY1 = data(rev_indx);
            [beta_hat1] = panel_ols(YY1, pos1, w, th_up);
            [a1, b1, c1] = deal(beta_hat1(3), beta_hat1(2), beta_hat1(1));
            p=[a1, b1, c1];
            % ��������ߵ�ֵ��Ϊ��error���ݵ��޲�
            pos_fix=pos(error_data_etmp);
            if polyval(p,pos_fix)>0 & polyval(p,pos_fix) <550
                cleaned_data(error_data_etmp)=polyval(p,pos_fix);
            end
        end
    end
    if indx==4
        error_data_indx_double=error_data{1,1};
        error_data_indx_p=error_data_indx_double(error_data_indx_double<zb0);
        error_data_indx_b=error_data_indx_double(error_data_indx_double>zb0);
        % ȷ��B���޲�����
        if length(error_data_indx_b)>=1
            raw_indx=zb0+1:min(zb0+length(error_data_indx_b),27);
            if (zb0+length(error_data_indx_b))>27
                raw_indx=zb0+1:27;
            end
            if length(raw_indx)<3
                raw_indx=[zb0 raw_indx];
            end
            rev_indx=raw_indx;
            error_data_etmp=error_data_indx_b;
            for j=1:length(error_data_indx_b)
                
                temp_con_indx=(rev_indx==error_data_etmp(j));
                rev_indx(temp_con_indx)=[];
            end
            if length(rev_indx)<3
                rev_indx=[zb0 rev_indx];
            end
            pos1=pos(rev_indx);
            YY1 = data(rev_indx);
            [beta_hat1] = panel_ols(YY1, pos1, w, th_up);
            [a2, b2, c2] = deal(beta_hat1(3), beta_hat1(2), beta_hat1(1));
            p=[a2, b2, c2];
            % ��������ߵ�ֵ��Ϊ��error���ݵ��޲�
            pos_fix=pos(error_data_etmp);
            if polyval(p,pos_fix)>0 &  polyval(p,pos_fix)<550
                cleaned_data(error_data_etmp)=polyval(p,pos_fix);
            end
        end
        % ȷ��P���޲�����
        if length(error_data_indx_p)>=1 & length(error_data_indx_p)<3
            
            raw_indx=-length(error_data_indx_p)+zb0-4:zb0-1;
            rev_indx=raw_indx;
            error_data_etmp=error_data_indx_p;
            %         temp_con_indx=find(raw_indx==error_data_etmp);
            rev_indx=setdiff(raw_indx,error_data_etmp);
            %         rev_indx([temp_con_indx-1 temp_con_indx])=[];
            pos1=pos(rev_indx);
            YY1 = data(rev_indx);
            [beta_hat1] = panel_ols(YY1, pos1, w, th_up);
            [a1, b1, c1] = deal(beta_hat1(3), beta_hat1(2), beta_hat1(1));
            p=[a1, b1, c1];
            pos_fix=pos(error_data_etmp);
            temp_fix_data=polyval(p,pos_fix);
            if temp_fix_data>0 & temp_fix_data<550
                cleaned_data(error_data_etmp)=temp_fix_data;
            end
        end
    end
    
end
end
