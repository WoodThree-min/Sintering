function [cleaned_data]=clean_the_data0(data,error_data,pre_clean_data_rev_indx,pos,indx)
% 清洗思路：依据需要修补位置的前后风箱数据，进行约束最小二乘，修补需要修补的位置
%% 变量说明
% 输出变量----------------------------------------------------------------
% cleaned_data：清洗后的数据，nrows * 27，nrows=sum(pre_clean_data_rev_indx)
% [a1 b1 c1]、[a2 b2 c2] H-LSE拟合所得参数
% 输入变量----------------------------------------------------------------
% data：源数据， nrows * 27 nrows=size(data,1)
% error_data：find_ideal_sinter脚本所得数据需要修补的位置，cell类型，第一列是位置，第二列是对应温度最高点风箱
% pre_clean_data_rev_indx： 需要修补的数据的索引，n个元素的向量,n=size(data,1)
% zb0：对应温度最高点风箱数据 ，n个元素的向量,n=size(data,1)
% pos:温度点对应现实位置
% indx:find_ideal_sinter脚本所得数据情况标记，n个元素的向量,n=size(data,1)

data=data(size(data,1),:);
zb0 = find(data == max(data),1);

cleaned_data=data;
w=0.75;th_up=0.1;
if pre_clean_data_rev_indx==1
    % 只左侧不满足单增约束
    if indx==1
        % 确定修补索引
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
            % 以拟合曲线的值作为对error数据的修补
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
            % 以拟合曲线的值作为对error数据的修补
            pos_fix=pos(error_data_etmp);
            if polyval(p,pos_fix)>0 & polyval(p,pos_fix)<550
            cleaned_data(error_data_etmp)=polyval(p,pos_fix);
            end
        end
    end
    % 只右侧不满足单减约束
    if indx==2
        % 确定修补索引
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
        % 以拟合曲线的值作为对error数据的修补
        pos_fix=pos(error_data_etmp);
        if polyval(b,pos_fix)>0 & polyval(b,pos_fix)<550
        cleaned_data(error_data_etmp)=polyval(b,pos_fix);
        end
    end
    % 双侧不满足单减约束
    if indx==3
        %先分别截取两侧的索引
        error_data_indx_double=error_data{1,1};
        error_data_indx_p=error_data_indx_double(error_data_indx_double<zb0);
        error_data_indx_b=error_data_indx_double(error_data_indx_double>zb0);
        % 确定B层修补索引
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
        % 以拟合曲线的值作为对error数据的修补
        pos_fix=pos(error_data_etmp);
        if polyval(p,pos_fix)>0 & polyval(p,pos_fix)<550
            cleaned_data(error_data_etmp)=polyval(p,pos_fix);
        end
        % 确定P层修补索引
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
            % 以拟合曲线的值作为对error数据的修补
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
            % 以拟合曲线的值作为对error数据的修补
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
        % 确定B层修补索引
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
            % 以拟合曲线的值作为对error数据的修补
            pos_fix=pos(error_data_etmp);
            if polyval(p,pos_fix)>0 &  polyval(p,pos_fix)<550
                cleaned_data(error_data_etmp)=polyval(p,pos_fix);
            end
        end
        % 确定P层修补索引
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
