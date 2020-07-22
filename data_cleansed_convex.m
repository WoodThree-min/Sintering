function data = data_cleansed_convex(data,pos)
w = 0.75; th_up = 0.1;
data = data(size(data,1), :);  % 取最后一个
[~, zb0] = max(data);

if zb0>=21
    % P层修复
    pos1 = pos([zb0-4,zb0-3,zb0-2,zb0]);
    YY1 = data([zb0-4,zb0-3,zb0-2,zb0]);
    [beta_hat1] = panel_ols(YY1, pos1, w, th_up);
    [a1, b1, c1] = deal(beta_hat1(3), beta_hat1(2), beta_hat1(1));
    p = [a1, b1, c1];
    pos_fix = pos(zb0-1);
    data(zb0-1) = polyval(p, pos_fix);


    % B层修复
    pos2 = pos(zb0:27);
    YY2 = data(zb0:27);
    [beta_hat1] = panel_ols(YY2, pos2, w, th_up);
    [a2, b2, c2] = deal(beta_hat1(3), beta_hat1(2), beta_hat1(1));
    b = [a2, b2, c2];
    if zb0<27
        pos_fix = pos(zb0+1);
        data(zb0+1) = polyval(b, pos_fix);
    end
end

