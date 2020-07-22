        function [beta_hat1,err2] = panel_ols(y_mat, pos, w, th_up)
[nry, ncy] = size(y_mat);

y_mat2 = y_mat.';
y_vec = y_mat2(:);
%z_mat = repmat([ones(length(pos),1), pos', pos'.^2], [nry,1]);
%W = diag(repelem(w.^(0:nry-1), ncy)*(1-w)/(1-w^nry));
z_mat2 = pos';
z_vec = z_mat2(:);
z_mat = [ones(length(z_vec),1),z_vec,z_vec.^2];
C = diag(repelem(1./sqrt(w.^(0:nry-1)), ncy)*(1-w^nry)/(1-w)); % 越靠近现在的数据权重越大

z_matc = C*z_mat; y_vecc = C*y_vec;
beta_hat = (z_matc.'*z_matc) \ (z_matc.'*y_vecc);

if beta_hat(3) > th_up 
    H = [0 0 1];
    CC = th_up;
    beta_hat1 = beta_hat + inv(z_matc.'*z_matc)* H.' * inv(H*inv(z_matc.'*z_matc)*H.')*(CC-H*beta_hat);
else
    beta_hat1 =beta_hat;
end

y_hat = z_mat * beta_hat1;

err2 = sum((y_vec - y_hat).^2)/length(y_vec);

end