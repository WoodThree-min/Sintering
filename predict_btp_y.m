function [z1, z2, sigma1, sigma2] = predict_btp_y(zb, p)

L1 = 5*12;  % 5 minutes
L2 = 10*12; % 10 minutes

if p>L2
    z1 = mean(zb(p-L1+1:p));
    sigma1 = std(zb(p-L1+1:p));
    z2 = mean(zb(p-L2+1:p));
    sigma2 = std(zb(p-L2+1:p));
elseif p>L1
    z1 = mean(zb(p-L1+1:p));
    sigma1 = std(zb(p-L1+1:p));
    z2 = mean(zb(1:p));
    sigma2 = std(zb(1:p));
else
    z1 = mean(zb(1:p));
    sigma1 = std(zb(1:p));
    z2 = mean(zb(1:p));
    sigma2 = sigma1;
end

