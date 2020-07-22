%=========================================================================
%|                                                                       |
%|           [ZB, TB, BET, Yh, R] = BTP_Search_Z(ZB0, ys)                |
%|                                                                       |
%|  This program performs the crude min search using one-dim points in   |
%| the neighbour of [ ZB0 ] w.r.t the target function Q = Q_form_Z.      |
%|                                                                       |
%|    ZB0      :   initial point in the one-dim interval;                |
%|    ys       :   temperature records at boxes: [21]---[28];            |
%|    Q(x,par) :   SSE of the candidate [ZB].                            |
%|                                                                       |
%|  The Outputs are                                                      |
%|     [ZB, TB] :  Optimized virtual BTP [ZB], and the temperature at it.|                                                                     |
%|        BET   :  [a1,b1,a2,b2]' LSE estimates of the parameters.       |
%|        Yh    :       Fitting values of [21]-[28] boxes' temperature;  |
%|        R     :       Sum of squares of the residuals.                 |
%=========================================================================

function [ZB, TB, R,BETA,zbs,tbs] = BTP_Search(ZB0, ys)

%% Parameter setting-------------------------------------------
ns   = [50,20];
di  =3;
az = max(ZB0 - di,57.1);
bz =  min(ZB0 + di,77.9);  
zm = 0.5 * (az + bz);
dz  =  (bz - az);
for ii = 1:2
    n = ns(ii);
    zbs = zm +  dz * ((1:n)/(n+1) - 0.5);
    q  = zeros(n,1);
    tbs = q;
    ind = (1:n);    
    for k = 1:n
        ZB = zbs(k);
        [R, TB,~] = nQ(ZB, ys);   % Evaluation Q at net points;
        q(k)  = R;
        tbs(k)  = TB;
    end
    
    j       =     ind(q == min(q));
    ZB      =     zbs(j);
    
    if ii == 1
        zm  =  mean(ZB);
        dz  = dz/n;
    end
end
ZB=mean(ZB);
[R, TB,BETA] = nQ(ZB, ys);
return
