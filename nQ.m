%=========================================================================
%|                                                                       |
%|           Q(x,par).m                                                  |
%|                                                                       |
%|    Target function to test the Newton-Raphson optimization method.    |
%|                                                                       |
%|    x    :   point in Euclidean R^n;                                   |
%|    par  :   parameters to be transmitted.                             |
%|                                                                       |
%=========================================================================

function [R,TB,BETA] = nQ(x,ys)

Z_B = x;

alfs   =     [46,49,52,55,58,61,64,67,70,72,74]'+2;
bets   =     [49,52,55,58,61,64,67,70,72,74,76]'+2;
delts  =     0.5 * ( bets - alfs );
z      =     alfs + delts;
Y      =    ys(18:28)';

ind     =  (1:11)';
i_max   =  max(ind(alfs <= Z_B) );

n1     =   i_max - 1;
n2     =   11 - i_max;

z1     =   z( 1 : n1 );             %  Seperating interpreating variables;
z2     =   z( (i_max+1) : 11);

alf    = alfs(i_max);            %  Interval of max flue gas temperature;
bet    = bets(i_max);
delt   = bet - alf;




%-----------------------   Design Matrix ------------------------------

X1  = [z1.^2 + ( (1/3) * delts(1:n1).^2 -Z_B^2 ),...
    z1 - Z_B * ones(n1,1)];
X2  = [z2.^2 + ( (1/3) * delts((i_max+1):11).^2 -Z_B^2 ) ,...
    z2 - Z_B * ones(n2,1)];

tmp  = (0.5 / delt);

x1   = tmp * (  (Z_B^3 - alf^3)/3 - Z_B^2 * (Z_B - alf ));
x2   = tmp * (  (Z_B^2 - alf^2)/2 - Z_B   * (Z_B - alf ));
x3   = tmp * (  (bet^3 - Z_B^3)/3 - Z_B^2 * (bet - Z_B ));
x4   = tmp * (  (bet^2 - Z_B^2)/2 - Z_B   * (bet - Z_B ));

X    = [X1,           zeros(n1,2);...
    x1,  x2,      x3,  x4;...
    zeros(n2,2),  X2            ];          % Design Matrix

%------------------------- W-LSE -------------------------------------
% Modified responses;
star=i_max-3;over=i_max+3;
Y=Y(max(1,star):min(over,11));
if n2>0 
nX = X;
nX=[nX,ones(size(nX,1),1)];nX=nX(max(1,star):min(over,11),:);
BET =  (nX' * nX) \ (nX' * Y);T_B1=BET(end);
a1=BET(1);b1=BET(2);a2=BET(3);b2=BET(4);
par1=[a1,b1,a2,b2,Z_B];
I_flag   =   Constrains_Test(par1);
    if I_flag == 0
        D = Y - nX*BET;
        Q1 = D.'*D;
        R=Q1;
        TB=T_B1;
        BETA=par1(1:4);
    else
        QQ=0;TT=0;BB=[0 0 0 0];
        %do 4:a1,a2,b2
        nX1  = (z1-Z_B).^2 +  (1/3) * delts(1:n1).^2;
        nx1   = tmp * (  (Z_B^3 - alf^3)/3 - Z_B^2 * (Z_B - alf ))...
            -2*Z_B*(tmp * (  (Z_B^2 - alf^2)/2 - Z_B   * (Z_B - alf )));
        nX    = [nX1,           zeros(n1,2);...
            nx1,           x3,  x4;...
            zeros(n2,1),  X2            ];
        nX=[nX,ones(size(nX,1),1)];nX=nX(max(1,star):min(over,11),:);
        BET =  (nX' * nX) \ (nX' * Y);T_B4=BET(end);
        a1=BET(1);a2=BET(2);b2=BET(3);b1=-2*a1*Z_B;
        
        par4=[a1,b1,a2,b2,Z_B];
        I_flag   =   Constrains_Test(par4);
        if I_flag == 0
            D = Y - nX*BET;
            Q4 = D.'*D;
            QQ=[QQ,Q4];
            TT=[TT,T_B4];
            BB=[BB;par4(1:4)];
        end
        
        %do 2:a1, b1, a2
        nX2 = (z2-Z_B).^2 +  (1/3) * delts((i_max+1):11).^2;
        nx3 = tmp * (  (bet^3 - Z_B^3)/3 - Z_B^2 * (bet - Z_B ))...
            - 2*Z_B*(tmp * (  (bet^2 - Z_B^2)/2 - Z_B   * (bet - Z_B )));
        nX    = [X1,           zeros(n1,1);...
            x1,  x2,      nx3 ;...
            zeros(n2,2),  nX2            ];
        nX=[nX,ones(size(nX,1),1)];nX=nX(max(1,star):min(over,11),:);
        BET =  (nX' * nX) \ (nX' * Y);T_B2=BET(end);
        a1=BET(1);b1=BET(2);a2=BET(3);b2=-2*a2*Z_B;
        par2=[a1,b1,a2,b2,Z_B];
        I_flag   =   Constrains_Test(par2);
        if I_flag == 0
            D = Y - nX*BET;
            Q2 = D.'*D;
            QQ=[QQ,Q2];
            TT=[TT,T_B2];
            BB=[BB;par2(1:4)];
        end
        
        %do 5:a1,a2
        nX1  = (z1-Z_B).^2 +  (1/3) * delts(1:n1).^2;
        nx1   = tmp * (  (Z_B^3 - alf^3)/3 - Z_B^2 * (Z_B - alf ))...
            -2*Z_B*(tmp * (  (Z_B^2 - alf^2)/2 - Z_B   * (Z_B - alf )));
        nX2 = (z2-Z_B).^2 +  (1/3) * delts((i_max+1):11).^2;
        nx3 = tmp * (  (bet^3 - Z_B^3)/3 - Z_B^2 * (bet - Z_B ))...
            - 2*Z_B*(tmp * (  (bet^2 - Z_B^2)/2 - Z_B   * (bet - Z_B )));
        nX    = [nX1,           zeros(n1,1);...
            nx1,           nx3;...
            zeros(n2,1),  nX2            ];
        nX=[nX,ones(size(nX,1),1)];nX=nX(max(1,star):min(over,11),:);
        BET =  (nX' * nX) \ (nX' * Y);T_B5=BET(end);
        a1=BET(1);a2=BET(2);b2=-2*a2*Z_B;b1=-2*a1*Z_B;
        
        par5=[a1,b1,a2,b2,Z_B];
        I_flag   =   Constrains_Test(par5);
        if I_flag == 0
            D = Y - nX*BET;
            Q5 = D.'*D;
            QQ=[QQ,Q5];
            TT=[TT,T_B5];
            BB=[BB;par5(1:4)];
        end
        
        %do 3:a1,b1,b2
        nX = X(:,[1,2,4]);
        nX=[nX,ones(size(nX,1),1)];nX=nX(max(1,star):min(over,11),:);
        BET =  (nX' * nX) \ (nX' * Y);T_B3=BET(end);
        a1=BET(1);b1=BET(2);b2=BET(3);a2=0;
        
        par3=[a1,b1,a2,b2,Z_B];
        I_flag   =   Constrains_Test(par3);
        if I_flag == 0
            D = Y - nX*BET;
            Q3 = D.'*D;
            QQ=[QQ,Q3];
            TT=[TT,T_B3];
            BB=[BB;par3(1:4)];
        end
        
        %do 7:b1,a2,b2
        nX = X(:, 2:4);
        nX=[nX,ones(size(nX,1),1)];nX=nX(max(1,star):min(over,11),:);
        BET =  (nX' * nX) \ (nX' * Y);T_B7=BET(end);
        b1=BET(1);a2=BET(2);b2=BET(3);a1=0;
        par7=[a1,b1,a2,b2,Z_B];
        I_flag   =   Constrains_Test(par7);
        if I_flag == 0
            D = Y - nX*BET;
            Q7 = D.'*D;
            QQ=[QQ,Q7];
            TT=[TT,T_B7];
            BB=[BB;par7(1:4)];
        end
        
        %do 8:b1,a2
        nX1 = z1 - Z_B * ones(n1,1);
        nX2 = (z2-Z_B).^2 +  (1/3) * delts((i_max+1):11).^2;
        nx3 = tmp * (  (bet^3 - Z_B^3)/3 - Z_B^2 * (bet - Z_B ))...
            - 2*Z_B*(tmp * (  (bet^2 - Z_B^2)/2 - Z_B   * (bet - Z_B )));
        nX    = [nX1,           zeros(n1,1);...
            x2,           nx3 ;...
            zeros(n2,1),  nX2            ];
        nX=[nX,ones(size(nX,1),1)];nX=nX(max(1,star):min(over,11),:);
        BET =  (nX' * nX) \ (nX' * Y);T_B8=BET(end);
        b1=BET(1);a2=BET(2);a1=0;b2=-2*a2*Z_B;
        
        par8=[a1,b1,a2,b2,Z_B];
        I_flag   =   Constrains_Test(par8);
        if I_flag == 0
            D = Y - nX*BET;
            Q8 = D.'*D;
            QQ=[QQ,Q8];
            TT=[TT,T_B8];
            BB=[BB;par8(1:4)];
        end
        
        %do 6:a1,b2
        nX1  = (z1-Z_B).^2 +  (1/3) * delts(1:n1).^2;
        nx1   = tmp * (  (Z_B^3 - alf^3)/3 - Z_B^2 * (Z_B - alf ))...
            -2*Z_B*(tmp * (  (Z_B^2 - alf^2)/2 - Z_B   * (Z_B - alf )));
        nX2 = z2 - Z_B * ones(n2,1);
        nX    = [nX1,           zeros(n1,1);...
            nx1,           x4;...
            zeros(n2,1),  nX2            ];
        nX=[nX,ones(size(nX,1),1)];nX=nX(max(1,star):min(over,11),:);
        BET =  (nX' * nX) \ (nX' * Y);T_B6=BET(end);
        a1=BET(1);b2=BET(2);a2=0;b1=-2*a1*Z_B;
        
        par6=[a1,b1,a2,b2,Z_B];
        I_flag   =   Constrains_Test(par6);
        if I_flag == 0
            D = Y - nX*BET;
            Q6 = D.'*D;
            QQ=[QQ,Q6];
            TT=[TT,T_B6];
            BB=[BB;par6(1:4)];
        end
        %do 9:b1,b2
        nX = X(:,[2,4]);
        nX=[nX,ones(size(nX,1),1)];nX=nX(max(1,star):min(over,11),:);
        BET =  (nX' * nX) \ (nX' * Y);T_B9=BET(end); 
        par9 =[0,BET(1),0,BET(2),Z_B];
        D = Y - nX*BET;
        Q9 = D.'*D;
        QQ=[QQ,Q9];
        TT=[TT,T_B9];
        BB=[BB;par9(1:4)];
        QQ=QQ(2:end);
        TT=TT(2:end);
        BB=BB(2:end,:);
        [R,k]=min(QQ);
        TB=TT(k);
        BETA=BB(k,:);
    end
else
    QQ=0;TT=0;BB=[0 0 0 0];
    %do 3:a1,b1,b2
    nX = X(:,[1,2,4]);
    nX=[nX,ones(size(nX,1),1)];nX=nX(max(1,star):min(over,11),:);
    BET =  (nX' * nX) \ (nX' * Y);T_B3=BET(end);
    a1=BET(1);b1=BET(2);b2=BET(3);a2=0;
    
    par3=[a1,b1,a2,b2,Z_B];
    I_flag   =   Constrains_Test(par3);
    if I_flag == 0
        D = Y - nX*BET;
        Q3 = D.'*D;
        QQ=[QQ,Q3];
        TT=[TT,T_B3];
        BB=[BB;par3(1:4)];
    end
    %do 6:a1,b2
    nX1  = (z1-Z_B).^2 +  (1/3) * delts(1:n1).^2;
    nx1   = tmp * (  (Z_B^3 - alf^3)/3 - Z_B^2 * (Z_B - alf ))...
        -2*Z_B*(tmp * (  (Z_B^2 - alf^2)/2 - Z_B   * (Z_B - alf )));
    nX2 = z2 - Z_B * ones(n2,1);
    nX    = [nX1,           zeros(n1,1);...
        nx1,           x4;...
        zeros(n2,1),  nX2            ];
    nX=[nX,ones(size(nX,1),1)];nX=nX(max(1,star):min(over,11),:);
    BET =  (nX' * nX) \ (nX' * Y);T_B6=BET(end);
    a1=BET(1);b2=BET(2);a2=0;b1=-2*a1*Z_B;
    
    par6=[a1,b1,a2,b2,Z_B];
    I_flag   =   Constrains_Test(par6);
    if I_flag == 0
        D = Y - nX*BET;
        Q6 = D.'*D;
        QQ=[QQ,Q6];
        TT=[TT,T_B6];
        BB=[BB;par6(1:4)];
    end
    %do 9:b1,b2
    nX = X(:,[2,4]);
    nX=[nX,ones(size(nX,1),1)];nX=nX(max(1,star):min(over,11),:);
    BET =  (nX' * nX) \ (nX' * Y);T_B9=BET(end);
    par9 =[0,BET(1),0,BET(2),Z_B];
    D = Y - nX*BET;
    Q9 = D.'*D;
    QQ=[QQ,Q9];
    TT=[TT,T_B9];
    BB=[BB;par9(1:4)];
    QQ=QQ(2:end);
    TT=TT(2:end);
    BB=BB(2:end,:);
    [R,k]=min(QQ);
    TB=TT(k);
    BETA=BB(k,:);
end
end




