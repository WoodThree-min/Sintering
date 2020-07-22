%=========================================================================
%|         FeatureFitting.m                                              |
%|                                                                       |
%|   [X,BET] = FeatureFiting(ys, Z_B, T_B)                               |
%|                                                                       |
%|   This function returns the design matrix and the W-LSE estimates of  |
%| the parobalas charatering the stable iron ore sintering process.      |
%|                                                                       |
%|   The input consists of :                                             |
%|        ys:   the last 8 flue gas temperatre records;                  |
%|        Z_B:  assumed virtual BTP;                                     |
%|        T_B:  assumed virtual background maximum flue gas temperatre.  |
%|                                                                       |
%|                       School of Math. Shandong University             |
%|                                2019.08.18.                            |
%|                                                                       |
%=========================================================================

function  [X,BET] = FeatureFitting(ys, Z_B, T_B)

alfs   =     [55,58,61,64,67,70,72,74]';
bets   =        [58,61,64,67,70,72,74,76]';
delts  =     0.5 * ( bets - alfs );
z      =     alfs + delts; 

Y      =    ys(21:28)';
 
ind     =  (1:8)';
i_max   =  max(ind(alfs <= Z_B) );


n1     =   i_max - 1;
n2     =   8 - i_max; 

z1     =   z( 1 : n1 );             %  Seperating interpreating variables;
z2     =   z( (i_max+1) : 8);

alf    = alfs(i_max);            %  Interval of max flue gas temperature;
bet    = bets(i_max);
delt   = bet - alf;



%-----------------------   Design Matrix ------------------------------

  X1  = [z1.^2 + ( (1/3) * delts(1:n1).^2 -Z_B^2 ) .* ones(n1,1),...
           z1 - Z_B * ones(n1,1)];
  X2  = [z2.^2 + ( (1/3) * delts((i_max+1):8).^2 -Z_B^2 ) .* ones(n2,1),...
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

       Y_m = Y - T_B;                          % Modified responses;
      BET =  (X' * X) \ (X' * Y_m);           % Parameter estimates;