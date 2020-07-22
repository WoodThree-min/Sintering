 
load ys.mat;                 % Single sample of flue gastemperatures;

figure (3); subplot(1,2,1);plot(ys(2:length(ys)),'.-');
 
alfs   =     [55,58,61,64,67,70,72,74]';
bets   =        [58,61,64,67,70,72,74,76]';
delts  =     0.5 * ( bets - alfs );
z      =     alfs + delts; 

Y      =    ys(21:28)';
 
ind     =  (1:8)';
i_max   =   min( ind( Y==max(Y) ) );    %  Pin-point the index of Y_B;
 
n1     =   i_max - 1;
n2     =   8 - i_max; 

z1     =   z( 1 : n1 );             %  Seperating interpreating variables;
z2     =   z( (i_max+1) : 8);

alf    = alfs(i_max);            %  Interval of max flue gas temperature;
bet    = bets(i_max);
delt   = bet - alf;

Z_B    =   z(i_max)  ;
T_B    =   Y(i_max)+5  ;

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
      BET = inv(X' * X) * (X' * Y_m);           % Parameter estimates;
       a1 = BET(1); b1 = BET(2);      
       a2 = BET(3); b2 = BET(4);
       
    Y_h   = X * BET + T_B;               % Response fitting;
    
    
    zf    = z(1) + (0:400)/400 * (z(8)-z(1));   % Back temperature curves;
     T_01    = a1 * (zf.^2 +  - Z_B^2)+...
              b1 * (zf - Z_B) + T_B;
    
    T_02    = a2 * (zf.^2 +  - Z_B^2)+...
              b2 * (zf - Z_B) + T_B;
    
    
    figure (3);   subplot(1,2,2);
        plot(z,Y,'.-','MarkerSize',30); hold on;
        plot(z,Y_h,'o','MarkerSize',10);
        line(zf,T_01,'Color','k');
        line(zf,T_02,'Color','k');
        grid on; grid minor; 
        xlim([55,76]);
        ylim([200,450]);
        hold off;
        
       