 
load ys.mat;                 % Single sample of flue gastemperatures;

figure (3); subplot(1,2,1);plot(ys(2:length(ys)),'.-',...
    'MarkerSize',35,'LineWidth',3);

alfs   =     [55,58,61,64,67,70,72,74]';
bets   =        [58,61,64,67,70,72,74,76]';
delts  =     0.5 * ( bets - alfs );
z      =     alfs + delts;  
Y      =    ys(21:28)';
 

Z_B    =    69.9412;% 68;
T_B    =  424.2394; % max(Y)  ;

[X,BET] = FeatureFitting(ys, Z_B, T_B);    % Performing W-LSE;

       a1 = BET(1); b1 = BET(2);      
       a2 = BET(3); b2 = BET(4);
       
    Y_h   = X * BET + T_B;               % Response fitting;
      
    zf    = z(1) + (0:400)/400 * (z(8)-z(1));   % Back temperature curves;
     T_01    = a1 * (zf.^2   - Z_B^2)+...
              b1 * (zf - Z_B) + T_B;
    
    T_02    = a2 * (zf.^2   - Z_B^2)+...
              b2 * (zf - Z_B) + T_B; 
    
figure (3);   subplot(1,2,2);
        plot(z,Y,'.-','MarkerSize',40,'LineWidth',2); hold on;
        plot(z,Y_h,'o','MarkerSize',20,'LineWidth',3);
        line(zf,T_01,'Color','r','LineWidth',2);
        line(zf,T_02,'Color','r','LineWidth',2);
         stem(Z_B,Y_h(5),'Color','m','LineWidth',2,'MarkerSize',20);
        grid on; grid minor; 
        xlim([55,76]);
        ylim([200,450]);
        hold off;
        
       