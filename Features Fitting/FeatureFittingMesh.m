 
load ys.mat;                 % Single sample of flue gastemperatures;

   
alfs   =     [55,58,61,64,67,70,72,74]';
bets   =        [58,61,64,67,70,72,74,76]';
delts  =     0.5 * ( bets - alfs );
z      =     alfs + delts;  
Y      =    ys(21:28)'; 

   M   =   40;
   N   =   50;
   J   =   zeros(M,N) + 250000;; 
   Z_BS =  alfs(5) + (1:N)/(N+1) * (alfs(6)-alfs(5));
      Y_max = max(Y);
        dd = 0.02;
   T_BS =  Y_max * (1-dd +(1:M)/(M+1)*2*dd );
   
   
   
   for j = 1:N,
         Z_B  = Z_BS(j);
    for i = 1:M
         T_B  = T_BS(i);

[X,BET] = FeatureFitting(ys, Z_B, T_B);    % Performing W-LSE;
         
      Y_h   = X * BET + T_B;               % Response fitting;
      
       D  =  Y - Y_h;
    J(i,j) =  D' * D;                       % Squared residuals;
    
         if J(i,j)==min(min(J))
            BETm = BET; Xm = X; im = i; jm =j; Jm = J(i,j); Y_hm =Y_h;
            T_Bm = T_B; Z_Bm = Z_B;
         end
         
     end
   end
      
       BET = BETm; X = Xm; T_B = T_Bm; Z_B = Z_Bm;
   
       Y_h   = X * BET + T_B;               % Response fitting;
      
           a1 = BET(1); b1 = BET(2);  a2 = BET(3); b2 = BET(4);
       
    zf    = z(1) + (0:400)/400 * (z(8)-z(1));   % Back temperature curves;
     T_01    = a1 * (zf.^2   - Z_B^2)+...
              b1 * (zf - Z_B) + T_B; 
     T_02    = a2 * (zf.^2   - Z_B^2)+...
              b2 * (zf - Z_B) + T_B; 
    
figure (4); subplot(2,2,1); 
  bar((2:28),ys(2:28),'b'); xlim([1,29]); grid on;

          
          subplot(2,2,3);
        plot(z,Y,'.-','MarkerSize',40,'LineWidth',2); hold on;
        plot(z,Y_h,'o','MarkerSize',20,'LineWidth',3);
        line(zf,T_01,'Color','k','LineWidth',2);
        line(zf,T_02,'Color','k','LineWidth',2);
        for i=2:8,
             line([alfs(i),alfs(i)],[200,450],'LineWidth',3,'Color','r');
        end
        stem(Z_BS(jm),Y_h(5),'Color','m','LineWidth',2,'MarkerSize',20);
         z_jq = -.5 * b2 / a2;  
        line([z_jq,z_jq],[200,450],'Color','b','LineWidth',2); 
        grid on; grid minor; 
        xlim([55,76]);
        ylim([200,450]);
        hold off;
  
       subplot(2,2,2); 
          surfc(Z_BS,T_BS,J); shading flat; alpha(0.8);
              xlabel('V-BTP');
              ylabel('V-T.max');
            xlim([min(Z_BS),max(Z_BS)]);
            ylim([min(T_BS),max(T_BS)]);
       
       subplot(2,2,4); 
          contour (Z_BS,T_BS,J,150) ; hold on;
          plot(Z_B, T_B,'om','MarkerSize',10); hold off;
            xlabel('V-BTP');
              ylabel('V-T.max');
       
       