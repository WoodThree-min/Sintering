%==========================================================================
%|                                                                        |
%|   T_flag = Constrains_Test ( par )                                     |
%|                                                                        |
%|   The function test if the constrains are satisfied or not.            |
%|  The returned value  is I_flag:                                        |
%|                                                                        |
%|            I_flag = 0,       no constrains violated;                   |
%|            I_flag = 1,       at least one of the constrains violated.  |
%|                                                                        |
%|     The constrains are :                                               |
%|                                                                        |
%|        {  a1 <=0,  a2  <= 0,  ZB <= M1,  ZB >= M2. }                   |
%|                                                                        |
%|     The inputs are [par] £º                                            |
%|                                                                        |
%|        par = [a1, b1, a2, b2, ZB]';                                    |  
%|                                                                        |
%==========================================================================


  function I_flag   =   Constrains_Test(par)
   
       a1 = par(1); b1 = par(2); a2 = par(3); b2 = par(4); ZB = par(5);
       
       te  = 0.00001;           % tolerance error avoiding spikes and kinks;
      
          if abs(a1) > te 
                M1 = -0.5 * b1 / a1;
          else
                M1 = Inf;
          end
          if abs(a2) > te
                M2 = -0.5 * b2 / a2;
          else
                M2 = -Inf;
          end
          
       if     a1 > 0
           I_flag =  1; return;
       elseif a2 > 0 
           I_flag = 1; return;
       elseif a1 < 0 && ZB > M1 +  te
           I_flag = 1;  return 
        elseif a2 < 0 && ZB < M2 - te
           I_flag = 1;   return 
       end
   
       I_flag = 0;
   return
