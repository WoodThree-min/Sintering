clear;
% load ys.mat
% load ys.mat;YY=ys;
% load nys.mat;
% load YY.mat;
% load YY1.mat;
% load('yy_4_250_3550.mat', 'yy_4_250_3550')
% C =yy_4_250_3550;
% load('temp_328300_333220.mat')
% C=temp_328300_333220;
% load data_cleansed.mat;
% C = (A+B) / 2;
% load data_cleansed5.mat;
% C = (A5+B5) / 2;
% load('A_side_cleaned_data_70450-78400.mat')
% C=A2;
% load('data5_A_88850-95000.mat')
% C=YY;
load('data5_A_88850-95000.mat')
C1=YY;
load('data5_B_88850-95000.mat')
C2=YY;
warning off;
nrC= size(C1,1);C1=[zeros(nrC,1),C1];
C2=[zeros(nrC,1),C2];
alfs   =     [55,58,61,64,67,70,72,74]+2;
bets   =     [58,61,64,67,70,72,74,76]+2;
delts  =     bets-alfs;
pos=[1 3 5 4.5+3*(1:22) 73 75 77];
Z=50:0.1:78;
[zb_hat1,TB_hat1,zb01,R1,zb_hat2,TB_hat2,zb02,R2] = deal(zeros(1,nrC));
[BETA1,BETA2]=deal(zeros(nrC,4));
m=12;k=1;[zb1,zb2]=deal(zeros(1,nrC-m+1));
%%
tic;
for i = 1:nrC
    [~,zb01(i)]=max(C1(i,:));
    [zb_hat1(i), TB_hat1(i), R1(i),BETA1(i,:),~,~] = BTP_Search((alfs(zb01(i)-20)+bets(zb01(i)-20))/2,C1(i,:));
    a11=BETA1(i,1);b11=BETA1(i,2);a21=BETA1(i,3);b21=BETA1(i,4);
    [~,zb02(i)]=max(C2(i,:));
    [zb_hat2(i), TB_hat2(i), R2(i),BETA2(i,:),~,~] = BTP_Search((alfs(zb02(i)-20)+bets(zb02(i)-20))/2,C2(i,:));
    a12=BETA2(i,1);b12=BETA2(i,2);a22=BETA2(i,3);b22=BETA2(i,4);    
    if i>=m
        c11 = TB_hat1(i) + a11 * ( - zb_hat1(i)^2) + b11 *(-zb_hat1(i));   % Ignore a * (delta ^2 /3);
        c21 = TB_hat1(i) + a21 * ( - zb_hat1(i)^2) + b21*(-zb_hat1(i));
        y11 = a11 * Z.^2 + b11 * Z + c11;
        y21 = a21 * Z.^2 + b21 * Z + c21;
        zb1(i-11) = mean(zb_hat1(i-m+1:i));
        c12 = TB_hat2(i) + a12 * ( - zb_hat2(i)^2) + b12 *(-zb_hat2(i));   % Ignore a * (delta ^2 /3);
        c22 = TB_hat2(i) + a22 * ( - zb_hat2(i)^2) + b22*(-zb_hat2(i));
        y12 = a12 * Z.^2 + b12 * Z + c12;
        y22 = a22 * Z.^2 + b22 * Z + c22;
        zb2(i-11) = mean(zb_hat2(i-m+1:i));
        figure(30);
        set(gcf,'position',[20   45   1500   600 ]);
        subplot(2,3,[1,2]);
        ylim([50,600]); xlim([50,78]);
        str2 = ['SSE = ', num2str(R1(i)) ];
        str3 = ['Total No. of iterations: [ ',num2str(nrC),'  ] '];
        str4 = ['The iteration No.:       [ ',  num2str(i), ' ]'];
        plot(pos(19:28),C1(i,19:28),'.-','MarkerSize',20);
        title('北侧风箱曲线拟合及残差');
        text(51,570,str2,'FontSize',10);
        text(51,520,str3,'FontSize',10);
        text(51,470,str4,'FontSize',10);
        hold on;
        plot(zb_hat1(i),TB_hat1(i),'o','MarkerSize',10);
        plot(Z,y11,'r',Z,y21,'r','linewidth',1.1);
        line([70,70]+2,[50,600],'color','y','linewidth',1.5);
        line([72,72]+2,[50,600],'color','y','linewidth',1.5);
        line([74,74]+2,[50,600],'color','y','linewidth',1.5);
        line([67,67]+2,[50,600],'color','y','linewidth',1.5);
        line([zb_hat1(i),zb_hat1(i)],[50,TB_hat1(i)],'color','m','linewidth',1.5);
        line([50,zb_hat1(i)],[TB_hat1(i),TB_hat1(i)],'color','m');        
        ylim([50,600]); xlim([50,78]);
        grid minor;
        hold off;
        subplot(2,3,[4,5]);
        ylim([50,600]); xlim([50,78]);
        str2 = ['SSE = ', num2str(R1(i)) ];
        str3 = ['Total No. of iterations: [ ',num2str(nrC),'  ] '];
        str4 = ['The iteration No.:       [ ',  num2str(i), ' ]'];
        plot(pos(19:28),C2(i,19:28),'.-','MarkerSize',20);
        title('南侧风箱曲线拟合及残差');
        text(51,570,str2,'FontSize',10);
        text(51,520,str3,'FontSize',10);
        text(51,470,str4,'FontSize',10);
        hold on;
        plot(zb_hat2(i),TB_hat2(i),'o','MarkerSize',10);
        plot(Z,y12,'r',Z,y22,'r','linewidth',1.1);
        line([70,70]+2,[50,600],'color','y','linewidth',1.5);
        line([72,72]+2,[50,600],'color','y','linewidth',1.5);
        line([74,74]+2,[50,600],'color','y','linewidth',1.5);
        line([67,67]+2,[50,600],'color','y','linewidth',1.5);
        line([zb_hat2(i),zb_hat2(i)],[50,TB_hat2(i)],'color','m','linewidth',1.5);
        line([50,zb_hat2(i)],[TB_hat2(i),TB_hat2(i)],'color','m');           
        ylim([50,600]); xlim([50,78]);
        grid minor;
        hold off;
        subplot(1,3,3);
        if i<23
            p1=plot(zb1(1:i-11),'.-','MarkerSize',20);
            hold on;
            p2=plot(zb2(1:i-11),'.-','MarkerSize',20);
            plot(k,zb1(i-11),'ro','MarkerSize',10);
            plot(k,zb2(i-11),'ro','MarkerSize',10);
            k=k+1;
        else
            p1=plot(zb1(i-m-10:i-11),'.-','MarkerSize',20);
            hold on;
            plot(12,zb1(i-11),'ro','MarkerSize',10);
            p2=plot(zb2(i-m-10:i-11),'.-','MarkerSize',20);
            plot(12,zb2(i-11),'ro','MarkerSize',10);
        end
        title('BTP实时在线估计（平滑后）');
        ylim([69,76]); xlim([0,13]);
        grid minor;
       % line([0,13],[69,69],'color','y');
       line([0,13],[72,72],'color','y','linewidth',1.5);
        line([0,13],[74,74],'color','y','linewidth',1.5);
       % line([0,13],[76,76],'color','y');
       legend([p1,p2],'北侧BTP','南侧BTP');
        hold off;
    end
end
toc;
% % nzb_hat=zeros(1,nrC-m+1);
% % for i = 1:(nrC-m+1)
% %     nzb_hat(i) = mean(zb_hat(i:i+m-1));
% % end
% figure(31);subplot(2,1,1);plot(1:nrC-m+1,zb1,'.-',1:nrC-m+1,zb2,'.-');
% ylim([69,76]);
% %ylim([min(zb_hat),max(zb_hat)]);
% subplot(2,1,2);plot(zb01);
% figure(22);subplot(2,1,1);histogram(zb1,10);subplot(2,1,2);histogram(zb2,10);