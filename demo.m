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
load('data5_A_267000-271950.mat')
C=YY;
% load('matlab2.mat','cleaned_data2')
% C=cleaned_data2(1:10000,:);
warning off;
nrC = size(C,1);C=[zeros(nrC,1),C];
alfs   =     [55,58,61,64,67,70,72,74]+2;
bets   =     [58,61,64,67,70,72,74,76]+2;
delts  =     bets-alfs;
pos=[1 3 5 4.5+3*(1:22) 73 75 77];
Z=50:0.1:78;
[zb_hat,TB_hat,zb0,R] = deal(zeros(1,nrC));
BETA=zeros(nrC,4);k=1;
m=12;zb1=zeros(1,nrC-m+1);
%%
tic;
k1=1;k2=nrC;
for i = k1:k2
    [yb0,zb0(i)]=max(C(i,:));
    [zb_hat(i), TB_hat(i), R(i),BETA(i,:),~,~] = BTP_Search((alfs(zb0(i)-20)+bets(zb0(i)-20))/2,C(i,:));
    a1=BETA(i,1);b1=BETA(i,2);a2=BETA(i,3);b2=BETA(i,4);
    if i>=m
        c1 = TB_hat(i) + a1 * ( - zb_hat(i)^2) + b1 *(-zb_hat(i));   % Ignore a * (delta ^2 /3);
        c2 = TB_hat(i) + a2 * ( - zb_hat(i)^2) + b2 *(-zb_hat(i));
        y1 = a1 * Z.^2 + b1 * Z + c1;
        y2 = a2 * Z.^2 + b2 * Z + c2;
        zb1(i-11) = mean(zb_hat(i-m+1:i));
        figure(20);
        
        subplot(1,2,1);
        ylim([50,600]); xlim([50,78]);
        str2 = ['  SSE   =    ', num2str(R(i)) ];
        str3 = ['Total No. of iterations: [ ',num2str(nrC),'  ] '];
        str4 = ['The iteration No.:         [ ',  num2str(i), ' ]'];
        plot(pos(19:28),C(i,19:28),'.-','MarkerSize',20);
        title('实时曲线拟合及残差');
        text(51,570,str2,'FontSize',15);
        text(51,520,str3,'FontSize',15);
        text(51,470,str4,'FontSize',15);
        hold on;
        plot(zb_hat(i),TB_hat(i),'o','markersize',10);
        plot(pos(zb0(i)),yb0,'r.','markersize',15);
        plot(Z,y1,'r',Z,y2,'r');
        line([70,70]+2,[50,600],'color','y','linewidth',1.5);
        line([72,72]+2,[50,600],'color','y','linewidth',1.5);
        line([74,74]+2,[50,600],'color','y','linewidth',1.5);
        line([76,76]+2,[50,600],'color','y','linewidth',1.5);
        line([67,67]+2,[50,600],'color','y','linewidth',1.5);
        line([zb_hat(i),zb_hat(i)],[50,TB_hat(i)],'color','m','linewidth',1.5);
        line([50,zb_hat(i)],[TB_hat(i),TB_hat(i)],'color','m');
        ylim([50,600]); xlim([50,78]);
        grid minor;
        hold off;
        subplot(1,2,2);
        zb3=zb1;
        if i<23
            
            plot(zb3(1:i-11),'.-','markersize',20);
            hold on;
            plot(k,zb3(i-11),'ro','MarkerSize',10);
            k=k+1;
        else
            plot(zb3(i-m-10:i-11),'.-','markersize',20);
            hold on;
            plot(12,zb3(i-11),'ro','MarkerSize',10);
        end
        title('BTP实时在线估计（平滑后）');
        ylim([69,78]); xlim([0,13]);
        grid minor;
        
       % line([0,13],[69,69],'color','y');
        line([0,13],[72,72],'color','y','linewidth',1.5);
        line([0,13],[74,74],'color','y','linewidth',1.5);
        line([0,13],[76,76],'color','y','linewidth',1.5);
        hold off;
        pause(0.5);
    end
    
end
toc;
% % nzb_hat=zeros(1,nrC-m+1);
% % for i = 1:(nrC-m+1)
% %     nzb_hat(i) = mean(zb_hat(i:i+m-1));
% % end
% figure(21);subplot(2,1,1);plot(zb1,'.');
% ylim([69,76]);
% %ylim([min(zb_hat),max(zb_hat)]);
% subplot(2,1,2);plot(zb0());
% figure(22);histogram(zb1,10);