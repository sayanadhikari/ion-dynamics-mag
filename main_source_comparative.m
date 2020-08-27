clc; clear;close all
global gamma_x gamma_z beta
alpha=89;
ep0=8.85E-12; AMU=1.67E-27;
mi=2*AMU;n0=1E16;e0=1.6E-19;
Te=1.5*e0;
cs=sqrt(Te/mi);
LD=sqrt(ep0*Te/(n0*e0^2));

B0=4;
Z =1E5;
gamma_x = sqrt(ep0/(n0*mi))*B0*cosd(alpha);
gamma_z = sqrt(ep0/(n0*mi))*B0*sind(alpha);
beta = sqrt((ep0*mi)/(n0*e0*e0))*Z;
%zspan=linspace(0,20,1000);
options=odeset('RelTol',1e-5);

[z1,y1]=ode45('diff_fun_source',[0 20],[0 0.01 1.0 0.01 0.01 1],options);
[z2,y2]=ode45('diff_fun_no_source',[0 40],[0 0.01 1 0.01 0.01 1],options);

%%%%%%%%%%%%% Exponetial Source%%%%%%%%%%%%%%%%%
N_e1 = exp(y1(:,1));
index1 = find(N_e1<=0.001); 
max_z1 = index1(1);
u1=y1((1:max_z1),4);
v1=y1((1:max_z1),5);
w1=y1((1:max_z1),6);
E1 = zeros(1,length(u1));
for i = 1:length(u1)
    E1(i)=sqrt((u1(i)^2+v1(i)^2+w1(i)^2)+(2*y1(i,1))-1);
end
% v_perp1 = sqrt(((u1.*cos(alpha)).^2)+((w1.*sin(alpha)).^2));
% v_para1 = sqrt(((u1.*sin(alpha)).^2)+((w1.*cos(alpha)).^2)+v1.^2);
% pitch1 = atand(v_perp1./v_para1);
  
% v_para1 = (u1.*cos(alpha))+(w1.*sin(alpha));
% v_perp1 = sqrt((((w1.*cos(alpha))-(u1.*sin(alpha))).^2)+v1.^2);
% pitch1 = atand(v_perp1./v_para1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v_para1 = sqrt((u1.^2)+(w1.^2));
    v_perp1 = sqrt(v1.^2);
    
    vel1=sqrt((u1.^2)+(v1.^2)+(w1.^2));
    r1=(mi.*v_perp1)/(e0.*B0);
    %T=(2*pi*r)./vel;
    T1=(2*pi*mi)./(e0*B0);
    pitch_length1=(T1*cs/LD).*v_para1;
    pitch1 = atand(v_perp1./v_para1);



%%%%%%%%%%%%%%% Heaviside Source %%%%%%%%%%%%%%%
N_e2 = exp(y2(:,1));
index2 = find(N_e2<=0.001); 
max_z2 = index2(1);
u2=y2((1:max_z2),4);
v2=y2((1:max_z2),5);
w2=y2((1:max_z2),6);
E2 = zeros(1,length(u2));
for i = 1:length(u2)
    E2(i)=sqrt((u2(i)^2+v2(i)^2+w2(i)^2)+(2*y2(i,1))-1);
end
% v_para2 = (u2.*cos(alpha))+(w2.*sin(alpha));
% v_perp2 = sqrt((((w2.*cos(alpha))-(u2.*sin(alpha))).^2)+v2.^2);
% pitch2 = atand(v_perp2./v_para2);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    v_para2 = sqrt((u2.^2)+(w2.^2));
    v_perp2 = sqrt(v2.^2);
    
    vel2=sqrt((u2.^2)+(v2.^2)+(w2.^2));
    r2=(mi.*v_perp2)/(e0.*B0);
    %T=(2*pi*r)./vel;
    T2=(2*pi*mi)./(e0*B0);
    pitch_length2=(T2*cs/LD).*v_para2;
    pitch2 = atand(v_perp2./v_para2);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% hold all
% figure(1);
%     subplot(431);
%     
%                     plot(z1(1:max_z1),y1((1:max_z1),1),'linewidth',2); hold on
%                     plot(z2(1:max_z2),y2((1:max_z2),1),'linewidth',2); hold on
%                     xlabel('z/\lambda_D'), ylabel('\eta'),grid on 
%                     legend('Exp.','Heaviside')
%                 
%     subplot(432); plot(z1(1:max_z1),y1((1:max_z1),2),'linewidth',2); hold on
%                   plot(z2(1:max_z2),y2((1:max_z2),2),'linewidth',2); hold on
%                 xlabel('z/\lambda_D'), ylabel('E'),grid on 
%                 legend('Exp.','Heaviside')
%                 
%                 title(['Plot for Mag B = ',num2str(B0),'T and Ang \alpha = ',num2str(alpha)])
% %                 title('Plot for Mag B = 5 T and Ang \alpha = 5');
%                 
%                 
%     subplot(433); plot(z1(1:max_z1),y1((1:max_z1),3),'linewidth',2); hold on
%                   plot(z2(1:max_z2),y2((1:max_z2),3),'linewidth',2); hold on
%                 xlabel('z/\lambda_D'), ylabel('N_i'),grid on
%                 plot(z1(1:max_z1),N_e1(1:max_z1),'linewidth',2);
%                 plot(z2(1:max_z2),N_e2(1:max_z2),'linewidth',2);
%                 xlabel('z/\lambda_D'), ylabel('N_e'),grid on  
%                 legend('N_i','N_i_{HS}','N_e','N_e_{HS}')
%     subplot(4,3,[4,5,6]); plot(z1(1:max_z1),y1((1:max_z1),4),'r','linewidth',2), grid on, hold on
%               plot(z1(1:max_z1),y1((1:max_z1),5),'g','linewidth',2), grid on, hold on
%               plot(z1(1:max_z1),y1((1:max_z1),6),'b','linewidth',2), grid on, hold on
%               legend('u','v','w')
%               title('Velocity Components for exponential ion source');
%               xlabel('z/\lambda_D'), ylabel('u,v,w')
%     subplot(4,3,[7,8,9]); plot(z2(1:max_z2),y2((1:max_z2),4),'r','linewidth',2), grid on, hold on
%               plot(z2(1:max_z2),y2((1:max_z2),5),'g','linewidth',2), grid on, hold on
%               plot(z2(1:max_z2),y2((1:max_z2),6),'b','linewidth',2), grid on, hold on
%               legend('u','v','w')
%               title('Velocity Components for Heaviside ion source');
%               xlabel('z/\lambda_D'), ylabel('u,v,w')        
%    sigma1 = y1(:,3)-exp(y1(:,1));
%    sigma2 = y2(:,3)-exp(y2(:,1));
%    subplot(4,3,[10,11,12]);plot(z1(1:max_z1),sigma1(1:max_z1),'linewidth',2);grid on, hold on
%                             plot(z2(1:max_z2),sigma2(1:max_z2),'linewidth',2);grid on, hold on
%                 xlabel('z/\lambda_D'), ylabel('Space charge (\sigma)')
%                 legend('Exp. ion source','Heaviside ion source')
%    
%   hold all
%    figure(2);
%    plot3(y1((1:max_z1),4),y1((1:max_z1),5),y1((1:max_z1),6),'o','linewidth',2);grid on;
%                 xlabel('u'), ylabel('v'),zlabel('w'),grid on, hold on
%    figure(2)
%    plot(z1(1:max_z1),E1,'linewidth',2), grid on, hold on
%     plot(z2(1:max_z2),E2,'linewidth',2), grid on, hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold all
figure(1);
    subplot(211);
    
                    plot(z1(1:max_z1),y1((1:max_z1),1),'linewidth',2); hold on
                    plot(z2(1:max_z2),y2((1:max_z2),1),'linewidth',2); hold on
                    xlabel('z/\lambda_D'), ylabel('\eta'),grid on 
                    legend('Exp.','Heaviside')
                     title(['Plot for Mag B = ',num2str(B0),'T and Ang \alpha = ',num2str(alpha),'\circ'])
                
    subplot(212); plot(z1(1:max_z1),y1((1:max_z1),2),'linewidth',2); hold on
                  plot(z2(1:max_z2),y2((1:max_z2),2),'linewidth',2); hold on
                xlabel('z/\lambda_D'), ylabel('E'),grid on 
                legend('Exp.','Heaviside')
                
               
%                 title('Plot for Mag B = 5 T and Ang \alpha = 5');
    hold all            
    figure(2);            
    plot(z1(1:max_z1),y1((1:max_z1),3),'linewidth',2); hold on
                  plot(z2(1:max_z2),y2((1:max_z2),3),'linewidth',2); hold on
                xlabel('z/\lambda_D'), ylabel('N_i'),grid on
                plot(z1(1:max_z1),N_e1(1:max_z1),'linewidth',2);
                plot(z2(1:max_z2),N_e2(1:max_z2),'linewidth',2);
                xlabel('z/\lambda_D'), ylabel('Total Density'),grid on  
                legend('N_i','N_i_{HS}','N_e','N_e_{HS}')
                title(['Plot for Mag B = ',num2str(B0),'T and Ang \alpha = ',num2str(alpha),'\circ'])

                
                
     hold all
   figure(3);
    subplot(2,1,1); plot(z1(1:max_z1),y1((1:max_z1),4),'r','linewidth',2), grid on, hold on
              plot(z1(1:max_z1),y1((1:max_z1),5),'g','linewidth',2), grid on, hold on
              plot(z1(1:max_z1),y1((1:max_z1),6),'b','linewidth',2), grid on, hold on
              legend('u','v','w')
              title('Velocities for exp. ion source');
              xlabel('z/\lambda_D'), ylabel('u,v,w')
    subplot(2,1,2); plot(z2(1:max_z2),y2((1:max_z2),4),'r','linewidth',2), grid on, hold on
              plot(z2(1:max_z2),y2((1:max_z2),5),'g','linewidth',2), grid on, hold on
              plot(z2(1:max_z2),y2((1:max_z2),6),'b','linewidth',2), grid on, hold on
              legend('u','v','w')
              title('Velocities for Heaviside ion source');
              xlabel('z/\lambda_D'), ylabel('u,v,w')        
%    sigma1 = y1(:,3)-exp(y1(:,1));
%    sigma2 = y2(:,3)-exp(y2(:,1));
%    subplot(4,3,[10,11,12]);plot(z1(1:max_z1),sigma1(1:max_z1),'linewidth',2);grid on, hold on
%                             plot(z2(1:max_z2),sigma2(1:max_z2),'linewidth',2);grid on, hold on
%                 xlabel('z/\lambda_D'), ylabel('Space charge (\sigma)')
%                 legend('Exp. ion source','Heaviside ion source')

   hold all
   figure(4);
       subplot(211); plot3(y1((1:max_z1),4),y1((1:max_z1),5),y1((1:max_z1),6),'o','linewidth',2);grid on;
                xlabel('u'), ylabel('v'),zlabel('w'),grid on, hold on
                plot3(y2((1:max_z2),4),y2((1:max_z2),5),y2((1:max_z2),6),'o','linewidth',2);grid on;
                xlabel('u'), ylabel('v'),zlabel('w'),grid on, hold on
                legend('Exp.','Heaviside','Location','northeast')
                 title('Ion flow velocity')
   
   subplot(2,1,2);plot(z1(1:max_z1),pitch1,'linewidth',2), grid on,  hold on
                    plot(z2(1:max_z2),pitch2,'linewidth',2), grid on,  hold on
                xlabel('z/\lambda_D'), ylabel('Pitch'),grid on, hold on
                legend('Exp.','Heaviside','Location','northeast')
 
%    subplot(2,1,2);plot(z1(1:max_z1),E1,'linewidth',2),  hold on
%                 plot(z2(1:max_z2),E2,'linewidth',2),  hold on
%    xlabel('z/\lambda_D'), ylabel('Energy'),grid on
%    legend('Exp.','Heaviside','Location','northeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold all            
    figure(5);            
    plot(z1(1:max_z1),y1((1:max_z1),3),'-- o','linewidth',2); hold on
                  plot(z2(1:max_z2),y2((1:max_z2),3),'linewidth',2); hold on
                xlabel('z/\lambda_D'), ylabel('Density'),grid on
                plot(z1(1:max_z1),N_e1(1:max_z1),'-- o','linewidth',2);
                plot(z2(1:max_z2),N_e2(1:max_z2),'linewidth',2);
                xlabel('z/\lambda_D'), ylabel('Total Density'),grid on  
                legend('N_i','N_i_{ZS}','N_e','N_e_{ZS}')
                title(['Plot for Mag B = ',num2str(B0),'T and Ang \alpha = ',num2str(alpha),'\circ'])
   
