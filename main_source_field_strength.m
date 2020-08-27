clc; clear;
global gamma_x gamma_z beta

B0=[1 2 3 4];
for k = 1:length(B0)
    
    ep0=8.85E-12; AMU=1.67E-27;
    mi=2*AMU;n0=1E16;e0=1.6E-19;
    alpha = 45;
    Z =1E5;
    Te=1.5*e0;
    cs=sqrt(Te/mi);
    LD=sqrt(ep0*Te/(n0*e0^2));
    gamma_x = sqrt(ep0/(n0*mi))*B0(k)*cosd(alpha);
    gamma_z = sqrt(ep0/(n0*mi))*B0(k)*sind(alpha);
    beta = sqrt((ep0*mi)/(n0*e0*e0))*Z;
    %zspan=linspace(0,20,1000);
    options=odeset('RelTol',1e-5);

    [z,y]=ode45('diff_fun_source',[0 20],[0 0.01 1 0.0 0.0 1],options);      
    % [z,y]=ode15s('diff_fun_source_heaviside',[0 20],[0 0.1 1 0.01 0.01 1],options);

    N_e = exp(y(:,1));
    index = find(N_e<=0.001); 
    max_z = index(1);
    u=y((1:max_z),4);
    v=y((1:max_z),5);
    w=y((1:max_z),6);
%     E = zeros(1,length(u));
%     for i = 1:length(u)
%         E(i)=(u(i)^2+v(i)^2+w(i)^2)+(2*y(i,1))-1;
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Energy   %%%%%%%%%%%%
    E=0.5*(u.^2+v.^2+w.^2);
    
    
    dis=z(1:max_z);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%% Lorentz force calculation %%%
   
   Fx=((LD*B0(k)*e0*cs)/(Te*sqrt(2))).*u;
   Fy=((LD*B0(k)*e0*cs)/(Te*sqrt(2))).*(w-u);
   Fz=-((LD*B0(k)*e0*cs)/(Te*sqrt(2)))*v;

   E1=-y((1:max_z),2);
   FE=E1;
   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     v_para = sqrt(((u.*cos(alpha)).^2)+((w.*sin(alpha)).^2));
%     v_perp = sqrt(((u.*sin(alpha)).^2)+((w.*cos(alpha)).^2)+v.^2);
    v_para = sqrt((u.^2)+(w.^2));
    v_perp = sqrt(v.^2);
    
    vel=sqrt((u.^2)+(v.^2)+(w.^2));
    r=(mi.*v_perp)/(e0.*B0(k));
    %T=(2*pi*r)./vel;
    T=(2*pi*mi)./(e0*B0(k));
    pitch_length=(T*cs/LD).*v_para;
    pitch = atand(v_perp./v_para);
    % hold all
    figure(1);
    plot(z(1:max_z),y((1:max_z),1),'linewidth',2); hold on
                    xlabel('z/\lambda_D'), ylabel('\eta'),grid on
                    legend(['B = ',num2str(B0(1)),'T'],['B = ',num2str(B0(2)),'T'],['B = ',num2str(B0(3)),'T'],['B = ',num2str(B0(4)),'T'])

    figure(2);
    plot(z(1:max_z),y((1:max_z),2),'linewidth',2); hold on
                    xlabel('z/\lambda_D'), ylabel('E'),grid on 
                    legend(['B = ',num2str(B0(1)),'T'],['B = ',num2str(B0(2)),'T'],['B = ',num2str(B0(3)),'T'],['B = ',num2str(B0(4)),'T'])
    figure(3);
    subplot(2,2,k);
    plot(z(1:max_z),y((1:max_z),3),'linewidth',2); hold on
                    xlabel('z/\lambda_D'), ylabel('N_i'),grid on
                    plot(z(1:max_z),N_e(1:max_z),'linewidth',2);
                    xlabel('z/\lambda_D'), ylabel('Density'),grid on  
                    legend('N_i','N_e','Orientation','horizontal')
%                     title(['B = ',num2str(B0(k)),'T,\alpha = ',num2str(alpha)])
                    title(['B = ',num2str(B0(k)),'T'])

    figure(4);
    subplot(4,1,k);
    plot(z(1:max_z),y((1:max_z),4),'r','linewidth',1.5), grid on, hold on
    plot(z(1:max_z),y((1:max_z),5),'g','linewidth',1.5), grid on, hold on
    plot(z(1:max_z),y((1:max_z),6),'b','linewidth',1.5), grid on, hold on
                  legend('u','v','w','Orientation','horizontal','Location','northwest')
                  xlabel('z/\lambda_D'), ylabel('u,v,w')
                  ylim([-1 4]);
                  title(['B = ',num2str(B0(k)),'T'])

    sigma = y(:,3)-exp(y(:,1));
    figure(5);
    plot(z(1:max_z),sigma(1:max_z),'linewidth',2);grid on, hold on
                    xlabel('z/\lambda_D'), ylabel('Space charge (\sigma)')
                    legend(['B = ',num2str(B0(1)),'T'],['B = ',num2str(B0(2)),'T'],['B = ',num2str(B0(3)),'T'],['B = ',num2str(B0(4)),'T'])

    %   hold all
    figure(6);
    subplot(2,2,k);
    plot3(y((1:max_z),4),y((1:max_z),5),y((1:max_z),6),'o','linewidth',2);grid on;
                    xlabel('u'), ylabel('v'),zlabel('w'),grid on, hold on
                    title(['B = ',num2str(B0(k)),'T'])

    %    hold all
    figure(7);
    subplot(4,1,k);
       plot(z(1:max_z),pitch,'linewidth',2), grid on
                    xlabel('z/\lambda_D'), ylabel('Pitch Angle'),grid on, hold on
                    title(['B = ',num2str(B0(k)),'T'])
                    
                    
                    
     figure(8);
    subplot(4,1,k);
       plot(z(1:max_z),pitch_length,'linewidth',2), grid on
                    xlabel('z/\lambda_D'), ylabel('Pitch Length'),grid on, hold on
                    title(['B = ',num2str(B0(k)),'T'])
                    
       figure(9)
   subplot(2,2,k);
   plot(dis,FE,'linewidth',2),grid on, hold on
   plot(dis,Fz,'linewidth',2),grid on, hold on
   
   
   plot(dis,FE+Fz,'linewidth',2),grid on, hold on
   plot(dis,Fx,'linewidth',2),grid on, hold on
   plot(dis,Fy,'linewidth',2),grid on, hold on

   xlabel('distance'),ylabel('Force')
   legend('FE','Fz','FE+Fz','Fx','Fy') 
   title(['B = ',num2str(B0(k)),'T'])
   
   figure(10)
   if k == 1
       plot(dis,FE,'linewidth',2),grid on, hold on
       plot(dis,Fz,'linewidth',2),grid on, hold on


       plot(dis,FE+Fz,'linewidth',2),grid on, hold on
       plot(dis,Fx,'linewidth',2),grid on, hold on
       plot(dis,Fy,'linewidth',2),grid on, hold on

       xlabel('distance'),ylabel('Force')
       legend('FE','Fz','FE+Fz','Fx','Fy')
   end
   
   figure(11)
   subplot(2,2,k)
   plot(dis,0.5*u.^2,'linewidth',2), grid on, hold on
   plot(dis,0.5*v.^2,'linewidth',2), grid on, hold on
   plot(dis,0.5*w.^2,'linewidth',2), grid on, hold on
   plot(dis,E,'linewidth',2), grid on, hold on
   xlabel('z/\lambda_D'), ylabel('Energy')
   %h=sprintf('B=%d',B0);
  title(['B = ',num2str(B0(k)),'T'])
   legend('u^2','v^2','w^2','Total energy')
    %    figure(4)
    %    plot(z(1:max_z),E,'linewidth',2), 
    %    xlabel('z/\lambda_D'), ylabel('Energy'),grid on
end
