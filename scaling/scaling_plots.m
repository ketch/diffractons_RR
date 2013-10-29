if exist('u1')==0
    load_data;
end
%Produce figures in section on scaling.
iy = 16;

a1=max(abs(u1(:,iy))); 
a2=max(abs(u2(:,iy)));
a3=max(abs(u3(:,iy)));
a4=max(abs(u4(:,iy)));
a5=max(abs(u5(:,iy)));
a6=max(abs(u6(:,iy)));

% x index of max point
[dummy,x_im1]=max(sig1(:,iy)); [dummy,x_im2]=max(sig2(:,iy)); [dummy,x_im3]=max(sig3(:,iy)); 
[dummy,x_im4]=max(sig4(:,iy)); [dummy,x_im5]=max(sig5(:,iy)); [dummy,x_im6]=max(sig6(:,iy));

xm1=x(x_im1); xm2=x(x_im2); xm3=x(x_im3);
xm4=x(x_im4); xm5=x(x_im5); xm6=x(x_im6);

% scale stegotons
afac = 1.;
%afac = 1.05;
sig1s = 1/a1^afac*sig1;
sig2s = 1/a2^afac*sig2;
sig3s = 1/a3^afac*sig3;
sig4s = 1/a4^afac*sig4;
sig5s = 1/a5^afac*sig5;
sig6s = 1/a6^afac*sig6;

u1s = 1/a1*u1;
u2s = 1/a2*u2;
u3s = 1/a3*u3;
u4s = 1/a4*u4;
u5s = 1/a5*u5;
u6s = 1/a6*u6;

afac = 1.4;
v1s = 1/a1^afac*v1;
v2s = 1/a2^afac*v2;
v3s = 1/a3^afac*v3;
v4s = 1/a4^afac*v4;
v5s = 1/a5^afac*v5;
v6s = 1/a6^afac*v6;

xx1s = sqrt(a1)*(xx-xm1);
xx2s = sqrt(a2)*(xx-xm2);
xx3s = sqrt(a3)*(xx-xm3);
xx4s = sqrt(a4)*(xx-xm4);
xx5s = sqrt(a5)*(xx-xm5);
xx6s = sqrt(a6)*(xx-xm6);

xx1sb = a1*(xx-xm1);
xx2sb = a2*(xx-xm2);
xx3sb = a3*(xx-xm3);
xx4sb = a4*(xx-xm4);
xx5sb = a5*(xx-xm5);
xx6sb = a6*(xx-xm6);

%% Plot of A(y)'s 
figure(4); clf; hold on
y=linspace(0,1,length(u1(x_im1,:)));
plot(y,abs(u1(x_im1,:)),'-b','linewidth',2)
plot(y,abs(u2(x_im2,:)),'-r','linewidth',2)
plot(y,abs(u3(x_im3,:)),'-c','linewidth',2)
plot(y,abs(u4(x_im4,:)),'-k','linewidth',2)
plot(y,abs(u5(x_im5,:)),'-g','linewidth',2)
plot(y,abs(u6(x_im6,:)),'-m','linewidth',2)
title('A(y)','fontsize',20)
set(gcf, 'Position', [50 50 800 300])
set(gca,'fontsize',20)

%% Plot for sig
figure(1); clf; hold on
plot(xx1s(:,iy),sig1s(:,iy),'b','linewidth',1)
plot(xx2s(:,iy),sig2s(:,iy),'r','linewidth',1)
plot(xx3s(:,iy),sig3s(:,iy),'c','linewidth',1)
plot(xx4s(:,iy),sig4s(:,iy),'k','linewidth',1)
plot(xx5s(:,iy),sig5s(:,iy),'g','linewidth',1)
plot(xx6s(:,iy),sig6s(:,iy),'m','linewidth',1)
yy = 1./a1*(sech(xx1s*1.5).^2);
plot(xx1s,yy,'--k','linewidth',2)
set(1, 'Position', [50 50 800 300])
set(gca,'fontsize',20)
xlim([-3,3])
xlabel('x'); ylabel('\sigma')
set(gcf, 'PaperPosition', [0 0 10 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [10 5]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'scaling_stress', 'pdf') %Save figure

%% Plot for u
figure(2); clf; hold on
plot(xx1s(:,iy),u1s(:,iy),'b','linewidth',1)
plot(xx2s(:,iy),u2s(:,iy),'r','linewidth',1)
plot(xx3s(:,iy),u3s(:,iy),'c','linewidth',1)
plot(xx4s(:,iy),u4s(:,iy),'k','linewidth',1)
plot(xx5s(:,iy),u5s(:,iy),'g','linewidth',1)
plot(xx6s(:,iy),u6s(:,iy),'m','linewidth',1)
set(2, 'Position', [50 50 800 300])
set(gca,'fontsize',20)
xlim([-3,3])
xlabel('x'); ylabel('u')
set(gcf, 'PaperPosition', [0 0 10 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [10 5]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'scaling_u', 'pdf') %Save figure

%% Plot for v
figure(3); clf; hold on
plot(xx1s(:,iy),v1s(:,iy),'b','linewidth',1)
plot(xx2s(:,iy),v2s(:,iy),'r','linewidth',1)
plot(xx3s(:,iy),v3s(:,iy),'c','linewidth',1)
plot(xx4s(:,iy),v4s(:,iy),'k','linewidth',1)
plot(xx5s(:,iy),v5s(:,iy),'g','linewidth',1)
plot(xx6s(:,iy),v6s(:,iy),'m','linewidth',1)
yy = -0.9*(sech(xx1s*1.5).^2.*tanh(xx1s*1.5));
plot(xx1s,yy,'--k','linewidth',2)
set(3, 'Position', [50 50 800 300])
set(gca,'fontsize',20)
xlim([-3,3])
xlabel('x'); ylabel('v')
set(gcf, 'PaperPosition', [0 0 10 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [10 5]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'scaling_v', 'pdf') %Save figure
