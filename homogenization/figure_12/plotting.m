clc

T=120;
nt=T/0.5+1; %dt=0.5

figure(1); clf; hold on
set(1, 'Position', [50 50 800 300])

% homogenized solution
load('solution.mat')
s=U(nt,:);
x=x-200;
xi=linspace(0,200,1000); si=spline(x,s,xi);
plot(xi,si,'-b','linewidth',2); 

% format plot
xlabel('x','fontsize',20);
ylabel('Stress','fontsize',20);
title(['t=' num2str(0.5*(nt-1))],'fontsize',20);
axis([150 200 -0.1 0.8])
set(gca,'fontsize',20)
