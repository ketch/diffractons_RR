function []=solver()
close('all'); clc; clear all; clf

%To pass parameters to Au
global alpha2  beta2 Kh ph

parameters(5/8,8/5,5/2,2/5)

dt=0.01;    %time step
tf=120;      %final time
td=0.5;     %time interval to display 

save_solution = 1; %flag to save solution
name_solution='solution.mat';

% physical domain
x_lower=0; x_upper=400;

mx=2^10; %Number of Fourier modes
Lx=x_upper-x_lower;
kx = (2*pi/Lx)*[0:(mx/2-1) (-mx/2):-1]; % Wavenumber vector in x

%discretized domain
dx = (x_upper-x_lower)/mx;
x = (0:(mx-1))*dx;

nit=floor(tf/dt); %number of iterations

%initial conditions
A=1;
x0=(x_upper-x_lower)/2;
varx=5;
s=A*exp(-(x-x0).^2/(2*varx)); %IC hom in y

u(1,:)=x.*0; %u
u(2,:)=s;   %sig

plot(x,s)
pause(2)

U(1,:)=s;
index=1;
for i=1:nit
    if(i*dt>=20 && i*dt<=21)
        u(:,1:end/2)=0;
    end
    disp('*********************************')
    disp('*********************************')
    disp(['Time step ' num2str(i) '. Time t=' num2str(i*dt)])
    % Four stages Runge-Kutta
    D1u=dt.*ps_discretization(u,kx);    
    D2u=dt.*ps_discretization(u+0.5*D1u,kx);
    D3u=dt.*ps_discretization(u+0.5*D2u,kx);
    D4u=dt.*ps_discretization(u+D3u,kx);
    u = u + (D1u+2*D2u+2*D3u+D4u)/6;
    if((i*dt-index*td)>=0)
        s=squeeze(u(2,:));
        plot(x,s)
        title(['t=' num2str(dt*i)]);
        U(index+1,:)=s;
        index=index+1;
        pause(0.1)        
    end
end
if save_solution == 1
    save(name_solution)
end
