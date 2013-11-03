function coeff
clc; clear all; clf
addpath '/Users/mquezada/Desktop/dispersion-relation/chebfun'

% compute material parameters
K1=5/8; p1=8/5;
K2=5/2; p2=2/5;
K=@(y)(K1+K2)/2+abs(K1-K2)/2*sin(2*pi*y);
p=@(y)1./K(y); 
%p=@(y)(p1+p2)/2+abs(p1-p2)/2*sin(2*pi*y);

pm=quad(p,0,1);
ph=quad(@(y)1./p(y),0,1)^-1
Kh=quad(@(y)1./K(y),0,1)^-1

disp('computing A')
Ap=@(y,A) Kh./K(y)-ph./p(y);
[y,A]=bvp(Ap);
plot(y,A); pause(1)

disp('computing H')
Hp=@(yd,H) ph^-1*p(yd).*get(y,A,yd);
[~,H]=bvp(Hp);
plot(y,H); pause(1)

% compute coefficients
disp('*********************************')
disp('*********************************')
disp('****** COMPUTE COEFFICIENTS *****')
disp('*********************************')
disp('*********************************')
disp('')

format long
% c-dispersion
alpha2=-Kh*average(y,K(y).^-1.*H)
beta2=ph*average(y,p(y).^-1.*H)

function [x,y] = bvp(fun)
    x=linspace(0,1,100);
    init.x=x; init.y=0*x;
    bc=@(l,r) l(1)-1;
    sol = bvp4c(fun,bc,init);
    mean=average(sol.x,sol.y(1,:));
    x=sol.x; y=sol.y(1,:)-mean;
end

function average = average(y,F)
    average=quad(@(yd)get(y,F,yd),0,1,eps);
end

function interp = get(y,F,yi)
    interp = interp1(y,F,yi,'cubic');
end


end