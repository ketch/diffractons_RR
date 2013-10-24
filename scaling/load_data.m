clf; clear all; clc

%create domain
x=linspace(0,300,9600);
y=linspace(0,1,128);
[yy,xx]=meshgrid(y,x);

disp('load coefficients')
K=load('K.txt');
rho=load('rho.txt');

disp('load sw1')
eps1=load('sw1_q0.txt');
sig1=exp(K.*eps1)-1;
u1=load('sw1_q1.txt');
v1=load('sw1_q2.txt');

disp('load sw2')
eps2=load('sw2_q0.txt');
sig2=exp(K.*eps2)-1;
u2=load('sw2_q1.txt');
v2=load('sw2_q2.txt');

disp('load sw3')
eps3=load('sw3_q0.txt');
sig3=exp(K.*eps3)-1;
u3=load('sw3_q1.txt');
v3=load('sw3_q2.txt');

disp('load sw4')
eps4=load('sw4_q0.txt');
sig4=exp(K.*eps4)-1;
u4=load('sw4_q1.txt');
v4=load('sw4_q2.txt');

disp('load sw5')
eps5=load('sw5_q0.txt');
sig5=exp(K.*eps5)-1;
u5=load('sw5_q1.txt');
v5=load('sw5_q2.txt');

disp('load sw6')
eps6=load('sw6_q0.txt');
sig6=exp(K.*eps6)-1;
u6=load('sw6_q1.txt');
v6=load('sw6_q2.txt');
