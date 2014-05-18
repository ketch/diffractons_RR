function get_speed
clc

% material parameters
K1=5./8;   rho1=8./5;                                                                                                                                                                                      
K2=5./2;   rho2=2./5;
Kh=2*K1*K2/(K1+K2);
rhoh=2*rho1*rho2/(rho1+rho2);
ceff=sqrt(Kh/rhoh);

data=load('sw1.txt');
tm1=data(:,1); xm1=data(:,2); sm1=mean(data(:,3)); um1=mean(data(:,4));
data=load('sw2.txt');
tm2=data(:,1); xm2=data(:,2); sm2=mean(data(:,3)); um2=mean(data(:,4));
data=load('sw3.txt');
tm3=data(:,1); xm3=data(:,2); sm3=mean(data(:,3)); um3=mean(data(:,4));
data=load('sw4.txt');
tm4=data(:,1); xm4=data(:,2); sm4=mean(data(:,3)); um4=mean(data(:,4));
data=load('sw5.txt');
tm5=data(:,1); xm5=data(:,2); sm5=mean(data(:,3)); um5=mean(data(:,4));
data=load('sw6.txt');
tm6=data(:,1); xm6=data(:,2); sm6=mean(data(:,3)); um6=mean(data(:,4));
data=load('sw7.txt');
tm7=data(:,1); xm7=data(:,2); sm7=mean(data(:,3)); um7=mean(data(:,4));
data=load('sw8.txt');
tm8=data(:,1); xm8=data(:,2); sm8=mean(data(:,3)); um8=mean(data(:,4));

xm1=shift_max_position(xm1);
xm2=shift_max_position(xm2);
xm3=shift_max_position(xm3);
xm4=shift_max_position(xm4);
xm5=shift_max_position(xm5);
xm6=shift_max_position(xm6);
xm7=shift_max_position(xm7);
xm8=shift_max_position(xm8);

% get speed
s1=mean(((xm1(3:end)-xm1(1:end-2))./(tm1(3:end)-tm1(1:end-2))));
s2=mean(((xm2(3:end)-xm2(1:end-2))./(tm2(3:end)-tm2(1:end-2))));
s3=mean(((xm3(3:end)-xm3(1:end-2))./(tm3(3:end)-tm3(1:end-2))));
s4=mean(((xm4(3:end)-xm4(1:end-2))./(tm4(3:end)-tm4(1:end-2))));
s5=mean(((xm5(3:end)-xm5(1:end-2))./(tm5(3:end)-tm5(1:end-2))));
s6=mean(((xm6(3:end)-xm6(1:end-2))./(tm6(3:end)-tm6(1:end-2))));
s7=mean(((xm7(3:end)-xm7(1:end-2))./(tm7(3:end)-tm7(1:end-2))));
s8=mean(((xm8(3:end)-xm8(1:end-2))./(tm8(3:end)-tm8(1:end-2))));


%% from original solitary waves
amplitudes_absu=abs([um8 um7 um6 um5 um4 um3 um2 um1]);
speeds=[s8 s7 s6 s5 s4 s3 s2 s1];
a=linspace(0,0.85,100);
s1_absu=get_fitting(a,speeds,amplitudes_absu,ceff,1);
s2_absu=get_fitting(a,speeds,amplitudes_absu,ceff,2);
s3_absu=get_fitting(a,speeds,amplitudes_absu,ceff,3);
figure(1)
clf; hold on
plot(a,s1_absu,'-r','linewidth',2)
plot(a,s2_absu,'--k','linewidth',2)
%plot(a,s3_absu,'--c','linewidth',1)
plot(amplitudes_absu,speeds,'s','linewidth',2,'markersize',20)
set(gca,'FontSize',20)
xlabel('Amplitude','fontsize',20); ylabel('Speed','fontsize',20)
%legend('linear','quadratic','cubic','location','best')
legend('linear fit','quadratic fit','location','best')
title('Speed-amplitude','fontsize',20)

%% from original sw and scaled sw
% read scaled sw
data=load('sw_1.25_speed.txt');
t_s1=data(:,1); x_s1=data(:,2); sig_s1=mean(data(:,3)); u_s1=mean(data(:,4));
data=load('sw_1.5_speed.txt');
t_s2=data(:,1); x_s2=data(:,2); sig_s2=mean(data(:,3)); u_s2=mean(data(:,4));
data=load('sw_1.75_speed.txt');
t_s3=data(:,1); x_s3=data(:,2); sig_s3=mean(data(:,3)); u_s3=mean(data(:,4));
data=load('sw_2_speed.txt');
t_s4=data(:,1); x_s4=data(:,2); sig_s4=mean(data(:,3)); u_s4=mean(data(:,4));

x_s1=shift_max_position(x_s1);
x_s2=shift_max_position(x_s2);
x_s3=shift_max_position(x_s3);
x_s4=shift_max_position(x_s4);

% get speed
s_s1=mean(((x_s1(3:end)-x_s1(1:end-2))./(t_s1(3:end)-t_s1(1:end-2))));
s_s2=mean(((x_s2(3:end)-x_s2(1:end-2))./(t_s2(3:end)-t_s2(1:end-2))));
s_s3=mean(((x_s3(3:end)-x_s3(1:end-2))./(t_s3(3:end)-t_s3(1:end-2))));
s_s4=mean(((x_s4(3:end)-x_s4(1:end-2))./(t_s4(3:end)-t_s4(1:end-2))));

amplitudes=[amplitudes_absu abs([u_s1 u_s2 u_s3 u_s4])];
speeds=[s8 s7 s6 s5 s4 s3 s2 s1 s_s1 s_s2 s_s3 s_s4];
sig = [sm8 sm7 sm6 sm5 sm4 sm3 sm2 sm1 sig_s1 sig_s2 sig_s3 sig_s4];

figure(2)
clf; hold on
plot(amplitudes_absu,[s8 s7 s6 s5 s4 s3 s2 s1],'s','linewidth',2,'markersize',20)
plot(abs([u_s1 u_s2 u_s3 u_s4]),[s_s1 s_s2 s_s3 s_s4],'o','linewidth',2,'markersize',20)
plot(amplitudes,sig./amplitudes,'xr','linewidth',2,'markersize',15)
set(gca,'FontSize',20)
xlabel('Amplitude','fontsize',20); ylabel('Speed','fontsize',20)
title('Speed-amplitude','fontsize',20)

function s = get_fitting(a, speeds, amplitudes, ceff, order)
    switch order
        case 1
            poly = amplitudes'\(speeds'-ceff*ones(length(speeds),1));
            s = poly(1)*a+ceff;            
        case 2
            poly = [amplitudes'.^2 amplitudes']\(speeds'-ceff*ones(length(speeds),1));
            s = poly(1)*a.^2+poly(2)*a+ceff;            
        case 3
            poly = [amplitudes'.^3 amplitudes'.^2 amplitudes']\(speeds'-ceff*ones(length(speeds),1));
            s = poly(1)*a.^3+poly(2)*a.^2+poly(3)*a+ceff;
        case 4
            poly = [amplitudes'.^4 amplitudes'.^3 amplitudes'.^2 amplitudes']\(speeds'-ceff*ones(length(speeds),1));  
            s = poly(1)*a.^4+poly(2)*a.^3+poly(3)*a.^2+poly(4)*a+ceff;
    end
end

function xm_shifted = shift_max_position(xm)
    is_shift = ( xm(2:end)-xm(1:end-1) < 0 );
    shifts_loc = find(is_shift)+1;
    for i=1:length(shifts_loc)
        xm(shifts_loc(i):end)=xm(shifts_loc(i):end)+300;
    end
    xm_shifted=xm;
end

end
