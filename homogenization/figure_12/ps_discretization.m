function Du = ps_discretization(u,k)

global alpha2  beta2 Kh ph

delta2=1;
NL=1;

uhat=fft(squeeze(u(1,:)));
shat=fft(squeeze(u(2,:)));
s=u(2,:);

dudx=real(ifft(1i*k.*uhat));
du2dx2=real(ifft(-k.^2.*uhat));
du3dx3=real(ifft(-1i*k.^3.*uhat));

dsdx=real(ifft(1i*k.*shat));
ds3dx3=real(ifft(-1i*k.^3.*shat));

Du(1,:)=ph^-1*    (dsdx  + delta2*beta2*ds3dx3);
Du(2,:)=Kh*       ((NL*s+1).*(dudx) ...
                            + delta2 *alpha2 *( (NL*s+1).*(du3dx3) + NL*2*dsdx.*(du2dx2) ));