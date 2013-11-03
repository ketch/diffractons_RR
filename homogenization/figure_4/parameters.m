function parameters(KA,pA,KB,pB)

global alpha2  beta2 Kh ph

cA=sqrt(KA/pA); cB=sqrt(KB/pB);

Km=(KA+KB)/2;
Kh=2*KA*KB/(KA+KB);
ph=2*pA*pB/(pA+pB);
pm=(pA+pB)/2;

alpha2=(KA-KB)/(192*Km^2)*(cA^2-cB^2)*pm;
beta2=-(pA-pB)/(192*Km)*(cA^2-cB^2);
