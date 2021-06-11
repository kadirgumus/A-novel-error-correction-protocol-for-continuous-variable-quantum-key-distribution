function [ SKR_theory_QPSK, SKR_theory_QPSK_recon, I_AB_het] = calcSKR( Vmod, neff, T, eps, vel, Beta, FER)
%CALCSKR Summary of this function goes here
%   Detailed explanation goes here

G = @(x) (x+1)*log2(x+1)-x*log2(x); %Given below Eq(5) in [1].

eps = max(eps,0);

n = 1;
V = Vmod+1;

xline = 1/T - 1 + eps;
xhet = (1+(1-neff)+2*vel)/neff;
xtot_het = xline+xhet/T;  
% VB_het = neff*T*(V+xtot_het)/2;
% VBA_het =  neff*T*(1+xtot_het)/2;




I_AB_het(n) = log2( (V+xtot_het) / (1 + xtot_het));


ZG = sqrt((Vmod^2+2*Vmod));

alpha = sqrt(Vmod/2);
e0 = 1/2*exp(-alpha^2)*(cosh(alpha^2) + cos(alpha^2));
e2 = 1/2*exp(-alpha^2)*(cosh(alpha^2) - cos(alpha^2));
e1 = 1/2*exp(-alpha^2)*(sinh(alpha^2) + sin(alpha^2));
e3 = 1/2*exp(-alpha^2)*(sinh(alpha^2) - sin(alpha^2));

Z4 = 2*alpha^2 * (e0^(3/2)*e1^(-1/2) + e1^(3/2)*e2^(-1/2) + e2^(3/2)*e3^(-1/2) + e3^(3/2)*e0^(-1/2));
A2 = V^2 + T^2*(V+xline)^2 - 2*T*Z4^2;
B2 = (T*V^2 + T*V*xline - T*Z4^2)^2;

C2 =  ( A2*xhet^2 + B2 + 1 + 2*xhet * (V *sqrt(B2) + T * (V+xline)) + 2*T*Z4^2 )  / (T*(V + xtot_het))^2;
D2 =  ((V + sqrt(B2)*xhet) / (T*(V + xtot_het)))^2;
lambda1_2 = sqrt( 1/2*(A2+sqrt(A2^2-4*B2)) ); %Equation (8) in [1]
lambda2_2 = sqrt( 1/2*(A2-sqrt(A2^2-4*B2)) ); %Equation (8) in [1]
lambda3_2 = sqrt( 1/2*(C2+sqrt(C2^2-4*D2)) ); %Equation (11) in [1]
lambda4_2 = sqrt( 1/2*(C2-sqrt(C2^2-4*D2)) ); %Equation (11) in [1]
H_BEhat2(n) = G((lambda1_2-1)/2)+G((lambda2_2-1)/2)-G((lambda3_2-1)/2)-G((lambda4_2-1)/2) - 0;

deltaI_Holevo_het_QPSK =  I_AB_het(n) - H_BEhat2(n);
SKR_theory_QPSK(n) = deltaI_Holevo_het_QPSK;
deltaI_Holevo_recon_het_QPSK =  (Beta*I_AB_het(n) - H_BEhat2(n))*(1-FER);
SKR_theory_QPSK_recon(n) = deltaI_Holevo_recon_het_QPSK;





