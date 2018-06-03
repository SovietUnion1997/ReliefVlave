% 为什么不直接解方程呢
clear
tic
syms p_C p_A p_3 p_4 p_5 x y F_EM 
syms m_x m_y
syms s
syms A_A A_C A_3 A_4 A_5
syms V_c V_a V_3 V_4 V_5
syms beta_e K_c K_q K_cy K_qy
syms R_1 R_2 R_3 R_4 R_L

eq1 = m_x*s^2*x == A_A*p_A - A_C*p_C
eq2 = (V_c/beta_e*s+1/R_1+1/R_2)*p_C - p_A/R_1- p_3/R_2 == A_C*s*x
eq3 = (V_a/beta_e*s+1/R_1+1/R_L + K_c)*p_A == p_C/R_1-(K_q+A_A*s)*x
eq4 = F_EM +A_5*p_5-A_3*p_3-A_4*p_4 == -m_y*s^2*y
eq5 = (V_4/beta_e*s+1/R_2+1/R_4+K_cy)*p_4-1/R_4*p_5+(A_4*s+K_qy)*y...
    == p_C*R_2
% eq6 = (V_3/beta_e*s+1/R_2+K_cy+1/R_3)*p_3 == p_C/R_2-(K_qy+A_3*s)*y+p_4/R_3
eq7 = (V_5/beta_e*s+1/R_4)*p_5 == p_4/R_4+A_5*s*y
sol1 = solve([eq1,eq2,eq3,eq4,eq5,eq6,eq7],[p_C,p_A,p_3,p_4,p_5,x,y])
% collect(sol1.p_C,s)
res = collect(sol1.p_A,s);
toc
% [num,den] = numden(collect(sol1.p_A,s))
% [c1,t1] = coeffs(den,s)
% roots(c1) 这个求根没法看
% collect(sol1.x,s)

%%

% clear
m_x = 0.06; 
m_y = 8e-4;

A_A = 1.7e-4;
A_C = 2e-4;
A_3 = 9.62e-6;
A_4 = 136.04e-6;
A_5 = 140.95e-6;

V_c = 1e-6;
V_a = 2e-3;
V_3 = 1.89e-7;
V_4 = 9.78e-8;
V_5 = 1.06e-7;

beta_e = 7e8;
K_c = 5e-10;
K_q = 3.3;
K_cy = K_c;
K_qy = K_q;

R_1 = 3e10;
R_2 = 2e11;
R_3 = R_2;
R_4 = R_3;
R_L = R_2;
F_EM = 1;

vpa(subs(res),3)
s_ = tf('s');

sys1 = -(1.64e40*R_4*s_^5 + (1.48e46*R_4 + 2.26e56)*s_^4 + (4.54e50*R_4 + 8.32e61)*s_^3 + (8.22e54*R_4 + 2.76e66)*s_^2 + (6.65e58*R_4 + 4.88e70)*s_ + 4.39e74)/(1.0*s_^3*(5.84e72*R_4 + 1.82e83) - s_^8*(1.76e28*R_4 + 1.24e38) - s_^6*(7.18e58*R_4 - 1.58e69) - 9.03e21*R_4*s_^9 - s_*(2.5e76*R_4 - 5.54e89) + 1.0*s_^4*(2.98e68*R_4 + 7.65e78) + 1.0*s_^7*(2.39e53*R_4 - 2.4e44) + 1.0*s_^5*(1.04e64*R_4 + 6.14e73) + 1.0*s_^2*(1.07e75*R_4 + 3.02e87) - 1.29e91);

sys2 = zpk(sys1)
% 


% sys1 = (R_L*m_x*p_3*beta_e^2*s^2 + A_A*A_C*R_1*R_L*p_3*beta_e^2*s + A_C*K_q*R_1*R_L*p_3*beta_e^2)/(R_1*R_2*R_L*V_a*V_c*m_x*s^4 + (R_1*R_2*V_c*beta_e*m_x + R_1*R_L*V_a*beta_e*m_x + R_2*R_L*V_a*beta_e*m_x + R_2*R_L*V_c*beta_e*m_x + K_c*R_1*R_2*R_L*V_c*beta_e*m_x)*s^3 + (R_1*beta_e^2*m_x + R_2*beta_e^2*m_x + R_L*beta_e^2*m_x + K_c*R_1*R_L*beta_e^2*m_x + K_c*R_2*R_L*beta_e^2*m_x + A_A^2*R_1*R_2*R_L*V_c*beta_e + A_C^2*R_1*R_2*R_L*V_a*beta_e)*s^2 + (A_C^2*R_1*R_2*beta_e^2 + A_A^2*R_1*R_L*beta_e^2 + A_A^2*R_2*R_L*beta_e^2 + A_C^2*R_2*R_L*beta_e^2 - 2*A_A*A_C*R_2*R_L*beta_e^2 + A_C^2*K_c*R_1*R_2*R_L*beta_e^2 + A_A*K_q*R_1*R_2*R_L*V_c*beta_e)*s + A_A*K_q*R_1*R_L*beta_e^2 + A_A*K_q*R_2*R_L*beta_e^2 - A_C*K_q*R_2*R_L*beta_e^2)
% syszpk = zpk(sys1)
% 
% step(sys1);
