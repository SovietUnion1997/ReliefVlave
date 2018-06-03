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
eq5 = (V_4/beta_e*s+1/R_3+1/R_4)*p_4-1/R_4*p_5+A_4*s*y == p_3/R_3
eq6 = (V_3/beta_e*s+1/R_2+K_cy+1/R_3)*p_3 == p_C/R_2-(K_qy+A_3*s)*y+p_4/R_3
eq7 = (V_5/beta_e*s+1/R_4)*p_5 == p_4/R_4+A_5*s*y
sol1 = solve([eq1,eq2,eq3,eq4,eq5,eq6,eq7],[p_C,p_A,p_3,p_4,p_5,x,y])
% collect(sol1.p_C,s)
% res = collect(sol1.p_A,s);
toc
% [num,den] = numden(collect(sol1.p_A,s))
% [c1,t1] = coeffs(den,s)
% roots(c1) 这个求根没法看
% collect(sol1.x,s)

%%
res = collect(sol1.p_A,s);
% clear
m_x = 60e-3; 
m_y = 25e-3;

A_A = 1.7e-4;
A_C = 2e-4;
A_3 = 9.62e-6;
A_4 = 136.04e-6;    % 影响极点
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
R_2 = R_1;  % 握草 他妈改了这个参数就稳定了
% 原来R_2按照路溢流阀的写的
% 改成R_1就能稳定了

rho = 870;
nu = 46e-6; % 单位都不对
l3 = 20e-3;
l4 = 4e-3;
d3 = 0.8e-3;
d4 = 0.5e-3;
R_3 = 1/(pi*d3^4/(128*rho*nu*l3));
R_4 = 1/(pi*d4^4/(128*rho*nu*l4));

% R_3 = R_2;
% R_4 = R_3;
R_L = R_2;
F_EM = 1;

res = vpa(subs(res),3);
[num,den] = numden(res);
[c,t] = coeffs(den)
roots(c)


%%
s = tf('s');
sys1 = (5.12e38*s^5 + 9.48e44*s^4 + 5.65e49*s^3 + 1.99e54*s^2 + 3.03e58*s + 3.29e62)/(1.32e21*s^9 + 3.02e27*s^8 + 7.43e32*s^7 + 7.26e37*s^6 + 3.51e42*s^5 + 7.71e46*s^4 + 5.94e50*s^3 + 1.27e53*s^2 + 7.37e55*s + 1.09e57)

sys2 = zpk(sys1)
step(sys2)
sys2.Z{1}
sys2.P{1}
% 


% sys1 = (R_L*m_x*p_3*beta_e^2*s^2 + A_A*A_C*R_1*R_L*p_3*beta_e^2*s + A_C*K_q*R_1*R_L*p_3*beta_e^2)/(R_1*R_2*R_L*V_a*V_c*m_x*s^4 + (R_1*R_2*V_c*beta_e*m_x + R_1*R_L*V_a*beta_e*m_x + R_2*R_L*V_a*beta_e*m_x + R_2*R_L*V_c*beta_e*m_x + K_c*R_1*R_2*R_L*V_c*beta_e*m_x)*s^3 + (R_1*beta_e^2*m_x + R_2*beta_e^2*m_x + R_L*beta_e^2*m_x + K_c*R_1*R_L*beta_e^2*m_x + K_c*R_2*R_L*beta_e^2*m_x + A_A^2*R_1*R_2*R_L*V_c*beta_e + A_C^2*R_1*R_2*R_L*V_a*beta_e)*s^2 + (A_C^2*R_1*R_2*beta_e^2 + A_A^2*R_1*R_L*beta_e^2 + A_A^2*R_2*R_L*beta_e^2 + A_C^2*R_2*R_L*beta_e^2 - 2*A_A*A_C*R_2*R_L*beta_e^2 + A_C^2*K_c*R_1*R_2*R_L*beta_e^2 + A_A*K_q*R_1*R_2*R_L*V_c*beta_e)*s + A_A*K_q*R_1*R_L*beta_e^2 + A_A*K_q*R_2*R_L*beta_e^2 - A_C*K_q*R_2*R_L*beta_e^2)
% syszpk = zpk(sys1)
% 
% step(sys1);
