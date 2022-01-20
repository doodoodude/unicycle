clear all
clc

syms t b2 ga2 g b1 ga1 b ga md mb mw Rw hb hd tau_b Iw Ib Id

%!!! hd and hb counts from wheels center!!

%% Model
%xw = Rw*b; zw = Rw; 
%xb = xw + hb*sin(ga); zb = zw + hb*cos(ga);
%xd = xw + hd*sin(ga); zd = zw + hd*cos(ga);

%xw1 = Rw*b1; zw1 = 0; 
%xb1 = Rw*b1 + hb*cos(ga)*ga1; zb1 = -hb*sin(ga)*ga1;
%xd1 = Rw*b1 + hd*cos(ga)*ga1; zd1 = -hd*sin(ga)*ga1;

%Tw = (mw*xw1*xw1+Iw*b1*b1)/2;
%Tb = (mb*(xb1*xb1 + zb1*zb1)+Ib*ga1*ga1)/2;
%Td = (md*(xd1*xd1 + zd1*zd1)+Id*ga1*ga1)/2;
%T = Tw + Tb + Td;
%V = mw*g*zw + mb*g*zb + md*g*zd;
%L = simplify(subs(T - V,[b1*t,t*ga1],[sym('b'),sym('ga')]))

%dLdga = diff(L,sym('ga'))
%dLdga1 = subs(diff(L,sym('ga1')),[sym('b1'),sym('ga1'),sym('b'),sym('ga')],[b2*t,t*ga2,b1*t,t*ga1])
%dLdga1dt = subs(diff(dLdga1,sym('t')),[b2*t,t*ga2,b1*t,t*ga1],[sym('b1'),sym('ga1'),sym('b'),sym('ga')])

%dLdb = diff(L,sym('b'))
%dLdb1 = subs(diff(L,sym('b1')),[sym('b1'),sym('ga1'),sym('b'),sym('ga')],[b2*t,t*ga2,b1*t,t*ga1])
%dLdb1dt = subs(diff(dLdb1,sym('t')),[b2*t,t*ga2,b1*t,t*ga1],[sym('b1'),sym('ga1'),sym('b'),sym('ga')])

%ga2_1 = subs(simplify(solve(dLdga1dt - dLdga==0,ga2)),[b2*t,t*ga2,b1*t,t*ga1],[sym('b1'),sym('ga1'),sym('b'),sym('ga')]);
%ga2_2 = subs(simplify(solve(dLdb1dt - dLdb==tau_b,ga2)),[b2*t,t*ga2,b1*t,t*ga1],[sym('b1'),sym('ga1'),sym('b'),sym('ga')]);

%b2 = simplify(solve(ga2_1==ga2_2,b2))
%ga2 = simplify(subs(ga2_1,sym('b2'),b2))

%% dLdga1dt - dLdga==0
%b2 =((tau_b + Rw*ga1^2*hb*mb*sin(ga) + Rw*ga1^2*hd*md*sin(ga))/(Rw*cos(ga)*(hb*mb + hd*md)) - (g*sin(ga)*(hb*mb + hd*md))/(mb*hb^2 + md*hd^2 + Ib + Id))/((Iw + Rw^2*mb + Rw^2*md + Rw^2*mw)/(Rw*cos(ga)*(hb*mb + hd*md)) - (Rw*cos(ga)*(hb*mb + hd*md))/(mb*hb^2 + md*hd^2 + Ib + Id));
%ga2 =((g*sin(ga) - (Rw*cos(ga)*((tau_b + Rw*ga1^2*hb*mb*sin(ga) + Rw*ga1^2*hd*md*sin(ga))/(Rw*cos(ga)*(hb*mb + hd*md)) - (g*sin(ga)*(hb*mb + hd*md))/(mb*hb^2 + md*hd^2 + Ib + Id)))/((Iw + Rw^2*mb + Rw^2*md + Rw^2*mw)/(Rw*cos(ga)*(hb*mb + hd*md)) - (Rw*cos(ga)*(hb*mb + hd*md))/(mb*hb^2 + md*hd^2 + Ib + Id)))*(hb*mb + hd*md))/(mb*hb^2 + md*hd^2 + Ib + Id);
 

%% dLdga1dt - dLdga== -tau_b  
%b2 = -((g*hb*mb*sin(ga) - tau_b + g*hd*md*sin(ga))/(mb*hb^2 + md*hd^2 + Ib + Id) - (tau_b + Rw*ga1^2*hb*mb*sin(ga) + Rw*ga1^2*hd*md*sin(ga))/(Rw*cos(ga)*(hb*mb + hd*md)))/((Iw + Rw^2*mb + Rw^2*md + Rw^2*mw)/(Rw*cos(ga)*(hb*mb + hd*md)) - (Rw*cos(ga)*(hb*mb + hd*md))/(mb*hb^2 + md*hd^2 + Ib + Id));
%ga2 = -((sin(2*ga)*Rw^2*ga1^2*hb^2*mb^2)/2 + sin(2*ga)*Rw^2*ga1^2*hb*hd*mb*md + (sin(2*ga)*Rw^2*ga1^2*hd^2*md^2)/2 - g*sin(ga)*Rw^2*hb*mb^2 - g*sin(ga)*Rw^2*hb*mb*md - g*mw*sin(ga)*Rw^2*hb*mb - g*sin(ga)*Rw^2*hd*mb*md - g*sin(ga)*Rw^2*hd*md^2 - g*mw*sin(ga)*Rw^2*hd*md + tau_b*Rw^2*mb + tau_b*Rw^2*md + mw*tau_b*Rw^2 + tau_b*cos(ga)*Rw*hb*mb + tau_b*cos(ga)*Rw*hd*md - Iw*g*sin(ga)*hb*mb - Iw*g*sin(ga)*hd*md + Iw*tau_b)/(Ib*Iw + Id*Iw + Rw^2*hb^2*mb^2 + Rw^2*hd^2*md^2 + Ib*Rw^2*mb + Ib*Rw^2*md + Id*Rw^2*mb + Id*Rw^2*md + Ib*Rw^2*mw + Id*Rw^2*mw + Iw*hb^2*mb + Iw*hd^2*md + Rw^2*hb^2*mb*md + Rw^2*hd^2*mb*md + Rw^2*hb^2*mb*mw + Rw^2*hd^2*md*mw - Rw^2*hb^2*mb^2*cos(ga)^2 - Rw^2*hd^2*md^2*cos(ga)^2 - 2*Rw^2*hb*hd*mb*md*cos(ga)^2);
 
%% проверка
%a2 = subs(a2,[OA,AB,mB,g],[0.5,0.25,10,9.821]);
%te2 = subs(te2,[OA,AB,mB,g],[0.5,0.25,10,9.821]);

%a2 = subs(a2,[M,te,te1,a,a1],[M,0,0,0,0])
%te2 = subs(te2,[M,te,te1,a,a1],[M,0,0,0,0])

%% State space implementaion
%!!! hd and hb counts from wheels center!!!!!!!!!!!!1
x_equi = [0,0,0,0,0.05,0.12,0.29,0.5,1.5,0.85,0.000625,0.00425,0.0013883,9.821];
%x_equi = [0,0,0,0,sym('OA'),sym('AB'),sym('mB'),sym('g')]; 

x_1 = sym('b');
x_2 = sym('b1');
x_3 = sym('ga');
x_4 = sym('ga1');
x = [x_1, x_2, x_3, x_4];
x1 = [x_2, b2, x_4, ga2]; 

%A = subs(jacobian(x1, x),[b,b1,ga,ga1,Rw,hb,hd,mw,mb,md,Iw,Ib,Id,g],x_equi)
%B = subs(jacobian(x1, sym('tau_b')),[b,b1,ga,ga1,Rw,hb,hd,mw,mb,md,Iw,Ib,Id,g],x_equi)
%C = [1 1 1 1];
%D = 0;

A = [0 1 0 0; 0 0 -287.8141 0; 0 0 0 1;0 0 104.5983 0];

%% dLdga1dt - dLdga==0
B = [0;  318.1032; 0; -68.7128];

%% dLdga1dt - dLdga== -tau_b   
%B = [0;  386.8159; 0; -93.6846];


%sys = ss(A,B,C,D);

%% Pole placement
%ctrb_matrix = ctrb(A,B)
%rank(ctrb_matrix)
%rank(A)

%eig(A)
%pole(sys)

eigs = [-7.5; -6; -8; -9];
eigs1 = [-229.07; -2.85; -2.86; -1.93];
%K = place(A,B,eigs)
%eig(A-B*K)

Q = [0.2 0 0 0; 
     0 0.003 0 0; 
     0 0 7 0;
     0 0 0 0.15];

R = 1.8;

K_lqr = lqr(A,B,Q,R)

eig(A-B*K_lqr)

%% discrete system

%Ts = 0.005; % 200 Hz
%dsys = c2d(sys,Ts);

%A_d = dsys.a;
%B_d = dsys.b;
%C_d = dsys.c;
%D_d = dsys.d;

%rank(ctrb(A_d,B_d));
%eig(A_d);
%eigs = [-3; -5; -7; -2];
%K_d = place(A_d,B_d,eigs)
%eig(A_d-B_d*K_d);








 
