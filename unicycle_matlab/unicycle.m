clear all
clc

syms t a2 te2 g a1 te1 a te md mb mw Rw hb hd tau_a Iw Ib Id

%% Model
%yw = Rw*sin(te); zw = Rw*cos(te); 
%yb = hb*sin(te); zb = hb*cos(te);
%yd = hd*sin(te); zd = hd*cos(te);

%yw1 = Rw*cos(te)*te1; zw1 = -Rw*sin(te)*te1; 
%yb1 = hb*cos(te)*te1; zb1 = -hb*sin(te)*te1;
%yd1 = hd*cos(te)*te1; zd1 = -hd*sin(te)*te1;

%Tw = (mw*(yw1*yw1 + zw1*zw1)+Iw*te1*te1)/2;
%Tb = (mb*(yb1*yb1 + zb1*zb1)+Ib*te1*te1)/2;
%Td = (md*(yd1*yd1 + zd1*zd1)+Id*(te1*te1+a1*a1))/2;
%T = Tw + Tb + Td;
%V = mw*g*zw + mb*g*zb + md*g*zd;
%L = simplify(subs(T - V,[a1*t,t*te1],[sym('a'),sym('te')]))

%L = (Id*a1^2)/2 + (Ib*te1^2)/2 + (Id*te1^2)/2 + (Iw*te1^2)/2 +(Rw^2*mw*te1^2)/2 + (hb^2*mb*te1^2)/2 + (hd^2*md*te1^2)/2 - Rw*g*mw*cos(te) - g*hb*mb*cos(te) - g*hd*md*cos(te);

%dLdte = diff(L,sym('te'))
%dLdte1 = subs(diff(L,sym('te1')),[sym('a1'),sym('te1'),sym('a'),sym('te')],[a2*t,t*te2,a1*t,t*te1])
%dLdte1dt = subs(diff(dLdte1,sym('t')),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')])

%dLda = diff(L,sym('a'))
%dLda1 = subs(diff(L,sym('a1')),[sym('a1'),sym('te1'),sym('a'),sym('te')],[a2*t,t*te2,a1*t,t*te1])
%dLda1dt = subs(diff(dLda1,sym('t')),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')])

%te2 = subs(simplify(solve(dLdte1dt - dLdte==-tau_a,te2)),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')])
%a2 = subs(simplify(solve(dLda1dt - dLda==tau_a,a2)),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')])

%a2 = tau_a/Id;
% te,a-->0, sin(a)=a, sin(te)=te, cos(a)=cos(te)=1, te1^2=a1^2=0
%te2 = (Rw*g*mw*sin(te) - tau_a + g*hb*mb*sin(te) + g*hd*md*sin(te))/(mw*Rw^2 + mb*hb^2 + md*hd^2 + Ib + Id + Iw);

%% проверка
%a2 = subs(a2,[OA,AB,mB,g],[0.5,0.25,10,9.821]);
%te2 = subs(te2,[OA,AB,mB,g],[0.5,0.25,10,9.821]);

%a2 = subs(a2,[M,te,te1,a,a1],[M,0,0,0,0])
%te2 = subs(te2,[M,te,te1,a,a1],[M,0,0,0,0])

%% State space implementaion
x_equi = [0,0,0,0,0.05,0.17,0.34,0.5,1.5,0.85,0.00033,0.00425,0.00272,9.821]; 
%x_equi = [0,0,0,0,sym('OA'),sym('AB'),sym('mB'),sym('g')]; 

x_1 = sym('a');
x_2 = sym('a1');
x_3 = sym('te');
x_4 = sym('te1');
x = [x_1, x_2, x_3, x_4];
x1 = [x_2, a2, x_4, te2]; 

%A = subs(jacobian(x1, x),[a,a1,te,te1,Rw,hb,hd,mw,mb,md,Iw,Ib,Id,g],x_equi)
%B = subs(jacobian(x1, sym('tau_a')),[a,a1,te,te1,Rw,hb,hd,mw,mb,md,Iw,Ib,Id,g],x_equi)
C = [1 1 1 1];
D = 0;

A = [0 1 0 0; 0 0 0 0; 0 0 0 1;0 0 5588149/150160 0];
B = [0;  6250/17; 0; -12500/1877];

%sys = ss(A,B,C,D);

%% Pole placement
%ctrb_matrix = ctrb(A,B);
%rank(ctrb_matrix)
%rank(A)

%eig(A)
%pole(sys)

%desired eigs
eigs = [-5; -5; -5; -5];

K = acker(A,B,eigs)

Q = [1 0 0 0; 
     0 1 0 0; 
     0 0 100 0;
     0 0 0 100];

%R = 0.1;

%K_lqr = lqr(A,B,Q,R)

%eig(A-B*K_lqr)

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








 
