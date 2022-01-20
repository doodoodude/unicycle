clear all
clc

syms OA AB mA mB t a2 te2 g a1 te1 a te M Iao Iab

%% Model
%a1 = a2*t;
%te1 = te2*t;
%a = a1*t;
%te = te1*t;

%Ioa = mA*OA*OA/3;
Iab = mB*AB*AB; 

%yA = OA*cos(te)/2; 
yB = OA*cos(te) + AB*cos(te+a); 

%vAx = te1*OA*cos(te)/2;
%vAy = -te1*OA*sin(te)/2;
vBx = te1*OA*cos(te) + (te1 + a1)*AB*cos(te+a);
vBy = -te1*OA*sin(te) - (te1 + a1)*AB*sin(te+a);

%Ek = expand((mB*(vBx*vBx + vBy*vBy) + Iab*a1*a1)/2);
Ek = expand(mB*(vBx*vBx + vBy*vBy)/2 + Iab*a1*a1/2);
Ep = expand(mB*g*yB);
L = simplify(subs(Ek - Ep,[a1*t,t*te1],[sym('a'),sym('te')]));

%L = (mB*(2*AB^2*a1^2 + 2*AB^2*a1*te1 + AB^2*te1^2 + 2*cos(a)*AB*OA*a1*te1 + 2*cos(a)*AB*OA*te1^2 - 2*g*cos(a + te)*AB + OA^2*te1^2 - 2*g*cos(te)*OA))/2;

dLdte = diff(L,sym('te'));
%dLdte= (mB*(2*AB*g*sin(a + te) + 2*OA*g*sin(te)))/2;
dLdte1 = subs(diff(L,sym('te1')),[sym('a1'),sym('te1'),sym('a'),sym('te')],[a2*t,t*te2,a1*t,t*te1]);
dLdte1dt = subs(diff(dLdte1,sym('t')),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')]);
%dLdte1dt = (mB*(2*AB^2*a2 + 2*AB^2*te2 + 2*OA^2*te2 + 2*AB*OA*a2*cos(a) + 4*AB*OA*te2*cos(a) - 2*AB*OA*a1^2*sin(a) - 4*AB*OA*a1*te1*sin(a)))/2;
   
dLda = diff(L,sym('a'));
%dLda = -(mB*(2*AB*OA*sin(a)*te1^2 + 2*AB*OA*a1*sin(a)*te1 - 2*AB*g*sin(a + te)))/2;
dLda1 = subs(diff(L,sym('a1')),[sym('a1'),sym('te1'),sym('a'),sym('te')],[a2*t,t*te2,a1*t,t*te1]);
dLda1dt = subs(diff(dLda1,sym('t')),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')]);
%dLda1dt = (mB*(4*AB^2*a2 + 2*AB^2*te2 + 2*AB*OA*te2*cos(a) - 2*AB*OA*a1*te1*sin(a)))/2;

te2_1 = subs(simplify(solve(dLdte1dt - dLdte==0,te2)),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')]);
te2_2 = subs(simplify(solve(dLda1dt - dLda==M,te2)),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')]);
a2 = simplify(solve(te2_1==te2_2,a2))
te2 = simplify(subs(te2_1,sym('a2'),a2))

%% диф.уравнения при dLdte1dt - dLdte==-M
%a2 = 
%te2 = 

%% диф.уравнения при dLdte1dt - dLdte==0
%a2 = -((AB*OA*sin(a)*a1^2 + 2*AB*OA*te1*sin(a)*a1 + AB*g*sin(a + te) + OA*g*sin(te))/(AB^2 + 2*cos(a)*AB*OA + OA^2) - (- AB*OA*mB*sin(a)*te1^2 + M + AB*g*mB*sin(a + te))/(AB*mB*(AB + OA*cos(a))))/((2*AB)/(AB + OA*cos(a)) - (AB*(AB + OA*cos(a)))/(AB^2 + 2*cos(a)*AB*OA + OA^2)); 
%te2 = (2*AB^2*g*mB*sin(a + te) - 2*M*OA*cos(a) - 2*AB*M + 4*AB^2*OA*a1^2*mB*sin(a) + 2*AB^2*OA*mB*te1^2*sin(a) + 3*AB*OA*g*mB*sin(te) + AB*OA^2*mB*te1^2*sin(2*a) - AB*OA*g*mB*sin(2*a + te) + 8*AB^2*OA*a1*mB*te1*sin(a))/(2*AB*mB*(AB^2 - OA^2*cos(a)^2 + 2*OA^2 + 2*AB*OA*cos(a)));

%% проверка
%a2 = subs(a2,[OA,AB,mB,g],[0.5,0.25,10,9.821]);
%te2 = subs(te2,[OA,AB,mB,g],[0.5,0.25,10,9.821]);

%a2 = subs(a2,[M,te,te1,a,a1],[M,0,0,0,0])
%te2 = subs(te2,[M,te,te1,a,a1],[M,0,0,0,0])

%% State space implementaion
x_equi = [0,0,0,0,0.5,0.25,10,9.821]; 
%x_equi = [0,0,0,0,sym('OA'),sym('AB'),sym('mB'),sym('g')]; 

x_1 = sym('te');
x_2 = sym('te1');
x_3 = sym('a');
x_4 = sym('a1');
x = [x_1, x_2, x_3, x_4];
x1 = [x_2, te2, x_4, a2]; 

A = subs(jacobian(x1, x),[te,te1,a,a1,OA,AB,mB,g],x_equi)
B = subs(jacobian(x1, sym('M')),[te,te1,a,a1,OA,AB,mB,g],x_equi)
%C = [1 1 1 1];
%D = 0;

%A = [0 1 0 0; 9821/750 0 -9821/2250 0; 0 0 0 1;0 0 9821/375 0];

%% для dLdte1dt - dLdte==-M
%B = [0;  -8/9; 0; 32/15];

%% для dLdte1dt - dLdte==0
%B = [0;  -8/15; 0; 8/5]; %НЕ управляема!

%sys = ss(A,B,C,D);

%% Pole placement
rank(ctrb(A,B))
%rank(A)

%eig(A)
%pole(sys)

%desired eigs
eigs = [-3; -5; -7; -2];

K = place(A,B,eigs)

eig(A-B*K)

%% discrete system

%Ts = 0.01; % 100 Hz
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








 
