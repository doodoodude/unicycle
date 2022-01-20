clear all
clc

syms OA AB mA mB t a2 te2 g a1 te1 a te M Iao Iab

%% Model
%a1 = a2*t;
%te1 = te2*t;
%a = a1*t;
%te = te1*t;

Ioa = mA*OA*OA/3;
Iab = mB*0.0765; %?????

%yA = OA*cos(te)/2; 
%yB = OA*cos(te) + AB*cos(te+a); 

%vAx = te1*OA*cos(te)/2;
%vAy = -te1*OA*sin(te)/2;
%vBx = te1*OA*cos(te) + (te1 + a1)*AB*cos(te+a);
%vBy = -te1*OA*sin(te) - (te1 + a1)*AB*sin(te+a);

%Ek = expand((mB*(vBx*vBx + vBy*vBy)/2) + Iab*a1*a1/2 + (mA*(vAx*vAx + vAy*vAy)/2) + Ioa*te1*te1/2);
%Ep = expand(mB*g*yB + mA*g*yA);
%L = simplify(subs(Ek - Ep,[a1*t,t*te1],[sym('a'),sym('te')]))

L = (153*a1^2*mB)/4000 + (AB^2*a1^2*mB)/2 + (AB^2*mB*te1^2)/2 + (7*OA^2*mA*te1^2)/24 + (OA^2*mB*te1^2)/2 - AB*g*mB*cos(a + te) + AB^2*a1*mB*te1 - (OA*g*mA*cos(te))/2 - OA*g*mB*cos(te) + AB*OA*mB*te1^2*cos(a) + AB*OA*a1*mB*te1*cos(a);

%dLdte = diff(L,sym('te'));
dLdte= AB*g*mB*sin(a + te) + (OA*g*mA*sin(te))/2 + OA*g*mB*sin(te);
%dLdte1 = subs(diff(L,sym('te1')),[sym('a1'),sym('te1'),sym('a'),sym('te')],[a2*t,t*te2,a1*t,t*te1]);
%dLdte1dt = subs(diff(dLdte1,sym('t')),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')]);
dLdte1dt = AB^2*a2*mB + AB^2*mB*te2 + (7*OA^2*mA*te2)/12 + OA^2*mB*te2 + AB*OA*a2*mB*cos(a) + 2*AB*OA*mB*te2*cos(a) - AB*OA*a1^2*mB*sin(a) - 2*AB*OA*a1*mB*te1*sin(a);
  
%dLda = diff(L,sym('a'));
dLda = - AB*OA*mB*sin(a)*te1^2 - AB*OA*a1*mB*sin(a)*te1 + AB*g*mB*sin(a + te);
%dLda1 = subs(diff(L,sym('a1')),[sym('a1'),sym('te1'),sym('a'),sym('te')],[a2*t,t*te2,a1*t,t*te1]);
%dLda1dt = subs(diff(dLda1,sym('t')),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')]);
dLda1dt = (153*a2*mB)/2000 + AB^2*a2*mB + AB^2*mB*te2 + AB*OA*mB*te2*cos(a) - AB*OA*a1*mB*te1*sin(a);

%te2_1 = subs(simplify(solve(dLdte1dt - dLdte==0,te2)),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')]);
%te2_2 = subs(simplify(solve(dLda1dt - dLda==M,te2)),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')]);
%a2 = simplify(solve(te2_1==te2_2,a2))
%te2 = simplify(subs(te2_1,sym('a2'),a2))

%% при dLdte1dt - dLdte==-M
%a2 = -((12*(AB*OA*mB*sin(a)*a1^2 + 2*AB*OA*mB*te1*sin(a)*a1 - M + AB*g*mB*sin(a + te) + (OA*g*mA*sin(te))/2 + OA*g*mB*sin(te)))/(12*AB^2*mB + 7*OA^2*mA + 12*OA^2*mB + 24*AB*OA*mB*cos(a)) - (- AB*OA*mB*sin(a)*te1^2 + M + AB*g*mB*sin(a + te))/(AB*mB*(AB + OA*cos(a))))/((2000*AB^2 + 153)/(2000*AB*(AB + OA*cos(a))) - (12*AB*mB*(AB + OA*cos(a)))/(12*AB^2*mB + 7*OA^2*mA + 12*OA^2*mB + 24*AB*OA*mB*cos(a)));
%te2 = (12*(306*AB*g*mB*sin(a + te) - 8000*AB^2*M - 4000*AB*M*OA*cos(a) - 306*M + 153*OA*g*mA*sin(te) + 306*OA*g*mB*sin(te) + 2000*AB^2*OA*g*mA*sin(te) + 2000*AB^2*OA*g*mB*sin(te) + 2000*AB^2*OA^2*mB*te1^2*sin(2*a) - 2000*AB^2*OA*g*mB*sin(2*a + te) + 4000*AB^3*OA*a1^2*mB*sin(a) + 4000*AB^3*OA*mB*te1^2*sin(a) + 306*AB*OA*a1^2*mB*sin(a) + 612*AB*OA*a1*mB*te1*sin(a) + 8000*AB^3*OA*a1*mB*te1*sin(a)))/(3672*AB^2*mB + 2142*OA^2*mA + 3672*OA^2*mB + 28000*AB^2*OA^2*mA + 48000*AB^2*OA^2*mB + 7344*AB*OA*mB*cos(a) - 48000*AB^2*OA^2*mB*cos(a)^2);

%% при dLdte1dt - dLdte==0
a2 = -((12*AB*OA*mB*sin(a)*a1^2 + 24*AB*OA*mB*te1*sin(a)*a1 + 12*AB*g*mB*sin(a + te) + 6*OA*g*mA*sin(te) + 12*OA*g*mB*sin(te))/(12*AB^2*mB + 7*OA^2*mA + 12*OA^2*mB + 24*AB*OA*mB*cos(a)) - (- AB*OA*mB*sin(a)*te1^2 + M + AB*g*mB*sin(a + te))/(AB*mB*(AB + OA*cos(a))))/((2000*AB^2 + 153)/(2000*AB*(AB + OA*cos(a))) - (12*AB*mB*(AB + OA*cos(a)))/(12*AB^2*mB + 7*OA^2*mA + 12*OA^2*mB + 24*AB*OA*mB*cos(a)));
te2 = (1836*AB*g*mB*sin(a + te) - 24000*AB*M*OA*cos(a) - 24000*AB^2*M + 918*OA*g*mA*sin(te) + 1836*OA*g*mB*sin(te) + 12000*AB^2*OA*g*mA*sin(te) + 12000*AB^2*OA*g*mB*sin(te) + 12000*AB^2*OA^2*mB*te1^2*sin(2*a) - 12000*AB^2*OA*g*mB*sin(2*a + te) + 24000*AB^3*OA*a1^2*mB*sin(a) + 24000*AB^3*OA*mB*te1^2*sin(a) + 1836*AB*OA*a1^2*mB*sin(a) + 3672*AB*OA*a1*mB*te1*sin(a) + 48000*AB^3*OA*a1*mB*te1*sin(a))/(1836*AB^2*mB + 1071*OA^2*mA + 1836*OA^2*mB + 14000*AB^2*OA^2*mA + 24000*AB^2*OA^2*mB + 3672*AB*OA*mB*cos(a) - 24000*AB^2*OA^2*mB*cos(a)^2);

%% проверка
%a2 = subs(a2,[OA,AB,mB,mA,g],[0.5,0.25,10,1,9.821])
%te2 = subs(te2,[OA,AB,mB,g],[0.5,0.25,10,1,9.821])

%a2 = subs(a2,[M,te,te1,a,a1],[M,0,0,0,0])
%te2 = subs(te2,[M,te,te1,a,a1],[M,0,0,0,0])

%% State space implementaion
x_equi = [0,0,0,0,0.5,0.25,10,0,9.821];  %mA = 0

x_1 = sym('te');
x_2 = sym('te1');
x_3 = sym('a');
x_4 = sym('a1');
x = [x_1, x_2, x_3, x_4];
x1 = [x_2, te2, x_4, a2]; 

%A = subs(jacobian(x1, x),[te,te1,a,a1,OA,AB,mB,mA,g],x_equi)
%B = subs(jacobian(x1, sym('M')),[te,te1,a,a1,OA,AB,mB,mA,g],x_equi)
%C = [1 0 0 0; 0 0 0 0; 0 1 0 0; 0 0 0 0]; %eye(4);
%D = 0;

A = [0 1 0 0; 35856471/2703500 0 -2857911/1081400 0; 0 0 0 1;-9821/43256 0 1836527/86512 0];

%% для dLdte1dt - dLdte==-M
%B = [0; -3918/5407; 0; 9175/5407]; 

%% для dLdte1dt - dLdte==0
B = [0; -2250/5407; 0; 6925/5407];

%sys = ss(A,B,C,D)

%% Pole placement
eig(A)
%pole(sys)

%desired eigs
eigs = [-1.5; -2; -1+0.5i; -1-0.5i];

K = place(A,B,eigs)

eig(A-B*K)














 
