clear all
clc

syms OA AB mA mB t a2 te2 g a1 te1 a te M

%% Model
%a1 = a2*t;
%te1 = te2*t;
%a = a1*t;
%te = te1*t;

%yB = OA*cos(te) + AB*cos(te-a); 

%vBx = te1*OA*cos(te) - (te1 - a1)*AB*cos(te-a);
%vBy = -te1*OA*sin(te) - (te1 - a1)*AB*sin(te-a);

%Ek = expand((mB*(vBx*vBx + vBy*vBy)/2) + mB*AB*AB*a1*a1/2);
%Ep = expand(mB*g*yB);
%L = simplify(subs(Ek - Ep,[a1*t,t*te1],[sym('a'),sym('te')]));

%dLdte = diff(L,sym('te'));
%dLdte1 = subs(diff(L,sym('te1')),[sym('a1'),sym('te1'),sym('a'),sym('te')],[a2*t,t*te2,a1*t,t*te1]);
%dLdte1dt = subs(diff(dLdte1,sym('t')),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')]);
 
%dLda = diff(L,sym('a'));
%dLda1 = subs(diff(L,sym('a1')),[sym('a1'),sym('te1'),sym('a'),sym('te')],[a2*t,t*te2,a1*t,t*te1]);
%dLda1dt = subs(diff(dLda1,sym('t')),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')]);

%te2_1 = subs(simplify(solve(dLdte1dt - dLdte==0,te2)),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')]);
%te2_2 = subs(simplify(solve(dLda1dt - dLda==M,te2)),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')]);
%a2 = simplify(solve(te2_1==te2_2,a2));
%te2 = simplify(subs(te2_1,sym('a2'),a2));

a2 = ((AB*OA*sin(a - 2*te)*a1^2 - 2*AB*OA*sin(a - 2*te)*a1*te1 + 2*AB*OA*sin(a - 2*te)*te1^2 - AB*g*sin(a - te) + OA*g*sin(te))/(AB^2 - 2*cos(a - 2*te)*AB*OA + OA^2) + (- AB*OA*mB*sin(a - 2*te)*te1^2 + M + AB*g*mB*sin(a - te))/(AB*mB*(AB - OA*cos(a - 2*te))))/((2*AB)/(AB - OA*cos(a - 2*te)) - (AB*(AB - OA*cos(a - 2*te)))/(AB^2 - 2*cos(a - 2*te)*AB*OA + OA^2));
te2 = -(2*AB*M - 2*M*OA*cos(a - 2*te) - 2*AB^2*g*mB*sin(a - te) - AB*OA*g*mB*sin(2*a - 3*te) + 4*AB^2*OA*a1^2*mB*sin(a - 2*te) + 6*AB^2*OA*mB*te1^2*sin(a - 2*te) + 3*AB*OA*g*mB*sin(te) + AB*OA^2*mB*te1^2*sin(2*a - 4*te) - 8*AB^2*OA*a1*mB*te1*sin(a - 2*te))/(AB*mB*(OA^2*cos(2*a - 4*te) - 2*AB^2 - 3*OA^2 + 4*AB*OA*cos(a - 2*te))); 

%% State space implementaion
x_equi = [0,0,0,0,0.5,0.25,10,9.821];

%a2 = subs(a2,[OA,AB,mB,g],[0.5,0.25,15,9.821])
%te2 = subs(te2,[OA,AB,mB,g],[0.5,0.25,15,9.821])

%a2 = subs(a2,[M,te,te1,a,a1],[M,1,1,1,1])
%te2 = subs(te2,[M,te,te1,a,a1],[M,1,1,1,1])

x_1 = sym('te');
x_2 = sym('te1');
x_3 = sym('a');
x_4 = sym('a1');
x = [x_1, x_2, x_3, x_4];
x1 = [x_2, te2, x_4, a2]; 

%A = subs(jacobian(x1, x),[te,te1,a,a1,OA,AB,mB,g],x_equi)
A = [0 1 0 0; 68747/250 0 -29463/250 0; 0 0 0 1;-19642/125 0 9821/125 0];
%B = subs(jacobian(x1, sym('M')),[te,te1,a,a1,OA,AB,mB,g],x_equi)
B = [0; -16/15; 0; 16/15];
C = [1 0 0 0; 0 0 0 0; 0 1 0 0; 0 0 0 0]; %eye(4);
D = 0;

%sys = ss(A,B,C,D)

%% Pole placement

eig(A)
%pole(sys)

%desired eigs
eigs = [-1; -1.1; -1.2; -1.3]

K = place(A,B,eigs)

eig(A-B*K)














 
