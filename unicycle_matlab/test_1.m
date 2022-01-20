clear
clc

syms OA AB mA mB t a2 te2 g a1 te1 a te M

%OA = 1;
%AB = 1;
%mA = 0.5;
%mB = 2;
%g = 10;

%a1 = a2*t;
%te1 = te2*t;
a = a1*t;
te = te1*t;

%xA = OA*sin(te);
%yA = OA*cos(te);
%xB = xA + AB*sin(te-a);
%yB = yA + AB*cos(te-a);

%vBx = diff(xB,t);
%vBy = diff(yB,t);
%vAx = diff(xA,t);
%vAy = diff(yA,t);

%Ek = expand(mB*(vBx*vBx + vBy*vBy)/2); % + mA*(vAx*vAx + vAy*vAy)/2;
%Ep = expand(mB*g*yB); % + mA*g*yA;
%L = simplify(subs(Ek - Ep,[a1*t,t*te1],[sym('a'),sym('te')]))
%old L = (mB*(AB^2*a1^2 - 2*AB^2*a1*te1 + AB^2*te1^2 - 2*cos(a)*AB*OA*a1*te1 + 2*cos(a)*AB*OA*te1^2 - 2*g*cos(a - te)*AB + OA^2*te1^2 - 2*g*cos(te)*OA))/2;

%dLdte = diff(L,sym('te'));
%dLdte1 = subs(diff(L,sym('te1')),[sym('a1'),sym('te1'),sym('a'),sym('te')],[a2*t,t*te2,a1*t,t*te1]);
%dLdte1dt = subs(diff(dLdte1,sym('t')),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')]);

%dLda = diff(L,sym('a'));
%dLda1 = subs(diff(L,sym('a1')),[sym('a1'),sym('te1'),sym('a'),sym('te')],[a2*t,t*te2,a1*t,t*te1]);
%dLda1dt = subs(diff(dLda1,sym('t')),[a2*t,t*te2,a1*t,t*te1],[sym('a1'),sym('te1'),sym('a'),sym('te')]);

%te2_1 = subs(simplify(solve(dLdte1dt - dLdte==0,te2)),[t*(2*a1 - 2*te1),t*(a1 - te1),a2*t,t*te2,a1*t,t*te1],[2*(sym('a')-sym('te')),sym('a')-sym('te'),sym('a1'),sym('te1'),sym('a'),sym('te')]);
%te2_2 = subs(simplify(solve(dLda1dt - dLda==M,te2)),[t*(2*a1 - 2*te1),t*(a1 - te1),a2*t,t*te2,a1*t,t*te1],[2*(sym('a')-sym('te')),sym('a')-sym('te'),sym('a1'),sym('te1'),sym('a'),sym('te')]);
%a2 = simplify(solve(te2_1==te2_2,a2));
%te2 = simplify(subs(te2_1,sym('a2'),a2));

%x_equi = [0,0,0,0,0];

%a2 = subs(a2,[M,te,te1,a,a1],x_equi)
%te2 = subs(te2,[M,te,te1,a,a1],x_equi)

%a = sym('a');
%te = sym('te');
%a1 = sym('a1');
%te1 = sym('te1');

%x_1 = sym('te');
%x_2 = sym('te1');
%x_3 = sym('a');
%x_4 = sym('a1');
%x = [x_1, x_2, x_3, x_4];
%x1 = [x_2, te2, x_4, a2]; 

%J = jacobian(x1, x)

%A = subs(J,[te,te1,a,a1],x_equi)
 

 

 
