clear all
clc

syms r m t a a1 a2 g M

I = m*r*r; 
y = r*cos(a);
vx = r*cos(a)*a1;
vy = -r*sin(a)*a1;

Ek = m*(vx*vx + vy*vy)/2 + I*a1*a1/2;
Ep = m*g*y;
L = Ek-Ep;

dLda = diff(L,sym('a'));
dLda1 = subs(diff(L,sym('a1')),[sym('a1'),sym('a')],[a2*t,a1*t]);
dLda1dt = subs(diff(dLda1,sym('t')),[a2*t,a1*t],[sym('a1'),sym('a')]);

a2 = subs(simplify(solve(dLda1dt - dLda==M,a2)),[a2*t,a1*t],[sym('a1'),sym('a')])
