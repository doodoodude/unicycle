clear all
clc

syms r tep ter x y z Sr Sp Cr Cp D 

Sr = sin(ter);
Sp = sin(tep);
D = sqrt(1-Sr*Sr-Sp*Sp);
Cr = cos(ter);
Cp = cos(tep);

x = r*Sp;
y = -r*Sr;
z = r*D;

J = [0 r*Cp Sp; -r*Cr 0 -Sr; -r*Cr*Sr/D -r*Cp*Sp/D D];
%J = jacobian([x,y,z],[ter,tep,r]);
%J_1 = subs(J_1,[sin(ter),sin(tep),sqrt(1-sin(ter)*sin(ter)-sin(tep)*sin(tep)),cos(ter),cos(tep)],[sym('Sr'), sym('Sp'), sym('D'), sym('Cr'), sym('Cp')])

