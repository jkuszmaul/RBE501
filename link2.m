syms a b c d e;
syms t1 t2 t1dot t2dot;

% tau = M * tddot + C * tdot + G
% tddot = -M^-1 * (C * tdot + G) + M^-1 * tau
% f(x) = -M^-1 * (C * tdot + G)
% g(x) = M^-1
M = [a + b + 2*c*cos(t2)  b + c * cos(t2);
     b + c * cos(t2)      b];
C = [-c * sin(t2) * t2dot -c * sin(t2) * (t1dot + t2dot);
     c * sin(t2) * t1dot  0];
G = [-d * sin(t1) - e * sin(t1 + t2);
     e * sin(t1 + t2)];
