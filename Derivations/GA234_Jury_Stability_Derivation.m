
clc; clear; close all;

syms un u1 dun du1 ddun ddu1 dddun dddu1 real
syms z xi om delta4 real
syms dt r positive

lambda = -xi + 1i*om;
alpha  = 1/(1 + r);

% ---------------------------------------------------
% Weighted linear combination of GA-2, GA-3 and GA-4
% ---------------------------------------------------

delta3 = ((r - 1)^2)/(2 * (r^2 - r + 1));

b2_0 = -(r - 3)/(2*(r + 1));
b2_1 =  (3*r - 1)/(2*(r + 1));

b3_0 =  (2*r^2 - 5*r + 11)/(6*(r + 1));
b3_1 = -(2*r^2 - 11*r + 5)/(6*(r + 1));
b3_2 = -(r^2 - r + 1)/(3*(r + 1));

b0_23 = (delta3 * b3_0 + (1 - delta3) * b2_0);
b1_23 = (delta3 * b3_1 + (1 - delta3) * b2_1);
b2_23 = (delta3 * b3_2);

b4_0 = (-3*r^3 + 7*r^2 - 13*r + 25) / (12*(r + 1));
b4_1 = ( 3*r^3 - 7*r^2 + 25*r - 13) / (12*(r + 1));
b4_2 = ( 3*r^3 - 7*r^2 +  7*r -  7) / (12*(r + 1));
b4_3 = ( r^3 - r^2 + r - 1) / (4*(r + 1)^2);

b0 = delta4*b4_0 + (1 - delta4)*b0_23;
b1 = delta4*b4_1 + (1 - delta4)*b1_23;
b2 = delta4*b4_2 + (1 - delta4)*b2_23;
b3 = delta4*b4_3;

% ---------------------------------------------------
% Amplification matrix and characteristic polynomial
% ---------------------------------------------------
A234 = amplification_GA234(b0, b1, b2, b3, alpha, r, dt, lambda, ...
                           un, u1, dun, du1, ddun, ddu1, dddun, dddu1);
p = charpoly(A234, z);

[pVec, ~] = coeffs(p, z);

a4 = pVec(1);
a3 = pVec(2);
a2 = pVec(3);
a1 = pVec(4);
a0 = pVec(5);

% ---------------------------------------------------
% Modified Jury table coefficients and conditions
% ---------------------------------------------------
k5 = a0 / conj(a4);

B3 = a4 - conj(a0)*k5;
B2 = a3 - conj(a1)*k5;
B1 = a2 - conj(a2)*k5;
B0 = a1 - conj(a3)*k5;

k4 = B0 / conj(B3);

C2 = B3 - conj(B0)*k4;
C1 = B2 - conj(B1)*k4;
C0 = B1 - conj(B2)*k4;

k3 = C0 / conj(C2);

D1 = C2 - conj(C0)*k3;
D0 = C1 - conj(C1)*k3;

k2 = D0 / conj(D1);

E0 = D1 - conj(D0)*k2;

con1 = local_real_rat(a4);
con2 = local_real_rat(B3);
con3 = local_real_rat(C2);
con4 = local_real_rat(D1);
con5 = local_real_rat(E0);

[Con2, ~, ~] = solve(con2 == 0, delta4, 'ReturnConditions', true);
Con2 = simplify(subs(Con2, dt, 0));

C3_dt0 = subs(con3, dt, 0);
[Con3, ~, ~] = solve(C3_dt0 == 0, delta4, 'ReturnConditions', true);
Con3 = simplify(Con3);

C4_dt0 = subs(con4, dt, 0);
[Con4, ~, ~] = solve(C4_dt0 == 0, delta4, 'ReturnConditions', true);
Con4 = simplify(Con4);

C5_sub = subs(con5, [om, xi], [1, 0]);
[Con5, ~, ~] = solve(C5_sub == 0, delta4, 'ReturnConditions', true);
Con5 = simplify(expand(subs(Con5, dt, 0)));

fprintf('From Condition (5) delta_4 = %s\n', char(Con5(2)));

% ============================================================
%  Local functions
% ============================================================

% ------
function A = amplification_GA234(b0, b1, b2, b3, alpha, r, dt, lambda, un, u1, dun, du1, ddun, ddu1, dddun, dddu1)

    g1 = 1 + r;
    g2 = 1 + r;
    g3 = 1 + r;

    eqs = [ ...
        0     == (b0*du1) + (b1*dun) + (b2*ddun*dt) + (b3*dddun*dt^2) - lambda*(alpha*u1 + (1-alpha)*un), ...
        du1   == g1*(u1 - un)/dt + (1 - g1)*dun, ...
        ddu1  == g2*(du1 - dun)/dt + (1 - g2)*ddun, ...
        dddu1 == g3*(ddu1 - ddun)/dt + (1 - g3)*dddun ];

    sol = solve(eqs, [u1, du1, ddu1, dddu1], 'ReturnConditions', true);

    scale = [1, dt, dt^2, dt^3];

    var1 = [sol.u1, sol.du1, sol.ddu1, sol.dddu1];
    varn = [un, dun, ddun, dddun];

    n = numel(var1);
    A = sym(zeros(n,n));

    for i = 1:n
        [c, v] = coeffs(var1(i), varn);
        row = sym(zeros(1,n));
        [~, ~, idx] = intersect(v, varn, 'stable');
        row(idx) = c(1:numel(idx));
        A(i,:) = row;
    end

    for i = 1:n
        A(i,:) = A(i,:) * scale(i);
        A(:,i) = A(:,i) / scale(i);
    end

    A = simplify(A);
end

% ------
function out = local_real_rat(expr)
      expr = simplifyFraction(expr);
    [N, D] = numden(expr);
    out = simplifyFraction((N*conj(D) + conj(N)*D) / (2*D*conj(D)));
end