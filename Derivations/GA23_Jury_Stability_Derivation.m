
clc; clear; close all;

syms un u1 dun du1 ddun ddu1 real
syms z xi om delta3 real
syms dt r positive

lambda = -xi + 1i*om;
alpha  = 1/(1 + r);

% ---------------------------------------------------
% Weighted linear combination of GA-2 and GA-3
% ---------------------------------------------------
b2_0 = -(r - 3)/(2*(r + 1));
b2_1 =  (3*r - 1)/(2*(r + 1));

b3_0 =  (2*r^2 - 5*r + 11)/(6*(r + 1));
b3_1 = -(2*r^2 - 11*r + 5)/(6*(r + 1));
b3_2 = -(r^2 - r + 1)/(3*(r + 1));

b0 = (delta3 * b3_0 + (1 - delta3) * b2_0);
b1 = (delta3 * b3_1 + (1 - delta3) * b2_1);
b2 = (delta3 * b3_2 + (1 - delta3) * 0);

% ---------------------------------------------------
% Amplification matrix and characteristic polynomial
% ---------------------------------------------------
A23 = amplification_GA23(b0, b1, b2, alpha, r, dt, lambda, un, u1, dun, du1, ddun, ddu1);

p = charpoly(A23, z);

[pVec, ~] = coeffs(p, z);

a3 = pVec(1);
a2 = pVec(2);
a1 = pVec(3);
a0 = pVec(4);

% ---------------------------------------------------
% Modified Jury table coefficients and conditions
% ---------------------------------------------------
k4 = a0 / conj(a3);

B2 = a3 - conj(a0)*k4;
B1 = a2 - conj(a1)*k4;
B0 = a1 - conj(a2)*k4;

k3 = B0 / conj(B2);

C1 = B2 - conj(B0)*k3;
C0 = B1 - conj(B1)*k3;

k2 = C0 / conj(C1);

D0 = C1 - conj(C0)*k2;

con1 = local_real_rat(a3);
con2 = local_real_rat(B2);
con3 = local_real_rat(C1);
con4 = local_real_rat(D0);

[Con2, ~, ~] = solve(con2 == 0, delta3, 'ReturnConditions', true);
Con2 = simplify(subs(Con2, dt, 0));

C3_dt0 = subs(con3, dt, 0);
[Con3, ~, ~] = solve(C3_dt0 == 0, delta3, 'ReturnConditions', true);
Con3 = simplify(Con3);

C4_sub = subs(con4, [om, xi], [1, 0]);
[Con4, ~, ~] = solve(C4_sub == 0, delta3, 'ReturnConditions', true);
Con4 = simplify(expand(subs(Con4, dt, 0)));

fprintf('Condition (1) is always satisfied\n');
fprintf('Condition (2) is satisfied if: delta_3 >= %s\n', char(Con2));
fprintf('Condition (3) is satisfied if: delta_3 >= %s\n', char(Con3));
fprintf('Condition (4) is satisfied if: delta_3 >= %s\n', char(Con4(2)));

% ============================================================
%  Local functions
% ============================================================

% ------
function A = amplification_GA23(b0, b1, b2, alpha, r, dt, lambda, un, u1, dun, du1, ddun, ddu1)

    g1 = 1 + r;
    g2 = 1 + r;

    eqs = [ ...
        0    == (b0*du1) + (b1*dun) + (b2*ddun*dt) - lambda*(alpha*u1 + (1-alpha)*un), ...
        du1  == g1*(u1 - un)/dt + (1 - g1)*dun, ...
        ddu1 == g2*(du1 - dun)/dt + (1 - g2)*ddun ];

    sol = solve(eqs, [u1, du1, ddu1], 'ReturnConditions', true);

    scale = [1, dt, dt^2];

    var1 = [sol.u1, sol.du1, sol.ddu1];
    varn = [un, dun, ddun];

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