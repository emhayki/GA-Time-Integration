
clc; clear; close all;

syms r dt lambda lam real
syms u1 un du1 dun ddu1 ddun dddu1 dddun real
syms b0 b1 b2 b3 real

taylorOrd = 10;
a = 1/(1 + r);

fprintf('GA2:\n');
[~, coeff2, betas2] = deriveGACoefficients(2, taylorOrd, ...
                                           r, a, dt, lambda, lam, ...
                                           u1, un, du1, dun, ddu1, ddun, dddu1, dddun, ...
                                           b0, b1, b2, b3);
printCoefficients(betas2, coeff2);


fprintf('\nGA3:\n');
[~, coeff3, betas3] = deriveGACoefficients(3, taylorOrd, ...
                                           r, a, dt, lambda, lam, ...
                                           u1, un, du1, dun, ddu1, ddun, dddu1, dddun, ...
                                           b0, b1, b2, b3);
printCoefficients(betas3, coeff3);


fprintf('\nGA4:\n');
[~, coeff4, betas4] = deriveGACoefficients(4, taylorOrd, ...
                                           r, a, dt, lambda, lam, ...
                                           u1, un, du1, dun, ddu1, ddun, dddu1, dddun, ...
                                           b0, b1, b2, b3);
printCoefficients(betas4, coeff4);



% ============================================================
%  Local functions
% ============================================================

% ------
function [A, coeff, betas] = deriveGACoefficients(Method, taylorOrd, ...
    r,a,dt,lambda,lam, u1,un,du1,dun,ddu1,ddun,dddu1,dddun, b0,b1,b2,b3)

    switch Method

        case 2
            g1 = 1 + r;
            coeff = [b0, b1];

            eqs = [   0 == (b0 * du1) + (b1 * dun) - lambda*(a * u1 + (1 - a) * un), ...
                    du1 == g1 * (u1 - un)/dt + (1 - g1) * dun ];

            sol  = solve(eqs, u1, du1);
            fact = [1, dt];

            A = amplification([sol.u1, sol.du1], [un, dun], fact);

        case 3
            g1 = 1 + r;
            g2 = 1 + r;
            coeff = [b0, b1, b2];

            eqs = [    0 == (b0 * du1) + (b1 * dun) + (b2 * ddun * dt) - lambda*(a * u1 + (1 - a) * un), ...
                     du1 == g1 * ( u1 - un)/dt  + (1 - g1) * dun, ...
                    ddu1 == g2 * (du1 - dun)/dt + (1 - g2) * ddun ];

            sol  = solve(eqs, u1, du1, ddu1);
            fact = [1, dt, dt^2];

            A = amplification([sol.u1, sol.du1, sol.ddu1], [un, dun, ddun], fact);

        case 4
            g1 = 1 + r;
            g2 = 1 + r;
            g3 = 1 + r;
            coeff = [b0, b1, b2, b3];

            eqs = [     0  == (b0 * du1) + (b1 * dun) + (b2 * ddun * dt) + (b3 * dddun * dt^2) - lambda * (a * u1 + (1 - a) * un), ...
                      du1  == g1 * (u1 - un)/dt + (1 - g1) * dun, ...
                     ddu1  == g2 * (du1 - dun)/dt + (1 - g2) * ddun, ...
                    dddu1  == g3 * (ddu1 - ddun)/dt + (1 - g3) * dddun ];

            sol  = solve(eqs, u1, du1, ddu1, dddu1);
            fact = [1, dt, dt^2, dt^3];

            A = amplification([sol.u1, sol.du1, sol.ddu1, sol.dddu1], [un, dun, ddun, dddun], fact);

        otherwise
            error('Method must be 2, 3, or 4.');
    end

    order = length(coeff);

    lamex = simplify(exp(lambda * dt));
    p     = charpoly(A, lam);
    err   = simplify(subs(p, lam, lamex));
    err   = taylor(err, dt, 'Order', taylorOrd);

    [co, ~] = coeffs(err, dt);

    equ = sym(zeros(1, order));
    for ii = 1:order
        equ(ii) = co(end - (ii - 1));
    end

    solb = solve(equ, coeff, 'ReturnConditions', false);

    betas = sym(zeros(1, order));
    for k = 1:order
        nm = char(coeff(k));
        betas(k) = simplify(solb.(nm));
    end
end

% ------
function printCoefficients(betas, coeff)
    for k = 1:length(betas)
        fprintf('  %s = %s\n', char(coeff(k)), char(betas(k)));
    end
end

% ------
function A = amplification(xnp1_expr_row, xn_row, fact_row)

    xnp1 = sym(xnp1_expr_row(:));
    xn   = sym(xn_row(:));
    fact = sym(fact_row(:));

    m = length(xnp1);
    n = length(xn);

    if length(fact) ~= m || length(fact) ~= n
        error('amplification: fact length must match state dimension.');
    end

    Xn = sym('Xn', [n,1], 'real');
    Xnp1_expr = fact .* subs(xnp1, xn, Xn./fact);

    A = sym(zeros(m,n));
    for i = 1:m
        for j = 1:n
            A(i,j) = coeffLinear(Xnp1_expr(i), Xn(j));
        end
    end

    A = simplify(A);
end

% ------
function c = coeffLinear(expr, var)
    [C,T] = coeffs(expr, var);
    c = sym(0);
    for k = 1:numel(T)
        if isequal(T(k), var)
            c = c + C(k);
        end
    end
end
