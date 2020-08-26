% INPUT: 3 diagonals (a,b,c) of a tridiagonal matrix and
% its constant vector d.
% OUTPUT: vector v which solves the system

function [x] = ThomasAlgorithm(a, b, c, d)
    N = length(b) - 1;
    for n=1:N
        b(n+1) = b(n+1) - a(n+1)/b(n)*c(n);
        d(n+1) = d(n+1) - a(n+1)/b(n)*d(n);
    end

    x(N+1) = d(N+1)/b(N+1);
    for n = N:-1:1
        x(n) = ( d(n) - c(n)*x(n+1) ) / b(n);
    end
end
