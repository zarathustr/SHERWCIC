function x = vec(X)
m = size(X, 1);
n = size(X, 2);
x = reshape(X, [m * n, 1]);
end