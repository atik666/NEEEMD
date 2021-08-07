function normalized = normFun(X)
X_min = min(X);
X_max = max(X);

A = X - X_min;
B= X_max - X_min;

for i= 1:size(A,2)
    normalized(:,i) = A(:,i)./B(:,i);
end
end