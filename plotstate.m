fileID = fopen('state.txt','r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);

n  = A(1);
ns = A(2);
nx = A(3);
N  = A(4);
d  = A(5);

A = A(6:length(A));

sol = zeros(d, N+1, n+1);

for k=1:n
    for j=1:ns
        for i=1:N+1
            for l=1:d
                sol(l, i, k) = A((k-1)*ns*nx + (j-1)*nx + (i-1)*d + l);
            end
        end
    end
end

A=A(n*ns*nx+1:n*ns*nx+nx);
length(A);
for i=1:N+1
    for l=1:d
        sol(l, i, n+1) = A((i-1)*d + l);
    end
end

for i=1:N+1
    plot(reshape(sol(1, i, :), 1, n+1), reshape(sol(2, i, :), 1, n+1));
    hold all
end