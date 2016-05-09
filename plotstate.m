fileID = fopen('state.txt','r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);

n  = A(1);
ns = A(2);
nx = A(3);
N  = A(4);
d  = A(5);
nMPC = A(6);


A = A(7:length(A));


sol = zeros(d, N+1, n+1, nMPC);

for z = 1:nMPC
    B = A((z-1)*(n*ns*nx+nx) + 1 : z*(n*ns*nx+nx));
    for k=1:n
        for j=1:ns
            for i=1:N+1
                for l=1:d
                    sol(l, i, k, z) = B( (k-1)*ns*nx + (j-1)*nx + (i-1)*d + l);
                end
            end
        end
    end

    C=B(n*ns*nx+1:n*ns*nx+nx);
    length(C);
    for i=1:N+1
        for l=1:d
            sol(l, i, n+1, z) = C((i-1)*d + l);
        end
    end

    for i=1:nMPC
        for j=1:N+1
            plot(reshape(sol(1, j, :, i), 1, n+1), reshape(sol(2, j, :, i), 1, n+1));
            hold all
        end
    end
end