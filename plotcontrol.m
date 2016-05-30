fileID = fopen('control.txt','r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);

n  = A(1);
ns = A(2);
nc = A(3);
N  = A(4);
d  = A(5);
nMPC = A(6);
t0 = A(7);
tf = A(8);


A = A(9:length(A));

U1 = zeros(1, n*ns);
U2 = zeros(1, n*ns);

k = 1;
for i = 1:2:n*ns
    U1(k) = A(i);
    U2(k) = A(i+1);
    k = k+1;
end

h = (tf-t0)/(n*ns-1);
T = t0:h:tf;

%figure
plot(T, U1);
hold all
plot(T, U2);