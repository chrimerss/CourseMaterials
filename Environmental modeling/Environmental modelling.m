c0 = 1;
k=20;
u=100;
A=10;
Q=A*u;
D=500;
delx=0.1;
n=500;

up = -ones(n-1,1)*D*A/delx;
low = -ones(n-1,1)*(D*A/delx+Q);
main = ones(n,1)*(2*D*A/delx+Q+k*delx*A);
main(1)=D*A/delx+k*delx*A+Q;
main(end)=D*A/delx+k*delx*A+Q;
matrix = diag(main,0)+diag(up,1)+diag(low,-1);
rhs = zeros(n,1);
rhs(1) = c0*Q;
c = inv(matrix)*rhs;
plot(c)
xlim([0,20])
hold on