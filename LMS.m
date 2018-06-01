N=1000; %信号的取样点数
M=10;%抽头的个数
iter=30;%迭代次数
n=1:N;
s=5*sin(2*pi*0.5*n/100);%期望信号
x=s+0.5*randn(1,N);%输入信号
d=s;
%初始化
w=zeros(1,M);%抽头系数
X_A=zeros(1,M);%输入信号的矩阵，以便求出Y矩阵
E=zeros(iter,N);%误差矩阵
e=zeros(1,N);%误差向量
u=0.001;%固定步长
y=zeros(1,N);
E=zeros(1,N);
%迭代的开始
for i=M:M+iter-1
    X_A=x(i:-1:i-M+1);
    y(i)=w*X_A';
    e(i)=d(i)-y(i);
    w=w+2*u*e(i)*X_A;
    E(i)=e(i)*e(i);
end
for i=iter+M:N
    X_A=x(i:-1:i-M+1);
    y(i)=w*X_A';
    e(i)=d(i)-y(i);
    E(i)=e(i)*e(i);
end
subplot(3,1,1)
plot(n,x)
subplot(3,1,2)
plot(n,y)
m=1:iter;
subplot(3,1,3)
plot(n,E)
