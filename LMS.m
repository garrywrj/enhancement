N=1000; %�źŵ�ȡ������
M=10;%��ͷ�ĸ���
iter=30;%��������
n=1:N;
s=5*sin(2*pi*0.5*n/100);%�����ź�
x=s+0.5*randn(1,N);%�����ź�
d=s;
%��ʼ��
w=zeros(1,M);%��ͷϵ��
X_A=zeros(1,M);%�����źŵľ����Ա����Y����
E=zeros(iter,N);%������
e=zeros(1,N);%�������
u=0.001;%�̶�����
y=zeros(1,N);
E=zeros(1,N);
%�����Ŀ�ʼ
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
