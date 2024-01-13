clear
clc

N = 100;
h = 2 / N;
x_i = linspace(-1, 1, N+1);

alpha1 = 0;
alpha2 = 1;
beta1 = 1;
beta2 = 0;

alpha = 0;
beta = 0;

K1 = zeros(1, N - 1);
K2 = zeros(1, N - 1);
K3 = zeros(1, N - 1);
F = zeros(1, N - 1);
for i = 1 : (N - 1)
    x_left = x_i(i);
    x_right = x_i(i+1);

    a_left = (x_left+4)/(x_left+5);
    a_right = (x_right+4)/(x_right+5);

    c_left = exp(x_left/4);
    c_right = exp(x_right/4);

    f_left = 2-x_left;
    f_right = 2-x_right;

    K1(i) = -a_left/h + c_left*h/6;
    K2(i) = (a_left + a_right)/h + (c_left + c_right)*h/3;
    K3(i) = -a_right/h + c_right*h/6;

    F(i) = (f_left + f_right)*h/2;
end

tic
S = zeros(1, N);
T = zeros(1, N);

S(1) = 0;
T(1) = alpha;

for i = 1 : (N - 1)
    S(i + 1) = -K3(i)/(K2(i) + K1(i)*S(i));
    T(i + 1) = -(K1(i)*T(i) - F(i))/(K2(i) + K1(i)*S(i));
end

Y1 = zeros(1, N + 1);
Y1(N + 1) = beta;

for i = N:-1:2
    Y1(i) = S(i)*Y1(i + 1) + T(i);
end
toc

plot(x_i, Y1)
hold on
grid on

x = 0.1:0.1:0.9;
x = x';
u = zeros(9, 1);
k1 = zeros(9, 1);
k2 = zeros(9, 1);
k3 = zeros(9, 1);
f = zeros(9, 1);
for i = 1:9
    u(i) = Y1(i*N/10);
    k1(i) = K1(i*N/10);
    k2(i) = K3(i*N/10);
    k3(i) = K3(i*N/10);
    f(i) = F(i*N/10);
end
T = table(x, u, k1, k2, k3, f)