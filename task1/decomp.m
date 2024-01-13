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

%определение параметра сигма
K = zeros(N - 1);
K(1, 1) = K2(1);
K(1, 2) = K3(1);
for i = 2:(N - 2)
    K(i, i - 1) = K1(i);
    K(i, i) = K2(i);
    K(i, i + 1) = K3(i);
end
K(N - 1, N - 2) = K1(N - 1);
K(N - 1, N - 1) = K2(N - 1);
evs = eig(diag(diag(K))^(-1)*K);
lmin = evs(1);
lmax = evs(N - 1);
sigma = 2/(lmax + lmin);

tic
while (err > 10^(-7))
    u_k1(1) = u_k(1) - sigma/K2(1)*(K2(1)*u_k(1) + K3(1)*u_k(2) - F(1));
    for i = 2:(N - 2)
        u_k1(i) = u_k(i) - sigma/K2(i)*(K1(i)*u_k(i - 1) + K2(i)*u_k(i) + K3(i)*u_k(i + 1) - F(i));
    end
    u_k1(N - 1) = u_k(N - 1) - sigma/K2(N - 1)*(K1(N - 1)*u_k(N - 2) + K2(N - 1)*u_k(N - 1) - F(N - 1));
    err = norm(u_k - u_k1);
    u_k = u_k1;
end
Y2 = zeros(1, N + 1);
Y2(2:N) = u_k;
toc

n = sqrt(N);
N = n * n;
h = 1 / N;
u_k = zeros(N + 1, 1);
u_k1 = zeros(N + 1, 1);
err = inf;
m = 0;
while (err > 10^(-7))
    m = m + 1;
    for l = 0:n:(N - 2*n)
        %tic
        d = zeros(2*n - 1, 1);
        if l == 0
            d(1) = sigma*(K2(l + 1)*u_k(l + 1) + K3(l + 1)*u_k(l + 2) - F(l + 1));
            for i = 2:(2*n - 1)
                d(i) = sigma*(K1(l + i)*u_k(l + i - 1) + K2(l + i)*u_k(l + i) + K3(l + i)*u_k(l + i + 1) - F(l + i));
            end
        elseif l == (N - 2*n)
            for i = 1:(2*n - 2)
                d(i) = sigma*(K1(l + i)*u_k(l + i - 1) + K2(l + i)*u_k(l + i) + K3(l + i)*u_k(l + i + 1) - F(l + i));
            end
            d(2*n - 1) = sigma*(K1(l + 2*n - 1)*u_k(l + 2*n - 2) + K2(l + 2*n - 1)*u_k(l + 2*n - 1) - F(l + 2*n - 1));
        else
            for i = 1:(2*n - 1)
                d(i) = sigma*(K1(l + i)*u_k(l + i - 1) + K2(l + i)*u_k(l + i) + K3(l + i)*u_k(l + i + 1) - F(l + i));
            end
        end
        S = zeros(1, 2 * n);
        T = zeros(1, 2 * n);
        S(1) = alpha2/(h*alpha1 + alpha2);
        T(1) = alpha*h/(h*alpha1 + alpha2);
        for i = 1:(2*n - 1)
            S(i + 1) = -K3(l + i)/(K2(l + i) + K1(l + i)*S(i));
            T(i + 1) = -(K1(l + i)*T(i) - d(i))/(K2(l + i) + K1(l + i)*S(i));
        end
        W = zeros(2*n + 1, 1);
        W(2*n + 1) = (beta2*T(2*n) + h*beta)/(h*beta1 + beta2 - beta2*S(2*n));
        for j = fliplr(1:(2*n))
            W(j) = S(j)*W(j + 1) + T(j);
        end
        u_k1((l + 1):(l + 2*n - 1)) = u_k1((l + 1):(l + 2*n - 1)) - W(2:2*n);
        %toc
    end
    err = norm(u_k - u_k1);
    u_k = u_k1;
end
plot(x_i, u_k1)