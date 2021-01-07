%% Initial Conditions

%In all the pictures presented in the present section, one uses the
%following values for the Poisson ratio and the parameters critical stress (Sigmac) and internal friction (k) of
%Drucker-Prager criterion:
PoissonRatio = 0.2;
CriticalStress = 0;
k = 0.2;

% On set:
% SigmaC := Rad(D1*E0)
% EpsilonC := Rad(D1/E0)
% Where  the material constant SigmaC fixes the scale of the stresses whereas
% The material constant EpsilonC fixes the scale of the strains.

% Modello paraetrizzato rispetto a E0 e D1, per semplicità si
% considereranno valori unitari
E0 = 1;
D1 = 1;

% The elastic moduli K0 and Mu0 are related to the young modulus E0 and the
% Poisson ratio of the sound material by:
% 3K0 = E0/(1 - 2*PoissonRatio)
% 2Mu0 = E0/(1 + PoissonRatio)
K0 = 1/3*E0/(1-2*PoissonRatio);

%The responses depend on the two exponents m and n, and on the two ratios R1/E0 and p0/Sigmac
m = 2;
n = 0.5;

% with:
% R1 > 0
% m > 1
% 0 < n < 1

% Considering the confining pressure equal to:
p0 = 0;

R1 = 0.1*E0;

alpha = 0:0.005:1;
R = zeros(size(alpha));
S = zeros(size(alpha));

%% Study of a family models
% Let us consider the following family on functions R(alpha):
% R(alpha) = R1 * ((1 - alpha)^ m)/(alpha ^ n)
% with:
% R1 > 0
% m > 1
% 0 < n < 1
% This choise of function is governed mainly by our intention to control
% the asymptotic behavior both at the initiation of damage (alpha~0) and at the
% end of the microcracking (alpha~1)

% When alpha~0 -> power law of n
% When alpha~1 -> power law of m

% Whe can observe that R(alpha) is strictly decreasing from +infinite to 0 when
% alpha grows from 0 to 1

% Because m > 1 > n > 0, R''(alpha) > 0 for every value of alpha.

for i = 1:length(alpha)
    
    R(i) = R1 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative(i) = - R1*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative(i)  = R1*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    
    %The compliance function S(alpha) and its first and second derivatives are
    %given by:
    S(i)= 1/R1* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative(i) = 1/R1 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative(i) = 1/R1*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    
    end

%% Graph R and S

%Function R(alpha)
%Whe can observe that R(alpha) is strictly decreasing from +infinite to 0 when
%alpha grows from 0 to 1
figure(1)
plot(alpha, R, 'g', 'LineWidth', 2)
grid on
grid minor
xlim([0  1])
ylim([0  1])
% Labels and Titles
xlabel('\alpha', 'FontSize', 15)
ylabel('R(\alpha)', 'FontSize', 12)
hold on
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 1)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 1)
hold on

%Function R'(alpha)
%The first derivative R'(alpha) is strictly increasing from -infinite to 0
%when alpha grows from 0 to 1. 
% Trp = sqrt((6*k^2*D1)/(-R'(alpha))) so one deduces that Trp grows from 0 to
% +infinite with the evolution of damage.
figure(2)
plot(alpha, RFirstDerivative, 'g', 'LineWidth', 2)
grid on
grid minor
xlim([0  1])
ylim([-100  0])
% Labels and Titles
xlabel('\alpha', 'FontSize', 15)
ylabel('R\prime(\alpha)', 'FontSize', 12)
hold on
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 1)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 1)
hold on

%Function R''(alpha)
%Because m > 1 > n > 0, R''(?) > 0 for every value of alpha.
figure(3)
%subplot (2,3,3)
plot(alpha, RSecondDerivative, 'g', 'LineWidth', 2)
grid on
grid minor
xlim([0  1])
ylim([0  100])
% Labels and Titles
xlabel('\alpha', 'FontSize', 15)
ylabel('R\prime\prime(\alpha)', 'FontSize', 12)
hold on
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 1)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 1)
hold on

%Function S(alpha)
figure (4)
%subplot (2,3,4)
plot(alpha, S, 'b', 'LineWidth', 2)
grid on
grid minor
xlim([0  1])
ylim([0  1000])
% Labels and Titles
xlabel('\alpha', 'FontSize', 15)
ylabel('S(\alpha)', 'FontSize', 12)
hold on
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 1)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 1)
hold on

%Function S'(alpha)
figure(5)
plot(alpha, SFirstDerivative, 'b', 'LineWidth', 2)
grid on
grid minor
xlim([0  1])
ylim([0  1000])
% Labels and Titles
xlabel('\alpha', 'FontSize', 15)
ylabel('S\prime(\alpha)', 'FontSize', 12)
hold on
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 1)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 1)
hold on

%Function S''(alpha)
figure(6)
plot(alpha, SSecondDerivative, 'b', 'LineWidth', 2)
grid on
grid minor
xlim([0  1])
ylim([-1000  10000])
% Labels and Titles
xlabel('\alpha', 'FontSize', 15)
ylabel('S\prime\prime(\alpha)', 'FontSize', 12)
hold on
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 1)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 1)
hold on

%The study of the sign of S''(alpha) shows that S''(alpha) is negative in the
%interval (0, alpha0), positive in the interval (alpha0, 1) with alpha0 given by:

alpha0 = 1/(m-n)*(sqrt((m*n)/(m-n+1))-n);

%Function S''(alpha)
figure(7)
title('S\prime\prime(\alpha)', 'FontSize', 15)
plot(alpha, SSecondDerivative, 'b', 'LineWidth', 2)
grid on
grid minor
xlim([0  1])
ylim([-1000  10000])
hold on
y = linspace(alpha0, alpha0, length(alpha));
plot(y, SSecondDerivative, 'k', 'LineWidth', 1, 'LineStyle', '--')
grid on
grid minor
xlim([0  1])
ylim([-1000  10000])
hold on
% Labels and Titles
xlabel('\alpha', 'FontSize', 15)
ylabel('S\prime\prime(\alpha)', 'FontSize', 12)
hold on
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 1)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 1)
hold on
grid minor
legend('S\prime\prime(\alpha)', '\alpha0')


