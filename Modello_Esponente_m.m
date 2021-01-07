clc, clear, close all

%% Initial Conditions

%In all the pictures presented in the present section, one uses the
%following values for the Poisson ratio and the parameters critical stress (Sigmac) and internal friction (k) of
%Drucker-Prager criterion:
PoissonRatio = 0.2;
CriticalStress = 0;
k = 0.2;

% On set:
% SigmaC := sqrt(D1*E0)
% EpsilonC := sqrt(D1/E0)
% Where  the material constant SigmaC fixes the scale of the stresses whereas
% The material constant EpsilonC fixes the scale of the strains.

% Modello paraetrizzato rispetto a E0 e D1, per semplicità si
% considereranno valori unitari
E0 = 1;
D1 = 1;

% The elastic moduli K0 and Mu0 are related to the young modulus E0 and the
% Poisson ratio ?0 of the sound material by:
% 3K0 = E0/(1 - 2*PoissonRatio)
% 2Mu0 = E0/(1 + PoissonRatio)
K0 = 1/3*E0/(1-2*PoissonRatio);

%The responses depend on the two exponents m and n, and on the two ratios R1/E0 and p0/?c

n = 0.5;
R1 = 0.1*E0;

%with:
% R1 > 0
% m > 1
% 0 < n < 1

% Considering the confining pressure equal to:
p0 = 0;

alpha = 0:0.005:1;
R = zeros(size(alpha));
S = zeros(size(alpha));


%% Calculate

m = 1;
for i = 1:length(alpha)
    R(i) = R1 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative(i) = - R1*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative(i)  = R1*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    S(i)= 1/R1* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative(i) = 1/R1 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative(i) = 1/R1*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    
    Fi = SSecondDerivative./RSecondDerivative.*((-RFirstDerivative).^(3/2))./((SFirstDerivative).^(3/2));
    
    Trp(i) = sqrt((6*k^2*D1)/(-RFirstDerivative(i)));
    
    SigmaZ1(i) = -3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative(i)))/(1-k);
    
    EpsilonZ1 (i) = SigmaZ1(i)/E0 - (1/k - 1)* Trp(i)/3 + 2*PoissonRatio*p0/E0;
    
    EpsilonV1 (i) = (SigmaZ1(i) - 2*p0)/(3*K0) + Trp (i);
end


m = 2;
for i = 1:length(alpha)
    R(i) = R1 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative(i) = - R1*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative(i)  = R1*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    S(i)= 1/R1* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative(i) = 1/R1 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative(i) = 1/R1*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    
    Fi = SSecondDerivative./RSecondDerivative.*((-RFirstDerivative).^(3/2))./((SFirstDerivative).^(3/2));
    
    Trp(i) = sqrt((6*k^2*D1)/(-RFirstDerivative(i)));
    
    SigmaZ2(i) = -3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative(i)))/(1-k);
    
    EpsilonZ2 (i) = SigmaZ2(i)/E0 - (1/k - 1)* Trp(i)/3 + 2*PoissonRatio*p0/E0;
    
    EpsilonV2 (i) = (SigmaZ2(i) - 2*p0)/(3*K0) + Trp (i);
end

m = 3;
for i = 1:length(alpha)
    R(i) = R1 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative(i) = - R1*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative(i)  = R1*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    S(i)= 1/R1* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative(i) = 1/R1 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative(i) = 1/R1*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    
    Fi = SSecondDerivative./RSecondDerivative.*((-RFirstDerivative).^(3/2))./((SFirstDerivative).^(3/2));
    
    Trp(i) = sqrt((6*k^2*D1)/(-RFirstDerivative(i)));
    
    SigmaZ3(i) = -3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative(i)))/(1-k);
    
    EpsilonZ3 (i) = SigmaZ3(i)/E0 - (1/k - 1)* Trp(i)/3 + 2*PoissonRatio*p0/E0;
    
    EpsilonV3 (i) = (SigmaZ3(i) - 2*p0)/(3*K0) + Trp (i);
end


m = 4;
for i = 1:length(alpha)
    R(i) = R1 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative(i) = - R1*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative(i)  = R1*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    S(i)= 1/R1* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative(i) = 1/R1 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative(i) = 1/R1*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    
    Fi = SSecondDerivative./RSecondDerivative.*((-RFirstDerivative).^(3/2))./((SFirstDerivative).^(3/2));
    
    Trp(i) = sqrt((6*k^2*D1)/(-RFirstDerivative(i)));
    
    SigmaZ4(i) = -3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative(i)))/(1-k);
    
    EpsilonZ4 (i) = SigmaZ4(i)/E0 - (1/k - 1)* Trp(i)/3 + 2*PoissonRatio*p0/E0;
    
    EpsilonV4 (i) = (SigmaZ4(i) - 2*p0)/(3*K0) + Trp (i);
end


m = 5;
for i = 1:length(alpha)
    R(i) = R1 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative(i) = - R1*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative(i)  = R1*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    S(i)= 1/R1* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative(i) = 1/R1 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative(i) = 1/R1*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    
    Fi = SSecondDerivative./RSecondDerivative.*((-RFirstDerivative).^(3/2))./((SFirstDerivative).^(3/2));
    
    Trp(i) = sqrt((6*k^2*D1)/(-RFirstDerivative(i)));
    
    SigmaZ5(i) = -3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative(i)))/(1-k);
    
    EpsilonZ5 (i) = SigmaZ5(i)/E0 - (1/k - 1)* Trp(i)/3 + 2*PoissonRatio*p0/E0;
    
    EpsilonV5 (i) = (SigmaZ5(i) - 2*p0)/(3*K0) + Trp (i);
end



%% Visualize the results - SigmaZ,EsilonZ

%The shape of the graphs giving the axial stress in term of the axial
%strain during an uniaxial compression test without confining pressure, is
%always the same whatever the values of m and n provided that they remain
%inside the admissible intervals m > 1 > n > 0.

%We can see the dependence of the response on m in the figure 1

figure (1)
plot(-EpsilonZ1, -SigmaZ1, 'r', 'LineWidth', 2, 'LineStyle', '-.' )
hold on
plot(-EpsilonZ2, -SigmaZ2, 'b', 'LineWidth', 2)
hold on
plot(-EpsilonZ3, -SigmaZ3, 'm', 'LineWidth', 2)
hold on
plot(-EpsilonZ4, -SigmaZ4, 'c', 'LineWidth', 2)
hold on
plot(-EpsilonZ5, -SigmaZ5, 'g', 'LineWidth', 2)
hold on
grid on
grid minor
xlim([0  4])
ylim([0  0.8])
% Labels and Titles
xlabel('-\epsilon_z/\epsilon_c', 'FontSize', 12)
ylabel('-\sigma_z/\sigma_c', 'FontSize', 12)
title('Axial stress \sigma_z - Axial strain \epsilon_z', 'FontSize', 14)
hold on
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 1)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 1)
hold on
legend({'m = 1', 'm = 2', 'm = 3', 'm = 4', 'm = 5'});

% We note taht as expected the m parameter influence mainly the post-pic
% stress behavior (for m = 1) the ultimate damage state achieved even for
% the finire value of epsilon Z) while n parameter alternate mainly the
% initial hardening state leaving the large epsilon Z asymptote intact.

%Remember that both the parameters (m and n) contribute to the overstress
%estimation

%% Visualize the results - EpsilonV,EsilonZ

%The shape of the graphs giving the axial stress in term of the axial
%strain during an uniaxial compression test without confining pressure, is
%always the same whatever the values of m and n provided that they remain
%inside the admissible intervals m > 1 > n > 0.

%We can see the dependence of the response on m in the figure 1

figure (5)
plot(-EpsilonZ1, EpsilonV1, 'r', 'LineWidth', 2, 'LineStyle', '-.' )
hold on
plot(-EpsilonZ2, EpsilonV2, 'b', 'LineWidth', 2)
hold on
plot(-EpsilonZ3, EpsilonV3, 'm', 'LineWidth', 2)
hold on
plot(-EpsilonZ4, EpsilonV4, 'c', 'LineWidth', 2)
hold on
plot(-EpsilonZ5, EpsilonV5, 'g', 'LineWidth', 2)
hold on
grid on
grid minor
xlim([0  2])
ylim([-0.25  1.5])
% Labels and Titles
xlabel('-\epsilon_z/\epsilon_c', 'FontSize', 12)
ylabel('-\epsilon_v/\epsilon_c', 'FontSize', 12)
title('Volumetric strain \epsilon_v - Axial strain \epsilon_z', 'FontSize', 14)
hold on
grid minor
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 1)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 1)
hold on
grid minor
legend({'m = 1', 'm = 2', 'm = 3', 'm = 4', 'm = 5'});



%% Overstress vs m
n = 0.5;
m = 0.5: 0.1 : 10;
alpha0 = zeros(size(m));

for i = 1:length(m)
    alpha0(i) = 1/(m(i)-n)*(sqrt((m(i)*n)/(m(i)-n+1))-n);
    SFirstDerivativeAlpha0(i) = 1/R1 * ((m(i) - n)* alpha0(i) + n)/(alpha0(i)^(1-n)*(1 - alpha0(i))^(m(i) + 1));
    DeltaSigmaZ(i) = sqrt((6*D1)/((1-k)^2*SFirstDerivativeAlpha0(i)));
end 

figure(2)
plot(m, DeltaSigmaZ, 'g', 'LineWidth', 2)
grid on
grid minor
xlim([0  10])
ylim([0  1])
xlabel('m', 'FontSize', 14)
ylabel('\Delta\sigma_z/\sigma_c', 'FontSize', 12)
title('Overstress \Delta\sigma_z - m', 'FontSize', 14)
hold on
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 1)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 1)
hold on
legend('n = 0.5')
