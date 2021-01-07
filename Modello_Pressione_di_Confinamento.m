clc, clear, close all
%% Initial Conditions

%In all the pictures presented in the present section, one uses the
%following values for the Poisson ratio and the parameters critical stress (?c) and internal friction (k) of
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

%The responses depend on the two exponents m and n, and on the two ratios R1/E0 and p0/sigmac
m = 2;
n = 0.5;
R1 = 0.1*E0;

%with:
% R1 > 0
% m > 1
% 0 < n < 1

% Considering the confining pressure equal to:
p00 = 0;
p01 = 0.2;
p02 = 0.4;
p03 = 0.6;
p04 = 0.8;
p05 = 1;

alpha = 0:0.005:1;
R = zeros(size(alpha));
S = zeros(size(alpha));

%% Calculate

for i = 1:length(alpha)
    R(i) = R1 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative(i) = - R1*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative(i)  = R1*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    S(i)= 1/R1* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative(i) = 1/R1 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative(i) = 1/R1*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    
    Fi = SSecondDerivative./RSecondDerivative.*((-RFirstDerivative).^(3/2))./((SFirstDerivative).^(3/2));
    
    Trp(i) = sqrt((6*k^2*D1)/(-RFirstDerivative(i)));
   
    SigmaZ00(i) = -3*CriticalStress/(1-k) - p00*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative(i)))/(1-k);
    SigmaZ01(i) = -3*CriticalStress/(1-k) - p01*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative(i)))/(1-k);
    SigmaZ02(i) = -3*CriticalStress/(1-k) - p02*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative(i)))/(1-k);
    SigmaZ03(i) = -3*CriticalStress/(1-k) - p03*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative(i)))/(1-k);
    SigmaZ04(i) = -3*CriticalStress/(1-k) - p04*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative(i)))/(1-k);
    SigmaZ05(i) = -3*CriticalStress/(1-k) - p05*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative(i)))/(1-k);
    
    EpsilonZ00 (i) = SigmaZ00(i)/E0 - (1/k - 1)* Trp(i)/3 + 2*PoissonRatio*p00/E0;
    EpsilonZ01 (i) = SigmaZ01(i)/E0 - (1/k - 1)* Trp(i)/3 + 2*PoissonRatio*p01/E0;
    EpsilonZ02 (i) = SigmaZ02(i)/E0 - (1/k - 1)* Trp(i)/3 + 2*PoissonRatio*p02/E0;
    EpsilonZ03 (i) = SigmaZ03(i)/E0 - (1/k - 1)* Trp(i)/3 + 2*PoissonRatio*p03/E0;
    EpsilonZ04 (i) = SigmaZ04(i)/E0 - (1/k - 1)* Trp(i)/3 + 2*PoissonRatio*p04/E0;
    EpsilonZ05 (i) = SigmaZ05(i)/E0 - (1/k - 1)* Trp(i)/3 + 2*PoissonRatio*p05/E0;
    
    EpsilonV00 (i) = (SigmaZ00(i) - 2*p00)/(3*K0) + Trp (i);
    EpsilonV01 (i) = (SigmaZ01(i) - 2*p01)/(3*K0) + Trp (i);
    EpsilonV02 (i) = (SigmaZ02(i) - 2*p02)/(3*K0) + Trp (i);
    EpsilonV03 (i) = (SigmaZ03(i) - 2*p03)/(3*K0) + Trp (i);
    EpsilonV04 (i) = (SigmaZ04(i) - 2*p04)/(3*K0) + Trp (i);
    EpsilonV05 (i) = (SigmaZ05(i) - 2*p05)/(3*K0) + Trp (i);
end

%% Graph

%The graph of the axial stress vs the axial strain for a given confining
%pressure during the damage.plasticity phase is simply a translation of
%that corresponding to the uniaxial compression test, the translation being
%given by the values of the acial strain and the axial stress at the end of
%the elastic stage.

figure (1)
plot(EpsilonZ00, SigmaZ00, 'm', 'LineWidth', 2)
hold on
plot(EpsilonZ01, SigmaZ01, 'g', 'LineWidth', 2)
hold on
plot(EpsilonZ02, SigmaZ02, 'b', 'LineWidth', 2)
hold on
plot(EpsilonZ03, SigmaZ03, 'r', 'LineWidth', 2)
hold on
plot(EpsilonZ04, SigmaZ04, 'c', 'LineWidth', 2)
hold on
plot(EpsilonZ05, SigmaZ05, 'k', 'LineWidth', 2)
hold on
grid on
grid minor
xlim([-5    0])
ylim([-3  0])
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
legend({'p0 = 0.0', 'p0 = 0.2', 'p0 = 0.4', 'p0 = 0.6', 'p0 = 0.8', 'p0 = 1.0'});


figure (2)
plot(EpsilonZ00, EpsilonV00, 'm', 'LineWidth', 2)
hold on
plot(EpsilonZ01, EpsilonV01, 'g', 'LineWidth', 2)
hold on
plot(EpsilonZ02, EpsilonV02, 'b', 'LineWidth', 2)
hold on
plot(EpsilonZ03, EpsilonV03, 'r', 'LineWidth', 2)
hold on
plot(EpsilonZ04, EpsilonV04, 'c', 'LineWidth', 2)
hold on
plot(EpsilonZ05, EpsilonV05, 'k', 'LineWidth', 2)
hold on
grid on
grid minor
xlim([-5  0])
ylim([-3  5])
% Labels and Titles
xlabel('-\epsilon_z/\epsilon_c', 'FontSize', 12)
ylabel('\epsilon_v/\epsilon_c', 'FontSize', 12)
title('Volumetric strain \epsilon_v - Axial strain \epsilon_z', 'FontSize', 14)
hold on
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 1)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 1)
hold on
legend({'p0 = 0.0', 'p0 = 0.2', 'p0 = 0.4', 'p0 = 0.6', 'p0 = 0.8', 'p0 = 1.0'});


q0 = SigmaZ00 - p00;
q1 = SigmaZ01 - p01;
q2 = SigmaZ02 - p02;
q3 = SigmaZ03 - p03;
q4 = SigmaZ04 - p04;
q5 = SigmaZ05 - p05;

figure (3)
plot(EpsilonZ00, q0, 'm', 'LineWidth', 2)
hold on
plot(EpsilonZ01, q1, 'g', 'LineWidth', 2)
hold on
plot(EpsilonZ02, q2, 'b', 'LineWidth', 2)
hold on
plot(EpsilonZ03, q3, 'r', 'LineWidth', 2)
hold on
plot(EpsilonZ04, q4, 'c', 'LineWidth', 2)
hold on
plot(EpsilonZ05, q5, 'k', 'LineWidth', 2)
hold on
grid on
grid minor
xlim([-5   0])
ylim([-3.5  0])
% Labels and Titles
xlabel('-\epsilon_z/\epsilon_c', 'FontSize', 12)
ylabel('-(\sigma_z-\sigma_x)/\sigma_c', 'FontSize', 12)
title('Deviatoric stress q - Axial strain \epsilon_z', 'FontSize', 14)
hold on
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 1)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 1)
hold on
legend({'p0 = 0.0', 'p0 = 0.2', 'p0 = 0.4', 'p0 = 0.6', 'p0 = 0.8', 'p0 = 1.0'});
