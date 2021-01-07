clc, clear, close all
%% Initial Conditions

% In all the pictures presented in the present section, one uses the
% following values for the Poisson ratio and the parameters critical stress (Sigmac) and internal friction (k) of
% Drucker-Prager criterion:
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

% The responses depend on the two exponents m and n, and on the two ratios R1/E0 and p0/?c
m = 2;
n = 0.5;
R1 = 0.05*E0;
R2 = 0.1*E0;
R3 = 0.15*E0;
R4 = 0.2*E0;
R5 = 0.25*E0;

% With:
% R1 > 0
% m > 1
% 0 < n < 1

% Considering the confining pressure equal to:
p0 = 0;

alpha = 0:0.005:1;
R = zeros(size(alpha));
S = zeros(size(alpha));

%% Calculate  

%R1 = 0.05E0
for i = 1:length(alpha)
    R005(i) = R1 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative1 (i) = - R1*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative1(i)  = R1*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    S1(i)= 1/R1* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative1(i) = 1/R1 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative1(i) = 1/R1*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    
    Fi1 = SSecondDerivative1./RSecondDerivative1.*((-RFirstDerivative1).^(3/2))./((SFirstDerivative1).^(3/2));
    
    Trp1(i) = sqrt((6*k^2*D1)/(-RFirstDerivative1(i)));
    
    SigmaZ1(i) = -3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative1(i)))/(1-k);
    
    EpsilonZ1 (i) = SigmaZ1(i)/E0 - (1/k - 1)* Trp1(i)/3 + 2*PoissonRatio*p0/E0;
    
    EpsilonV1 (i) = (SigmaZ1(i) - 2*p0)/(3*K0) + Trp1 (i);
end

%R1 = 0.10E0
for i = 1:length(alpha)
    R01(i) = R2 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative2 (i) = - R2*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative2(i)  = R2*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    S2(i)= 1/R2* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative2(i) = 1/R2 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative2(i) = 1/R2*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    
    Fi2 = SSecondDerivative2./RSecondDerivative2.*((-RFirstDerivative2).^(3/2))./((SFirstDerivative2).^(3/2));
    
    Trp2(i) = sqrt((6*k^2*D1)/(-RFirstDerivative2(i)));
    
    SigmaZ2(i) = -3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative2(i)))/(1-k);
    
    EpsilonZ2 (i) = SigmaZ2(i)/E0 - (1/k - 1)* Trp2(i)/3 + 2*PoissonRatio*p0/E0;
    
    EpsilonV2 (i) = (SigmaZ2(i) - 2*p0)/(3*K0) + Trp2 (i);
end

%R1 = 0.15E0
for i = 1:length(alpha)
    R015(i) = R3 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative3 (i) = - R3*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative3(i)  = R3*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    S3(i)= 1/R3* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative3(i) = 1/R3 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative3(i) = 1/R3*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    
    Fi3 = SSecondDerivative3./RSecondDerivative3.*((-RFirstDerivative3).^(3/2))./((SFirstDerivative3).^(3/2));
    
    Trp3(i) = sqrt((6*k^2*D1)/(-RFirstDerivative3(i)));
    
    SigmaZ3(i) = -3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative3(i)))/(1-k);
    
    EpsilonZ3 (i) = SigmaZ3(i)/E0 - (1/k - 1)* Trp3(i)/3 + 2*PoissonRatio*p0/E0;
    
    EpsilonV3 (i) = (SigmaZ3(i) - 2*p0)/(3*K0) + Trp3 (i);
end


%R1 = 0.2E0
for i = 1:length(alpha)
    R02(i) = R4 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative4 (i) = - R4*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative4(i)  = R4*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    S4(i)= 1/R4* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative4(i) = 1/R4 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative4(i) = 1/R4*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    
    Fi4 = SSecondDerivative4./RSecondDerivative4.*((-RFirstDerivative4).^(3/2))./((SFirstDerivative4).^(3/2));
    
    Trp4(i) = sqrt((6*k^2*D1)/(-RFirstDerivative4(i)));
    
    SigmaZ4(i) = -3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative4(i)))/(1-k);
    
    EpsilonZ4 (i) = SigmaZ4(i)/E0 - (1/k - 1)* Trp4(i)/3 + 2*PoissonRatio*p0/E0;
    
    EpsilonV4 (i) = (SigmaZ4(i) - 2*p0)/(3*K0) + Trp4 (i);
end


%R1 = 0.25E0
for i = 1:length(alpha)
    R025(i) = R5 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative5 (i) = - R5*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative5(i)  = R5*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    S5(i)= 1/R5* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative5(i) = 1/R5 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative5(i) = 1/R5*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    
    Fi5 = SSecondDerivative5./RSecondDerivative5.*((-RFirstDerivative5).^(3/2))./((SFirstDerivative5).^(3/2));
    
    Trp5(i) = sqrt((6*k^2*D1)/(-RFirstDerivative5(i)));
    
    SigmaZ5(i) = -3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative5(i)))/(1-k);
    
    EpsilonZ5 (i) = SigmaZ5(i)/E0 - (1/k - 1)* Trp5(i)/3 + 2*PoissonRatio*p0/E0;
    
    EpsilonV5 (i) = (SigmaZ5(i) - 2*p0)/(3*K0) + Trp5 (i);
end

%% Visualize the results - SigmaZ,EsilonZ

%We are plotted some graphs for different values of R1/E0. One can note
%that the larger R1, the larger the overstress but also more rapid is the
%decrease of the axial stress after the peak.

%A similar behavior can be seen for the volumetric strain: the larger R1,
%the greater the maximal contraction, but also more rapid the growing of
%the dilatancy after the peak.

figure (1)
plot(-EpsilonZ1, -SigmaZ1, 'm', 'LineWidth', 2)
hold on
plot(-EpsilonZ2, -SigmaZ2, 'g', 'LineWidth', 2)
hold on
plot(-EpsilonZ3, -SigmaZ3, 'b', 'LineWidth', 2)
hold on
plot(-EpsilonZ4, -SigmaZ4, 'r', 'LineWidth', 2)
hold on
plot(-EpsilonZ5, -SigmaZ5, 'c', 'LineWidth', 2)
hold on
grid on
grid minor
xlim([0  4])
ylim([0  1])
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
legend({'R1 = 0.05', 'R1 = 0.10', 'R1 = 0.15', 'R1 = 0.20', 'R1 = 0.25'});
%Figure 1 shown the influence of R1/E0 on the graph of the axial stress
%versus the axial strain when there is no confining pressure (uniaxial
%compression test)
%% Visualize the results - EpsilonV,EsilonZ

%A similar behavior can be seen for the volumetric strain: the larger R1,
%the greater the maximal contraction, but also more rapid the growing of
%the dilatancy after the peak.

figure(2)
plot(-EpsilonZ1, EpsilonV1, 'm', 'LineWidth', 2)
hold on
plot(-EpsilonZ2, EpsilonV2, 'g', 'LineWidth', 2)
hold on
plot(-EpsilonZ3, EpsilonV3, 'b', 'LineWidth', 2)
hold on
plot(-EpsilonZ4, EpsilonV4, 'r', 'LineWidth', 2)
hold on
plot(-EpsilonZ5, EpsilonV5, 'c', 'LineWidth', 2)
hold on
grid on
grid minor
xlim([0   3])
ylim([-0.5  2.2])
xlabel('\epsilon_z/\epsilon_c', 'FontSize', 12)
ylabel('\epsilon_v/\epsilon_c', 'FontSize', 12)
title('Volumetric strain \epsilon_v - Axial strain \epsilon_z', 'FontSize', 14)
hold on
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 1)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 1)
hold on
legend({'R1 = 0.05', 'R1 = 0.10', 'R1 = 0.15', 'R1 = 0.20', 'R1 = 0.25'});
%Figure 1 shown the influence of R1/E0 on the graph of the volumetric
%strain versus the axial strain when there is no confining pressure (uniaxial
%compression test)

%% Snap-back
%When R1 is grater than 0.275E0, the graph of SigmaZ contains a snap-back
%which induces, during a test where the axial strain is controlled and
%continuously decreasing, a discontinuity of the evolution of the axial
%stress and the damage

R6 = 0.4*E0;

%R1 = 0.4E0
for i = 1:length(alpha)
    R04(i) = R6 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative6 (i) = - R6*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative6(i)  = R6*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    S6(i)= 1/R6* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative6(i) = 1/R6 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative6(i) = 1/R6*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    
    Fi6 = SSecondDerivative6./RSecondDerivative6.*((-RFirstDerivative6).^(3/2))./((SFirstDerivative6).^(3/2));
    
    Trp6(i) = sqrt((6*k^2*D1)/(-RFirstDerivative6(i)));
    
    SigmaZ6(i) = -3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative6(i)))/(1-k);
           
    EpsilonZ6 (i) = SigmaZ6(i)/E0 - (1/k - 1)* Trp6(i)/3 + 2*PoissonRatio*p0/E0;
       
    EpsilonV6 (i) = (SigmaZ6(i) - 2*p0)/(3*K0) + Trp6 (i);
end

EpsilonZ6_1 = EpsilonZ6(1 : 35);
EpsilonZ6_2 = EpsilonZ6(36 : 144);
EpsilonZ6_3 = EpsilonZ6(145 : 201);

SigmaZ6_1 = SigmaZ6(1 : 35);
SigmaZ6_2 = SigmaZ6(36 : 144);
SigmaZ6_3 = SigmaZ6(145 : 201);

figure(3)
plot(-EpsilonZ6_1, -SigmaZ6_1, 'b', 'LineWidth', 2 )
hold on
plot(-EpsilonZ6_3, -SigmaZ6_3, 'b', 'LineWidth', 2 )
hold on
plot(-EpsilonZ6_2, -SigmaZ6_2, 'r', 'LineWidth', 1, 'LineStyle', '--')
grid on
grid minor
xlim([0  2])
ylim([0  1.2])
%Labels and Titles
xlabel('-\epsilon_z/\epsilon_c', 'FontSize', 12)
ylabel('-\sigma_z/\sigma_c', 'FontSize', 12)
title('Snap-back - Axial stress \sigma_z VS the Axial strain \epsilon_z', 'FontSize', 12)
hold on
%Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 1)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 1)
hold on
legend({'R1 = 0.40'});


