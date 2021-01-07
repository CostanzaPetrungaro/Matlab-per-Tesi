clc, clear, close all

%% Initial Conditions

%In all the pictures presented in the present section, one uses the
%following values for the Poisson ratio (?0) and the parameters critical stress (?c) and internal friction (k) of
%Drucker-Prager criterion:
PoissonRatio = 0.2;
CriticalStress = 0;
k = 0.2;

% On set:
% SigmaC := ?(D1*E0)
% EpsilonC := ?(D1/E0)
%Where  the material constant SigmaC fixes the scale of the stresses whereas
%the material constant EpsilonC fixes the scale of the strains.

% Modello paraetrizzato rispetto a E0 e D1, per semplicità si
% considereranno valori unitari
E0 = 1;
D1 = 1;

% The elastic moduli K0 and Mu0 are related to the young modulus E0 and the
% Poisson ratio ?0 of the sound material by:
%3K0 = E0/(1 - 2*PoissonRatio)
%2Mu0 = E0/(1 + PoissonRatio)
K0 = 1/3*E0/(1-2*PoissonRatio);

%The responses depend on the two exponents m and n, and on the two ratios R1/E0 and p0/?c
m = 2;
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
for i = 1:length(alpha)
    R(i) = R1 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative(i) = - R1*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative(i)  = R1*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    S(i)= 1/R1* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative(i) = 1/R1 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative(i) = 1/R1*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    
    Fi = SSecondDerivative./RSecondDerivative.*((-RFirstDerivative).^(3/2))./((SFirstDerivative).^(3/2));
    FiApprossimata(i) = (m+1)/(m-1)*(1-alpha(i))^m;
    Fi0(i) = -(1+n)/(1-n)*alpha(i)^(-n);
    
    Trp(i)= sqrt((6*k^2*D1)/(-RFirstDerivative(i)));
    
    SigmaZ(i) = -3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative(i)))/(1-k);
    
    EpsilonZ (i) = SigmaZ(i)/E0 - (1/k - 1)* Trp(i)/3 + 2*PoissonRatio*p0/E0;
    
    EpsilonV (i) = (SigmaZ(i) - 2*p0)/(3*K0) + Trp (i);
end

%% Snap-back condition
%When there is snap-back, the control of the axial strain (?z) necessarily
%implies a discontinuous evolution of the damage and that can lead to a
%brutal rupture of the volume element with a sudden jump of ? to 1 when the
%axial strain reaches its limit value.

%The snap-back exists id and only is the following condiion is satisfied:
% (s''(alpha)*(-R'(alpha)^1.5)/ (R''(alpha)*(S'(alpha)^1.5)) <= (1 - k)^2 * E0/3 
% The left hand side above can read as:
%(s''(alpha)*(-R'(alpha)^1.5)/ (R''(alpha)*(S'(alpha)^1.5)) = R1*Fi(alpha) where the function
%Fi(alpha) depends only the exponents m and n of the model.

% The non snap-back condition can as :
% R1 <= (1-k)^2*E0/(3*max(Fi(alpha)))
% the function Fi(alpha) reaches its upper bound in the interval (alpha0, 1) so the
% maximum of the function will be for alpha in (alpha0, 1)
FiMax = max(Fi);

%So the snap-back condition required that:
R1Limit = (1 - k)^2*E0/(3 * FiMax);

%There isn't snap-back if the harnening modulus R1 is small enough.
% R1 = 0.1E0  (hence there isn't snap-back)

%% Overstress 
alpha0 = 1/(m-n)*(sqrt((m*n)/(m-n+1))-n);
SFirstDerivativeAlpha0 = 1/R1 * ((m - n)* alpha0 + n)/(alpha0^(1-n)*(1 - alpha0)^(m + 1));
DeltaSigmaZ = sqrt((6*D1)/((1-k)^2*SFirstDerivativeAlpha0));

%% function Fi

% figure(1)
% plot(alpha, Fi, 'g', 'LineWidth', 2)
% grid on
% grid minor
% xlim([0  1])
% ylim([-0.5  0.1])
% hold on
% grid minor
% y = linspace(alpha0, alpha0, length(alpha));
% plot(y, Fi, 'k', 'LineWidth', 1, 'LineStyle', '--')
% grid on
% grid minor
% xlim([0  1])
% ylim([-0.5  0.1])
% hold on
% xlabel('\alpha', 'FontSize', 15)
% ylabel('\phi(\alpha)', 'FontSize', 12)
% title('\phi(\alpha)', 'FontSize', 14)
% hold on
% % Black line x = 0
% plot([0 0], ylim, 'k', 'LineWidth', 1)
% % Black line y = 0
% plot(xlim, [0 0], 'k', 'LineWidth', 1)
% hold on
% legend('\phi(\alpha)', '\alpha0')
% 
% % One observes a contractance during the phase where  Fi(alpha) is negative,
% % that is only during the stress-hardening phases. During the stress
% % softening.phases, since Fi(alpha)>0, one observes a dilatancy.
% 
% %Therefore a contractance can happen only at the beginninng of the test, as
% %long as alpha < alpha0 (indeed the snap-back can happen only for alpha > alpha0)
% 
% %% Fi approssimata
% 
% figure(2)
% plot(alpha, Fi, 'g', 'LineWidth', 2)
% grid on
% grid minor
% xlim([0  1])
% ylim([-0.5  2])
% hold on
% grid minor
% plot(alpha, FiApprossimata, 'b', 'LineWidth', 2)
% grid on
% grid minor
% xlim([0  1])
% ylim([-0.5  2])
% hold on
% xlabel('\alpha', 'FontSize', 15)
% ylabel('\phi(\alpha)', 'FontSize', 12)
% title('\phi(\alpha)', 'FontSize', 14)
% hold on
% grid minor
% % Black line x = 0
% plot([0 0], ylim, 'k', 'LineWidth', 1)
% % Black line y = 0
% plot(xlim, [0 0], 'k', 'LineWidth', 1)
% hold on
% grid minor
% legend('\phi(\alpha)','\phi Approssimata')
% %% Visualize the results - SigmaZ,EsilonZ
% 
% % Let us study the relation between the axial strain and the axial stress
% % in the case when condition of non snap-back is satisfied.
% % R1 = 0.1E0  (hence there isn't snap-back)
% 
% figure (3)
% plot(-EpsilonZ, -SigmaZ, 'b', 'LineWidth', 2)
% grid on
% grid minor
% hold on
% xlim([0  2])
% ylim([0 0.7])
% % Labels and Titles
% xlabel('-\epsilon_z/\epsilon_c', 'FontSize', 12)
% ylabel('-\sigma_z/\sigma_c', 'FontSize', 12)
% title('Axial Stress \sigma_z - Axial Strain \epsilon_z', 'FontSize', 12)
% hold on
% % Black line x = 0
% plot([0 0], ylim, 'k', 'LineWidth', 1)
% % Black line y = 0
% plot(xlim, [0 0], 'k', 'LineWidth', 1)
% hold on
% 
% %% Visualize the results - EpsilonV,EsilonZ
% 
% figure(4)
% plot(-EpsilonZ, EpsilonV, 'g', 'LineWidth', 2)
% grid on
% grid minor
% xlim([0  2.3])
% ylim([-0.3 1.4])
% xlabel('\epsilon_z/\epsilon_c', 'FontSize', 12)
% ylabel('\epsilon_v/\epsilon_c', 'FontSize', 12)
% title('Volumetric Strain \epsilon_v - Axial Strain \epsilon_z', 'FontSize', 12)
% hold on
% % Black line x = 0
% plot([0 0], ylim, 'k', 'LineWidth', 1)
% % Black line y = 0
% plot(xlim, [0 0], 'k', 'LineWidth', 1)
% hold on
% 
% %% Trp
% 
% figure(5)
% plot(alpha, Trp, 'r', 'LineWidth', 2)
% grid on
% grid minor
% xlim([0  1])
% ylim([0  20])
% hold on
% grid on
% grid minor
% % Labels and Titles
% xlabel('\alpha', 'FontSize', 15)
% ylabel('Trp', 'FontSize', 12)
% title('TrP(\alpha)', 'FontSize', 12)
% hold on
% grid minor
% % Black line x = 0
% plot([0 0], ylim, 'k', 'LineWidth', 1)
% % Black line y = 0
% plot(xlim, [0 0], 'k', 'LineWidth', 1)
% hold on
% 
% %% function FiApprossimata: alpha=0
% 
% figure(6)
% plot(alpha, Fi, 'g', 'LineWidth', 2)
% grid on
% grid minor
% xlim([0  1])
% ylim([-25  2])
% hold on
% grid minor
% plot(alpha, Fi0, 'm', 'LineWidth', 2)
% grid on
% grid minor
% xlim([0  1])
% ylim([-25  2])
% hold on
% xlabel('\alpha', 'FontSize', 15)
% ylabel('\phi(\alpha)', 'FontSize', 12)
% title('\phi(\alpha)', 'FontSize', 14)
% hold on
% grid minor
% % Black line x = 0
% plot([0 0], ylim, 'k', 'LineWidth', 1)
% % Black line y = 0
% plot(xlim, [0 0], 'k', 'LineWidth', 1)
% hold on
% grid minor
% legend('\phi(\alpha)','\phi(\alpha)Approssimata')


%% SigmaZ VS alpha

%Therefore, the absolute value of the axialstress increases and one
%observes a stress-hardening behavior when alpha grows from 0 to alpha0 and
%decreases and one observes a stress-softening behavior when alpha grows from
%alpha0 and 1.

%Stress behavior with alpha
figure(7)
plot(alpha, -SigmaZ, 'g', 'LineWidth', 2)
grid on
grid minor
xlim([0  1])
ylim([0  1])
hold on
y = linspace(alpha0, alpha0, length(alpha));
plot(y, alpha, 'k', 'LineWidth', 1, 'LineStyle', '--')
grid on
grid minor
xlim([0  1])
ylim([0  1])
% Labels and Titles
xlabel('\alpha', 'FontSize', 15)
ylabel('-\sigma_z/\sigma_c', 'FontSize', 12)
title('Axial stress \sigma_z - Damage \alpha', 'FontSize', 13)
hold on
grid minor 
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 1)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 1)
hold on
legend('\sigma_z', '\alpha0')


%% Overstress
%Moreover, since S'(0) = S'(1) = +infinite, SigmaZ starts from the value
%corresponding at the end of the elastic stage, then is increasing (in
%absolute value) until its maximum value corresponding to the time at which
%the damage reaches the value alpha0, and finally is decreasing (in absolute
%value) to come back to the initial value reached at the end of the elastic
%stage.

%The absolute value of the overstress (DeltaSigmaz) associated with the
%stress-hardening phase is given by:

SFirstDerivativeAlpha0 = 1/R1 * ((m - n)* alpha0 + n)/(alpha0^(1-n)*(1 - alpha0)^(m + 1));

DeltaSigmaZ = sqrt((6*D1)/((1-k)^2*SFirstDerivativeAlpha0));

