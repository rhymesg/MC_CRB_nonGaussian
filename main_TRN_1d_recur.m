% Author: Youngjoo Kim
% Y. Kim and H. Bang, "Monte-Carlo Calculation of Cramer-Rao Bound for non-Gaussian Recursive
% Filtering", APISAT, 2017.

%%% Recursive filtering

clear all;
%%

LENG = 40;
T_LENG = 100;
Dist = 1; % 0 for Gaussian, 1 for Laplacian

x = zeros(1,LENG);
sig_Z = 2;
sig = 5;
sig_p = 1;

L = 200; % Number of state realization

%%% generate trajectory
for k = 1:1:LENG;
    x(k) = (k-1)*T_LENG/LENG;
end


%% particle filter
Np = 200; % number of particles
N_MC = 500; % number of MC
x_est = zeros(N_MC, LENG);
x_err = zeros(N_MC, LENG);
init_err = sig*randn(1,N_MC);
for i = 1:1:N_MC;
    rng('shuffle');
    particle = x(1) + init_err(i)*ones(1,Np) + sig*randn(1,Np);
    x_est(i,1) = x(1) + init_err(i); % MLSE estimate 
    x_err(i,1) = x_est(i,1) - x(1);
    p = ones(1,Np);
for k = 2:1:LENG;
    
    particle = particle + (x(k) - x(k-1)) + sig_p*randn*ones(1,Np);
    
    particle = particle + sig_p*randn(1,Np);

    if Dist == 0
        Z_mea = meas_TRN_1d(x(k)) + sig_Z*randn;
    else
        Z_mea = meas_TRN_1d(x(k)) + laprnd(1,1,0,sig_Z);
    end
        
    for j = 1:1:Np;
        Z_est = meas_TRN_1d(particle(j));
       
        if Dist == 0
            var = sig_Z^2;
            p(j) = p(j) * (1/sqrt(2*pi*var)) * exp( -(Z_mea - Z_est)^2/(2*var) );
        else
            b = sig_Z / sqrt(2);
            p(j) = p(j) * (1/(2*b)) * exp( -abs(Z_mea - Z_est)/b);
        end
    end
    
%     [tmp, ind] = max(p); % MAP estimate
%     x_est(i,k) = particle(ind);
    
    p = p/sum(p);
    x_est(i,k) = dot(p, particle); % MLSE estimate 
    x_err(i,k) = x_est(i,k) - x(k);
    
    % resampling 
    Neff = 1/sum(p.^2);
    if (Neff < 0.3*Np)
        w_cdf = cumsum(p);
        for j = 1:1:Np;
            ind = find(rand <= w_cdf, 1);
            if (isempty(ind)) 
                [a, ind] = max(p); 
            end
            particle_tmp(j) = particle(ind);
        end
        particle = particle_tmp;
        p = ones(1,Np);
    end
end
end

x_err_RMS = sqrt(mean(x_err.^2,1));
disp(['PF completed']);

%%
%%% FIM estimation
J1 = 1/sig^2;
J2 = 1/sig^2;
LB1_res(1) = sqrt(inv(J1));
LB2_res(1) = sqrt(inv(J2));
for k = 2:1:LENG; % for each time step
    rng('shuffle');
    s = 0.001;
    prior = inv(J1)*randn(1,L);
    if (L == 1)
        prior = 0;
    end
    FIM_lin = 0;
    for j = 1:1:L;
        stat = x(k) + prior(j);
        %%% linearized estimation
        z_x1 = meas_TRN_1d(stat - s);
        z_x2 = meas_TRN_1d(stat + s);

        H = (z_x2 - z_x1)/(2*s);
        R = sig_Z^2;

        FIM_lin_t = H'*inv(R)*H;
        
        FIM_lin = (((j-1)/j)*FIM_lin) + ((1/j)*FIM_lin_t);
    end
    
    FIM_lin_res(k) = FIM_lin;
    J1 = inv(inv(J1) + sig_p^2) + FIM_lin;
    LB1_res(k) = sqrt(inv(J1));
    
    %%% Monte-Carlo estimation
    % setting for FIM estimation
    N = 200; %0.5e6; % Number of pseudodata vector, i.e., number of Hessian estimates.
    c = 0.001;
    c_tlde = 0.0011;
    
    HhatLbar = 0;
    FIM_mon = 0;
    prior = inv(J2)*randn(1,L);
    if (L == 1)
        prior = 0;
    end
    for j = 1:1:L
        thet = x(k) + prior(j);
        height = meas_TRN_1d(thet);
        for iN = 1:1:N
            %%% Generate Zpseudo
            if Dist == 0
                Zpseudo = height + sig_Z*randn;  % Gaussian
            else
                Zpseudo = height + laprnd(1,1,0,sig_Z); % Laplacian
            end
            
            Delk=2*round(rand())-1;
            Delk_tlde=2*round(rand())-1;

            %
            %%%% Estimation of FIM based on log-likelihood measurements %%%%
            thetpp = thet + (c_tlde*Delk_tlde) + (c*Delk);
            thetp = thet + (c*Delk);
            Z_pp = meas_TRN_1d(thetpp);
            Z_p = meas_TRN_1d(thetp);
            G1plus = (1/c_tlde)*(loglikelihood_TRN(Z_pp,Zpseudo,sig_Z,Dist)- ...
                loglikelihood_TRN(Z_p,Zpseudo,sig_Z,Dist))*(1./Delk_tlde);
            %
            thetpm = thet + (c_tlde*Delk_tlde) - (c*Delk);
            thetm = thet - (c*Delk);
            Z_pm = meas_TRN_1d(thetpm);
            Z_m = meas_TRN_1d(thetm);
            G1minus = (1/c_tlde)*(loglikelihood_TRN(Z_pm,Zpseudo,sig_Z,Dist)- ...
                loglikelihood_TRN(Z_m,Zpseudo,sig_Z,Dist))*(1./Delk_tlde);
            %
            JhatL = (1/(2*c))*(G1plus - G1minus)*(1./Delk)'; % Jacobian estimate
                            % based on log-likelihood measurements by using Spall 2005
            %
            Hhat0L = (1/2)*(JhatL+JhatL'); % Hessian estimate based on 
                            % log-likelihood measurements by using Spall 2005.
            %
            HhatL = Hhat0L;
            HhatLbar = (((iN-1)/iN)*HhatLbar) + ((1/iN)*HhatL);
            FhatL = -HhatLbar;
            if min(eig(FhatL)) < 0; FhatL = sqrtm(FhatL*FhatL);end
        end
        FIM_mon = (((j-1)/j)*FIM_mon) + ((1/j)*FhatL);
    end

    FIM_mon_res(k) = FIM_mon;
    J2 = inv(inv(J2) +sig_p^2) + FIM_mon;
    LB2_res(k) = sqrt(inv(J2));
end


%%

for k = 1:1:1000;
    height(k) = meas_TRN_1d(k*T_LENG/1000);
end

% figure;
% plot(T_LENG/1000:T_LENG/1000:T_LENG, height, 'k', 'linewidth', 1.5);
% axis equal
% xlabel('x Position');
% ylabel('Terrain Height');
% grid on;

figure
plot(LENG/1000:LENG/1000:LENG, height/40-1.0, 'k', 'linewidth', 1.5); hold on;
plot(1:1:LENG, x_err_RMS, 'k:', 'linewidth', 1.5); hold on;
plot(1:1:LENG, LB1_res, 'b', 'linewidth', 1.5); hold on;
plot(1:1:LENG, LB2_res, 'r--', 'linewidth', 1.5);
grid on;
xlabel('Time Step');
ylabel('RMSE');
legend('Terrain', 'PF RMSE', 'CRB 1', 'CRB 2');
axis([0 40 0 3.5]);

% figure
% plot(x, FIM_lin_res, 'b'); hold on;
% plot(x, FIM_mon_res, 'r');

%% distributions
% 
% for k = 1:1:800;
%     z = -4 + 8*(k-1)/800;
%     Gau(k) = exp(loglikelihood_TRN(z*3,0,3,0));
%     Lap(k) = exp(loglikelihood_TRN(z*3,0,3,1));
% end
% 
% figure;
% plot(Gau, 'b', 'linewidth', 2); hold on;
% plot(Lap, 'r', 'linewidth', 2); hold on;
% ylabel('Probability');
% grid on;

