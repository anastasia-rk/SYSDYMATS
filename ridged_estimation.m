my_init
%% Load estimated thetas
dataset = 'C'; % 'D'; % 
metaFileName = ['Meta_',dataset];
load(metaFileName);
n_y = 0;
n_u = 4;
T = 2000;
folder = 'Results';                                                         % specify category where to save files
names = {'set','ny','nu'};                                                  % names used to define results folder name (no more than 3).
folderName = make_folder(folder,names,dataset,n_y,n_u);                     % create results folder
fileName = [folderName,'/OLS_results_T_',num2str(T),'.mat'];
load(fileName);
addpath('LASSO_Shmidt')                                                     % for lasso library
%% Select independent variables for the regression
load('External_parameters');
L_cut_all = [values{1}(:, 9);values{2}(:, 9)];
D_rlx_all = [values{1}(:,11);values{2}(:,11)];
A_imp_all = [values{1}(:, 6);values{2}(:, 6)];
V_imp_all = [values{1}(:, 7);values{2}(:, 7)];

index = find(Files <=10);
Files_sub = Files(index);
x = L_cut_all(Files_sub,1);
y = D_rlx_all(Files_sub,1);
%% Form the regression matrix
id = ones(size(x));                                                         % create unit vector for constants
A = [id x y x.*y x.^2 y.^2];                                               % create the matrix for ls
R_aa = A'*A;
names = sym('beta_0')
names = [names sym('beta_',[1 size(A,2)-1])];
if length(x) <= 4
    A = A(:,1:4);
end
M_d = size(A,2);
clear beta
B = Theta(:,Files_sub)';
betas_ols = A\B;
B_hat = A*betas_ols;
%% Collinearity analysis
% Diagnose collinearity based on condition number of the scaled regression
% matrix. Small eigenvalues of it's square also indicate how many columns
% have near-linear dependency

% Normalise data matrix
A_2 = A(:,2:end)
for j=1:2 %size(A_2,2)
    m(j) = mean(A_2(:,j));
    c(j) = sqrt(sum((A_2(:,j) - m(j)).^2));
end
c(3) = c(1)*c(2);
m(3) = m(1)*m(2);
c(4) = c(1)*c(1);
m(4) = m(1)*m(1);
c(5) = c(2)*c(2);
m(5) = m(2)*m(2);
for j=1:size(A_2,2)
    A_norm(:,j) = (A_2(:,j)- m(j))/c(j);
    A_n(:,j) = (A_2(:,j))/c(j);
end
n = size(A,1)
R_a2 = A_norm'*A_norm;
R_22 = R_a2(2:end,2:end);
O = zeros(size(A_norm,2),1);
for j = 2:size(A_norm,2)
r(j-1,1) = A_norm(:,1)'*A_norm(:,j);
end
R_new = [1 r'; r R_22];
R_norm = [n O'; O R_new];
A_s = [id A_norm];
R_ss = A_s'*A_s

A_sn = [id A_n];
R_n = A_sn'*A_sn
eig_r = eig(R_n);
kappa = cond(A_sn);
thr = 100;
if kappa > thr
    fprintf(['Condition number of the scaled matrix is ', num2str(kappa)]);
    for j = 1:size(A_norm,2)
        a_star{j} = A_norm(:,j);                                            % Regress this column on other covariates
        a_mean = mean(a_star{j});
        X_star = A_norm;
        X_star(:,j) = [];
        Hat = X_star*inv(X_star'*X_star)*X_star';                           % Hat matrix
        Eye_n  = eye(size(Hat));
        VIF(j) = norm((Eye_n - Hat)*a_star{j})^(-2);                        % =  1 - R^2_j                  
%         % Check is determination coefficients computed correctly 
%         % 1. Via norm
         Rs_1(j) = 1 - norm((Eye_n - Hat)*a_star{j})^2;
%         b = X_star\a_star{j};
%         a_hat{j} = X_star*b;
%         % 2. Via direct formula using estimates
%         Rs_2(j) = r_squared(a_star{j},a_hat{j});
%         % 3. Via RSS
%         RSS{j} = 1 - a_star{j}'*Hat*a_star{j};                              % Residual sum of squares
%         Rs_3(j) = 1 - RSS{j};                                               % For normalised data denominator is 1
        
    end
end

%% Ridged estimation from normalised centred data (in closed form) 
nPlots = round(size(A,2)/2);
log_max = 0;
log_min = -6;
vec = [0:1/50:1];
coeffs = sort(10.^(log_min + (log_max-log_min)*vec));
I_k = eye(size(R_ss));
for iTheta = 1:finalTerm
        ik = 0;
        for k = coeffs
            ik = ik + 1;
            th_norm = normalize(B(:,iTheta));
            gain                   = inv(R_ss + k*I_k)*A_s';
            beta{iTheta}(:,ik)     = gain*th_norm;
            Hat                    = A_s*gain;
            F_inv                  = inv(R_ss + k*I_k);
            sing_F                 = svd(F_inv);
            sigma_b{iTheta,ik}(:,:)= gain*A_s*inv(R_ss + k*I_k);
            mean{iTheta}(:,ik)     = A_s*beta{iTheta}(:,ik);
            sigma_t{iTheta,ik}(:,:) = diag((th_norm - mean{iTheta}(:,ik)).^2);
            Lik                    = mvnpdf(th_norm,mean{iTheta}(:,ik),sigma_t{iTheta,ik}(:,:));
            Cov_cmplx{iTheta}(ik)  = (M_d/2)*log((1/M_d)*sum(sing_F)*prod(sing_F)^(1/M_d));
            AIC{iTheta}(:,ik)      = 2*trace(Hat) - 2*log(Lik);
            ICOMP{iTheta}(:,ik)    = AIC{iTheta}(:,ik) + Cov_cmplx{iTheta}(ik);
            % Raw data
            gain                   = inv(R_aa + k*I_k)*A';
            beta_raw{iTheta}(:,ik) = gain*B(:,iTheta);
            Hat                    = A*gain;
            F_inv                  = inv(R_aa + k*I_k);
            sing_F                 = svd(F_inv);
            sigma_b_raw{iTheta,ik}(:,:)= gain*A*inv(R_aa + k*I_k);
            mean_th                = A*beta_raw{iTheta}(:,ik);
            sigma_th               = diag((th_norm - mean_th).^2);
            Lik                    = mvnpdf(th_norm,mean_th,sigma_th);
            Cov_cmplx_raw{iTheta}(ik)  = (M_d/2)*log((1/M_d)*sum(sing_F)*prod(sing_F)^(1/M_d));
            AIC_raw{iTheta}(:,ik)  = 2*trace(Hat) - 2*log(Lik);
            ICOMP_raw{iTheta}(:,ik)= AIC_raw{iTheta}(:,ik) + Cov_cmplx_raw{iTheta}(ik);
            betas_lasso{iTheta}(:,ik) =  LassoShooting(A_s,th_norm,k,'verbose',0);
            betas_lasso_raw{iTheta}(:,ik) = LassoShooting(A,B(:,iTheta),k,'verbose',0);
        end
        [AIC_min(iTheta),ind(iTheta)] = min(AIC{iTheta});
        [AIC_raw_min(iTheta),ind_raw(iTheta)] = min(AIC_raw{iTheta});
        [ICOMP_min(iTheta),ind_comp(iTheta)] = min(ICOMP{iTheta});
        [ICOMP_raw_min(iTheta),ind_comp_raw(iTheta)] = min(ICOMP_raw{iTheta});
%         [betas_lasso{iTheta},fit_lasso{iTheta}] =  Lasso(A_s,B(:,iTheta));
%         [betas_lasso_raw{iTheta},fit_lasso_raw{iTheta}] =  Lasso(A,B(:,iTheta));
end

%% Plot ridge traces from normalised and raw data
for iBeta=1:size(A,2)
      temp = arrayfun(@char, names(iBeta), 'uniform', 0);
    if length(temp) > 0
        str = temp{1};
        for iString=2:length(temp)
            str = [str,temp{iString}];
        end
    end
    Legends{iBeta} = strcat('$\',str,'$');
    clear temp
end
for iTheta = 1:finalTerm
        figName = char(symb_term{S(iTheta)});
        fig =  figure('Name',figName,'NumberTitle','off');
        subplot(2,2,1);
        for ib = 1:size(A,2)
            h1 = semilogx(coeffs,beta{iTheta}(ib,:),'-o','MarkerSize',5); hold on;
            set(h1, 'markerfacecolor', get(h1, 'color')); 
        end
        for ib = 1:size(A,2)
            h2 = semilogx(coeffs(ind(iTheta)),beta{iTheta}(ib,ind(iTheta)),'-k'); hold on;
            set(h2, 'markerfacecolor', get(h2, 'color')); 
        end
        title('Tikhonov normalised')
        xlabel('$\gamma$');
        ylabel('Standardised estimate')
        subplot(2,2,2);
        for ib = 1:size(A,2)
            h1 = semilogx(coeffs,beta_raw{iTheta}(ib,:),'-o','MarkerSize',5); hold on;
            set(h1, 'markerfacecolor', get(h1, 'color')); 
        end
        for ib = 1:size(A,2)
            h2 = semilogx(coeffs(ind_raw(iTheta)),beta{iTheta}(ib,ind_raw(iTheta)),'-k'); hold on;
            set(h2, 'markerfacecolor', get(h2, 'color')); 
        end
        title('Tikhonov raw')
        xlabel('$\gamma$');
        ylabel('Raw data estimate')
        legend(Legends);
        subplot(2,2,3);
        for ib = 1:size(A,2)
%             h1 = semilogx(fit_lasso{iTheta}.Lambda,betas_lasso{iTheta}(ib,:),'-o','MarkerSize',5); hold on;
            h1 = semilogx(coeffs,betas_lasso{iTheta}(ib,:),'-o','MarkerSize',5); hold on;
            set(h1, 'markerfacecolor', get(h1, 'color')); 
        end
        title('LASSO normalised')
        xlabel('$\gamma$');
        ylabel('Standardised estimate')
        subplot(2,2,4);
        for ib = 1:size(A,2)
%             h1 = semilogx(fit_lasso_raw{iTheta}.Lambda,betas_lasso_raw{iTheta}(ib,:),'-o','MarkerSize',5); hold on;
            h1 = semilogx(coeffs,betas_lasso_raw{iTheta}(ib,:),'-o','MarkerSize',5); hold on;
            set(h1, 'markerfacecolor', get(h1, 'color')); 
        end
        title('LASSO raw')
        xlabel('$\gamma$');
        ylabel('Raw data estimate')
        legend(Legends);
        
%          tikzName = [folderName,'/ridges_th_',num2str(iTheta),'.tikz'];
%         cleanfigure;
%         matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
%         false, 'height', '12cm', 'width','15cm','checkForUpdates',false);
        pngName = [folderName,'/ridges_th_',num2str(iTheta),'.png'];
        saveas(fig,pngName);
end
%% Criterion 1 = Akaike information criterion
figName = 'AIC_norm';
figure('Name',figName,'NumberTitle','off');
subplot(1,2,1);
for iTheta = 1:finalTerm
    h1 = semilogx(coeffs,AIC{iTheta},'-o','MarkerSize',5); hold on;
    set(h1, 'markerfacecolor', get(h1, 'color')); 
    Leg_theta{iTheta} = strcat('$\theta_',num2str(iTheta),'$');
    h2 = semilogx(coeffs(ind(iTheta)),AIC_min(iTheta),'ok','MarkerSize',10); hold on;
    set(h2, 'markerfacecolor', get(h2, 'color')); 
end
xlabel('$\lambda$');
ylabel('AIC')
subplot(1,2,2);
for iTheta = 1:finalTerm
    h1 = semilogx(coeffs,ICOMP{iTheta},'-o','MarkerSize',5); hold on;
    set(h1, 'markerfacecolor', get(h1, 'color')); 
    Leg_theta{iTheta} = strcat('$\theta_',num2str(iTheta),'$');
    h2 = semilogx(coeffs(ind_comp(iTheta)),ICOMP_min(iTheta),'ok','MarkerSize',10); hold on;
    set(h2, 'markerfacecolor', get(h2, 'color')); 
end
xlabel('$\lambda$');
ylabel('ICOMP')
tikzName = [folderName,'/',figName,'.tikz'];
        cleanfigure;
        matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '7cm', 'width','12cm','checkForUpdates',false);
        
figName = 'AIC_raw';
figure('Name',figName,'NumberTitle','off');
subplot(1,2,1);
for iTheta = 1:finalTerm
    h1 = semilogx(coeffs,AIC_raw{iTheta},'-o','MarkerSize',5); hold on;
    set(h1, 'markerfacecolor', get(h1, 'color')); 
    Leg_theta{iTheta} = strcat('$\theta_',num2str(iTheta),'$');
    h2 = semilogx(coeffs(ind_raw(iTheta)),AIC_raw_min(iTheta),'ok','MarkerSize',10); hold on;
    set(h2, 'markerfacecolor', get(h2, 'color')); 
end
xlabel('$\lambda$');
ylabel('AIC')
subplot(1,2,2);
for iTheta = 1:finalTerm
    h1 = semilogx(coeffs,ICOMP_raw{iTheta},'-o','MarkerSize',5); hold on;
    set(h1, 'markerfacecolor', get(h1, 'color')); 
    Leg_theta{iTheta} = strcat('$\theta_',num2str(iTheta),'$');
    h2 = semilogx(coeffs(ind_comp(iTheta)),ICOMP_min(iTheta),'ok','MarkerSize',10); hold on;
    set(h2, 'markerfacecolor', get(h2, 'color')); 
end
xlabel('$\lambda$');
ylabel('ICOMP')

tikzName = [folderName,'/',figName,'.tikz'];
        cleanfigure;
        matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '7cm', 'width','12cm','checkForUpdates',false);
        
%% Covariance interndependence of variables
figName = 'Covariance compplexity';
figure('Name',figName,'NumberTitle','off');
subplot(1,2,1);
for iTheta = 1:finalTerm
    h1 = semilogx(coeffs,Cov_cmplx{iTheta},'-o','MarkerSize',5); hold on;
    set(h1, 'markerfacecolor', get(h1, 'color')); 
    Leg_theta{iTheta} = strcat('$\theta_',num2str(iTheta),'$');
end
xlabel('$\lambda$');
ylabel('Covariance complexity')
title('Normalised data')
subplot(1,2,2);
for iTheta = 1:finalTerm
    h1 = semilogx(coeffs,Cov_cmplx_raw{iTheta},'-o','MarkerSize',5); hold on;
    set(h1, 'markerfacecolor', get(h1, 'color')); 
    Leg_theta{iTheta} = strcat('$\theta_',num2str(iTheta),'$');
end
xlabel('$\lambda$');
ylabel('Covariance complexity')
title('Raw data')

