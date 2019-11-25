my_init
%% Load estimated thetas
dataset = 'C'; % 'D'; % 
metaFileName = ['Meta_',dataset];
load(metaFileName);
T = 2000;
folder = 'Results';                                                         % specify category where to save files
names = {'set','ny','nu'};                                                  % names used to define results folder name (no more than 3).
folderName = make_folder(folder,names,dataset,n_y,n_u);                     % create results folder
fileName = [folderName,'/OLS_results_T_',num2str(T),'.mat'];
load(fileName);
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
A = [id x y x.*y x.^2 y.^2];                                                % create the matrix for ls
names = sym('beta_',[1 size(A,2)]);
if length(x) <= 4
    A = A(:,1:4);
end
R_aa = A'*A;
clear beta
B = Theta(:,Files_sub)';
betas_ols = A\B;
%% Collinearity analysis
% Diagnose collinearity based on condition number of the scaled regression
% matrix. Small eigenvalues of it's square also indicate how many columns
% have near-linear dependency
% % check
for j=1:size(A,2)
    if j > 1
        m = mean(A(:,j));
    else
        m = 0;
    end
    c = sqrt(sum(A(:,j).^2));
    A_norm(:,j) = (A(:,j)- m)/c;
end
R_norm = A_norm'*A_norm;


A_n = normalize(A,'norm')
R_n = A_n'*A_n
eig_r = eig(R_n);
kappa = cond(A_n);
thr = 1000;
% if kappa > thr
%     print('Collinearity diagnosed');
%     print('Identifying variance inflation')
%     [R_A,P_A] = corrcoef(A_n);                                              % corr coeffs and p-values. If P_A elmement < 0.05 correlation is significant
%     for j=1:size(P_A,2)
%         corr_ind(:,j) = find(P_A(:,j) < 0.05);
%     end
%     for iTheta = 1:size(B,2)
%         M = [A B(:,iTheta)];
%         R{iTheta} = corrcoef(M);    
%     end
% end
% % check
% for j=1:size(A,2)
%     c = sqrt(sum(A(:,j).^2));
%     A_norm(:,j) = A(:,j)/c;
% end

%% Plot ridge traces
n = 18;
nPlots = round(size(A,2)/2);
log_max = 0;
log_min = -6;
vec = [0:1/50:1];
coeffs = sort(10.^(log_min + (log_max-log_min)*vec));
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

%% Ridged estimation (in closed form) 
I_k = eye(size(R_n));
for iTheta = 1:finalTerm
        ik = 0;
        for k = coeffs
            ik = ik + 1;
            beta{iTheta}(:,ik) = inv(R_norm + k*I_k)*A_norm'*B(:,iTheta);
        end
end
%% Plot ridge traces 
for iTheta = 1:finalTerm
        figName = char(symb_term{S(iTheta)});
        figure('Name',figName,'NumberTitle','off');
        for ib = 1:size(A,2)
            h1 = semilogx(coeffs,beta{iTheta}(ib,:),'-o','MarkerSize',5); hold on;
            set(h1, 'markerfacecolor', get(h1, 'color')); 
        end
        xlabel('$\gamma$');
        ylabel('Standarised estimate')
        legend(Legends);
end
    
%% Estimate with selected value of k

