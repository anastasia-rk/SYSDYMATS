my_init
%% parameters from foam report
% Dictionary
params = {'D_init','D_imposed','L_init','L_imposed','R_imposed',...
         'A_imposed','V_imposed','D_relax','L_cut','W_relax','Dens_relax'};
L = numel(params);
values{1}(:,1)  = [30,30,30,30,30]';
values{2}(:,1)  = [50,50,50,50,50]';
values{1}(:,2)  = [19,19,19,19,19]';
values{2}(:,2)  = [19,19,19,19,19]';
values{1}(:,3)  = [170,150,130,110,90]';
values{2}(:,3)  = [160,140,120,100,90]';
values{1}(:,4)  = [80,80,80,80,80]';
values{2}(:,4)  = [75,75,75,75,75]';
values{1}(:,5)  = [1.58,1.58,1.58,1.58,1.58]';
values{2}(:,5)  = [2.63,2.63,2.63,2.63,2.63]';
values{1}(:,6)  = [2.13,1.88,1.63,1.38,1.13]';
values{2}(:,6)  = [2.13,1.87,1.60,1.33,1.20]';
values{1}(:,7)  = [5.30,4.67,4.05,3.43,2.80]';
values{2}(:,7)  = [14.8,12.9,11.1,9.2,8.3]';
values{1}(:,8)  = [20,20,20,20,20]';
values{2}(:,8)  = [20,20,20,20,20]';
values{1}(:,9)  = [74,72,69,61,56]';
values{2}(:,9)  = [67,66,64,62,57]';
values{1}(:,10) = [2.87,2.56,2.26,1.83,1.64]';
values{2}(:,10) = [5.8,5.28,4.48,4.2,3.45]';
values{1}(:,11) = [0.123,0.113,0.104,0.095,0.093]';
values{2}(:,11) = [0.276,0.255,0.230,0.209,0.193]';
%% total number of permutations
M       = permn(1:L,2);                                                     % get all permutations with repetition
ind     = find(M(:,2)>=M(:,1));                                             % sort out only increasing indeces
indeces = M(ind,:);                                                         % updated set
K       = length(indeces)                                                   % total number of sample pairs
clear ind
%% Check correlations of sample pairs
check_var = 0/0;
for i = 1:2
    % for all pairs
    for k = 1:K   
        % perform correlation test (possible types: Kendall, Soearman, Pearson)
        [rho{i}(k),pval{i}(k)]= corr(values{i}(:,indeces(k,1)),values{i}(:,indeces(k,2)));
    end
    % Find all pairs that don't reject non-correlation hypothesis
    ind             = find(~isnan(pval{i}) & pval{i}>0.0005);
    % update arrays
    rho{i}          = rho{i}(ind);
    pval{i}         = pval{i}(ind);    
    indeces_new{i}  = indeces(ind,:);
    [rho_max,k_max] = min(rho{i});
    param_select{i} = params(indeces_new{i}(k_max,:));
end

%% Load estimated thetas
dataset = 'C'; % 'D';
n_y = 0;
n_u = 4;
fileName = ['OLS_results_',dataset,'_ny_',num2str(n_y),'_nu_',num2str(n_u),'.mat'];
load(fileName);

%% Create custom fit function and estimate coefficients
% Files = [1 2 4 5]; % smaller dateset size
L_cut_all = [values{1}(:, 9);values{2}(:, 9)];
L_cut = L_cut_all(Files,1);
D_rlx_all = [values{1}(:,11);values{2}(:,11)];
D_rlx = D_rlx_all(Files,1);
A_imp_all = [values{1}(:, 6);values{2}(:, 6)];
A_imp = A_imp_all(Files,1);
V_imp_all = [values{1}(:, 7);values{2}(:, 7)];
V_imp = V_imp_all(Files,1);
clear x y
fo = fitoptions('Method','NonlinearLeastSquares');
g = fittype(@(b0,b1,b2,b3,b4,b5,x,y) b0 + b1*x + b2*y + b3*x.^2 + b4*y.^2 + b5*x.*y, 'independent',{'x','y'},'dependent','z','options',fo);
for iTerm = 1:finalTerm
z = Theta(iTerm,Files)';
[ft{iTerm},gof{iTerm},outp{iTerm}]= fit([A_imp,V_imp],z,g);
cfs(iTerm,:) = coeffvalues(ft{iTerm});
end
%% Save to table
Tab = table(Terms);
for iCoeff = 1:size(cfs,2)
    Parameters = round(cfs(:,iCoeff),3);
    varName = ['$b_',num2str(iCoeff-1),'$'];
    Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
end
Table_coeffs = Tab

tableName = ['Betas_',dataset,'_ny_',num2str(n_y),'_nu_',num2str(n_u),'_size_',num2str(T)];
table2latex(Table_coeffs,tableName);

%% Compute parameters for validation
iFile = 3;
A_test = A_imp_all(iFile,1);
V_test = V_imp_all(iFile,1);
for iTerm = 1:finalTerm
    theta_test(iTerm,1) = g(cfs(iTerm,1),cfs(iTerm,2),cfs(iTerm,3),cfs(iTerm,4),cfs(iTerm,5),cfs(iTerm,6),A_test,V_test);
end
fileName = ['Dict_',dataset,num2str(iFile)];
File = matfile(fileName,'Writable',true);


%%

% for iTerm = 1:finalTerm
% z = Theta(iTerm,Files)';
% [ft{iTerm},gof{iTerm},outp{iTerm}]= fit([L_cut,D_rlx],z,g);
% end