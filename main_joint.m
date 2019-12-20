my_init;

foamset = questdlg('Select data folder', ...
    'Data to process',...
	'foam_2010','foam_2019','');
switch foamset
    case 'foam_2010'
        dataset = questdlg('Select data set', ...
        'Choice of set',...
        'C','D','');
    case 'foam_2019'
        dataset = questdlg('Select data set', ...
        'Choice of set',...
        'S','Y','Z','');
endmetaFileName = ['Meta_',dataset];
load(metaFileName);
d           = n_y + n_u;                                                    % size of input vector x
dict_set    = ['dict_',dataset];
fileNames   = sym(dict_set,[1 K]);                                          % vector of filenames
if normC == 1
    folder = 'Results';                                                     % specify category where to save files
else
    folder = 'Results_norm';
end
T = 2000;
names = {'set','ny','nu'};                                                  % names used to define results folder name (no more than 3).
folderName = make_folder(folder,names,dataset,n_y,n_u);                     % create results folder
%% Select first significant basis vector for all datasets
index = (1:T);                                                              % length of the sample
Files =  [1 2 4 5 6 7 9 10]; % 1:K; %                                       % id of the sample
K = length(Files);
iTerm = 1;                                                                  % the first significant term
AEER{iTerm} = zeros(nTerms,1);                                              % placeholder for AERR criteria
for iFile=Files                                                             % over all datasets
    fName = [dictFolder,'/',char(fileNames(iFile))];
    File  = matfile(fName,'Writable',true);
    residual_init{iFile} =  File.y_narx(index,1);                           % initial residual
    for jTerm = dict_terms                                                  % over all polynomial terms in the dictionary
        term0 = File.term(index,jTerm);
        cf(iFile,jTerm) = cor_sqr(residual_init{iFile},term0);              % squared correlation coefficient for the dataset and the polynomial term
        AEER{iTerm}(jTerm) = AEER{iTerm}(jTerm) + cf(iFile,jTerm);          % Average error reduction ration over all datasets
        clear term0
    end
    clear File
end
AEER{iTerm}(:,:) = AEER{iTerm}(:,:)/K;
[AEER_m,iMax] = max(AEER{iTerm});                                           % find the index of the term with the highest criterion across all datasets
AEER_mm(iTerm,1) = AEER_m;
S(iTerm) = iMax;                                                            % save index of the term
dict_terms(iMax) = [];                                                      % reduce the dictionary of available terms
AMDL_sum = 0;
for iFile=Files                                                             % over all datasets
    fName = [dictFolder,'/',char(fileNames(iFile))];
    File  = matfile(fName,'Writable',true);
    alpha{iFile}(:,iTerm)    = File.term(index,iMax);                       % the corresponding basis candidate term    
    phi  {iFile}(:,iTerm)    = File.term(index,iMax);                       % the corresponding basis vector 
    residual{iFile}(:,iTerm) = residual_update(residual_init{iFile},...     % the corresponding model residual
                                               phi{iFile}(:,iTerm));                                                        
    AMDL_sum = AMDL_sum + AMDL(residual{iFile}(:,iTerm),nNarx,iTerm);       % AMDL for the iFile dataset 
    clear File
end
AAMDL_all(iTerm) = AMDL_sum/K;                                              % average AMDL over all sets
significant_term{iTerm} =  symb_term{S(iTerm)};
disp(['Significant term ', num2str(iTerm),':'])
significant_term{iTerm}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main loop - forward selection

if length(dict_terms) + 1 < 30                                              % Maximum significant terms (if the algorithm is not terminated by the criterion)
    maxSign = length(dict_terms) + 1;
else
    maxSign     = 30;                                                           
end
converged   = false;
iTerm       = 2;
while(iTerm <= maxSign) && ~converged                                       % loop over the number of significant terms
    AEER{iTerm} = zeros(nTerms,1);                                          % placeholder for AERR criteria
    for iFile=Files                                                         % over all datasets
        fName = [dictFolder,'/',char(fileNames(iFile))];
        File  = matfile(fName,'Writable',true);
        for jTerm = dict_terms                                              % over all polynomial terms in the dictionary
            p{iTerm,iFile}(:,jTerm) = orthogonalise(File.term(index,jTerm),...
                                                    phi{iFile},iTerm);      % orthogonalise basis
            cf(iFile,jTerm)         = cor_sqr(residual_init{iFile},...
                                              p{iTerm,iFile}(:,jTerm));     % squared correlation coefficient for the dataset and the polynomial term
            AEER{iTerm}(jTerm) = AEER{iTerm}(jTerm) + cf(iFile,jTerm);      % Average error reduction ration over all datasets
        end
        clear File
    end
    AEER{iTerm}(:,:) = AEER{iTerm}(:,:)/K;
    [AEER_m,iMax] = max(AEER{iTerm});                                       % Find the index of the term with the highest criterion across all datasets
    AEER_mm(iTerm,1)   = AEER_m;
    S(iTerm) = iMax;                                                        % Save index of the term  
    ind = find(dict_terms == iMax);
    dict_terms(ind) = [];                                                   % Reduce the dictionary of available terms
    AMDL_sum = 0;
    for iFile=Files
        fName = [dictFolder,'/',char(fileNames(iFile))];
        File = matfile(fName,'Writable',true);
        alpha{iFile}(:,iTerm) = File.term(index,S(iTerm));                  % the corresponding basis candidate term    
        phi{iFile}(:,iTerm)   = p{iTerm,iFile}(index,S(iTerm));             % the corresponding basis vector 
        residual{iFile}(:,iTerm) = residual_update(residual{iFile}(:,iTerm-1),...
                                                   phi{iFile}(:,iTerm));    % the corresponding model residual                                 
        AMDL_sum = AMDL_sum + AMDL(residual{iFile}(:,iTerm),nNarx,iTerm);   % AMDL for the iFile dataset
        clear File
    end
    significant_term{iTerm} = symb_term{S(iTerm)};
    disp(['Significant term ', num2str(iTerm),':'])
    significant_term{iTerm}
    AAMDL_all(iTerm) = AMDL_sum/K;                                          % average AMDL over all sets
    converged = (AAMDL_all(iTerm) < -10);                                   % check convergence PLACEHOLDER
    iTerm = iTerm + 1;                                                      % increase the number of significant terms
end
clear AEER
[min_aamdl,i_min] = min(AAMDL_all);
finalTerm = i_min;

%% Direct estimation of polynomial coefficients
times = [1:T];
Phi_bar = [];
Y_bar   = [];
for iFile = Files
    fName = [dictFolder,'/',char(fileNames(iFile))];
    File  = matfile(fName,'Writable',true);
    indSign   = S(1:finalTerm);                                             % select the indeces of significant terms from the ordered set
    Phi_file  = File.term(times,:);                                         % extract all terms into a vector - cannot reorder directly in files
    y_file    = File.y_narx(times,:);                                       % extract output
    Phi       = Phi_file(times,indSign);                                    % select only signficant terms
    Phi_bar   = blkdiag(Phi_bar,Phi);                                       % diagonal block, size TKxNK
    Y_bar     = [Y_bar; y_file];                                            % block vector, size TKx1
end
% Crate the matrix of external design parameters
load('External_parameters');
L_cut_all = [values{1}(:, 9);values{2}(:, 9)];
D_rlx_all = [values{1}(:,11);values{2}(:,11)];
A_imp_all = [values{1}(:, 6);values{2}(:, 6)];
V_imp_all = [values{1}(:, 7);values{2}(:, 7)];
x = L_cut_all(Files,1);
y = D_rlx_all(Files,1);
id = ones(size(x));
% Matrix A should have dimension KxL where L < N
A = [id x y x.*y x.^2 y.^2];                                                % rectangular matrix, size KxL
I = eye(finalTerm);                                                         % unit matrix, size NxN
% \Theta(NxK) = B(NxL) A'(LxK)), vectorise using Kronecker product
Kr = kron(A,I);
M  = Phi_bar*Kr;                                                            % LS matrix - increased dimension does not guarantee increase in rank
B_bar = M\Y_bar;
Y_hat = M*B_bar;
R_2 = r_squared(Y_bar,Y_hat);
L = size(A,2);
Betas = reshape(B_bar,[finalTerm,L]);
%% Create column of names
finalTerm = i_min;
for iTerm=1:finalTerm
     temp = arrayfun(@char, significant_term{iTerm}, 'uniform', 0);
    if length(temp) > 0
        str = temp{1};
        for iString=2:length(temp)
            str = [str,',',temp{iString}];
        end
    end
    Terms{iTerm,1} = strcat('$',str,'$');
    clear temp
end
%% Store to table
Step = [1:finalTerm]';
Tab = table(Step,Terms);
for iBeta=1:L
    Parameters = round(Betas(:,iBeta),2);
    varName = ['$\beta_{',num2str(iBeta-1),'}$'];
    Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
end
Table_all = Tab
tableName = [folderName,'/Betas_direct_T_',num2str(T)];
table2latex(Table_all,tableName);
%% Validate the model
testFiles = [3 8];
for iFile = testFiles
x  = L_cut_all(iFile,1);
y  = D_rlx_all(iFile,1);
id = ones(size(x));
% Matrix A should have dimension KxL where L < N
A_k         = [id x y x.*y x.^2 y.^2]';                                             % vector
% fileName  = ['Dict_',dataset,num2str(iFile)];
% File      = matfile(fileName,'Writable',true);
fName       = [dictFolder,'/',char(fileNames(iFile))];
File        = matfile(fName,'Writable',true);
indSign     = S(1:finalTerm);                                                 % select the indeces of significant terms from the ordered set
Phi_file    = File.term(times,:);                                             % extract all terms into a vector - cannot reorder directly in files
Phi         = Phi_file(times,indSign);                                        % select only signficant terms
y_file      = File.y_narx(times,:);                                           % extract output
y_model     = Phi*Betas*A_k;                                                  % model NARMAX output
RMSE(iFile) = sqrt(mean((y_file - y_model).^2));                            % Root Mean Squared Error
%% Compare outputs
indPlot = [1:500];
figure('Name','Modelled output','NumberTitle','off');
colormap(my_map);
plot(indPlot+File.t_0,File.y_narx(indPlot,1)); hold on;
plot(indPlot+File.t_0,y_model(indPlot,1),'--'); hold on;
legend('True output','Generated output');

tikzName = [folderName,'/',dataset,num2str(iFile),'_gen_y_direct_T_',num2str(T),'.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '6cm', 'width','12cm','checkForUpdates',false);
        
clear File Phi_all Phi y_model
end

%% Ridge estimation
% Normalise M
% Normalise data matrix
A_2 = [x y];
for j=1:2 %size(A_2,2)
    m(j) = mean(A_2(:,j));
    c(j) = sqrt(sum((A_2(:,j) - m(j)).^2));
    d_norm(:,j) = (A_2(:,j)- m(j))/c(j);
    d_scale(:,j) = (A_2(:,j))/c(j);
end
A_n = [id d_norm d_norm(:,1).*d_norm(:,2) d_norm(:,1).^2 d_norm(:,2).^2];                                               % create the matrix for ls
A_s = [id d_scale d_scale(:,1).*d_scale(:,2) d_scale(:,1).^2 d_scale(:,2).^2]; 
A_norm = [d_norm d_norm(:,1).*d_norm(:,2) d_norm(:,1).^2 d_norm(:,2).^2];  
n = size(A,1)
R_a2 = A_norm'*A_norm;
R_22 = R_a2(2:end,2:end);
O = zeros(size(A_norm,2),1);
for j = 2:size(A_norm,2)
r(j-1,1) = A_norm(:,1)'*A_norm(:,j);
end
R_new = [1 r'; r R_22];
R_norm = [n O'; O R_new];
R_nn = A_n'*A_n;
Kr_norm = kron(A,I);
M_norm  = Phi_bar*Kr_norm;                                                            % LS matrix - increased dimension does not guarantee increase in rank
R_mm    = M_norm'*M_norm;
log_max = 0;
log_min = -6;
vec = [0:1/50:1];
coeffs = sort(10.^(log_min + (log_max-log_min)*vec));
I_k = eye(size(R_mm));
ik = 0;
for k = coeffs
    ik = ik + 1;
%     betas_rls(:,ik)     = pinv(R_norm + k*I_k)*M_norm'*Y_bar;
    betas_rls_raw(:,ik) = inv(R_mm + k*I_k)*M'*Y_bar;
%             betas_lasso{iTheta}(:,ik) =  LassoShooting(A_s,B(:,iTheta),k,'verbose',0);
%             betas_lasso_raw{iTheta}(:,ik) = LassoShooting(A,B(:,iTheta),k,'verbose',0);
end
% Betas_rls = reshape(betas_rls_raw,[finalTerm,L]);


figure;
for ib = 1:6
            h1 = semilogx(coeffs,betas_rls_raw(ib,:),'-o','MarkerSize',5); hold on;
            set(h1, 'markerfacecolor', get(h1, 'color')); 
end
title('Tikhonov normalised')
xlabel('$\gamma$');
ylabel('Standardised estimate')

