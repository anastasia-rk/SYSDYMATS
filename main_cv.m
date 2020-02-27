local_init;
%%
foamset = questdlg('Select data folder', ...
    'Data to process',...
	'foam_2010','foam_2019','vdpo','');
switch foamset
    case 'foam_2010'
        dataset = questdlg('Select data set', ...
        'Choice of set',...
        'C','D','');
    case 'foam_2019'
        dataset = questdlg('Select data set', ...
        'Choice of set',...
        'S','Y','Z','');
    case 'vdpo'
        dataset = 'V';   
end
folder = 'results_cv';                                                     % specify category where to save files
dFolder = 'dictionaries';
regressors = questdlg('Select the domain of regressors', ...
    'Domain choice',...
	'Shift','Delta','');
switch regressors
    case 'Shift'
        metaFileName = ['Meta_',dataset];
        load(metaFileName);
        names = {'set','ny','nu'};                                          % names used to define results folder name (no more than 3).
        if normC ~= 1
            folder = [folder,'_norm'];
            dFolder = [dFolder,'_norm'];
        end
        folderName = make_folder(folder,names,dataset,n_y,n_u);             % create results folder
        dictFolder = make_folder(dFolder,names,dataset,n_y,n_u);            % create results folder
        d       = n_y + n_u;                                                % size of input vector x
    case 'Delta'
        metaFileName = ['Meta_delta_',dataset];
        load(metaFileName);
         folder = ['delta_',folder];
         dFolder = ['delta_',dFolder];
         if normC ~= 1
            folder = [folder,'_norm'];
            dFolder = [dFolder,'_norm'];
         end
        direction = questdlg('Type of delta operator', ...
        'Causality',...
        'Forward','Backward','');
        switch direction
            case 'Backward'
                folder = [folder,'_b'];
                dFolder = [dFolder,'_b'];
            case 'Forward'
                folder = [folder,'_f'];
                dFolder = [dFolder,'_f'];
        end
         names = {'set','lambda'};
         folderName = make_folder(folder,names,dataset,lambda);             % create results folder
         dictFolder = make_folder(dFolder,names,dataset,lambda);            % create results folder
         n_u = 1+lambda;
         n_y = 1;
         d = lambda*2;
end
dict_set = ['dict_',dataset];                                   
fileNames = sym(dict_set,[1 K]);                                            % vector of filenames
Files =  1:K;                                                               % ids of the sample files
% Set maximum number of covariates
if length(dict_terms) + 1 < 30                                              % Maximum significant terms (if the algorithm is not terminated by the criterion)
    maxSign = length(dict_terms) + 1;
else
    maxSign = 30;                                                           
end
dict_terms_all = dict_terms;
%% Block CV with K-folds
nFolds = 10;
times = 1:nNarx;
cvpart = cvpartition(nNarx,'kFold',nFolds);
for iFold = 1:nFolds
    timesTrain{iFold} = times(cvpart.training(iFold));
    timesTest{iFold}  = times(cvpart.test(iFold));
    % block holdout - for time series
%     timesTrain{iFold} = [1:iFold*blockLength];
%     timesTest{iFold} = [iFold*blockLength+1:(iFold+1)*blockLength];
end
    
%% Dataset separation plot
figure;
for iFold=1:nFolds
    plot(timesTrain{iFold},iFold*ones(size(timesTrain{iFold})),'o','Color',my_map(110,:),'Linewidth',1); hold on;
    plot( timesTest{iFold},iFold* ones(size(timesTest{iFold})),'o','Color',my_map(250,:),'Linewidth',2); hold on;
end
colormap(my_map);
xlabel('sample index'); ylabel('CV block index');
tikzName = [folderName,'/K_folds.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '6cm', 'width','6cm','checkForUpdates',false);

%% Create the regression matrix based on the dataset (does not depend on CV parameters)
load(['External_parameters_',dataset]);
x = values(Files,1);
% x_m = mean(x);
% x_s = sqrt((x-x_m)'*(x-x_m));
% x   = x/x_s;
if size(values,2) > 1
    y = values(Files,2);
%     y_m = mean(y);
%     y_s = sqrt((y-y_m)'*(y-y_m));
%     y   = y/y_s;
else 
    y = [];
end
A = ones(size(x));                                                          % create unit vector for constants 
A_symb{1} = '1'
if ~isempty(y)                                                              % unknown mapping is a surface
   powers = permn(0:2,2);                                                   % permuntations of all 
   powers = powers(2:end,:);    
   nCols = min(size(powers,1),K);                                           % number of terms in the model shouldn't be higher then K
   for iCol = 1:nCols
       xCol = x.^powers(iCol,1);
       yCol = y.^powers(iCol,2);
       A = [A xCol.*yCol];
       A_symb{iCol+1} = ['$x^',num2str(powers(iCol,1)),'$ $y^',num2str(powers(iCol,2)),'$']; 
   end
else                                                                        % unknown mapping is a curve
    nCols =  min(3,K);                                                       % number of terms in the model shouldn't be higher then K
    for iCol = 1:nCols                                                      % limit order of the model by the number of experimants
       A = [A x.^(iCol)];
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Structure identification with CV model selection
    index       = times;                                                    % structure id over all times
    dict_terms  = dict_terms_all;                                           % refill the dictionary
    %% Select first significant basis vector for all datasets
    iTerm = 1;                                                              % the first significant term
    AERR{iTerm} = zeros(nTerms,1);                                          % placeholder for AERR criteria
    for iFile = Files                                                       % over all datasets
        fName   = [dictFolder,'/',char(fileNames(iFile))];
        File    = matfile(fName,'Writable',true);
        y_n     = File.y_narx(:,1);
        residual_init{iFile} =  y_n(index,1);                               % initial residual
        for jTerm = dict_terms                                              % over all polynomial terms in the dictionary
            term_all    = File.term(:,jTerm);
            term0       = term_all(index,:);
            cf(iFile,jTerm)     = cor_sqr(residual_init{iFile},term0);      % squared correlation coefficient for the dataset and the polynomial term
            AERR{iTerm}(jTerm)  = AERR{iTerm}(jTerm) + cf(iFile,jTerm);     % Average error reduction ration over all datasets
            clear term0 term_all
        end
        clear File y_n 
    end
    AERR{iTerm}(:,:)    = AERR{iTerm}(:,:)/K;
    [AERR_m,iMax]       = max(AERR{iTerm});                                 % find the index of the term with the highest criterion across all datasets
    AERR_mm(iTerm,1)    = AERR_m;
    S(iTerm)            = iMax;                                             % save index of the term
    dict_terms(iMax)    = [];                                               % reduce the dictionary of available term
    BIC_sum             = 0;
    for iFile = Files                                                       % over all datasets
        fName = [dictFolder,'/',char(fileNames(iFile))];
        File  = matfile(fName,'Writable',true);
        term_all    = File.term(:,iMax);
        alpha{iFile}(:,iTerm)    = term_all(index,:);                       % the corresponding basis candidate term    
        phi  {iFile}(:,iTerm)    = term_all(index,:);                       % the corresponding basis vector 
        residual{iFile}(:,iTerm) = residual_update(residual_init{iFile},... % the corresponding model residual
                                               phi{iFile}(:,iTerm));                                                        
        BIC_sum  = BIC_sum  +  BIC(residual{iFile}(:,iTerm),nNarx,iTerm);   % BIC for the iFile dataset
        clear File term_all
    end
    BIC_all(iFold,iTerm)            = BIC_sum/K;                            % average AMDL over all sets
    significant_term{iFold,iTerm}   = symb_term{S(iTerm)};
%% Main loop   
    converged   = false;
    bics        = [];
    while(iTerm < maxSign) %&& ~converged                                   % loop over the number of significant terms
        iTerm = iTerm + 1;                                                  % increase the number of significant terms
        AERR{iTerm} = zeros(nTerms,1);                                      % placeholder for AERR criteria
        for iFile = Files                                                   % over all datasets
            fName   = [dictFolder,'/',char(fileNames(iFile))];
            File    = matfile(fName,'Writable',true);
            for jTerm = dict_terms                                          % over all polynomial terms in the dictionary
                term_all    = File.term(:,jTerm);
                p{iTerm,iFile}(:,jTerm) = orthogonalise(term_all(index,:),...
                                                    phi{iFile},iTerm);      % orthogonalise basis
                cf(iFile,jTerm)         = cor_sqr(residual_init{iFile},...
                                              p{iTerm,iFile}(:,jTerm));     % squared correlation coefficient for the dataset and the polynomial term
                AERR{iTerm}(jTerm) = AERR{iTerm}(jTerm) + cf(iFile,jTerm);  % average error reduction ration over all datasets
            end
            clear File
        end
        AERR{iTerm}(:,:)    = AERR{iTerm}(:,:)/K;
        [AERR_m,iMax]       = max(AERR{iTerm});                             % find the index of the term with the highest criterion across all datasets
        AERR_mm(iTerm,1)    = AERR_m;
        S(iTerm)            = iMax;                                         % save index of the term  
        ind = find(dict_terms == iMax);
        dict_terms(ind) = [];                                               % Reduce the dictionary of available terms
        BIC_sum         = 0;
        for iFile = Files
            fName   = [dictFolder,'/',char(fileNames(iFile))];
            File    = matfile(fName,'Writable',true);
            alpha{iFile}(:,iTerm) = File.term(index,S(iTerm));              % the corresponding basis candidate term    
            phi{iFile}(:,iTerm)   = p{iTerm,iFile}(:,S(iTerm));             % the corresponding basis vector 
            residual{iFile}(:,iTerm) = residual_update(residual{iFile}(:,iTerm-1),...
                                                   phi{iFile}(:,iTerm));    % the corresponding model residual                                 
            BIC_sum  = BIC_sum  +  BIC(residual{iFile}(:,iTerm),nNarx,iTerm); % BIC for the iFile dataset
            clear File x_n
        end
        significant_term{iFold,iTerm} = symb_term{S(iTerm)};
        BIC_all(iTerm) = BIC_sum/K;                                         % average AMDL over all sets
        converged_BIC = (abs((BIC_all(iTerm) - BIC_all(iTerm-1))/BIC_all(iTerm)) < 0.002); % check convergence
        if converged_BIC
            bics  = [bics,iTerm];
        end
    end
if isempty(bics)
    finalTerm = 20;
else
    finalTerm = bics(1);
end
    BIC_trunc = BIC_all(1:finalTerm);
figure;
plot([1:maxSign],BIC_all); hold on;
plot(finalTerm,BIC_trunc(end),'*');
xlabel('Terms');ylabel('BIC');

tikzName = [folderName,'/BIC_all_folds.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '6cm', 'width','6cm','checkForUpdates',false);
    
    %% Create column of names
    for iTerm=1:finalTerm
        temp = arrayfun(@char, significant_term{iFold,iTerm}, 'uniform', 0);
        if length(temp) > 0
            str = temp{1};
            for iString=2:length(temp)
                str = [str,temp{iString}];
            end
        end
        Terms{iTerm,1} = strcat('$',str,'$');
        clear temp
    end
    Step = [1:finalTerm]';
    Tab = table(Step,Terms);
    clear AERR
    AERR  = round(AERR_mm(1:finalTerm,1)*100,3);
%% Parameter estimation
    for iFile=Files
        U{iFile} = zeros(finalTerm,finalTerm);                              % placeholder for upper-trig unit matrix
        iTerm = 1;                                                          % for the first term
        g{iFile}(iTerm) = (residual_init{iFile}'*phi{iFile}(:,iTerm))/...
                          (phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));       % top raw of the rotation matrix
        for jTerm =iTerm:finalTerm
            U{iFile}(iTerm,jTerm) = alpha{iFile}(:,jTerm)'*phi{iFile}(:,iTerm)/...
                                (phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));
        end
        for iTerm = 2:finalTerm                                             % loop over significant terms
            g{iFile}(iTerm,1) = (residual{iFile}(:,iTerm-1)'*phi{iFile}...
                    (:,iTerm))/(phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));  % righthand side (normalised)
            for jTerm =iTerm:finalTerm
                U{iFile}(iTerm,jTerm) = alpha{iFile}(:,jTerm)'*phi{iFile}...
                     (:,iTerm)/(phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));  % upper triangular unit matrix
            end
        end
        Dets = det(U{iFile});
        Theta(:,iFile) = linsolve(U{iFile},g{iFile},struct('UT', true));    % solve upper triangular system via backward substitution
        Parameters = round(Theta(:,iFile),2);
        varName = [dataset,num2str(iFile)];
        Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
    end
%% Saving internal parameters to table
    Tab   = addvars(Tab,AERR,'NewVariableNames',{'AERR($\%$)'});
    internalParams = addvars(Tab,BIC_trunc,'NewVariableNames',{'BIC'})
    tableName = [folderName,'/Thetas_overall'];
    table2latex(internalParams,tableName);
    clear Tab tableName
    clear AERR alpha U phi residual p Theta g

%% Regularised LS (joint) + CV across all datasets
log_max = 0; log_min = -6; vec = [0:1/20:1];
Coeffs = sort(10.^(log_min + (log_max-log_min)*vec));
% vectorisation for joint estimation
I = eye(finalTerm);                                                         % unit matrix, size NxN
Kr = kron(A,I);
L = size(A,2);
iLambda = 0;
for lambda = Coeffs                                                         % across regiularisation coeffs
    iLambda   = iLambda + 1;
    for iFold = 1:nFolds                                                    % across folds
        index       = timesTrain{iFold};
        index_test  = timesTest{iFold};
%% Estimation procedure
        Phi_bar = [];
        Y_bar   = [];
        for iFile = Files
            fName = [dictFolder,'/',char(fileNames(iFile))];
            File  = matfile(fName,'Writable',true);
            indSign   = S(1:finalTerm);                                     % select the indeces of significant terms from the ordered set
            Phi_file  = File.term(:,:);                                     % extract all terms into a vector - cannot reorder directly in files
            y_file    = File.y_narx(:,:);                                   % extract output
            Phi       = Phi_file(index,indSign);                            % select only signficant terms
            Phi_bar   = blkdiag(Phi_bar,Phi);                               % diagonal block, size TKxNK
            y_n       = y_file(index,:);
            Y_bar     = [Y_bar; y_n];                                       % block vector, size TKx1
        end
        M  = Phi_bar*Kr;                                                    % LS matrix - increased dimension does not guarantee increase in rank
        R_mm  = M'*M;
        gain  = inv(R_mm + lambda*eye(size(R_mm)))*M';                      % RLS gain
        B_bar = gain*Y_bar;
        B_lasso =  LassoShooting(M,Y_bar,lambda,'verbose',0);
        Betas{iLambda,iFold} = reshape(B_bar,[finalTerm,L]);
        Betas_lasso{iLambda,iFold} = reshape(B_lasso,[finalTerm,L]);
        clear M R_mm gain 
%% Validation procedure  
        Theta_test{iFold}   = Betas{iLambda,iFold}*A';
        Theta_lasso{iFold}  = Betas_lasso{iLambda,iFold}*A';
%         iRSS    = 0;
        Phi_bar = [];
        Y_bar   = [];
        for iFile = Files
            fName = [dictFolder,'/',char(fileNames(iFile))];
            File  = matfile(fName,'Writable',true);
            indSign   = S(1:finalTerm);                                     % select the indeces of significant terms from the ordered set
            Phi_file  = File.term(:,:);                                     % extract all terms into a vector - cannot reorder directly in files
            Phi       = Phi_file(index_test,indSign);                       % select only signficant terms
            Phi_bar   = blkdiag(Phi_bar,Phi);                               % diagonal block, size TKxNK
            y_file    = File.y_narx(:,:);
            y_n       = y_file(index_test,:);
            Y_bar     = [Y_bar; y_n];                                       % block vector, size TKx1
        end
       M        = Phi_bar*Kr;
       Y_hat    = M*B_bar;
       Y_lasso  = M*B_lasso;
       R_mm     = M'*M;
       Hat_tikh = M*inv(R_mm + lambda*eye(size(R_mm)))*M';
%        Hat      = M*inv(R_mm)*M';
       nData(iFold)             = length(index_test);
       PE        = Y_bar - Y_hat;
       PE_lasso  = Y_bar - Y_lasso;
       RSS(iLambda,iFold)       = PE'*PE/nData(iFold);
       RSS_lasso(iLambda,iFold) = PE_lasso'*PE_lasso/nData(iFold);
       BIC_tikh(iLambda,iFold)  = BIC(PE,nData(iFold),length(B_bar));
       BIC_lasso(iLambda,iFold) = BIC(PE_lasso,nData(iFold),length(B_lasso));
       clear PE PE_lasso
%% Compute betas for each file separately
%         for iFile = Files
%             fName   = [dictFolder,'/Dict_',dataset,num2str(iFile)];
%             File    = matfile(fName,'Writable',true);
%             indSign = S(1:finalTerm);                                       % select the indeces of significant terms from the ordered set
%             Phi_all = File.term;                                            % extract all terms into a vector
%             Y_all   = File.y_narx;
%             Phi     = Phi_all(index_test,indSign);                          % select only signficant terms
%             y_model = Phi*Theta_test{iFold}(:,iFile);                       % model NARMAX output
%             y_lasso = Phi*Theta_lasso{iFold}(:,iFile);
%             iRSS   = iRSS + 1;
%             RSS{iLambda}(iFold,iRSS) = (Y_all(index_test,1) - y_model)'*...
%                                        (Y_all(index_test,1) - y_model);     % Root Mean Squared Error
%             RSS_lasso{iLambda}(iFold,iRSS) = (Y_all(index_test,1) - y_lasso)'*...
%                                         (Y_all(index_test,1) - y_lasso);    % Root Mean Squared Error           
%             clear Phi_all Phi Y_all
%         end
%          SRSS(iLambda,iFold) = sum(RSS{iLambda}(iFold,:));
%          SRSS_lasso(iLambda,iFold) = sum(RSS_lasso{iLambda}(iFold,:));
          
         
         % for BIC
%

%          nTerms  = length(B_bar);
%          nTerms_lasso  = length(B_lasso);
%          SRSS(iLambda,iFold) = sum(RSS{iLambda}(iFold,:));
%          BIC_CV(iLambda,iFold) = log(nData)*nTerms + nData*log(SRSS(iLambda,iFold)/nData);
%          SRSS_lasso(iLambda,iFold) = sum(RSS_lasso{iLambda}(iFold,:));
%          BIC_lasso(iLambda,iFold) = log(nData)*nTerms + nData*log(SRSS_lasso(iLambda,iFold)/nData);
    end
%     nAll    = sum(nData);
%     
%     BIC(iLambda)            = BIC(PE(iLambda,:),nAll,length(B_bar));
%     BIC_lasso(iLambda)      = BIC(PE_lasso(iLambda,:),nAll,length(B_lasso));
    PRESS(iLambda)          = mean(RSS(iLambda,:));                                    % %mean(MRMSE(iLambda,:));
    PRESS_lasso(iLambda)    = mean(RSS_lasso(iLambda,:));
end
Folds = [1:nFolds];
%% Plots
[plotFolds,plotCoeffs] = meshgrid(Folds,Coeffs);                              % create meshgrid
figure('Name','RSS all','NumberTitle','off');
plot3(plotCoeffs,plotFolds,RSS);
set(gca,'XScale','log')
xlabel('$\lambda$');ylabel('Folds');zlabel('BIC');
%% individual plots
% figure;
% subplot(1,2,1);
% for iFold = Folds
% plot(Coeffs,RSS(:,iFold));hold on;
% end
% set(gca,'XScale','log')
% xlabel('$\lambda$');ylabel('RSS');
% title('Tikhonov regularisation');
% subplot(1,2,2);
% for iFold = Folds
% plot(Coeffs,RSS_lasso(:,iFold));hold on;
% end
% set(gca,'XScale','log')
% xlabel('$\lambda$');ylabel('RSS');
% title('LASSO regularisation');
% tikzName = [folderName,'/RSS_all_folds.tikz'];
% cleanfigure;
% matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
%             false, 'height', '6cm', 'width','6cm','checkForUpdates',false);

[PRESS_min,i_min] = min(PRESS);
[PRESS_min_l,i_min_l] = min(PRESS_lasso);
figure('Name','PRESS','NumberTitle','off');
subplot(1,2,1);
plot(Coeffs,PRESS,'linewidth',2); hold on;
plot(Coeffs(i_min), PRESS_min, '*','linewidth',3);
set(gca,'XScale','log')
xlabel('$\lambda$');ylabel('PRESS');
title('Tikhonov regularisation');
subplot(1,2,2);
plot(Coeffs,PRESS_lasso,'linewidth',2); hold on;
plot(Coeffs(i_min_l), PRESS_min_l, '*','linewidth',3);
set(gca,'XScale','log')
xlabel('$\lambda$');ylabel('PRESS');
title('LASSO regularisation');
xlabel('$\lambda$');ylabel('PRESS');
tikzName = [folderName,'/PRESS_min.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '5cm', 'width','12cm','checkForUpdates',false);
%% Estimating with optimal regularisation coefficient
lambda_opt = Coeffs(i_min);
lambda_lasso_opt = Coeffs(i_min_l);
index       = times;
Phi_bar = [];
Y_bar   = [];
for iFile = Files
     fName = [dictFolder,'/',char(fileNames(iFile))];
     File  = matfile(fName,'Writable',true);
     indSign   = S(1:finalTerm);                                            % select the indeces of significant terms from the ordered set
     Phi_file  = File.term(:,:);                                            % extract all terms into a vector - cannot reorder directly in files
     y_file    = File.y_narx(:,:);                                          % extract output
     Phi       = Phi_file(index,indSign);                                   % select only signficant terms
     Phi_bar   = blkdiag(Phi_bar,Phi);                                      % diagonal block, size TKxNK
     y_n       = y_file(index,:);
     Y_bar     = [Y_bar; y_n];                                              % block vector, size TKx1
end
M       = Phi_bar*Kr;                                                            % LS matrix - increased dimension does not guarantee increase in rank
R_mm    = M'*M;
gain    = inv(R_mm + lambda_opt*eye(size(R_mm)))*M';                          % RLS gain
B_tikh  = gain*Y_bar;
B_bar   = M\Y_bar;
B_lasso =  LassoShooting(M,Y_bar,lambda_lasso_opt,'verbose',0);
L       = size(A,2);
Betas_opt       = reshape(B_bar,[finalTerm,L])
Betas_tikh_opt  = reshape(B_tikh,[finalTerm,L])
Betas_lasso_opt = reshape(B_lasso,[finalTerm,L])
%% Saving external parameters to table
Tab = table(Step,Terms);
for iBeta=1:L
    Parameters = round(Betas_tikh_opt(:,iBeta),2);
    varName = ['$\beta_{',num2str(iBeta-1),'}$'];
    Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
end
tableName = [folderName,'/Betas_tikhonov_',num2str(iFold)];
table2latex(Tab,tableName);
 clear Tab
Tab = table(Step,Terms);
for iBeta=1:L
    Parameters = round(Betas_tikh_opt(:,iBeta),2);
    varName = ['$\beta_{',num2str(iBeta-1),'}$'];
    Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
end
tableName = [folderName,'/Betas_ls_',num2str(iFold)];
table2latex(Tab,tableName);
%% Validate Tikhonov reg
testFiles   = Files;
Theta_test  = Betas_opt*A';
iRMSE       = 0;
figure('Name','Outputs','NumberTitle','off');
L2 = 4;
index_test  = 1000:3000;
for iFile = testFiles
    fName   = [dictFolder,'/Dict_',dataset,num2str(iFile)];
    File    = matfile(fName,'Writable',true);
    indSign = S(1:finalTerm);                                               % select the indeces of significant terms from the ordered set
    Phi_all = File.term(index_test,:);                                      % extract all terms into a vector
    Phi     = Phi_all(:,indSign);                                           % select only signficant terms
    y_model = Phi*Theta_test(:,iFile);                                      % model NARMAX output
    iRMSE   = iRMSE + 1;
    RMSE(iRMSE) = sqrt(mean((File.y_narx(index_test,1) - y_model).^2));     % Root Mean Squared Error
% Compare outputs
    subplot(L2,2,iFile);
    plot(index_test(1:500)+File.t_0,File.y_narx(index_test(1:500),1)); hold on;
    plot(index_test(1:500)+File.t_0,y_model(index_test(1:500),1),'--'); hold on;
    legend('True output','Generated output');
    xlabel('Sample index'); ylabel(['$',y_str,'$']);
    title(['File ',num2str(iFile),', RMSE = ',num2str(RMSE(iRMSE))]);
       
 clear File Phi_all Phi y_model
end
tikzName = [folderName,'/','Tikhonov_validation.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '4cm', 'width','15cm','checkForUpdates',false);
 
%% Validate lasso reg
Theta_test  = Betas_lasso_opt*A';
iRMSE       = 0;
figure('Name','Outputs','NumberTitle','off');
L = 4;
for iFile = testFiles
    fName   = [dictFolder,'/Dict_',dataset,num2str(iFile)];
    File    = matfile(fName,'Writable',true);
    indSign = S(1:finalTerm);                                               % select the indeces of significant terms from the ordered set
    Phi_all = File.term(index_test,:);                                      % extract all terms into a vector
    Phi     = Phi_all(:,indSign);                                           % select only signficant terms
    y_model = Phi*Theta_test(:,iFile);                                      % model NARMAX output
    iRMSE   = iRMSE + 1;
    RMSE(iRMSE) = sqrt(mean((File.y_narx(index_test,1) - y_model).^2));     % Root Mean Squared Error
% Compare outputs
    subplot(L2,2,iFile);
    plot(index_test(1:500),File.y_narx(index_test(1:500),1)); hold on;
    plot(index_test(1:500),y_model(index_test(1:500),1),'--'); hold on;
    legend('True output','Generated output');
    xlabel('Sample index'); ylabel(['$',y_str,'$']);
    title(['File ',num2str(iFile),', RMSE = ',num2str(RMSE(iRMSE))]);
 clear File Phi_all Phi y_model
end
tikzName = [folderName,'/','LASSO_validation.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '4cm', 'width','15cm','checkForUpdates',false);

