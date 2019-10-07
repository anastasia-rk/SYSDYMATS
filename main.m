my_init;

dataset =  'C'; %  'D'; %                                                   % name of dataset
metaFileName = ['Meta_',dataset];
load(metaFileName);
d       = n_y + n_u;                                                        % size of input vector x
dict_set = ['dict_',dataset];
fileNames = sym(dict_set,[1 K]);                                            % vector of filenames
T = 4000;
%% Select first significant basis vector for all datasets
index = (1:T);                                                              % length of the sample
Files = [1 2 4 5 6 7 9 10]; % 1:K                                           % id of the sample
iTerm = 1;                                                                  % the first significant term
AEER{iTerm} = zeros(nTerms,1);                                              % placeholder for AERR criteria
for iFile=Files                                                             % over all datasets
    File = matfile(char(fileNames(iFile)),'Writable',true);
    residual_init{iFile} =  File.y_narx(index,1);                           % initial residual
    for jTerm = dict_terms                                                  % over all polynomial terms in the dictionary
        term0 = File.term(index,jTerm);
        cf(iFile,jTerm) = cor_sqr(residual_init{iFile},term0);              % squared correlation coefficient for the dataset and the polynomial term
        AEER{iTerm}(jTerm) = AEER{iTerm}(jTerm) + cf(iFile,jTerm);          % Average error reduction ration over all datasets
        clear term0
    end
    clear File
end
[AEER_m,iMax] = max(AEER{iTerm});                                         % find the index of the term with the highest criterion across all datasets
AEER_mm(iTerm,1) = AEER_m;
S(iTerm) = iMax;                                                            % save index of the term
dict_terms(iMax) = [];                                                      % reduce the dictionary of available terms
AMDL_sum = 0;
for iFile=Files                                                             % over all datasets
    File = matfile(char(fileNames(iFile)),'Writable',true);
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
        File = matfile(char(fileNames(iFile)),'Writable',true);
        for jTerm = dict_terms                                              % over all polynomial terms in the dictionary
            p{iTerm,iFile}(:,jTerm) = orthogonalise(File.term(index,jTerm),...
                                                    phi{iFile},iTerm);      % orthogonalise basis
            cf(iFile,jTerm)         = cor_sqr(residual_init{iFile},...
                                              p{iTerm,iFile}(:,jTerm));     % squared correlation coefficient for the dataset and the polynomial term
            AEER{iTerm}(jTerm) = AEER{iTerm}(jTerm) + cf(iFile,jTerm);      % Average error reduction ration over all datasets
        end
        clear File
    end
    [AEER_m,iMax] = max(AEER{iTerm});                                       % Find the index of the term with the highest criterion across all datasets
    AEER_mm(iTerm,1)   = AEER_m;
    S(iTerm) = iMax;                                                        % Save index of the term  
    ind = find(dict_terms == iMax);
    dict_terms(ind) = [];                                                   % Reduce the dictionary of available terms
    AMDL_sum = 0;
    for iFile=Files
        File = matfile(char(fileNames(iFile)),'Writable',true);
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

%% Select optimal number of terms
[min_aamdl,i_min] = min(AAMDL_all);
figure('Name','AAMDL','NumberTitle','off');
plot(AAMDL_all,'o'); hold on;
plot(i_min,min_aamdl,'*','LineWidth',5);
xlim([1 maxSign])
xlabel('Number of terms');
ylabel('AAMDL')
%%
tikzName = ['AAMDL_',dataset,'_ny_',num2str(n_y),'_nu_',num2str(n_u),'_size_',num2str(T),'.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '6cm', 'width','8cm','checkForUpdates',false);
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
Step = [1:finalTerm]';
Tab = table(Step,Terms);
%% Parameter estimation
for iFile=Files
    U{iFile} = zeros(finalTerm,finalTerm);                                  % placeholder for upper-trig unit matrix
    iTerm = 1;                                                              % for the first term
    g{iFile}(iTerm) = (residual_init{iFile}'*phi{iFile}(:,iTerm))/...
                      (phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));           % top raw of the rotation matrix
    for jTerm =iTerm:finalTerm
        U{iFile}(iTerm,jTerm) = alpha{iFile}(:,jTerm)'*phi{iFile}(:,iTerm)/...
                                (phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));
    end
    for iTerm = 2:finalTerm                                                 % loop over significant terms
        g{iFile}(iTerm,1) = (residual{iFile}(:,iTerm-1)'*phi{iFile}...
                    (:,iTerm))/(phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));  % righthand side (normalised)
        for jTerm =iTerm:finalTerm
            U{iFile}(iTerm,jTerm) = alpha{iFile}(:,jTerm)'*phi{iFile}...
                     (:,iTerm)/(phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));  % upper triangular unit matrix
        end
    end
    Dets = det(U{iFile});
    Theta(:,iFile) = linsolve(U{iFile},g{iFile},struct('UT', true));        % solve upper triangular system via backward substitution
    Parameters = round(Theta(:,iFile),2);
    varName = [dataset,num2str(iFile)];
    Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
end
%% Store to table
temp      = AAMDL_all(1:finalTerm)';
AAMDL     = round(temp,3);
AEER  = round(AEER_mm(1:finalTerm,1),5);
Table_all = addvars(Tab,AEER,AAMDL,'NewVariableNames',{'AEER','AAMDL'})
tableName = ['Thetas_',dataset,'_ny_',num2str(n_y),'_nu_',num2str(n_u),'_size_',num2str(T)];
table2latex(Table_all,tableName);
%% 

%% Store data in table
workspaceName = ['OLS_results_',dataset,'_ny_',num2str(n_y),'_nu_',num2str(n_u),'.mat'];
save(workspaceName);