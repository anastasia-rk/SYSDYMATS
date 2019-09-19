my_init;

dataset =  'D'; %  'C'; %                                                   % name of dataset
metaFileName = ['Meta_',dataset];
load(metaFileName);
d       = n_y + n_u;                                                        % size of input vector x
dict_set = ['dict_',dataset];
fileNames = sym(dict_set,[1 K]);                                            % vector of filenames

%% Select first significant basis vector for all datasets
iFile = 1;                                                                  % id of the sample
iTerm = 1;                                                                  % the first significant term
index = (1:2000);
AEER{iTerm} = zeros(nTerms,1);                                              % placeholder for AERR criteria
for iFile=1:K                                                               % over all datasets
    File = matfile(char(fileNames(iFile)),'Writable',true);
    residual_init{iFile} =  File.y_narx(index,1);                                    % initial residual
    for jTerm = dict_terms                                                  % over all polynomial terms in the dictionary
        cf(iFile,jTerm) = cor_sqr(residual_init{iFile},File.term(index,jTerm)); % squared correlation coefficient for the dataset and the polynomial term
        AEER{iTerm}(jTerm) = AEER{iTerm}(jTerm) + cf(iFile,jTerm);          % Average error reduction ration over all datasets
    end
    clear File
end
[AEER_max,iMax] = max(AEER{iTerm});                                         % Find the index of the term with the highest criterion across all datasets
S(iTerm) = iMax;                                                            % Save index of the term
dict_terms(iMax) = [];                                                      % Reduce the dictionary of available terms
for iFile=1:K                                                               % over all datasets
    File = matfile(char(fileNames(iFile)),'Writable',true);
    alpha{iFile}(:,iTerm) = File.term(index,iMax);                              % the corresponding basis candidate term    
    phi{iFile}(:,iTerm)   = File.term(index,iMax);                              % The corresponding basis vector 
    residual{iFile}(:,iTerm) = residual_update(residual_init{iFile},...     % the corresponding model residual
                                   phi{iFile}(:,iTerm));
    clear File
end
significant_term{iTerm} = File.symb_term{S(iTerm)}
disp(['Significant term ', num2str(iTerm),':'])
significant_term{iTerm}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main loop - forward selection
maxSign     = 30;                                                           % Maximum significant terms (if the algorithm is not terminated by the criterion)
converged   = false;
iTerm       = 2;
while(iTerm < maxSign) && ~converged                                        % loop over the number of significant terms
    AEER{iTerm} = zeros(nTerms,1);                                          % placeholder for AERR criteria
    for iFile=1:K                                                           % over all datasets
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
    [AEER_max,iMax] = max(AEER{iTerm});                                     % Find the index of the term with the highest criterion across all datasets
    S(iTerm) = iMax;                                                        % Save index of the term  
    index = find(dict_terms == iMax);
    dict_terms(index) = [];                                                 % Reduce the dictionary of available terms
    AMDL_sum = 0;
    for iFile=1:K
        File = matfile(char(fileNames(iFile)),'Writable',true);
        alpha{iFile}(:,iTerm) = File.term(index,S(iTerm));                      % the corresponding basis candidate term    
        phi{iFile}(:,iTerm)   = p{iTerm,iFile}(index,S(iTerm));                 % the corresponding basis vector 
        residual{iFile}(:,iTerm) = residual_update(residual{iFile}(:,iTerm-1),...
                                                   phi{iFile}(:,iTerm));    % the corresponding model residual                                 
        AMDL_sum = AMDL_sum + AMDL(residual{iFile}(:,iTerm),nNarx,iTerm);   % AMDL for the iFile dataset
        clear File
    end
    significant_term{iTerm} = symb_term{S(iTerm)};
    disp(['Significant term ', num2str(iTerm),':'])
    significant_term{iTerm}
    AAMDL(iTerm) = AMDL_sum/K;                                              % average AMDL over all sets
    converged = (AAMDL(iTerm) < -10);                                       % check convergence PLACEHOLDER
    iTerm = iTerm + 1;                                                      % increase the number of significant terms
end
%%
figure('Name','AAMDL','NumberTitle','off'); 
plot(AAMDL(2:end),'o');

%% Parameter estimation
[min_aamdl,i_min] = min(AAMDL);
finalTerm = i_min;

for iFile=1:K
    U{iFile} = zeros(finalTerm,finalTerm);                                  % placeholder for upper-trig unit matrix
    iTerm = 1;                                                              % for the first term
    g{iFile}(iTerm) = (residual_init{iFile}'*phi{iFile}(:,iTerm))/...
                       norm(phi{iFile}(:,iTerm));                           % Top raw of the  
    for jTerm =iTerm:finalTerm
        U{iFile}(iTerm,jTerm) = alpha{iFile}(:,jTerm)'*phi{iFile}(:,iTerm)/(phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));
    end
    for iTerm = 2:finalTerm                                                 % loop over significant terms
        g{iFile}(iTerm,1) = (residual{iFile}(:,iTerm-1)'*phi{iFile}(:,iTerm))/(phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));
        for jTerm =iTerm:finalTerm
            U{iFile}(iTerm,jTerm) = alpha{iFile}(:,jTerm)'*phi{iFile}(:,iTerm)/norm(phi{iFile}(:,iTerm));
        end
    end
end