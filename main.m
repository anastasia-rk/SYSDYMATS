my_init;

dataset = 'C'; % 'S'; % 'Z'; % 'Y'; % 'S'; %  'D'; %                        % name of dataset
metaFileName = ['Meta_',dataset];
load(metaFileName);
d       = n_y + n_u;                                                        % size of input vector x
T = 2000;
dict_set = ['dict_',dataset];                                   % 
fileNames = sym(dict_set,[1 K]);                                            % vector of filenames
folder = 'Results';                                                         % specify category where to save files
names = {'set','ny','nu'};                                                  % names used to define results folder name (no more than 3).
folderName = make_folder(folder,names,dataset,n_y,n_u);                     % create results folder
index = (1:T);                                                              % length of the sample
Files =  1:K; %                                                             % id of the sample
Files(3) = [];Files(8) = [];
K = length(Files);
%% Show correlation between autoregressive and input terms
if n_y ~= 0
figure;
L2 = round(K/2);
iPlot = 0;
for iFile=Files                                                             % over all datasets
        fName = [dictFolder,'/',char(fileNames(iFile))];
        File = matfile(fName,'Writable',true);
        ind_y= n_y;
        iPlot = iPlot + 1;
        tht = File.term(index,ind_y)\File.y_narx(index,1);
        y_regressed = File.term(index,ind_y)*tht;
        Rr = round(corrcoef(File.term(index,ind_y),File.y_narx(index,1)),4);
        subplot(L2,2,iPlot)
        scatter(File.term(index,ind_y),File.y_narx(index,1),'filled'); hold on;
        plot(File.term(index,ind_y),y_regressed,'Linewidth',2); 
        title([dataset,num2str(iFile),', R = ',num2str(Rr(2,1))]); xlabel(char(symb_term{ind_y})); ylabel('y(t)');
        clear File tht Rr
end
    tikzName = [folderName,'/Regressions_yt_T_',num2str(T),'.tikz'];
    cleanfigure;
    matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '18cm', 'width','12cm','checkForUpdates',false);
end
%%     
figure;
L2 = round(length(Files)/2);
iPlot = 0;
for iFile=Files                                                             % over all datasets
        fName = [dictFolder,'/',char(fileNames(iFile))];
        File = matfile(fName,'Writable',true);
        ind_u= n_y + n_u;
        iPlot = iPlot + 1;
        tht = File.term(index,ind_u)\File.y_narx(index,1);
        y_regressed = File.term(index,ind_u)*tht;
        Rr = round(corrcoef(File.term(index,ind_u),File.y_narx(index,1)),4);
        subplot(L2,2,iPlot)
        scatter(File.term(index,ind_u),File.y_narx(index,1),'filled'); hold on;
        plot(File.term(index,ind_u),y_regressed,'Linewidth',2); 
        title([dataset,num2str(iFile),', R = ',num2str(Rr(2,1))]); xlabel(char(symb_term{ind_u})); ylabel('y(t)');
        clear File tht
end
tikzName = [folderName,'/Regressions_ut_T_',num2str(T),'.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '18cm', 'width','12cm','checkForUpdates',false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select first significant basis vector for all datasets
iTerm = 1;                                                                  % the first significant term
AEER{iTerm} = zeros(nTerms,1);                                              % placeholder for AERR criteria
for iFile=Files                                                             % over all datasets
    fName = [dictFolder,'/',char(fileNames(iFile))];
    File = matfile(fName,'Writable',true);
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
    File = matfile(fName,'Writable',true);
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
%% Plot AEER across terms
for j=1:nTerms                                                 % assign ticks
    x_ticklabels{j} = char(symb_term{j});
end
figure;
plot(AEER{iTerm}*100,'-o'); hold on;
ylabel('AEER, $\%$'); xlabel('Terms')
set(gca,'XTick',[1:nTerms]);
set(gca,'XTickLabel',x_ticklabels);
xtickangle(90)
tikzName = [folderName,'/AEER_first_',num2str(T),'.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '5cm', 'width','12cm','checkForUpdates',false);

%% Plot correlations
figure;
L2 = round(length(Files)/2);
iPlot = 0;
for iFile=Files                                                             % over all datasets
        fName = [dictFolder,'/',char(fileNames(iFile))];
        File = matfile(fName,'Writable',true);
        ind_u= S(iTerm);
        iPlot = iPlot + 1;
        tht = File.term(index,ind_u)\File.y_narx(index,1);
        y_regressed = File.term(index,ind_u)*tht;
        Rr = round(corrcoef(File.term(index,ind_u),File.y_narx(index,1)),4);
        subplot(L2,2,iPlot)
        scatter(File.term(index,ind_u),File.y_narx(index,1),'filled'); hold on;
        plot(File.term(index,ind_u),y_regressed,'Linewidth',2); 
        title([dataset,num2str(iFile),', R = ',num2str(Rr(2,1))]); xlabel(char(symb_term{ind_u})); ylabel('y(t)');
        clear File tht
end
tikzName = [folderName,'/Regressions_significant_T_',num2str(T),'.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '18cm', 'width','12cm','checkForUpdates',false);

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
        File = matfile(fName,'Writable',true);
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
%% Select optimal number of terms
[min_aamdl,i_min] = min(AAMDL_all);
figure('Name','AAMDL','NumberTitle','off');
plot(AAMDL_all,'o'); hold on;
plot(i_min,min_aamdl,'*','LineWidth',5);
xlim([1 maxSign])
xlabel('Number of terms');
ylabel('AAMDL')
tikzName = [folderName,'/AAMDL_T_',num2str(T),'.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '4cm', 'width','6cm','checkForUpdates',false);
%% Create column of names
finalTerm = i_min;
for iTerm=1:finalTerm
      temp = arrayfun(@char, significant_term{iTerm}, 'uniform', 0);
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
%% Plot AEER surface
figure;
for iTerm=2:finalTerm
plot3(iTerm*ones(nTerms,1),[1:nTerms],AEER{iTerm}*100,'-o'); hold on;
zlabel('AEER, $\%$'); ylabel('Terms'); xlabel('Iteration')
end
set(gca,'YTick',[1:nTerms]);
set(gca,'YTickLabel',x_ticklabels);
ytickangle(45)
tikzName = [folderName,'/All_AEER_T_',num2str(T),'.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '4cm', 'width','12cm','checkForUpdates',false);

clear AEER
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
AEER  = round(AEER_mm(1:finalTerm,1)*100,3);
Table_all = addvars(Tab,AEER,'NewVariableNames',{'AEER($\%$)'})
tableName = [folderName,'/Thetas_T_',num2str(T)];
table2latex(Table_all,tableName);
%% 
L2 = round(finalTerm/2);
figure('Name','Internal parameters','NumberTitle','off');
for iFig=1:finalTerm
    subplot(2,L2,iFig)
    plot(Files,Theta(iFig,Files),'o','LineWidth',2); hold on;
    ylabel(Terms{iFig,1});
    xlabel('Dataset index');
end
tikzName = [folderName,'/Estimates_T_',num2str(T),'.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '6cm', 'width','12cm','checkForUpdates',false);
%% Store data in table
workspaceName = [folderName,'/OLS_results_T_',num2str(T),'.mat'];
save(workspaceName,'Theta','Terms','Files','finalTerm','T','n_y','n_u','S');