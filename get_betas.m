my_init
%% Load estimated thetas
dataset = 'C'; % 'D';
n_y = 0;
n_u = 4;
T = 2000;
fileName = ['OLS_results_',dataset,'_ny_',num2str(n_y),'_nu_',num2str(n_u),'_size_',num2str(T),'.mat'];
load(fileName);
%% Create matrices for fitting
load('External_parameters');
L_cut_all = [values{1}(:, 9);values{2}(:, 9)];
D_rlx_all = [values{1}(:,11);values{2}(:,11)];
A_imp_all = [values{1}(:, 6);values{2}(:, 6)];
V_imp_all = [values{1}(:, 7);values{2}(:, 7)];

index = find(Files <=10);
Files_sub = Files(index);
x = L_cut_all(Files_sub,1);
y = D_rlx_all(Files_sub,1);
%% Estimate coefficients via curve fitting
fo = fitoptions('Method','NonlinearLeastSquares');
if length(x) >= 6
g = fittype(@(b0,b1,b2,b3,b4,b5,x,y) b0 + b1*x + b2*y + b3*x.*y + b4*x.^2 + b5*y.^2, 'independent',{'x','y'},'dependent','z','options',fo); % ); % 
else
g = fittype(@(b0,b1,b2,b3,x,y) b0 + b1*x + b2*y + b3*x.*y, 'independent',{'x','y'},'dependent','z','options',fo); % ); %   
end
clear cfs
for iTerm = 1:finalTerm
z = Theta(iTerm,Files_sub)';
[ft{iTerm},gof{iTerm},outp{iTerm}]= fit([x,y],z,g);
cfs(iTerm,:) = coeffvalues(ft{iTerm});
end
Tab = table(Terms);
for iCoeff = 1:size(cfs,2)
    Parameters = round(cfs(:,iCoeff),2);
    varName = ['$\beta_',num2str(iCoeff-1),'$'];
    Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
end
Table_coeffs_nls = Tab
%% Estimate coefficients with LS
id = ones(size(x));                                                         % create unit vector for constants
A = [id x y x.*y x.^2 y.^2];                                                % create the matrix for ls
if length(x) <= 4
    A = A(:,1:4);
end
clear beta
B = Theta(:,Files_sub)';
beta = A\B; 
beta = beta';

Tab = table(Terms);
for iCoeff = 1:size(beta,2)
    Parameters = round(beta(:,iCoeff),2);
    varName = ['$\beta_',num2str(iCoeff-1),'$'];
    Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
end
Table_coeffs_ls = Tab
%% Save table to tex
tableName = ['Betas_',dataset,'_ny_',num2str(n_y),'_nu_',num2str(n_u),'_size_',num2str(T)];
% table2latex(Table_coeffs_ls,tableName);

%% Compute parameters for validation
testFiles = [3 8];
for iFile = testFiles
L_test = L_cut_all(iFile,1);
D_test = D_rlx_all(iFile,1);
for iTerm = 1:finalTerm
    theta_test{iFile}(iTerm,1) = g(cfs(iTerm,1),cfs(iTerm,2),cfs(iTerm,3),cfs(iTerm,4),cfs(iTerm,5),cfs(iTerm,6),L_test,D_test);
end
fileName = ['Dict_',dataset,num2str(iFile)];
File = matfile(fileName,'Writable',true);
indSign = S(1:finalTerm);                                                   % select the indeces of significant terms from the ordered set
Phi_all = File.term;                                                        % extract all terms into a vector
Phi     = Phi_all(:,indSign);                                               % select only signficant terms
y_model = Phi*theta_test{iFile};                                            % model NARMAX output
RMSE(iFile) = sqrt(mean((File.y_narx - y_model).^2));                                        % Root Mean Squared Error
%% Compare outputs
indPlot = [1:500];
figure('Name','Modelled output','NumberTitle','off');
colormap(my_map);
plot(indPlot+File.t_0,File.y_narx(indPlot,1)); hold on;
plot(indPlot+File.t_0,y_model(indPlot,1),'--'); hold on;
legend('True output','Generated output');

tikzName = ['Generated_output_',dataset,num2str(iFile),'_size_',num2str(T),'.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '6cm', 'width','12cm','checkForUpdates',false);
        
clear File Phi_all Phi y_model
end
        %%
% finalTerm = 8;
index1 = find(Files <=5);
index2 = find(Files > 5);
az = -140;
el =   50;
figure('Name','Parameter surfaces','NumberTitle','off');
colormap(my_map);
L2 = round(finalTerm/2);
for iTerm=1:finalTerm
subplot(L2,2,iTerm);
z = Theta(iTerm,Files)'; 
scatter3(L_cut_all(Files(index1)),D_rlx_all(Files(index1)),z(index1),'filled','LineWidth',5); hold on;
scatter3(L_cut_all(Files(index2)),D_rlx_all(Files(index2)),z(index2),'filled','LineWidth',5); hold on;
for iFile = testFiles
    scatter3(L_cut_all(iFile,1),D_rlx_all(iFile,1), theta_test{iFile}(iTerm,1),'filled','k','LineWidth',5); hold on;
end
srf = plot(ft{iTerm}); 
alpha(srf,0.7);
shading interp
grid on;
xlabel('$L_{cut}$');
ylabel('$D_{rlx}$');
zlabel(Terms{iTerm});
view(az,el)
end
% legend('Foam 1','Foam 2');
%%
% tikzName = ['Theta_surfaces_',dataset,'_ny_',num2str(n_y),'_nu_',num2str(n_u),'_size_',num2str(T),'.tikz'];
% cleanfigure;
% matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
%             false, 'height', '20cm', 'width','12cm','checkForUpdates',false);