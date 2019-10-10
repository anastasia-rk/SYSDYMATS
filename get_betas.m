my_init
load('External_parameters');
%% Load estimated thetas
dataset = 'C'; % 'D';
n_y = 0;
n_u = 4;
fileName = ['OLS_results_',dataset,'_ny_',num2str(n_y),'_nu_',num2str(n_u),'.mat'];
load(fileName);
%% Create matrices for fitting
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
    Parameters = round(cfs(:,iCoeff),3);
    varName = ['$b_',num2str(iCoeff-1),'$'];
    Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
end
Table_coeffs_nls = Tab
%% Estimate coefficients with LS
id = ones(size(x));
A = [id x y x.*y x.^2 y.^2];
if length(x) <= 5
    A = A(:,1:4);
end
clear beta
for iTerm = 1:finalTerm
B = Theta(iTerm,Files_sub)';
beta(:,iTerm) = A\B;
end 
beta = beta';

Tab = table(Terms);
for iCoeff = 1:size(beta,2)
    Parameters = round(beta(:,iCoeff),3);
    varName = ['$b_',num2str(iCoeff-1),'$'];
    Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
end
Table_coeffs_ls = Tab
%%
% finalTerm = 8;
index1 = find(Files <=5);
index2 = find(Files > 5);
az = -140;
el =   50;
figure;
colormap(my_map);
L2 = round(finalTerm/2);
for iTerm=1:finalTerm
subplot(2,L2,iTerm);
z = Theta(iTerm,Files)'; 
plot3(L_cut_all(Files(index1)),D_rlx_all(Files(index1)),z(index1),'*','LineWidth',5); hold on;
plot3(L_cut_all(Files(index2)),D_rlx_all(Files(index2)),z(index2),'*','LineWidth',5); hold on;
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
%% Save table to tex
tableName = ['Betas_',dataset,'_ny_',num2str(n_y),'_nu_',num2str(n_u),'_size_',num2str(T)];
table2latex(Table_coeffs_nls,tableName);

%% Compute parameters for validation
iFile = 3;
A_test = A_imp_all(iFile,1);
V_test = V_imp_all(iFile,1);
for iTerm = 1:finalTerm
    theta_test(iTerm,1) = g(cfs(iTerm,1),cfs(iTerm,2),cfs(iTerm,3),cfs(iTerm,4),cfs(iTerm,5),cfs(iTerm,6),A_test,V_test);
end
fileName = ['Dict_',dataset,num2str(iFile)];
File = matfile(fileName,'Writable',true);