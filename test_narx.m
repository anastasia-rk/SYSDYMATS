my_init;
%% Set up and simulate the model

narx = @(u,y,n_u,n_y,Theta) Theta(1:n_u)*u(:,end-n_u+1:end)' + Theta(n_u+1:n_u+n_y)*y(:,end-n_y+1:end)';

theta = ones(1,10);
u = [1:10];
y = [1:10];
n_u = 3;
n_y = 2;
%%
narx(u,y,n_u,n_y,theta)

