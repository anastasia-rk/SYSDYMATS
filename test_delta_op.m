%%
my_init;
dT = 0.1;
x = 1:dT:100;
y = sin(x);
y1 = cos(x);
y3 = -sin(x);
%%
for t=2:length(x)-1;
    y2(t-1) = delta_operator(1,y,t,dT,'backward');
end

for t=3:length(x)-2;
    y4(t-2) = delta_operator(2,y,t,dT,'backward');
end
%%
figure;
subplot(2,1,1)
plot(y,'b'); hold on;
plot(y1,'k','LineWidth',2); hold on;
plot(y2,'r--','LineWidth',2); hold on;
legend('function','analytical','$\delta$-operator')
title('$\lambda$ = 1');
subplot(2,1,2)
plot(y,'b'); hold on;
plot(y3,'k','LineWidth',2); hold on;
plot(y4,'r--','LineWidth',2); hold on;
legend('function','analytical','$\delta$-operator')
title('$\lambda$ = 2');
