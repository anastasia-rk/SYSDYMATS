function dydt = vdpo(t,y,d,w,phase,ak,Mu);
% Van der Pol oscillator with an input in ODE format
% generate input signal
for k=1:d
        fun(k) = ak*cos(2*pi*w(k)*t + phase(k)); 
end
u = sum(fun);
% right hand side of ODE    
dydt = [y(2); Mu*(1-y(1)^2)*y(2)-y(1)+u];