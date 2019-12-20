function[d_op] = delta_operator(order,f,t,dT,direction)
switch direction
    case 'forward'
        tq = @(x,shift) x + (order-shift);
    case 'backward'
         tq = @(x,shift) x - shift;
end
sum = 0;
for l=1:order+1
    k = l-1;
      sum = sum + (-1)^(k)*nchoosek(order,k)*f(tq(t,k));
end
d_op = sum/(dT)^order;

        
        