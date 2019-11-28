function[R_sq] = r_squared(Y_bar,Y_hat)
R_sq = 1 - sum((Y_bar - Y_hat).^2)/sum((Y_bar - mean(Y_bar)).^2);