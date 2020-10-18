function jump = adaptive_independence_sir(win, corr_coef_th)
n_tot = size(win,1);
k = 1;
while( k <= n_tot-1 )
    idxes = [1:k:n_tot];
    correlation_coef = [];
    total_jump = length(idxes)-1;
    for j = 1: total_jump 
        temp_corr = corrcoef(win(idxes(j),:),win(idxes(j+1),:));
        correlation_coef = [ correlation_coef, temp_corr(1,2) ];
    end
    if mean(correlation_coef) < corr_coef_th
        break
    else
        k= k+1;
    end
end
if k == 1
    jump = k;
else
    jump = k-1;
end
