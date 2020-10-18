function jump = adaptive_independence(win, corr_coef_th)
i=1;
jump_adap = [];
n_tot = size(win,1);
while (1)    
    for j=i:n_tot
        corr = corrcoef(win(i,:),win(j,:));
        if corr(1,2) < corr_coef_th
            break;
        end
    end
    jump_adap= [jump_adap ,j-i];
    i=j;
    if i==n_tot
        break;
    end
end 
jump= mean(jump_adap);
