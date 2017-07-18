function [] = Convergence(x)

[T,din,N] = size(x);
tt = 1:20:T;
nn = tt*N;
R_stat = nan(length(tt),din);
for t = 1:length(tt)
    t_end = tt(t);
    t_start = max(1,floor(0.5*t_end));
    R_stat(t,:) = Gelman(x(t_start:t_end,:,:));
end
figure('units','normalized','outerposition',[0 0 1 1],'Color',[1 1 1]);
plot(nn,R_stat); hold on
xlabel('Number of model evaluations','fontsize',14,'fontweight','bold','fontname','Times');
ylabel('R_{stat}','fontsize',14,'fontweight','bold','fontname','Times');
plot([0 nn(end)],[1.2 1.2],'k--','linewidth',2);
axis([0 nn(end) 0.8 5]);
title('Convergence of sampled chains','fontsize',14,'fontweight','bold','fontname','Times');

end

function R_stat = Gelman(chain)

[t,din,N] = size(chain);

if (t < 10),
    R_stat = NaN(1,din);
else
    % Determine the _chainuence means
    mean_chain = mean(chain); mean_chain = reshape(mean_chain(:),din,N)';
    
    % Determine the variance between the _chainuence means 
    B = t * var(mean_chain); %#ok<UDIM>  
    
    % Compute the variance of the various chain
    for zz = 1:N,
        var_chain(zz,:) = var(chain(:,:,zz));
    end;   
    
    % Calculate the average of the within _chainuence variances
    W = mean(var_chain);  
    
    % Estimate the target variance
    sigma2 = ((t - 1)/t) * W + (1/t) * B;  
    
    % Compute the R-statistic
    R_stat = sqrt((N + 1)/N * sigma2 ./ W - (t-1)/N/t);
    
end;
end