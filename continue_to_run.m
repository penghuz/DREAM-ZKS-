function [x,p_x,Z,y] = continue_to_run(N,T,din,dout,Obs,sd,range,alpha)

load DREAM_ZKS
t_last = t;
for t = t_last+1:T
    disp([num2str(t/T*100,'%4.2f'),'% finished'])                     % Show the process
    Xp = nan(N,din); X_ind = zeros(N,2); p_acc = nan(N,1);            % X_ind records the use of KF updater
    dX = zeros(N,din);                                                % Set N jump vectors to zero
    lambda = unifrnd(-c,c,N,1);                                       % Draw N lambda values
    std_X = std(X);                                                   % Compute std of each dimension
    R = randsample(1:m,N*n_d,'false'); R = reshape(R,n_d,N);          % Sample N*n_d integer values from [1..m]
    method = randsample({'parallel' 'snooker'},1,'true',[1-p_s p_s]); % Mix of parallel and snooker updates
    for i = 1:N
        D = randsample(1:delta,1,'true');                             % Select delta (equal select. prob.)
        a = R(1:D,i); b = R(D+1:2*D,i); c_sn = R(2*D+1:3*D,i);        % Define a and b (parallel) + c (snooker)
        if strcmp(method,'parallel'),                                 % PARALLEL DIRECTION update
            id(i) = randsample(1:n_CR,1,'true',p_CR);                 % Select index of crossover value
            z = rand(1,din);                                          % Draw d values from U[0,1]
            A = find(z < CR(id(i)));                                  % Derive subset A selected dimensions
            d_star = numel(A);                                        % How many dimensions sampled?
            if d_star == 0, [~,A] = min(z); d_star = 1; end           % A must contain at least one value
            gamma_d = 2.38/sqrt(2*D*d_star);                          % Calculate jump rate
            g = randsample([gamma_d 1],1,'true',[1-p_g p_g]);         % Select gamma: 80/20 mix [default 1]
            dX(i,A) = c_star*randn(1,d_star) + ...
                (1+lambda(i))*g*sum(Z(a,A)-Z(b,A),1);                 % Compute ith jump diff. evol.
        elseif strcmp(method,'snooker'),                              % SNOOKER update
            id(i) = n_CR;                                             % Full crossover
            F = X(i,1:din)-Z(a,1:din); D_e = max(F*F',1e-300);        % Define projection X(i,1:d) - Z(a,1:d)
            zP = F*( sum((Z(b,1:din)-Z(c_sn,1:din)).*F ) / D_e );     % Orthogonally project zR1 and zR2 onto F
            g = 1.2 + rand;                                           % Determine jump rate
            dX(i,1:din) = c_star*randn(1,din) + (1+lambda(i))*g*zP;   % Calculate ith jump snooker update
        end
        if rand < alpha && t > N1 + 2                          % Kalman filter-based update
            X_ind(i,1) = 1;                                    % Record the use of KF update
            a1 = union(randsample(1:t-2,N1),t-N1:t-2);         % Use a mixture of global and local samples
            b1 = setdiff(1:N,i); c1 = 2;
            d1 = randsample(b1,c1);                            % Choose c1 chain pairs without replacement
            xf = [GenParSet(x(a1,:,d1)); X(i,:)]';             % Ensemble of the model parameters
            yf = [GenParSet(y(a1,:,d1)); Y(i,:)]';             % Ensemble of the model responses
            Ne = size(xf,2);                                   % Number of samples in the ensemble
            Cd = eye(dout);                                    % Covariance of the measurement errors
            for ii = 1:dout
                Cd(ii,ii) = sd(ii)^2;
            end
            meanxf = repmat(mean(xf,2),1,Ne);
            meanyf = repmat(mean(yf,2),1,Ne);
            Cxy =  (xf-meanxf)*(yf-meanyf)'/(Ne-1);            % Cross-covariance between xf and yf
            Cyy =  (yf-meanyf)*(yf-meanyf)'/(Ne-1);            % Auto-covariance of yf
            K{i,1} = Cxy/(Cyy+Cd);                             % Kalman gain
            K{i,2} = K{i,1}*(Obs-Y(i,:)');
            dX(i,1:din) = K{i,1}*(randn(dout,1).*sd) + K{i,2}; % Calculate the jump distance
        end
        if X_ind_p(i,1) == 1                                   % If the previous state is from the KF update
            dX(i,1:din) = -K{i,1}*(randn(dout,1).*sd) - K{i,2};% Backward jump of the KF update
        end
        Xp(i,:) = X(i,1:din) + dX(i,1:din);                    % Compute the ith proposal
        Xp(i,:) = Boundary_handling(Xp(i,:),range,'reflect');  % Boundary handling
    end
    
    [p_Xp,Yp] = pdf(Xp,Obs,sd);                                       % Calculate the log pdf and the model responses 
    if strcmp(method,'snooker'),                                      % Snooker correction: non-symmetry
        alfa_sn = (sum((Xp - Z(R(1,1:N),1:din)).^2,2)./sum((X - ...   % Proposal distribution
            Z(R(1,1:N),1:din)).^2,2)).^((din-1)/2);
    else                                                              % Otherwise no correction needed
        alfa_sn = ones(N,1);
    end
    for i = 1:N                                                       % Accept/reject proposals (can use "parfor")
        if X_ind_p(i,1) == 1 || X_ind(i,1) == 1;
            alfa_sn(i) = 1;                                           % If use the forward or backward KF jump
        end
        p_acc(i) = min(1,alfa_sn(i)*exp(p_Xp(i,1)-p_X(i,1)));         % Compute acceptance probability
        if p_acc(i) > rand,                                           % p_acc(i) larger than U[0,1]?
            X(i,:) = Xp(i,:); p_X(i,1) = p_Xp(i,1);                   % True: Accept proposal
            Y(i,:) = Yp(i,:); X_ind(i,2) = 1;
        else
            dX(i,1:din) = 0;                                          % Set jump back to zero for p_CR
        end
        J(id(i)) = J(id(i)) + sum((dX(i,1:din)./std_X).^2);           % Update jump distance crossover idx
        n_id(id(i)) = n_id(id(i)) + 1;                                % How many times idx crossover used
    end
    x(t,1:din,1:N) = reshape(X',1,din,N); p_x(t,1:N) = p_X';          % Append current X and density
    y(t,1:dout,1:N) = reshape(Y',1,dout,N);                           % Append model responses;
    X_ind_p = X_ind(:,1).*X_ind(:,2);                                 % Store the record of KF update
    
    if (mod(t,k) == 0),                                               % Check whether to append X to archive Z
        Z(m+1:m+N,1:din+1) = [X p_X]; m = m + N;                      % Append current values to Z after k generations
        if t<T/10,
            p_CR = J./n_id; p_CR = p_CR/sum(p_CR);                    % Update selection prob. crossover
        end
    end
    if mod(t,100) == 0
        save DREAM_ZKS.mat                                            % Save intermediate results
    end
end

end

function [x] = Boundary_handling(x,range,flag)
% Function to check whether parameter values remain within prior bounds

% First determine the size of x
[m,~] = size(x);
% Now replicate min and max
min_d = repmat(range(:,1)',m,1); max_d = repmat(range(:,2)',m,1);
% Now find which elements of x are smaller than their respective bound
[ii_low] = find(x < min_d);
% Now find which elements of x are larger than their respective bound
[ii_up] = find(x > max_d);

switch flag
    case 'reflect'  % reflection
        % reflect in min
        x(ii_low)= 2 * min_d(ii_low) - x(ii_low);
        % reflect in max
        x(ii_up)= 2 * max_d(ii_up) - x(ii_up);
    case 'bound'    % set to bound
        % set lower values to min
        x(ii_low)= min_d(ii_low);
        % set upper values to max
        x(ii_up)= max_d(ii_up);
    case 'fold'     % folding
        % Fold parameter space lower values
        x(ii_low) = max_d(ii_low) - ( min_d(ii_low) - x(ii_low) );
        % Fold parameter space upper values
        x(ii_up) = min_d(ii_up) + ( x(ii_up) - max_d(ii_up) );
        % Now double check in case elements are still out of bound -- this is
        % theoretically possible if values are very small or large
    otherwise
        disp('do not know this boundary handling option! - treat as unbounded parameter space')
end

% Now double check if all elements are within bounds
% If parameter so far out of space then possible that reflection or
% folding still creates point that is out of bound - this is very
% rare but needs explicit consideration
if strcmp(flag,'reflect') || strcmp(flag,'fold')
    % Lower bound
    [ii_low] = find(x < min_d); x(ii_low) = min_d(ii_low) + rand(size(ii_low)).* ( max_d(ii_low) - min_d(ii_low) );
    % Upper bound
    [ii_up]  = find(x > max_d); x(ii_up)  = min_d(ii_up)  + rand(size(ii_up)).*  ( max_d(ii_up)  - min_d(ii_up)  );
end
end
