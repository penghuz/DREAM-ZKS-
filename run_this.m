clear all;clc

N = 3;               % Number of parallel chains
T = 1500;            % Number of iterations
din = 5;             % Dimension of the model parameters
dout = 731;          % Dimension of the model responses

xmin = [1.0 0.10 0.10 0.00 0.10];      % Lower bound of the prior model parameters
xmax = [500 2.00 0.99 0.10 0.99];      % Upper bound of the prior model parameters
range = [xmin' xmax'];                 % Range of the prior model parameters

xreal = prior(1,din);                  % The reference parameters drawn from the prior distribution
yreal = forwardmodel(xreal);
sd = yreal*0.1;                        % Standard deviation of the measurement errors
Obs = yreal + sd.*randn(size(yreal));  % The measurements

alpha = 0.1;                           % The probability of using the Kalman filter-based updater
[x,~,~,fx]  = dream_zks(N,T,din,dout,Obs,sd,range,alpha);
chain_es = GenParSet(x);
drawplot(chain_es,xreal,range,2,3)     % Show the trace plot of the model parameters

alpha = 0;                             % Using the original dream_zs
[x,~,~,fx]  = dream_zks(N,T,din,dout,Obs,sd,range,alpha);
chain = GenParSet(x);
drawplot(chain,xreal,range,2,3)        % Show the trace plot of the model parameters

save results
delete DREAM_ZKS.mat                   % The intermediate results saved when running dream_zks
beep on; beep;