function f = cost_tractor_mincon(U,z0,Nu,Ns,parameters, constr_param, MODE)
% COST_TRACTOR_MINCON retrieves the cost function of the non linear problem
% according to notation of the fmincon matlab solver
%   INPUTS:
%       - U                 = optimization variable vector 
%       - z0                = initial guess for the optimization variables
%       - parameters        = model parameters
%       - constr_param      = struct containing all required parameters to describe constraints
%       - MODE              = selector for the model to be used between
%       tractor-only and tractor and implement
%   OUTPUTS:
%       - h                 = nonlinear inequality constraints vector
%       - g                 = nonlinear equality constraints vector

%% Retrieve variables and parameters
zf          =   constr_param.zf;

Np          =   ceil((Ns+1)/Nu);

u_in        =   [U(1:Np,1)';
                U(Np+1:2*Np,1)'];

Ts          =   U(end,1);

%% Run simulation with RK2

n_mode      =   size(z0,1);             % Indicates how many states are involved according to the selected model

z_sim       =   zeros(n_mode,Ns+1);
z_sim(:,1)  =   z0;                     % Initialized states vector with initial states z0

p           =   ones(size(zf));         % Weighting vector p to be used in the cost function
p(3)        =   5;

if strcmp(MODE,'01')
    p(7)        = p(3);
end

% Retrieve the model function according to the selected mode  
Tractor_model_used = str2func(['Tractor_',MODE, '_trail_model']);

% Start numerical integration algorithm including downsampled input 
for ind=2:Ns+1
    
    u                   =  u_in(:,ceil(ind/Nu));
    % RK2 method
    z_prime_rk2         =   z_sim(:,ind-1)+Ts/2*Tractor_model_used(z_sim(:,ind-1),u,parameters);
    z_sim(:,ind)        =   z_sim(:,ind-1)+Ts*Tractor_model_used(z_prime_rk2,u,parameters);

end 

%% Cost Function

f1          =   p'*abs(z_sim(:,end)-zf);    % Final states accuracy term 
f2          =   Ns*Ts;                      % Total execution time
gamma       =   0.15;                       % Trade off parameter

f           =   gamma*f1 + (1-gamma)*f2;    % Multi-objective cost function

f           =   f*50;                       % Scaling factor 
end