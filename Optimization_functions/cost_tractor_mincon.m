function f = cost_tractor_mincon(U,z0,Nu,Ns,parameters, constr_param, MODE)
% COST_TRACTOR_MINCON retrieves the cost function of the non linear problem
% in the form  to be given as input to the fmincon matlab solver
%   INPUTS:
%       - U                 =
%       - z0                =
%       - parameters        = 
%       - optimization_opt  = 
%       - constr_param      =
%       - MODE              = 
%   OUTPUTS:
%       - h                 = nonlinear inequality constraints
%       - g                 = nonlinear equality constraints

zf = constr_param.zf;

Np=ceil((Ns+1)/Nu);

u_in        =   [U(1:Np,1)';
                U(Np+1:2*Np,1)'];

Ts=     U(end,1);

%% Run simulation with FFD

n_mode      = size(z0,1);

z_prime_rk2=0;
z_sim       =   zeros(n_mode,Ns+1);
z_sim(:,1)  =   z0;
f           =   0;
f1          =   0;
f2          =   0;

p           =   ones(size(zf));
p(3)        = 5;
if strcmp(MODE,'01')
    p(7)        = p(3);
end

Tractor_model_used = str2func(['Tractor_',MODE, '_trail_model']);

for ind=2:Ns+1
    
    u               =  u_in(:,ceil(ind/Nu));
    % FFD
    %zdot                =   Tractor_model_used(z_sim(:,ind-1),u,parameters);
    %z_sim(:,ind)        =   z_sim(:,ind-1)+Ts*zdot;

    % RK2
    z_prime_rk2         =   z_sim(:,ind-1)+Ts/2*Tractor_model_used(z_sim(:,ind-1),u,parameters);
    z_sim(:,ind)        =   z_sim(:,ind-1)+Ts*Tractor_model_used(z_prime_rk2,u,parameters);

    
end 

f1=p'*abs(z_sim(:,end)-zf);
f2 = Ns*Ts;

gamma =0.15;
f = gamma*f1 + (1-gamma)*f2; 
%disp(["f1 = ", num2str(f1),"f2= ",num2str(f2)]);

f=f*50;      %questo serve per scalare la funzione. Serve perchè fmincon non può settare i valori di linsearch e quindi con questo riusciamo a cambiarli (credo).
                %se è più alto la ricerca è più lenta ma più precisa es(50
                %o 100)sembrano funzionare bene
end