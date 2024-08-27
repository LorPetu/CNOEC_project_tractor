function [h, g] = constr_tractor_mincon(U,z0,Nu,Ns,parameters, constr_param,MODE)
% CONSTR_TRACTOR_MINCON retrieves the constraints of the non linear problem
% in the form g(x) = 0 and h(x)>=0 to be given as input to the fmincon matlab solver
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

vsat = constr_param.vsat;
delta_psi_sat = constr_param.delta_psi_sat;

m_up= constr_param.m(1);
m_down= constr_param.m(2);
q_up = constr_param.q(1);
q_down = constr_param.q(2);

zf = constr_param.zf;

c_vel = constr_param.c_vel;

tol_f = constr_param.tol_f;


Np=ceil((Ns+1)/Nu);

u_in        =   [U(1:Np,1)';
                U(Np+1:2*Np,1)'];

Ts=     U(end,1);

%% Run simulation with FFD
n_mode      = size(z0,1);

zdot        =   zeros(n_mode,1);
z_sim       =   zeros(n_mode,Ns+1);
z_sim(:,1)  =   z0;

Tractor_model_used = str2func(['Tractor_',MODE, '_trail_model']);
for ind=2:Ns+1
    u                   =   u_in(:,ceil(ind/Nu));
    % FFD
    %zdot                =   Tractor_model_used(z_sim(:,ind-1),u,parameters);
    %z_sim(:,ind)        =   z_sim(:,ind-1)+Ts*zdot;

    % RK2
    z_prime_rk2         =   z_sim(:,ind-1)+Ts/2*Tractor_model_used(z_sim(:,ind-1),u,parameters);
    z_sim(:,ind)        =   z_sim(:,ind-1)+Ts*Tractor_model_used(z_prime_rk2,u,parameters);
    

end 


%% equality constraints 
g=[];
    

%% Inequality constraints h(x)

tol_f = tol_f(1:n_mode);

% ## fmincon info ##
% non linear inequalities constraints are specified in the form h(x)<=0

% Constraints for mode 00
if strcmp(MODE,'00')
    
    h = [
        % Velocity saturation
        (-z_sim(4,:)-c_vel*vsat*ones(1,Ns+1))'; %c_vel can be 0 or 1, is used to select trajectories with vel>0 only USER DEFINED
        (+z_sim(4,:)-vsat*ones(1,Ns+1))';

        % Parametrized boundaries constraints
        (z_sim(2,:)-m_up*z_sim(1,:)-(q_up)*ones(1,Ns+1))';      % y < m*x+q --> y - m*x - q  <=  0 
        (-z_sim(2,:)+m_down*z_sim(1,:)+q_down*ones(1,Ns+1))';   % y > m*x+q --> -y + m*x + q <=  0 

        % Final states tolerances
        (abs(z_sim(:,end)-zf)-tol_f)
        ];


% Constraints for mode 01
else  

    delta_psi=abs(z_sim(3,:)-z_sim(7,:));

    h = [
        % Velocity saturation
        (-z_sim(4,:)-c_vel*vsat*ones(1,Ns+1))'; %c_vel can be 0 or 1, is used to select trajectories with vel>0 only USER DEFINED
        (+z_sim(4,:)-vsat*ones(1,Ns+1))';
  
        % Parametrized boundaries constrained
        (z_sim(2,:)-m_up*z_sim(1,:)-q_up*ones(1,Ns+1))'; % y < m*x+q --> y - m*x - q  <=  0 
        (z_sim(6,:)-m_up*z_sim(5,:)-q_up*ones(1,Ns+1))'; % y < m*x+q --> y - m*x - q  <=  0 
        (-z_sim(6,:)+m_down*z_sim(5,:)+q_down*ones(1,Ns+1))'; % y > m*x+q --> -y + m*x + q <=  0 ONLY FOR THE IMPLEMENT
        
        % Final states tolerances
        (abs(z_sim(:,end)-zf)-tol_f);
        (delta_psi-delta_psi_sat)'  %difference between the orientation of the 2 must be lower than a threshold
        ];     

end






end