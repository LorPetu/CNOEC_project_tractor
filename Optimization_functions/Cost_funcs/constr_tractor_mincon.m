function [h, g] = constr_tractor_mincon(U,z0,parameters,Optimization_opt, constr_param,MODE)
% Function that computes the trajectory of the tractor exiting from a row
% and the constraint

vsat = Optimization_opt.vsat;
delta_psi_sat = Optimization_opt.delta_psi_sat;

Ns = Optimization_opt.Ns;
Nu=Optimization_opt.Nu;

m_up= constr_param.m(1);
m_down= constr_param.m(2);
q_up = constr_param.q(1);
q_down = constr_param.q(2);

zf = constr_param.zf;

lb_vel = constr_param.lb_vel;


Np=ceil((Ns+1)/Nu);

u_in        =   [U(1:Np,1)';
                U(Np+1:2*Np,1)'];

Ts=     U(end,1);

%% Run simulation with FFD
n_mode      = size(z0,1);

zdot        =   zeros(n_mode,1);
z_sim       =   zeros(n_mode,Ns+1);
z_sim(:,1) =   z0;

Tractor_model_used = str2func(['Tractor_',MODE, '_trail_model']);
for ind=2:Ns+1
    u               =  u_in(:,ceil(ind/Nu));
    zdot               =   Tractor_model_used(z_sim(:,ind-1),u,parameters);
    z_sim(:,ind)       =   z_sim(:,ind-1)+Ts*zdot;

end 


%% equality constraints 
g=[];
    

%% Inequality constraints h(x)

lim = [0.05,0.05,5*pi/180,0.5/3.6,0.05,0.05,5*pi/180,0.5/3.6]';

lim = lim(1:n_mode);

% Constraints for mode 00
if strcmp(MODE,'00')
    
    h = [(-z_sim(4,:)-lb_vel*vsat*ones(1,Ns+1))'; %lb_vel can be 0 or 1, is used to select trajectories with vel>0 only USER DEFINED
        (+z_sim(4,:)-vsat*ones(1,Ns+1))';
        (z_sim(2,:)-m_up*z_sim(1,:)-(q_up)*ones(1,Ns+1))'; % -y + m*x + q > 0 y < m*x+q
        (-z_sim(2,:)+m_down*z_sim(1,:)+q_down*ones(1,Ns+1))'; % y - m*x - q > 0
        (abs(z_sim(:,end)-zf)-lim)];


% Constraints for mode 01
else  

    delta_psi=abs(z_sim(3,:)-z_sim(7,:));

    h = [(-z_sim(4,:)-lb_vel*vsat*ones(1,Ns+1))'; %lb_vel can be 0 or 1, is used to select trajectories with vel>0 only USER DEFINED
    (+z_sim(4,:)-vsat*ones(1,Ns+1))';
  
% Parametrized constrained
    (z_sim(2,:)-m_up*z_sim(1,:)-(q_up)*ones(1,Ns+1))'; % -y + m*x + q > 0 y < m*x+q
    (z_sim(6,:)-m_up*z_sim(5,:)-(q_up)*ones(1,Ns+1))';
    (-z_sim(6,:)+m_down*z_sim(5,:)+q_down*ones(1,Ns+1))'; % y - m*x - q > 0
% Final position tolerances
    (abs(z_sim(:,end)-zf)-lim);
    (delta_psi-delta_psi_sat)'];      %difference between the orientation of the 2 must be lower than a threshold

end






end