function [h, g] = constr_tractor_mincon(U,z0,parameters,Optimization_opt, constr_param)
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

s =    U(2*Np+1:end-1,1);

Ts=     U(end,1);

%% Run simulation with FFD

zdot        =   zeros(8,1);
z_sim      =   zeros(8,Ns+1);
z_sim(:,1) =   z0;
for ind=2:Ns+1
    u               =  u_in(:,ceil(ind/Nu));
    zdot               =   Tractor_01_trail_model(z_sim(:,ind-1),u,parameters);
    z_sim(:,ind)       =   z_sim(:,ind-1)+Ts*zdot;

end 

delta_psi=abs(z_sim(3,:)-z_sim(7,:));
%% equality constraints 


g= [(abs(z_sim(1, end)-zf(1))-s(1,1));                %uguaglianza con slask variables sulla fine
    (abs(z_sim(2, end)-zf(2))-s(2,1));
    (abs(z_sim(3, end)-zf(3))-s(3,1));
    (abs(z_sim(4, end)-zf(4))-s(4,1));
    (abs(z_sim(5, end)-zf(5))-s(5,1));
    (abs(z_sim(6, end)-zf(6))-s(6,1));
    (abs(z_sim(7, end)-zf(7))-s(7,1));
    (abs(z_sim(8, end)-zf(8))-s(8,1))];
    

%% Inequality constraints h(x)

h = [(-z_sim(4,:)-lb_vel*vsat*ones(1,Ns+1))'; %lb_vel can be 0 or 1, is used to select trajectories with vel>0 only USER DEFINED
    (+z_sim(4,:)-vsat*ones(1,Ns+1))';
  
% Parametrized constrained
   (z_sim(2,:)-m_up*z_sim(1,:)-(q_up)*ones(1,Ns+1))'; % -y + m*x + q > 0 y < m*x+q
   (z_sim(6,:)-m_up*z_sim(5,:)-(q_up)*ones(1,Ns+1))';
   (-z_sim(6,:)+m_down*z_sim(5,:)+q_down*ones(1,Ns+1))'; % y - m*x - q > 0
   (delta_psi-delta_psi_sat)'];      %difference between the orientation of the 2 must be lower than a threshold



end