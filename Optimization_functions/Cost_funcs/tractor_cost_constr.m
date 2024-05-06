function v = tractor_cost_constr(U,Ts,Np,th)
% Function that computes the trajectory of the tractor exiting from a row
% and the constraint

%% Build vector of inputs
t_in        =   [0:Ts:(Np-1)*Ts]';
z0         =   [0;x(1:2,1);0;x(3,1);0];
u_in        =   [U(1:Np,1)';
                U(Np+1:end,1)'];


%% Run simulation with FFD
time_FFD    =   [0:0.01:(Np-1)*Ts];
Nblock      =   Ts/0.01;
Nsim_FFD    =   length(time_FFD);

z_sim      =   zeros(6,Nsim_FFD);
z_sim(:,1) =   z0;
for ind=2:Nsim_FFD
    u                   =   u_in(:,1+floor(time_FFD(ind)/Ts));
    zdot               =   tractor_model(0,z_sim(:,ind-1),u,0,th);
    z_sim(:,ind)       =   z_sim(:,ind-1)+Ts/Nblock*zdot;
end

X_sim       =   z_sim(1,1:Nblock:end)';
Y_sim       =   z_sim(2,1:Nblock:end)';

%% Compute path constraints h(x)
h           =   [Y_sim-(tanh((X_sim-100)/2e1)*10+5);
                -Y_sim+(tanh((X_sim-75)/2e1)*10+15)];

%% Constraints
% Questi vanno riscritti con X_sim e Y_sim 
h =         vsat+z_sim(4,ind) 
            vsat-z_sim(4,ind)
            deltasat+u(1,1)
            deltasat-u(1,1)
            asat+u(2,1)
            asat-u(2,1)

%% Cost function
    f   =   

%% Compute cost function f(x)
% Come funziona questa?
delta_diff  =   (x(Np+5:end,1)-x(Np+4:end-1,1));
Td_diff     =   (x(5:Np+3,1)-x(4:Np+2,1));
f           =   -X_sim(end,1)+1e3*z_sim(5,end)^2+1e4*(delta_diff'*delta_diff)+1e-1*(Td_diff'*Td_diff);

%% Stack cost and constraints
v           =   [f;h];

end