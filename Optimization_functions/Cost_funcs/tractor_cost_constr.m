function v = tractor_cost_constr(U,Ts,Np,th)
% Function that computes the trajectory of the tractor exiting from a row
% and the constraint


vsat        =   10;                     % Input saturation
asat        =   3;                      % Cart position limits
deltasat    =   30*pi/180;

%% initial states
xt        =      0;                 % inertial X position (m)
yt        =      0;                 % inertial Y position (m)
psit      =      pi/2;              % yaw angle (rad)
vt        =      4;                     % body x velocity (m/s) 

z0=[xt;yt;psit;vt];
%% final state
xf      =       3;
yf      =       0;
psif    =       -pi/2;
vf      =       4;

zf      =       [xf;yf;psif;vf];

%% Build vector of inputs
t_in        =   [0:Ts:(Np-1)*Ts]';
u_in        =   [U(1:Np,1)';
                U(Np+1:end,1)'];


%% Run simulation with FFD
time_FFD    =   [0:0.01:(Np-1)*Ts];
Nblock      =   Ts/0.01;
Nsim_FFD    =   length(time_FFD);

z_sim      =   zeros(4,Nsim_FFD);
z_sim(:,1) =   z0;
f=0;
for ind=2:Nsim_FFD
    u                  =   u_in(:,1+floor(time_FFD(ind)/Ts));
    zdot               =   tractor_model(z_sim(:,ind-1),u,th);
    z_sim(:,ind)       =   z_sim(:,ind-1)+Ts/Nblock*zdot;

    f= f +  1/(z_sim(4,ind)*cos(z_sim(3,ind)))*Ts/Nblock*zdot(1);
end

X_sim       =   z_sim(1,1:end)';
Y_sim       =   z_sim(2,1:end)';


%% Equality constraints g(x)
g = [z_sim(:,end)-zf(:,1)];

%% Inequality constraints h(x)

h = [(z_sim(4,:)+vsat*ones(1,Nsim_FFD))';
    (-z_sim(4,:)+vsat*ones(1,Nsim_FFD))'];
%% Compute cost function f(x)
% Come funziona questa?
% delta_diff  =   (x(Np+5:end,1)-x(Np+4:end-1,1));
% Td_diff     =   (x(5:Np+3,1)-x(4:Np+2,1));

%% Stack cost and constraints
v           =   [f;g;h];

end