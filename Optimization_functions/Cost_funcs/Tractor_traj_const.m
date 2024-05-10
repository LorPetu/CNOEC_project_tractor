function [z_sim] = Tractor_traj_const(U,z0,zf,Np,Ts,parameters)
%% Build vector of inputs

u_in        =   [   U(1:Np,1)';
                    U(Np+1:end,1)'];

%% Run simulation with FFD
time_FFD    =   [0:0.01:(Np-1)*Ts];
Nblock      =   Ts/0.01;
Nsim_FFD    =   length(time_FFD);

ztemp=z0;
z_sim      =   zeros(4,Nsim_FFD);
z_sim(:,1) =   z0;

for ind=2:Nsim_FFD
    u                  =   u_in(:,1+floor(time_FFD(ind)/Ts));
    zdot               =   tractor_model(z_sim(:,ind-1),u,parameters);
    z_sim(:,ind)       =   z_sim(:,ind-1)+Ts/Nblock*zdot;
    

  end

%% Post-processing the results

plx     = z_sim(1,:);
ply     = z_sim(2,:);
ang     = z_sim(3,:);
vel     = z_sim(4,:);


vsat        =   10;                     % Input saturation
asat        =   3;                      % Cart position limits
deltasat    =   30*pi/180;
maxvsat = vsat*ones(Nsim_FFD);
minvsat = -vsat*ones(Nsim_FFD);
maxdeltasat = deltasat*ones(Np);
mindeltasat = -deltasat*ones(Np);
maxasat = asat*ones(Np);
minasat = -asat*ones(Np);

del     = u_in(1,:);
acc     = u_in(2,:);


plxf = zf(1,1);
plyf = zf(2,1);

time_s=time_FFD;%linspace(0,Ts_s,Nsim_FFD);
time_p=linspace(0,Ts,Np);

figure(1);
subplot(5,2,1:6),plot(plx,ply,plxf,plyf,"xr"),daspect([1,1,1]),grid on,axis([-10 10 -10 10])
xlabel('x'), ylabel('y'),title('traiettoria')
subplot(5,2,7),plot(time_s,ang),xlabel('Time (s)'),ylabel('psi'),
subplot(5,2,8),plot(time_s,vel,time_s,maxvsat,time_s,minvsat),xlabel('Time (s)'),ylabel('velocità')
subplot(5,2,9);plot(time_p,del,time_p,maxdeltasat,time_p,mindeltasat),xlabel('Time (s)'),ylabel('delta')
subplot(5,2,10);plot(time_p,acc,time_p,maxasat,time_p,minasat),xlabel('Time (s)'),ylabel('acc');



end

