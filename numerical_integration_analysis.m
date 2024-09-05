clear all
close all
clc

addpath(genpath('Models/'));
addpath(genpath('Optimization_functions/'));
%% Model Parameters

Lt        =   3;                   % Wheelbase [m]
Li        =   2;                   % Wheelbase of implements [m]
d         =   4.0;                 % Row width [m]

parameters=[Lt;Li;d];

%% initial states

psit      =    pi/2;               % tractor yaw angle [rad]
vt        =    15/3.6;             % tractor velocity [m/s] 
xt        =    0;                  % tractor COG X position [m]
yt        =    0;                  % tractor COG Y position [m]

z0=[xt;yt;psit;vt];

%% Input variables


Tend=10;                           % total time of the simulation

j=[0.01,0.025,0.05,0.1,0.25,0.5];  % vector of sampling time Ts to use in the simulations

n=length(j);

time=zeros(n,5);                   % initialization of the array of the time required to perform the simulations for each sampling time and numerical integration method 
max_error=zeros(n,4);              % initialization of the array of the maximum error with respect to ODE45 obtained in the simulations for each sampling time and numerical integration method 
mean_error=zeros(n,4);             % initialization of the array of the mean error with respect to ODE45 obtained in the simulations for each sampling time and numerical integration method 

%% Execution of the simulations
for k=1:length(j)                  % for cicle to run the simulation for every sampling time in j
       
    Ts=j(k);                       
    disp( num2str(Ts))             % display which Ts is used in the actual cicle
    Ns=Tend/Ts;                    % number of simulation steps 

    delta=-30*pi/180*ones(1,Ns);   % steering angle applied to the system (maximum allowed value)
    acc=zeros(1,Ns);               % acceleration applied to the system (zero since already at maximum allowed speed)
    u=[delta; acc];                % creation of input vector
    
    
    errore_ffd=zeros(1,Ns-1);      % creation errror vector for FFD
    errore_rk2=zeros(1,Ns-1);      % creation errror vector for RK2
    errore_rk3=zeros(1,Ns-1);      % creation errror vector for RK3
    errore_rk4=zeros(1,Ns-1);      % creation errror vector for RK4

    zdot=zeros(4,1);               % initialization of vectors used for FFD
    z_sim_ffd=zeros(4,Ns);
    z_sim_ffd(:,1)=z0;

    z_prime_rk2=zeros(4,1);        % initialization of vectors used for RK2
    z_sim_rk2=zeros(4,Ns);
    z_sim_rk2(:,1)=z0;

    z_prime_rk3_1  =   zeros(4,Ns);% initialization of vectors used for RK3
    z_prime_rk3_2 =   zeros(4,Ns);
    z_sim_rk3    =   zeros(4,Ns);
    z_sim_rk3(:,1) =   z0;


    z_sim_rk4    =   zeros(4,Ns);  % initialization of vectors used for RK4
    z_sim_rk4(:,1) =   z0;

    z_sim_ode=zeros(4,Ns);         % initialization of vectors used for ODE45
    z_sim_ode(:,1)=z0;
 %% simulation using ode45

    tic

    % simulation of the system using ODE45

    for ind=2:Ns                   
        
        z_temp= ode45(@(t,z)Tractor_00_trail_model(z, u(:, ind-1), parameters),[0,Ts], z_sim_ode(:, ind-1));
        z_sim_ode(:, ind) = z_temp.y(:,end);

    end 
    time(k,1)=toc;                 % saving time of computation


    
    %% simulation using rk2

    tic 

    % simulation of the system using RK2

    for ind=2:Ns                                
    
        z_prime_rk2       =   z_sim_rk2(:,ind-1)+Ts/2*Tractor_00_trail_model(z_sim_rk2(:,ind-1),u(:,ind-1),parameters);
        z_sim_rk2(:,ind)  =   z_sim_rk2(:,ind-1)+Ts*Tractor_00_trail_model(z_prime_rk2,u(:,ind-1),parameters);
     
         errore_rk2(k,ind-1)=sqrt((z_sim_rk2(1,ind)-z_sim_ode(1,ind))^2+(z_sim_rk2(2,ind)-z_sim_ode(2,ind))^2);
    end 
    time(k, 2)=toc;                             % saving time of computation
    max_error(k, 1)= max(errore_rk2(k,:));      % saving maximum error
    mean_error(k, 1)= mean(errore_rk2(k,:));    % saving mean error
   
    
    %% RK3
    
    tic

    % simulation of the system using RK3

    for ind = 2:Ns                                
        
        k1 = Tractor_00_trail_model(z_sim_rk3(:,ind-1), u(:,ind-1), parameters);
        z_prime_rk3_1 = z_sim_rk3(:,ind-1) + Ts/2 * k1;
        k2 = Tractor_00_trail_model(z_prime_rk3_1, u(:,ind-1), parameters);
        z_prime_rk3_2 = z_sim_rk3(:,ind-1) + Ts * k2;
        k3 = Tractor_00_trail_model(z_prime_rk3_2, u(:,ind-1), parameters);
        z_sim_rk3(:,ind) = z_sim_rk3(:,ind-1) + Ts/6 * (k1 + 4*k2 + k3);

        errore_rk3(k,ind-1)=sqrt((z_sim_rk3(1,ind)-z_sim_ode(1,ind))^2+(z_sim_rk3(2,ind)-z_sim_ode(2,ind))^2);
    end
    time(k,3)=toc;                                % saving time of computation
     max_error(k, 2)= max(errore_rk3(k,:));       % saving maximum error
     mean_error(k, 2)= mean(errore_rk3(k,:));     % saving mean error

    %% RK4
        
    tic 
    % simulation of the system using RK4

    for ind = 2:Ns                                
        
        k1 = Tractor_00_trail_model(z_sim_rk4(:,ind-1), u(:,ind-1), parameters);
        k2 = Tractor_00_trail_model(z_sim_rk4(:,ind-1) + Ts/2 * k1, u(:,ind-1), parameters);
        k3 = Tractor_00_trail_model(z_sim_rk4(:,ind-1) + Ts/2 * k2, u(:,ind-1), parameters);
        k4 = Tractor_00_trail_model(z_sim_rk4(:,ind-1) + Ts * k3, u(:,ind-1), parameters);
        z_sim_rk4(:,ind) = z_sim_rk4(:,ind-1) + Ts/6 * (k1 + 2*k2 + 2*k3 + k4);
        
        errore_rk4(k,ind-1)=sqrt((z_sim_rk4(1,ind)-z_sim_ode(1,ind))^2+(z_sim_rk4(2,ind)-z_sim_ode(2,ind))^2);
    end
    time(k,4)=toc;                                % saving time of computation
     max_error(k, 3)= max(errore_rk4(k,:));       % saving maximum error
     mean_error(k, 3)= mean(errore_rk4(k,:));     % saving mean error
   
    %% simulation using FFD

    tic
    
    % simulation of the system using FFD

    for ind=2:Ns                                 
        
        zdot               =   Tractor_00_trail_model(z_sim_ffd(:,ind-1),u(:,ind-1),parameters);
        z_sim_ffd(:,ind)       =   z_sim_ffd(:,ind-1)+Ts*zdot;
        
        errore_ffd(k,ind-1)=sqrt((z_sim_ffd(1,ind)-z_sim_ode(1,ind))^2+(z_sim_ffd(2,ind)-z_sim_ode(2,ind))^2);
    end 
    time(k,5)=toc;                                % saving time of computation
    max_error(k, 4)= max(errore_ffd(k,:));        % saving maximum error
    mean_error(k, 4)= mean(errore_ffd(k,:));      % saving mean error


    % plot of the error for each numerical integration method and for the
    % actual considered sampling time

    figure()                                      
    plot(linspace(0, Tend, Ns-1), errore_rk2(k,:)); hold on;
    plot(linspace(0, Tend, Ns-1), errore_rk3(k,:))
    plot(linspace(0, Tend, Ns-1), errore_rk4(k,:))
    plot(linspace(0, Tend, Ns-1), errore_ffd(k,:))
    xlabel('Time [s]'),ylabel('Error [m]'),grid on;
    legend('RK2', 'RK3', 'RK4','FFD');
    title('Errors with respect to ODE45', ['Ts=' num2str(j(k))])

    
    % plot of the trajectory of the system obtained for each numerical integration method and for the
    % actual considered sampling time

    figure()
    plot(z_sim_rk2(1,:), z_sim_rk2(2,:));hold on;
    plot(z_sim_rk3(1,:), z_sim_rk3(2,:));hold on;
    plot(z_sim_rk4(1,:), z_sim_rk4(2,:));hold on;
    plot(z_sim_ffd(1,:), z_sim_ffd(2,:));hold on;
    plot(z_sim_ode(1,:), z_sim_ode(2,:));hold on;
    daspect([1,1,1]),xlabel('X [m]'),ylabel('Y [m]'),grid on;
    legend('RK2', 'RK3', 'RK4','FFD', 'ODE45');
    title('Comparison between trajectories', ['Ts=' num2str(j(k))])
end

%% Print of the table of the results

rownames= string(j);                     % name of the raws (sampling times analized)
colnames = {'RK2', 'RK3', 'RK4','FFD'};  % name of the columns (integration method)

% creation of the maximum error array using the scientific notation

max_error = arrayfun(@(x) sprintf('%.2e', x), max_error, 'UniformOutput', false);

% print of the table of maximum error

for i = 1:size(max_error, 4)
    fprintf('\nTable of maximum error for each integration method and each sampling time tested:\n', i);
    T = array2table(max_error(:,:,i), 'VariableNames', colnames,'RowNames', rownames); 
    disp(T); 
end

% creation of the mean error array using the scientific notation

mean_error = arrayfun(@(x) sprintf('%.2e', x), mean_error, 'UniformOutput', false); 

% print of the table of mean error

for i = 1:size(mean_error, 4)
    fprintf('\nTable of mean error for each integration method and each sampling time tested:\n', i);
    T = array2table(mean_error(:,:,i), 'VariableNames', colnames,'RowNames', rownames); 
    disp(T); 
end

colnames = {'ODE45','RK2', 'RK3', 'RK4','FFD'};

% creation of the simulation time array using the scientific notation

time = arrayfun(@(x) sprintf('%.2e', x), time, 'UniformOutput', false); 

% print of the table of simulation time

for i = 1:size(time, 5)
    fprintf('\nTable of computation time for each integration method and each sampling time tested:\n', i);
    T = array2table(time(:,:,i), 'VariableNames', colnames,'RowNames', rownames); 
    disp(T); 
end


