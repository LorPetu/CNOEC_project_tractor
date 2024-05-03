function [F,zsim] = modified_tractor_cost_function_F(u, theta, z0, zf, Ts, Q, Qf, Qdot, R, alpha_barrier, beta_barrier)


asat    =   theta(3,1);                 % maximum accelleration (m/s^2)   
deltamax=   theta(4,1);                 % maximum steering angle (rad)
deltamin=   theta(5,1);                 % minimum steering angle (rad)
vsat    =   theta(6,1);                 % maximum velocity (m/s)


N=length(u)/2;

zsim        =   zeros(4,N+1);
zsim(:,1)   =   z0;
ztemp       =   z0;


%inizializzazione funzione di costo gauss- newton

F           =   [   zeros(4*N,1); ...                                          % Tracking error weight 
                    zeros(4*N,1); ...                                       % Derivative of the states
                    zeros(2*N,1)]; ...                                       % Actuator effort weight
                    % zeros(2*N,1); ...                                       % Velocity saturation Barrier function
                    % zeros(2*N,1); ...                                       % Steering angle saturation Barrier function
                    % zeros(2*N,1)];                                           % Accelleration saturation Barrier function
F_temp = zeros(10+0,N);

for ind=2:N+1
    % Update the states with Foreward Euler

    u_k                     = [u(ind-1,1); u(N+ind-1,1)];
    zdot                    = Tractor_00_trail_model(ztemp,u_k, theta);
    ztemp                   = ztemp + Ts*zdot;
    zsim(:,ind)             = ztemp;
    e                       = zf - ztemp;
   
   % Update the cost function 
    F_temp(1:4,ind-1)                       = Q*e(:,1);
    F_temp(5:8,ind-1)                       = Qdot*zdot;
    F_temp(9:10,ind-1)                      = R*u_k;
    % F_temp(11,ind-1)                        = alpha_barrier*exp(-beta_barrier*(vsat+ztemp(4,1)));
    % F_temp(12,ind-1)                        = alpha_barrier*exp(-beta_barrier*(vsat-ztemp(4,1)));
    % F_temp(13,ind-1)                        = alpha_barrier*exp(-beta_barrier*(deltamin+u(2*(ind-2)+1,1)));
    % F_temp(14,ind-1)                        = alpha_barrier*exp(-beta_barrier*(deltamax-u(2*(ind-2)+1,1)));
    % F_temp(15,ind-1)                        = alpha_barrier*exp(-beta_barrier*(asat+u(2*(ind-1),1)));
    % F_temp(16,ind-1)                        = alpha_barrier*exp(-beta_barrier*(asat-u(2*(ind-1),1)));
    % 
    %F((ind-2)*16+1:(ind-1)*16) = F_temp;

end

% Terminal Cost Term

    %F((ind-2)*16+1:(ind-2)*16+4) = Qf*e(:,1);

    F_temp(1:4,N) = Qf*e(:,1);

    F(1:N,1) = F_temp(1,:)';

    [n_terms, ~] = size(F_temp);

for ind = 2:n_terms

    F((ind-1)*N+1:ind*N,1) = F_temp(ind,:)';

end


