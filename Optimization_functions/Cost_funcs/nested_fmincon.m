function [x,fval,eflag,outpt] = nested_fmincon(z0,parameters,Optimization_opt,constr_param,x0,A,b,Aeq,beq,lb,ub,opts)

if nargin == 1 % No options supplied
    opts = [];
end

xLast = []; % Last place computeall was called
myf = []; % Use for objective at xLast
myc = []; % Use for nonlinear inequality constraint
myceq = []; % Use for nonlinear equality constraint

fun = @objfun; % The objective function, nested below
cfun = @constr; % The constraint function, nested below

% Call fmincon
[x,fval,eflag,outpt] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,cfun,opts);

    function y = objfun(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,myc] = tractor_cost_constr(x,z0,parameters,Optimization_opt,constr_param);
            xLast = x;
        end
        % Now compute objective function
        y = myf; %+ 20*(x(3) - x(4)^2)^2 + 5*(1 - x(4))^2;
    end

    function [c,ceq] = constr(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,myc,myceq] = tractor_cost_constr(x,z0,parameters,Optimization_opt,constr_param);
            xLast = x;
        end
        % Now compute constraint function
        c = myc; % In this case, the computation is trivial
        ceq = myceq;
    end

end