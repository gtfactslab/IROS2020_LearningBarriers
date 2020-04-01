function [h, u_bwdist] = barrier_bwdist(x,bwd,u_des)
    
    gamma = 1;  
    H = eye(2);
    
    h = bwd.evaluate(x);                        % Barrier Value
    grad_h = bwd.gradient(x);                   % Gradient of Barrier
    
    A = [grad_h(1,:), grad_h(2,:)];
      
    B = gamma*h^5;
    
    opts = optimoptions(@quadprog, 'Display', 'off');
    u_bwdist = quadprog(H, -u_des, A, B, [], [], [], [], [], opts);  % QP Solver
    
end
