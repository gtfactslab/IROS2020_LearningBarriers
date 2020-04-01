function [h, u] = nav_domain(x, dfx,dfy, g, u_des)
    
    gamma = 0.3;
    H = eye(2);
    
    [h, grad_x, grad_y] = gradient_compute(x, dfx,dfy, g);
    
    
    A1 = [grad_x, grad_y];
      
    B1 = gamma*((-h));   
    
    opts = optimoptions(@quadprog, 'Display', 'off');
    u = quadprog(H, -u_des, A1, B1, [], [], [], [], [], opts);  

end
