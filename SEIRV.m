function dy = SEIRV(t,y,param)

    lambda = param(1);
    k      = param(2);
    delta  = param(3);
    sigma  = param(4);

    beta_E = param(5);
    beta_I = param(6);
    beta_R = param(7);

    alpha  = param(8);
    gamma  = param(9);

    dy = zeros(6,1);
    
    S = y(1);  
    E = y(2);      
    I = y(3);
    R = y(4);

    dy(1) = -lambda*S*I;
    dy(2) = lambda*S*I - k*E;                               
    dy(3) = k*E - delta*I;
    dy(4) = delta*I - sigma*R;
    dy(5) = alpha*(1-gamma)*(beta_E*E + beta_I*I + beta_R*R);
    dy(6) = delta*I; 
    dy(7) = alpha*(1-gamma)*(beta_E*E);
    dy(8) = alpha*(1-gamma)*(beta_I*I);
    dy(9) = alpha*(1-gamma)*(beta_R*R);
end