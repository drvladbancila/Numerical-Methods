% Script for solving ODE with backward Euler (not the best method out
% there!)
% The ODE is given in the form a*x' + b*x = f(t).
% x(t) is then evaluated between 0 and Tmax using dt for the linear
% spacing. The function then outputs the coefficients of the n-degree
% polynomial that interpolates x(t) and the absolute error.
function c = backward_euler(a, b, f, x0, Tmax, dt, n)
    t = 0:dt:Tmax;
    
    % initialization
    x = zeros(1, length(t));
    x(1) = x0;
    
    % backward Euler
    % x' = (x(t_k) - x(t_(k-1)))/dt
    % a*((x(t_k) - x(t_(k-1)))/dt) + b*x(t_k) = f(t_k)
    % (a + b*dt)*x(t_k) = x(t_(k-1)) + dt*f(t_k)
    % x(t_k) = (x(t_(k-1)) + dt*f(t_k))/(a + b*dt)
    for k=2:length(t)
        x(k) = (a*x(k-1) + dt*f(t(k)))/(a + b*dt);
    end
    
    c = polyfit(t, x, n);
    z = polyval(c, t);
    err = norm(x - z, inf)/norm(x, inf)
    
    figure;
    hold all;
    title("ODE: ax'+bx=f(t)");
    plot(t, x);
    plot(t, z, '--', 'linewidth', 3);
    legend('Euler', 'Poly Interpolation');
    xlabel('t');
    ylabel('x(t)');
end