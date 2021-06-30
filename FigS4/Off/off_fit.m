function [error,y_t,y] = off_fit(t,data,SA,V,frac,Conc,S,k_off_p,k_on_l,k_off_l,k_off_d)
    tspan = [0 5];
    y0 = 0;
    %fun = @(t,y) [k_on_l * (S-y) * (1-y*SA/(Conc*602.2e6*V*frac)) * Conc - k_off_l*y];
    fun = @(t,y) [k_on_l * (S-y) * (1-y/Conc/602.2e6*3/7) * Conc - k_off_l*y];
    [~,y] = ode15s(fun, tspan, y0);

    y(end);
    y0 = [y(end),0];
    fun = @(t,y) [k_on_l * (S-y(1)-y(2)) * (1-(y(1)+y(2))*SA/(Conc*602.2e6*V*frac)) * Conc * exp(-k_off_p*t) - k_off_p * y(1) - k_off_l * y(1);k_off_p*y(1) - k_off_d*y(2)];
    sol = ode45(fun, [t(1),t(end)], y0);
    y = deval(sol,t);
    y = sum(y);
    y_t = y ./ max(y);
    d_t = data' ./ max(data);
    error = sum((d_t - y_t).^2);
end