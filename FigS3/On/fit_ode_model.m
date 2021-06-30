function [error,y,y0] = fit_ode_model(t,data,k_on_l,S,KD,SA,V,C0)
    k_off_l = 0.1e-6*k_on_l;
    y0 = S*C0'./(C0'+KD*1e-6);
    %fun = @(t,y) [k_on_l * (S-y) .* (1-(y-y0)./C0/602.2e6*SA/V).*C0 - k_off_l*y]';
    %@(t,y) ydot_fun(t,y,y0,k_on_l,k_off_l,S,KD,SA,V,C0)
    sol = ode15s(@(t,y) ydot_fun(t,y,y0,k_on_l,k_off_l,S,KD,SA,V,C0), [t(1),t(end)], y0);
    y = deval(sol,t);
    error = sum(sum((data' - (y-y0)).^2));
end

