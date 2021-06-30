function [ydot] = ydot_fun(t,y,y0,k_on_l,k_off_l,S,KD,SA,V,C0)
    ydot = [k_on_l * (S-y) .* (1-(y-y0)./C0'/602.2e6.*SA./V).*C0' - k_off_l*y];
end