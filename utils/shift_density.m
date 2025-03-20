function [t0, x0, u0] = shift_density(T, t0, x0, u, f, rho)
st = x0;
con = u(1,:)'/rho(1);


f_value = f(st,con);
st = st + (T*f_value);

x0 = full(st);
t0 = t0 + T;
u0 = [u(2:size(u,1),:); u(size(u,1),:)];
end