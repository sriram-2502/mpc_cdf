function b = obstacle_parabola(states)
% Inputs 
% states    : [x; y; x^2] -> casadi variables;
% outputs
% b(x)       : cbf function for circular obstalce

x1 = states(1);
x2 = states(2);

b = -(x2-x1^2-x1+1);

end
