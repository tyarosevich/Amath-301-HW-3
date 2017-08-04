function [x] = forward_euler(xt0, A, dt)

x = (eye(2, 2) + dt * A) * xt0;

end
