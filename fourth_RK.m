function x_k = fourth_RK(f, t, x0, h)

f1 = f(t, x0);
f2 = f((t+h/2), x0 + (h/2) * f1);
f3 = f((t + h/2), x0 + (h/2) * f2);
f4 = f((t+h), x0 + h*f3);
x_k = x0 + (h/6) * (f1 + 2*f2 + 2*f3 + f4);
end