function dy = pend(t,y)
dy(1,1) = y(2);
dy(2,1) = (-10/8)*y(1) - 3*y(2);
end
