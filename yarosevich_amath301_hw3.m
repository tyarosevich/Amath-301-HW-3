%% Homework 3 
% Exercise 1
clear all;close all;clc;

t = transpose(-90:10:110);
N = transpose([7.24 9.64 12.87 17.07 23.19 31.44 38.56 50.19 62.98 76.21 92.23 106.02 123.2 132.16 151.33 179.32 203.3 226.54 248.71 281.42 307.75]);
dt = 10;
A1 = zeros(21, 1);
for i = 2:20
    A1(i) = (N(i+1) - N(i-1))/(2*dt);
end
A1(1) = (1/(2*dt)) * ( -3*N(1) + 4*N(2) - N(3));
A1(21) = (1/(2*dt)) * (3*N(21) - 4*N(20) + N(19));
% test = zeros(21,1);
% for i = 2:21
%     test(i) = N(i-1) + A1(i-1)
% end
% plot(t, N, 'g')
% hold on
% plot(t(2:20), test(2:20), 'b')
save A1.dat A1 -ASCII
%% Exercise 2
clear all;clc;close all;
% As far as I can tell, THIS is the composite simpsons. It goes lockstep,
% accounting for repeats the whole time, and in that sense is different
% from the basic simpson's rule:
% integ a to b f(x) dx = h/3 * (f(x_o) + ( 2* sum j = 1 to n/2-1 f(x_2*j) ) + ( 4* sum j = 1 to n/2 f(x_2*j - 1) )  + f(x_n)
% note I dunno why it is written thsi way but the f(x_n) at the end is
% actually f(x_last term). Also the wikipedia article includes excellent
% python code at : https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson.27s_rule

r = transpose(.308:.017:.478);
T = transpose([640 794 885 943 1034 1064 1114 1152 1204 1222 1239]);

T_r_num = (r.*T) * 0.7051;
s_num = T_r_num(1) + T_r_num(end);
h = .017;
for i = 2:2:10
    s_num = s_num + 4*(T_r_num(i));
end
for i = 3:2:9
    s_num = s_num + 2*(T_r_num(i));
end
numerator = s_num * (h/3);

T_r_den = r*0.7051;
s_den = T_r_den(1) + T_r_den(end);
for i = 1:2:9
    s_den = s_den + 4*(T_r_den(i));
end
for i = 2:2:8
    s_den = s_den + 2*(T_r_den(i));
end
denominator = s_den * (h/3);

A2 = numerator/denominator;
save A2.dat A2 -ASCII


%% Exercise 2.b.)
trap_numerator = trapz(r, T_r_num);
trap_denominator = trapz(r, T_r_den);
A3 = trap_numerator/trap_denominator;

save A3.dat A3 -ASCII

%% Exercise 3.a.)
clc; clear all; close all;

% This looks something like x_k+1 = x_k + h*f(x_k)
% And dfor this problem, f(x) is the differential equation we have, h is
% .1, x0 = x_k and x0 = 0 incrementing in h to 1, to provide us with a
% trajectory of X values in relation to t, i.e. X(t).

f = @(t,x) -2 * (x / (.25 + x)) * x + 1.5 * sin(pi*t);
X_list = [];
X_list(1) = 1;
h = .1;
counter = 1;
for i = 0:.1:(1-h)
    X_list(counter+1) = X_list(counter) + h * f(i, (X_list(counter)));
    counter = counter + 1;
end
A4 = transpose(X_list);
save A4.dat A4 -ASCII

%% Exercise 3.b.)

X_list_RK = [];
X_list_RK(1) = 1;
counter =1;
for i = 0:.1:(1-h)
    X_list_RK(counter + 1) = fourth_RK(f, i, X_list_RK(counter), h);
    counter = counter + 1;
end
A5 = transpose(X_list_RK);
save A5.dat A5 -ASCII


%% Exercise 3.c.)

X_list2 = [];
X_list2(1) = 1;
h = .01;
counter = 1;
for i = 0:.01:(1-h)
    X_list2(counter+1) = X_list2(counter) + h * f(i, (X_list2(counter)));
    counter = counter + 1;
end

X_list_RK2 = [];
X_list_RK2(1) = 1;
h = .01;
counter = 0;
for i = 0:.01:(1-h)
    counter = counter + 1;
    X_list_RK2(counter + 1) = fourth_RK(f, i, X_list_RK2(counter), h);
end
A6 = [transpose(X_list2) transpose(X_list_RK2)];
save A6.dat A6 -ASCII


%% Exercise 3.d.)
[t, Xode] = ode45(f, [0; 1], 1);
% plot(t, Xode, 'b')
% hold on
% plot(0:.1:1, X_list, 'g')
% plot(0:.1:1, X_list_RK, 'k')
% plot(0:.01:1, X_list2, 'm')
% plot(0:.01:1, X_list_RK2, 'y')
A7 = [t Xode];
save A7.dat A7 -ASCII

%% Exercise 4.a.)
clc;clear all;close all;
A = [0 1;
    -5/4 -3];
save A8.dat A -ASCII

%% Exercise 4.b.)

dt_max = 4/5;
save A9.dat dt_max -ASCII

%% Exercise 4.c.)

dt1 = .8*dt_max;
x_t_List = zeros(2, length(0:dt1:50));
x_t_List(:,1) = [1; 0];
counter = 1;
for n = dt1:dt1:50
    x_t_List(:, (counter+1)) = forward_euler((x_t_List(:,(counter))),A,dt1);
    counter = counter + 1;
end
A10 = x_t_List(1,:);
save A10.dat A10 -ASCII

% plot(0:dt1:50, x_t_List(1,:), 'g')
% hold on
% plot(0:dt1:50, x_t_List(2,:), 'b')

%% Exercise 4.d.)
dt1 = 1.05*dt_max;
x_t_List2 = zeros(2, length(0:dt1:50));
x_t_List2(:,1) = [1; 0];
counter = 1;
for n = dt1:dt1:50
    x_t_List2(:, (counter+1)) = forward_euler((x_t_List2(:,(counter))),A,dt1);
    counter = counter + 1;
end
A11 = x_t_List2(1,:);
save A11.dat A11 -ASCII

% plot(0:dt1:50, x_t_List2(1,:), 'g')
% hold on
% plot(0:dt1:50, x_t_List2(2,:), 'b')

%% Exercise 4.e.)
t = 1:50;
y0 = [1;0];
[t,y] = ode45(@(t,y)pend(t,y), t, y0);

A12 = y(:,1);
A13 = y(:,2);
save A12.dat A12 -ASCII
save A13.dat A13 -ASCII


    
