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

%% Exercise 2.b.)
trap_numerator = trapz(r, T_r_num);
trap_denominator = trapz(r, T_r_den);
A3 = trap_numerator/trap_denominator;
