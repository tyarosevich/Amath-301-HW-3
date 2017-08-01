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

