function [Q_test,t_test,v_test,engage_angle,V_poly,test_vol_average,test_degrees] = Measured_roller_flow_solver(N,h)

% load('flow_solver_variables.mat');
% load('Test_data.mat');
% load('test_deg_vals.mat');
load('Measured_roller_volume_data.mat');

% volume = fliplr(volume);
% n = numel(volume);
% range = linspace(0,engage_angle,n);
% range2 = linspace(0,Degrees(end),n);

for i=1:numel(Test1)
    av_MX(i,:) = [Test1(i),Test2(i),Test3(i),Test4(i),Test5(i),Test6(i),Test7(i)];
    test_average(i) = mean(av_MX(i,:));
end
test_vol_average = test_average;

Test_deg_vals = open('test_deg_vals.mat');
Degrees = Test_deg_vals.test_deg_vals;

% Degrees(end) = [];
% test_average(1) = [];
% Test1(1) = [];
% Test2(1) = [];
% Test3(1) = [];
% Test4(1) = [];
% Test5(1) = [];
% Test6(1) = [];
% Test7(1) = [];

% test_average(end) = [];
% Degrees(end) = [];
% Test1(end) = [];
% Test2(end) = [];
% Test3(end) = [];
% Test4(end) = [];
% Test5(end) = [];
% Test6(end) = [];
% Test7(end) = [];


V_fun1 = polyfit(Degrees,test_average, 3);
V_fun2 = polyfit(Degrees,test_average, 4);
V_fun3 = polyfit(Degrees,test_average, 5);


omega = (N/60)*(360); 
angle_test = 0:0.1:Degrees(end);
t_max = Degrees(end)/omega; %with omega in degrees per second
rot_time = angle_test/omega;

t_check = 0:h:t_max;

syms  time

V_tf_1(time) = V_fun1(1).*(omega.*time).^3 ...
    +V_fun1(2).*(omega*time).^2 ...
    +V_fun1(3).*(omega*time) ...
    +V_fun1(4);

V_tf_2(time) = V_fun2(1).*(omega.*time).^4 ...
    +V_fun2(2).*(omega*time).^3 ...
    +V_fun2(3).*(omega*time).^2 ...
    +V_fun2(4).*(omega*time) ...
    +V_fun2(5);

V_tf_3(time) = V_fun3(1).*(omega.*time).^5 ...
    +V_fun3(2).*(omega*time).^4 ...
    +V_fun3(3).*(omega*time).^3 ...
    +V_fun3(4).*(omega*time).^2 ...
    +V_fun3(5).*(omega*time) ...
    +V_fun3(6);


Q_t_i1(time) = diff(V_tf_1, time);
Q_t_i2(time) = diff(V_tf_2, time);
Q_t_i3(time) = diff(V_tf_3, time);

%Induced Flow rate check variables
Q_ed1 = double(Q_t_i1(t_check));
Q_ed2 = double(Q_t_i2(t_check));
Q_ed3 = double(Q_t_i3(t_check));

%volume1 check variables
V_c1 = double(V_tf_1(t_check));
V_c2 = double(V_tf_2(t_check));
V_c3 = double(V_tf_3(t_check));


test_degrees = Degrees;
Q_test = Q_ed3;
t_test = t_check;
v_test = V_c3;
engage_angle = Degrees(end);
V_poly = V_fun3; %Polynomial as a function of rotational angle
end