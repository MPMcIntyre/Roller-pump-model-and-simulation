%M.P.MCINTYRE
%MPMCINTYRE@YAHOO.COM
%SIMULINK VARIABLE GENERATOR
%2019-09-06

%{
This Matlab file generates the variables needed to run the simulink
simulation of a peristaltic pump (found in Lumped_parameter_sim.slx).

This file calls the flow_solver.m function in order to generate
the volume and induced flow of the roller coming into contact with the
process tube.

The first two sections of this script are USER INPUTS to describe the
design properties of the pump followed by the working conditions of the
pump.
%}

%Clear workspace (keep it tidy)
clear all
clc
tic                             %timer

%% User input values
%If you would like figures, enable this feature
plot_enabled = false;

%If you want to use the measured roller volume data, enable this feature
measured_roller_values = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%  SIMULATION VALUES   %%%%%%%%%%%%%%%%%%%%%%%%
%IMPORTANT: The simulink simulation time must have the same value as s_time
s_time = 1;                     %Total simulated time
h = 0.001;                      %Time step (increment)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DESIGN VALUES   %%%%%%%%%%%%%%%%%%%%%%%%%%
N = 170;                        %RPM
r_offset = 39.5;                %mm (exp = 39.5)
r_roller = 19.87;               %mm (exp = 19.87)
d_o = 14;                       %outer tube diameter in mm (exp = 14)
d_i = 10;                       %inner tube diameter in mm (exp = 10)
disengage_angle = 30;           %degrees (exp = 30)
r_b = 62.76;                    %designed radius of casing (exp = 62.76)
NU = 3;                         %Number of rollers on the pump

%%%%%%%%%%%%%%%%%%%%%   PARAMETER VALUES    %%%%%%%%%%%%%%%%%%%%%%%%%%
P_res = 88.6372;                %Reservoir pressure
R_in = 0.1108;                  %Inlet resistance
R_out = 0.1108;                 %Outlet resistance
C_in = 0.036108499;             %Inlet compliance
C_out = 0.036108499;            %Outlet compliance
L_in = 0.0042;                  %Inlet inertia
L_out = 0.0042;                 %Outlet inertia


%%                      ANGLE AND TIME VALUE GENERATION
%{
In this section the engagements and disengagement angles and  corresponding
times are created in order to generate the induced flow variables in the
following section. The angles an times are first generated for each roller
individually, then sorted from smallest to largest in a single row vector
that is then used for referencing in the Q_EDI and Q_ed_out loops.
%}

%Occlusion variable
n_occ = 1-((r_b - r_roller - r_offset-(d_o-d_i))/(d_i));
n_occ(n_occ>1) = 1;

%d_tube = 2*d_i;
omega = (N/60)*2*pi;      %rad/s
omega_deg = (N/60)*360;   %rad2deg(omega);


%Calling the Flow_solver function to obtain the induced flow values for
%incremental time value equal to that of h with the provided design
%criterea
[sim_engage_angle,A,txt_recomend,poly,V_function_f,cyl_vol,r_dynamic,aprx_V,t_sim] = Flow_solver(N,r_offset,r_roller,d_o,d_i,disengage_angle,r_b,h,NU,plot_enabled);
%Calling the measured roller volume function to determine the roller
%induced flow with the measured values
[Q_test,t_test,v_test,test_engage_angle,V_poly,test_vol_average,test_degrees] = Measured_roller_flow_solver(N,h);

%Obtaining the range of values of the flow_solver/approximation functions
inflating_angle = A(2,:);
v_sim =  fliplr(A(5,:));
Q_sim = -fliplr(A(4,:));
normal_angle = inflating_angle(A(3,:)<10);
lambda = A(2,end);

%Approximation: Cyllindrical method of approximation
X = r_roller*sind(lambda);
A_cyl = X*d_i;

%Approximation: Relational method of approximation
A_roller = (lambda/360)*pi*r_roller^2 - 0.5*X*(r_roller - d_i);
V_approx =  (A_roller/A_cyl)*pi*((d_i/2)^2)*2*X/1000;

%If the measured roller volume is used, then this part is used to determine
%the roller induced flow, else the approximated roller volume displacement
%is used
if measured_roller_values == true;
    volume = v_test;
    Q_flow = Q_test;
    engage_angle = test_engage_angle;
else
    volume = v_sim;
    Q_flow = Q_sim;
    engage_angle = sim_engage_angle;
end

%Cancel negative values
Q_flow(Q_flow <0) =0;

%The working area of the process tube for flow
a = pi*(d_i/2)^2;

%Q_nom is the nominal flow induced by the rotation of the motor
Q_nom = omega*((r_b- d_o/2))*a*n_occ /1000;
t_rot = 360/omega_deg;

%The angles at which the rollers engage and dissengage are dependant on the
%amount of rollers located around the rotating axis, assuming that they are
%evenly spaced
f_n = 360/NU;

%Let the pump start with the roller fully in contact with the tube,
%indicating that the first pulsation will be that of the disengaging roller
disengage_intervals = 0:f_n:360;

%The engage_angles variable is a container for how far the motor must
%rotate for each roller to come into contact with the tube within one rotation

%Counter for the amount of rollers
ticker = 0:NU-1;
clear engaging_angles
%Determination of the rotation required for the roller to come into contact
engaging_angles(:,1:NU) = 180-(2*disengage_angle + engage_angle) + ticker*360/NU;

%If an angle is larger than 360 degrees, which the last value should be, it
%implies that the last roller to come out of contact with the tube is the
%first roller to engage the tube (or could possible already be engaged)

%This is fixed by subtracting one rotation from the angle's value.
engaging_angles(engaging_angles>360-engage_angle) = engaging_angles(engaging_angles>360-engage_angle) - 360;
engaging_angles = sort(engaging_angles);

%Rollers that are already in contact will display these values by leaving
%the second row of the engageing_angles values larger than 360. These
%values will repeat as a time variable from 0:h:s_time.

%Resetting any values larger than one rotation of the disengage_angles
%variables
disengage_intervals(disengage_intervals>=360) = [];
disengaging_angles = disengage_intervals;

%Calculating the amount of rotations in the simulation period with two
%additional rotation for computaional purposes
N_sim_rotations = (N/60)*s_time+2;
%The N_rotations is rounded down ALWAYS to account for the extra rotation
remainder = rem(N_sim_rotations,floor(N_sim_rotations));
N_sim_rotations = N_sim_rotations - remainder;

%-1 value accounts for the first value already being calculated and
%automatically rounds
rotations = 1:N_sim_rotations-1;

%Reference angle for setting up a matrix with all of the engaging and
%disengaging angles
ref_angle_en = engaging_angles;
ref_angle_dis = disengaging_angles;

%The engaging and disengaging angle matrix is made by adding multiples of
%one rotation (360 degrees) onto the calculated engaging and disengaging
%angles
for i = 1:NU
    engaging_angles(2:N_sim_rotations,i) = ref_angle_en(i) + rotations*360';
    disengaging_angles(2:N_sim_rotations,i) = ref_angle_dis(i) + rotations*360';
end

%The engaging and disengaging angles break the loop statement for creating
%Q_EDI and Q_ed_out if the values are negative (first value) and therefor
%are set to zero

engaging_angles(engaging_angles<0) = 0;
disengaging_angles(disengaging_angles<0) = 0;

%The values of where the roller is either fully in contact or not in
%contact at all can be calculated by simply adding the engage_angle
engaged_angles = engaging_angles+engage_angle;
disengaged_angles= disengaging_angles+engage_angle;

%Converting the angles to time constants
engaging_times = engaging_angles/omega_deg;
disengaging_times = disengaging_angles/omega_deg;
engaged_times = engaged_angles/omega_deg;
disengaged_times = disengaged_angles/omega_deg;

%Resizing the matrix as a 1xn vector for use in the for loop
%Engaging values (important)
eng_t_vec = sort(reshape(engaging_times,1,numel(engaging_times)));
%Engaged values (additional user info)
engd_t_vec = sort(reshape(engaged_times,1,numel(engaged_times)));

%Disengaging values (important)
deng_t_vec = sort(reshape(disengaging_times,1,numel(disengaging_times)));
%Disengaged values (additional user info)
dengd_t_vec = sort(reshape(disengaged_times,1,numel(disengaged_times)));

%Engaging and disengaging angle vectors for additional user info
eng_ang_vec = sort(reshape(engaging_angles,1,numel(engaging_angles)));
deng_ang_vec = sort(reshape(disengaging_angles,1,numel(disengaging_angles)));


%%                      INDUCED FLOW VARIABLE GENERATION
%{
In this section the induced flow variables (Q_ed_in & Q_ed_out) caused by the
rollers coming into and out of contact are created over the full length of
the simulation period (s_time). The variables are initially made the length
of the simulation time devided into increment value h and given the value
of 0. A for loop is then used to find and allocate coresponding flow values
to the corresponding time and angle values within Q_ed_in and Q_ed_out
variables.

The for loop increases in increments of 1, from 1 to the the value of the
simulation time devided by the simulation increment h. This implies that
every iteration of the loop represents one of the incremental time values
of h on top of the previous iterations.
%}

%Time it takes for the simulation devided into increments of h
%Time vector to create the flow vector
t_vec = 0:h:s_time;                    %The t_vec is used to plot later on
Q_ed_in(1,1:numel(t_vec)) = 0;          %Induced flow at pump inlet
Q_ed_out(1,1:numel(t_vec)) = 0;         %Induced flow at pump outlet

%The period from peak to peak
period = (360/NU)/omega_deg;

%For the outlet, the induced flow is both reversed and negative (inversed)
Q_flow2 = fliplr(Q_flow);
Q_flow2 = -Q_flow2;

%A for loop is created in order to give the variables Q_ed_in and Q_ed_out
%values of Q_flow and Q_flow2 depending on the corresponing times and
%corresponding opening
counter1 = 1;               %Counter for inlet roller and rotation
counter2 = 1;               %Counter for outlet roller and rotation

for i=1:round(s_time/h)+1
    
    %Induced flow at inlet:
    if i == round(eng_t_vec(counter1)/h)+1
        Q_ed_in = Q_ED_I_assign(i,Q_flow,Q_ed_in,s_time,h);
        %When the if statement is triggered ans it assigns values to
        %Q_ed_in the counter increments
        counter1 = counter1 + 1;
    end
    
    %Induced flow at outlet:
    if i == round(deng_t_vec(counter2)/h)+1
        Q_ed_out = Q_ED_II_assign(i,Q_flow2,Q_ed_out,s_time,h);
        %When the if statement is triggered ans it assigns values to
        %Q_ed_out the counter increments
        counter2 = counter2 + 1;
    end
    
end

%Only one of the induced flows are subtracted as it represents the volume
%lost for one roller for a full rotation
Q_inline_flow = Q_nom - Q_ed_in;
V_gain(1) = 0;
for i = 2:numel(Q_inline_flow)
    V_gain(i) = V_gain(i-1) + h*Q_inline_flow(i);
end

Q_avg_flow = mean(Q_inline_flow);

%Angle between each roller during full contact in radians
contact_angle = 2*pi/NU;
%Length of the contact arc
contact_length = contact_angle*(r_b - (d_o/2));
%Volume of tube for contact angle in mL minus the halves of the roller
%volume present on either side of the section
contact_vol = (pi*((d_i/2)^2)*contact_length/1000 -max(volume))*NU;

%Average flow rate via numerical calculations
Q_avg = contact_vol*(N/60);

%%                  PLOTTING AND DISPLAY
if plot_enabled ==true
    %Plotting the Q_ed_in and Q_ed_out values
    disp('Plotting pretty images...')
    
    figure
    plot(t_vec, Q_ed_in,'k--')
    hold on
    plot(t_vec, Q_ed_out, 'k')
    %axis([period 3*period min(Q_flow2)*1.1 max(Q_flow)*1.1])
    xlabel('Time [s]')
    ylabel('Flow rate [mL/s]')
    legend('$Q_{ed\_in}$', '$Q_{ed\_out}$', 'Interpreter','latex','Location','southeast')
    hold off
    
    
    figure
    plot(test_degrees(1:end-2), test_vol_average(2:end-1),'k')
    hold on
    plot(t_sim*omega_deg,v_sim,'k--')
    plot([test_degrees(1) test_degrees(end-2)],[aprx_V aprx_V],'k:')
    plot([test_degrees(1) test_degrees(end-2)],[cyl_vol cyl_vol],'k-.')
    legend('Test average', 'Continuous approximation','interpreter','latex')
    xlabel('Degrees [^o]')
    ylabel('Volume [mL]')
    %axis([0 ceil(max([t_sim(end)*omega_deg test_degrees(end-2)])) 0 max([test_vol_average(end) v_sim(end)])])
    axis([0 ceil(max([t_sim(end)*omega_deg test_degrees(end-2)])) 0 3])

    hold off
    
end

%Presenting the data to the user of the program in an understandable manner
if measured_roller_values == true;
    txt = sprintf('Simulink_variable_generator complete at %.4f seconds',toc);
    disp(txt)
    valtxt = sprintf('Values obtained: \n');
    txt1 = sprintf('Simulation time - %g seconds, increment h = %.4f', s_time,h);
    txt2 = sprintf('Number of rotations = %.2f',N_sim_rotations-2+remainder);
    txt3 = sprintf('Engage angle = %.2f degrees', engage_angle);
    txt4 = sprintf('Occlusion = %.2f', n_occ);
    txt5 = sprintf('Integrated volume approximation = %.2f mL', max(volume));
    txt6 = sprintf('Max induced flow value = %.2f mL/s', max(Q_flow));
    txt7 = sprintf('Min induced flow value = %.2f mL/s', min(Q_flow));
    txt8 = sprintf('Average flow through pump = %.2f mL/s', Q_avg_flow);
    txt9 = sprintf('\nTotal program run time is %.4f seconds',toc);
    f = msgbox({valtxt,txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8,txt9},'Success','help');
    disp('Ready for Simulink simulation...')
else
    txt = sprintf('Simulink_variable_generator complete at %.4f seconds',toc);
    disp(txt)
    valtxt = sprintf('Values obtained: \n');
    txt1 = sprintf('Simulation time - %g seconds, increment h = %.4f', s_time,h);
    txt2 = sprintf('Number of rotations = %.2f',N_sim_rotations-2+remainder);
    txt3 = sprintf('Engage angle = %.2f degrees', engage_angle);
    txt4 = sprintf('Occlusion = %.2f', n_occ);
    txt5_1 = sprintf('Integrated volume approximation = %.2f mL', max(volume));
    txt5_12 = sprintf('Area volume approximation = %.2f mL', aprx_V);
    txt5_2 = sprintf('Renference cylindrical volume displaced = %.2f mL', cyl_vol);
    txt6 = sprintf('Max induced flow value = %.2f mL/s', max(Q_flow));
    txt7 = sprintf('Min induced flow value = %.2f mL/s', min(Q_flow));
    txt8 = sprintf('Average flow through pump = %.2f mL/s', Q_avg_flow);
    txt9 = sprintf('Recomended curve - %s, polynomial order %g', txt_recomend, poly);
    txt10 = sprintf('\nTotal program run time is %.4f seconds',toc);
    f = msgbox({valtxt,txt1,txt2,txt3,txt4,txt5_1, txt5_12,txt5_2,txt6,txt7,txt8,txt9,txt10},'Success','help');
    disp('Ready for Simulink simulation...')
end
