%MICHAEL.P.MCINTYRE
%MPMCINTYRE@YAHOO.COM
%volume1 AND Flow rate APPROXIMATION
%2019-09-06

%{
This script is intended to utelise the design parameters of a peristaltic
pump in order to approximate the volume displaced inside the process tube
by the roller. This script is utelised as a function of the
Main.m script in order to obtain values needed for
the simulink pressure pulsation simulation of the peristaltic pump. 
%}
%%
function  [engage_angle,A,txt_recomend,poly,V_function_f,cyl_vol,r_dynamic,aprx_V,t_sim] = Flow_solver(rpm,r_inducer,r_roller,d_o,d_i,disengage_angle,r_casing,h,NU,plot_enabled)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Inputs can be uncommented and adjusted independantly during script editing
or simply using this script independantly from the
Main.m file. 
%}



% clear all
% clc
% tic
% 
% plot_enabled = false;
% 
% 
% h=0.001;
% rpm = 300;                         %RPM
% r_inducer = 39.5;                 %mm
% r_roller = 19.87;                  %mm
% d_o = 14;                       %outer tube diameter in mm
% d_i = 10;                       %inner tube diameter in mm
% disengage_angle = 15;           %degrees
% motor_torque = 3;               %Nm
% r_ring = 7.5;                   %Vertical radius of roller in mm
% r_casing = 62.76;                  %designed radius of casing
% NU = 2;           


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   ERROR CHECK   %%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Occlusion is calculated based on the outer casing,wall thickness,
roller radius, and inducer radius.
If the occlusion is less than zero, it implies that the roller does not
come into contact with the process tube, the occlusion can not be more than
one as 1 implies that the tube is fully occluded
%}
occlusion = 1-((r_casing - r_roller - r_inducer-(d_o-d_i))/(d_i));
occlusion(occlusion>1) = 1;


if (occlusion <=0)
    error('Design parameters are faulty, occlusion does not occur')
end

%Error if the end point of the roller exceeds the outer casing radius
if (r_casing<(r_inducer+r_roller))
    error('Casing radius cannot be smaller than endpoint of roller')
end

%{
This script allows for full compression of the process tube (seen by the
end point error calculation), however a warning is given to the user when
the compression seems unreasonable (comp_ratio>=50%) where the ratio is
determined by the tube wall thickness and the distance between the end
point of the roller and the casing
%}

%t is the wall thickness
t = (d_o-d_i)/2;

%determine the amount the tube's walls are compressed in mm
comp = 2*t - (r_casing - r_inducer - r_roller);

%Warn user if compression seems unreasonable (>50%)
comp_ratio = (comp/(2*t))*100;
percent_str = '%';

if (comp_ratio>=50)
    warning('Compression seems unreasonable  ---> %s %s',num2str(comp_ratio),num2str(percent_str));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   CALCULATIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%
%{
This section calculates the trigonometric values of the pump design
parameters, most importantly the engage_angle and the distance from the
roller's edge to the the outer casing's edge as it comes out of
contact with the process tube (d1)
%}

%Side length calculation
a = r_casing - d_o - r_roller;
b = r_inducer;
c = sqrt((b^2-a^2));

%Determine inner angle
alpha = asind(a/b);

%Determine working angle (engageing angle)
engage_angle = 90-alpha;

%Angle variable for calculations during engagement
%NOTE! - The angle variable is from 0->engage_angle , implying that this
%flow is for the roller disengaging the tube
angle = 0:0.1:engage_angle;

%The dynamic radius defines the radius from the leading point of the roller
%(the point closest to the outer casing wall) to the central axis of
%rotation and is used for validation

r_dynamic = sqrt((r_roller+ r_inducer.*sind(90-angle)).^2+(r_inducer.*cosd(90-angle)).^2);

%The distance d1 defines, from the point of contact of the roller and the
%process tube, the current distance between the leading edge of the roller
%and the minimum between the distance the roller achivieves from the outer
%radius of the casing

delta = (r_casing - (r_roller+r_inducer.*sind(90-angle)));
d1 = d_o - delta;


%{
If compression of the process tube occurs, the distance will have negative
values equal to that of the compression in mm. This negative distance does
not cause the roller to take up more volume1 and can thus be negated to
reduce computational error by making the negative values equal to zero.
%}

%%
%%%%%%%%%%%%%%%%%%%%%%%%   volume1 DETERMINATION  %%%%%%%%%%%%%%%%%%%%%%%%
%{
This section approximates the volume1 displaced by the roller and finds
different polynomial curve to fit the aproximation. The best approximation
should be the higher order polynomial, however this is not true for all
cases, thus the polynomials are compared to give a visual indication in the
plotting section.
%}

%Inner diameter of the process tube
r_ID = d_i/2;

%Inflating angle refers to the central contact point of the roller
%to the point where the roller edge no longer makes contact with the tube
inflating_angle = acosd((r_roller-(delta-(2*t-comp)))/r_roller);

%The A matrix is a control matrix where the variables generated can quickly
%be viewed and compared to other variables

n = numel(angle);

%Elpise variables:
%{
    The circumference is assumed to stay constant during deformation, this
    implies that the enlongated radius (r_long) can be determined with the
    elliptical circumference approximation p = 2*pi*sqrt((a^2+b^2)/2)
    where a is the short radius and b the long radius. The short radius is
    varied from r_ID to 0 based on the distance d1 (now recalculated with
    the inflating angle as d2) while the circumference is constant.
%}

% % Cylindrical and relational approximation:
%The roller angle defines the angle from the point of maximum distance
%(d1) to the point where the roller no longer makes contact with the tube
roller_angle = acosd((r_roller-d1)/r_roller);
max_alpha = max(roller_angle);
lin_length = r_roller*sind(roller_angle);
max_lin_length = max(lin_length);

%Cylindrical volume1 over area (for reference value)
cyl_vol = pi*(r_ID^2)*max_lin_length*2/1000;

A_c  = max_lin_length*d_i;
A_r_max = (max_alpha/360)*pi*r_roller^2 - 0.5*max_lin_length*(r_roller - d_i);
aprx_V = cyl_vol*(A_r_max/A_c);

A_r = (roller_angle./360).*pi.*r_roller.^2 - 0.5.*lin_length.*(r_roller - d_i);
V_aprx_1 = cyl_vol*(A_r/A_c);


%Setting up initial values to create reference matrices
dx = max_lin_length/n;

%The vector xinterp ranges vertically depending on the distance covered by
%the roller (determined by d1) and gives the linear length value (horisontal
%length, where d1 is the vertical length)
xinterp = [0:dx:max_lin_length];

%The complimentary angle of each horizontal value is found
comp_dx_angle = (asind((xinterp)/r_roller));
%The complimentary length describes the inverse theoretical distance of the
%roller over the engaging angle and roller angle
comp_length = r_roller-r_roller*cosd(comp_dx_angle);

%When the roller comes into contact, the roller bends the tube acros
%roughly the same distance as the maximum horisontal length. Thus a second
%length can be created to accompany an attempt to a more realistic model
ref_length = comp_length(:,1);

for o = 1:n
    for i = 1:n
    length(o,i) = d1(o) -(r_roller - r_roller*cosd(comp_dx_angle(i)));
    end
end

length(length<0) = 0;
length(length>d_i) = d_i;

%The values that are larger than the inner diameter of the tube is where
%the tube is compressed between the roller and outer casing, and do not add
%additional volume1 as the radius of the tube is fixed
comp_length(comp_length>d_i) = d_i;


%Constant circumference of the collapsing tube
circumference = 2*pi*r_ID;

%Calculating the long and short radii of the elipse
r_short = r_ID -length/2;
r_long= sqrt(2*(circumference/(2*pi))^2 - r_short.^2);


%Area of the ellipse 
Area = pi*r_ID^2 -  pi.*r_short.*r_long;

%Sectional volume displaced by roller
vol =  (2*Area.*dx)/1000;
%Total volume displaced by roller
for i = 1:n
    volume(i) = sum(vol(i,:));
end


%%
%%%%%%%%%%%%%%%%%%%%%%   POLYNIMIAL DETERMINATION  %%%%%%%%%%%%%%%%%%%%%%
%Creating polynomial fitting functions for the volume1 over the angle of
%engagement
V_fun1 = polyfit(angle,volume, 3);
V_fun2 = polyfit(angle,volume, 4);
V_fun3 = polyfit(angle,volume, 5);

%Polynomial volume1 function
syms time

omega = (rpm/60)*(360);                 %RPM in degrees per second

V_function1(time) = V_fun1(1).*(omega.*time).^3 ...
    +V_fun1(2).*(omega*time).^2 ...
    +V_fun1(3).*(omega*time) ...
    +V_fun1(4);

V_function2(time) = V_fun2(1).*(omega.*time).^4 ...
    +V_fun2(2).*(omega*time).^3 ...
    +V_fun2(3).*(omega*time).^2 ...
    +V_fun2(4).*(omega*time) ...
    +V_fun2(5);

V_function3(time) = V_fun3(1).*(omega.*time).^5 ...
    +V_fun3(2).*(omega*time).^4 ...
    +V_fun3(3).*(omega*time).^3 ...
    +V_fun3(4).*(omega*time).^2 ...
    +V_fun3(5).*(omega*time) ...
    +V_fun3(6);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   INDUCED FLOW   %%%%%%%%%%%%%%%%%%%%%%%%%%
%{
This section creates vectors containing values for the polynomials in order
to determine whether or not the values seem accurate. These values are
plotted against one another in the plot section.
%}

%Diffential values for the volume1 functions created
Q_induced1(time) = diff(V_function1, time);
Q_induced2(time) = diff(V_function2, time);
Q_induced3(time) = diff(V_function3, time);



%Time for the roller to rotate the engage_angle distance
t_max = engage_angle/omega; %with omega in degrees per second
rot_time = angle/omega;

t_check = 0:h:t_max;

%Induced Flow rate check variables
Q_ed1 = double(Q_induced1(t_check));
Q_ed2 = double(Q_induced2(t_check));
Q_ed3 = double(Q_induced3(t_check));

%volume1 check variables
V_check1 = double(V_function1(t_check));
V_check2 = double(V_function2(t_check));
V_check3 = double(V_function3(t_check));

%Save important variables in a container for external use
save('flow_solver_variables.mat', 'engage_angle', 'Q_ed3','V_function3', 'volume','t_check')

t_sim = t_check;

%The optimal curve is selected by the error between the test volume1 values
%created by the polynomial fitting and the actual volume1 data
vec_length_test = 1:numel(V_check1);
vec_length_req = linspace(1,numel(V_check1),n);

v_test1 = volume - interp1(vec_length_test,V_check1,vec_length_req) ;
v_test2 = volume -interp1(vec_length_test,V_check2,vec_length_req) ;
v_test3 = volume -interp1(vec_length_test,V_check3,vec_length_req) ;


% 
% diff1 = diff(Q_induced1,time);
% diff2 = diff(Q_induced2,time);
% diff3 = diff(Q_induced3,time);

%Selecting choice
ymin_vec = abs(double([max(v_test1) max(v_test2) max(v_test3)]));

%Colour and width of line for indicative plotting
width1 = 0.2;
width2 = 0.2;
width3 = 0.2;
c1 = 'r--';
c2 = 'b--';
c3 = 'g--';

%Conditional plotting to indicate the best suited induced Flow rate curve and
%polynomial degree. The curve with the final value closest to zero
%generally is the best suited curve due to the nature of the expected curve
%(the roller should not cause negative Flow rate).
if (ymin_vec(1)<ymin_vec(2)) && (ymin_vec(1)<ymin_vec(3))
    txt_recomend = 'Q1';
    poly = numel(V_fun1) - 1;
    c1 = 'k';
    width1 = 2;
    V_function_f = V_function1;
    Q_flow = Q_ed1;
    return
else
    if (ymin_vec(2)<ymin_vec(1)) &&(ymin_vec(2)<ymin_vec(3))
        txt_recomend = 'Q2';
        poly = numel(V_fun2) - 1;
        c2 = 'k';
        width2 = 2;
        V_function_f = V_function2;
        Q_flow = Q_ed2;
    else
        txt_recomend = 'Q3';
        poly = numel(V_fun3) - 1;
        c3 = 'k';
        width3 = 2;
        V_function_f = V_function3;
        Q_flow = Q_ed3;
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   PLOTTING   %%%%%%%%%%%%%%%%%%%%%%%%%%
%{
This section creates plots based off of the design parameters and the
values assosiated with the design. These plots are made to give an
indication to the user if there are any discrepancies within the design
parameter or in the script.
%}

if plot_enabled == true;
%figure 1 - Physical dimentions:
%Outer circumference of tube
theta = 90-disengage_angle:270+disengage_angle;
x_total_distance = r_casing*sind(theta);
y_total_distance = r_casing*cosd(theta);

%Backplate indication
x_backplate = 1.01*r_casing*sind(theta);
y_backplate = 1.01*r_casing*cosd(theta);

%Inner circumference of
x_tube = (r_casing-d_o)*sind(theta);
y_tube = (r_casing-d_o)*cosd(theta);

%Pipe exit - straight
r_t = r_casing-d_o;
x_cont1 = [-(r_t)*cosd(disengage_angle),-(r_t)*cosd(disengage_angle)+(r_t)*cosd(90-disengage_angle)];
y_cont1 = [(r_t)*sind(disengage_angle),(r_t)*sind(disengage_angle)+(r_t)*sind(90-disengage_angle)];
x_cont2 = [-(r_casing)*cosd(disengage_angle),-(r_casing)*cosd(disengage_angle)+(r_t)*cosd(90-disengage_angle)];
y_cont2 = [(r_casing)*sind(disengage_angle),(r_casing)*sind(disengage_angle)+(r_t)*sind(90-disengage_angle)];

%uncovered distance of rotation
theta2 = 0:360;
x_uncovered = r_casing*sind(theta2);
y_uncovered = r_casing*cosd(theta2);

%triroller arm
x_tri = [-r_inducer*cosd(disengage_angle+engage_angle),0];
y_tri = [r_inducer*sind(disengage_angle+engage_angle),0];

%Roller position (disengaged)
tri_x = -r_inducer*cosd(disengage_angle+engage_angle);
tri_y = r_inducer*sind(disengage_angle+engage_angle);
x_roller = tri_x +(r_roller)*cosd(theta2);
y_roller = tri_y +(r_roller)*sind(theta2);

tri_x2 = -r_inducer*cosd(disengage_angle);
tri_y2 = r_inducer*sind(disengage_angle);
x_roller2 = tri_x2 +(r_roller)*cosd(theta2);
y_roller2 = tri_y2 +(r_roller)*sind(theta2);

% Roller position (engaged)
for i = 1:NU
tri_x2(i) = -r_inducer*cosd(disengage_angle+((360/NU)*(i-1)));
tri_y2(i) = r_inducer*sind(disengage_angle+((360/NU)*(i-1)));
x_roller2(i,:) = tri_x2(i) +(r_roller)*cosd(theta2);
y_roller2(i,:) = tri_y2(i) +(r_roller)*sind(theta2);
end

% % Roller position (engaged)
% for i = 1:NU
% tri_x2(i) = -r_inducer*cosd(disengage_angle+((360/NU)*(i-1)));
% tri_y2(i) = r_inducer*sind(disengage_angle+((360/NU)*(i-1)));
% x_roller2(i,:) = tri_x2(i) +(r_roller)*cosd(theta2);
% y_roller2(i,:) = tri_y2(i) +(r_roller)*sind(theta2);
% end



%Engaging angle reference
x_ref1 = [-70*cosd(disengage_angle), 0];
y_ref1 = [70*sind(disengage_angle), 0];

%Engaging angle reference 2
x_ref2 = [70*cosd(disengage_angle), 0];
y_ref2 = [70*sind(disengage_angle), 0];

%theta angle display and trianlge
thet_disp = disengage_angle:disengage_angle+engage_angle;
x_t = -10*cosd(thet_disp);
y_t = 10*sind(thet_disp);
x_t2 = [tri_x,(r_casing-r_roller-d_o)*cosd(180-disengage_angle)];
y_t2 = [tri_y,(r_casing-r_roller-d_o)*sind(180-disengage_angle)];
x_t3 = [(r_casing-r_roller-d_o)*cosd(180-disengage_angle),0];
y_t3 = [(r_casing-r_roller-d_o)*sind(180-disengage_angle),0];

%zero-degree reference
xx = [-30,0];
yy = [0,0];
thet_xx = 0:disengage_angle;
xx_angle = -10*cosd(thet_xx);
yy_angle = 10*sind(thet_xx);

%Co-angle indication
x_i = -r_inducer*cosd(engage_angle+disengage_angle);
y_i = r_inducer*sind(engage_angle+disengage_angle);
x_ii = x_i-r_roller*cosd(disengage_angle);
y_ii = y_i+r_roller*sind(disengage_angle);
x_co_angle = [x_ii,x_ii-(r_roller+a)*cosd(180 - disengage_angle)];
y_co_angle = [y_ii,y_ii-(r_roller+a)*sind(180 - disengage_angle)];
xx_co_angle = [0,c*cosd(90-disengage_angle)];
yy_co_angle = [0,c*sind(90-disengage_angle)];

dynamix_radius_x = [0,x_i-r_roller*cosd(engage_angle)];
dynamix_radius_y = [0,y_i+r_roller*sind(engage_angle)];

%Tube wall plotting
dfrac = (r_casing-(d_o-d_i)/2)/r_casing;
dfrac2 = (r_casing-(d_o)+(d_o-d_i)/2)/r_casing;

%Full length sectional circle
theta_xx = 270+disengage_angle-15:270+disengage_angle+15;
x_total_sectional = (r_casing - (d_o-d_i))*sind(theta_xx);
y_total_sectional = (r_casing - (d_o-d_i))*cosd(theta_xx);


figure
hold on
plot(x_total_distance,y_total_distance,'k')
plot(x_backplate,y_backplate,'k','LineWidth',3)
plot(x_total_distance*dfrac,y_total_distance*dfrac,'k--')
plot(x_total_distance*dfrac2,y_total_distance*dfrac2,'k--')
plot(x_total_sectional,y_total_sectional,'r--')
plot(x_ref1,y_ref1, '-.b')
plot(x_ref2,y_ref2, '-.b')
plot(x_uncovered,y_uncovered,':k')
plot(x_tube,y_tube,'-.k')
plot(x_cont1,y_cont1,'-.k')
plot(x_cont2,y_cont2,'-.k')
% plot(x_tri,y_tri,'-k')
% plot(x_roller,y_roller,':k')
%plot(x_roller2,y_roller2,'k')
for i = 1:NU
plot(x_roller2(i,:),y_roller2(i,:),'k')
end
% plot(x_t,y_t, 'b')
% plot(x_t2,y_t2, 'b')
% plot(x_t3,y_t3, 'b')
% plot(xx,yy,'-.b');
% plot(xx_angle,yy_angle,'-.');
% plot(x_co_angle,y_co_angle,'-.r');
% plot(xx_co_angle,yy_co_angle,'-.r');
hold off

xticks([])
yticks([])% 


%figure 2 - inflating angle depiction:
%inflation angle max
max_diff_occ = (r_casing - r_roller - r_inducer-(d_o-d_i));

figure
hold on
%Roller
plot(r_roller*cosd(theta2),r_roller*sind(theta2)+r_roller+delta(1),'k-', 'LineWidth', 1.5)

%Tube walls
plot([-r_roller*1.5 r_roller*1.5],[0 0],'r-', 'LineWidth', 1.5)
plot([-r_roller*1.5 r_roller*1.5],[d_o d_o],'r-', 'LineWidth', 1.5)
plot([-r_roller*1.5 r_roller*1.5],[t t],'r--', 'LineWidth', 1.5)
plot([-r_roller*1.5 r_roller*1.5],[t+d_i t+d_i],'r--', 'LineWidth', 1.5)
plot([-r_roller*0.5 r_roller],[(d_o-d_i) (d_o-d_i)],'r--', 'LineWidth', 1.5)


%plotting increments
for i = 1:31:n
    plot([max_lin_length-xinterp(1,i) max_lin_length-xinterp(1,i)],[d_o-length(1,n-i) d_o],'k:')
    plot([(max_lin_length-xinterp(1,i)) (max_lin_length-xinterp(1,i))],[2*t (d_i+2*t)-length(1,n-i)],'b--')
    plot([0 max_lin_length-xinterp(1,i)],[r_roller+(2*t-comp) (d_i+2*t)-length(1,n-i)] ,'k--')
end

%Inflating angle
plot([0 0],[r_roller+(2*t-comp) (2*t-comp)],'b')
plot_angle = 270- inflating_angle(end);
x_point = -r_roller*cosd(plot_angle);
y_point = r_roller+(2*t-comp) + r_roller*sind(plot_angle);
plot([0 x_point],[r_roller+(2*t-comp) y_point],'b')
axis([-r_roller*1.5 r_roller*1.5 -5 r_roller*2.5])
hold off


%figure 3 - Tube colapse indication (ellipsation):
%Inner walls
%Half collapsed ellipse sketch
a_ellipse = r_short(1,round(n/2));
b_ellipse = r_long(1,round(n/2));
%Uncollapsed ellipse sketch
a1_ellipse = r_short(1,1);
b1_ellipse = r_long(1,1);
%Fully collapsed ellipse sketch
a2_ellipse = r_short(1,n);
b2_ellipse = r_long(1,n);

%Outer wall
%Half collapsed ellipse sketch
aa_ellipse = r_short(1,round(n/2))+t;
bb_ellipse = r_long(1,round(n/2))+t;
%Uncollapsed ellipse sketch
aa1_ellipse = r_short(1,1)+(t-comp/2);
bb1_ellipse = r_long(1,1)+(t-comp/2);
%Fully collapsed ellipse sketch
aa2_ellipse = r_short(1,n)+t;
bb2_ellipse = r_long(1,n)+t;


figure
subplot(1,3,1)
hold on
plot(bb2_ellipse*cosd(theta2),aa2_ellipse*sind(theta2),'k')
plot(b2_ellipse*cosd(theta2),a2_ellipse*sind(theta2),'k')
axis([-r_ID*2 r_ID*2 -r_ID*2 r_ID*2])
hold off

subplot(1,3,2)
hold on
plot(bb_ellipse*cosd(theta2), aa_ellipse*sind(theta2),'k')
plot(b_ellipse*cosd(theta2), a_ellipse*sind(theta2),'k')
axis([-r_ID*2 r_ID*2 -r_ID*2 r_ID*2])
hold off

subplot(1,3,3)
hold on
plot(bb1_ellipse*cosd(theta2),aa1_ellipse*sind(theta2),'k')
plot(b1_ellipse*cosd(theta2),a1_ellipse*sind(theta2),'k')
axis([-r_ID*2 r_ID*2 -r_ID*2 r_ID*2])
hold off

%
%figure 5 - volume1 and induced Flow rate co-plots for comparison
figure
subplot(1,3,1)
plot(t_check,Q_ed1)
ylabel('Flow rate [$\frac{m\ell}{s}$]','Interpreter','latex')
yyaxis right
hold on
plot(t_check,V_check1)
plot(rot_time, volume, 'k:')
hold off
xlabel('Time [s]','Interpreter','latex')
ylabel('volume1 [$m\ell$]','Interpreter','latex')
title('3rd degree polynomial')

subplot(1,3,2)
plot(t_check,Q_ed2)
ylabel('Flow rate [$\frac{m\ell}{s}$]','Interpreter','latex')
yyaxis right
hold on
plot(t_check,V_check2)
plot(rot_time, volume, 'k:')
hold off
xlabel('Time [s]','Interpreter','latex')
ylabel('volume1 [$m\ell$]','Interpreter','latex')
title('4th degree polynomial')

subplot(1,3,3)
plot(t_check,Q_ed3)
ylabel('Flow rate [$\frac{m\ell}{s}$]','Interpreter','latex')
yyaxis right
hold on
plot(t_check,V_check3)
plot(rot_time, volume, 'k:')
hold off
xlabel('Time [s]','Interpreter','latex')
ylabel('volume1 [$m\ell$]','Interpreter','latex')
title('5th degree polynomial')

figure
hold on
plot(t_check,V_check1,'k:')
plot(t_check,V_check2, 'k -o')
plot(t_check,V_check3, 'k--')
plot(rot_time, volume, 'k')
hold off
xlabel('Time [s]')
ylabel('Volume [mL]')
legend('3rd degree polynomial','4th degree polynomial','5th degree polynomial','Integration volume')


% figure 5 - Polynomial flow rate values and indication:


%plotting figure 5
figure
hold on
% plot(t_check,Q_ed1,c1,'LineWidth', width1)
% plot(t_check,Q_ed2,c2,'LineWidth', width2)
% plot(t_check,Q_ed3,c3,'LineWidth', width3)
plot(t_check,Q_ed1,'k--')
plot(t_check,Q_ed2,'k-.')
plot(t_check,Q_ed3,'k')
hold off
legend('3rd degree polynomial','4rth degree polynomial','5th degree polynomial')
title('Induced flow rate curve comparison')
% xlabel('Time [s]','Interpreter','latex')
% ylabel('Flow rate [$\frac{m\ell}{s}$]','Interpreter','latex')
xlabel('Time [s]')
ylabel('Flow rate [mL/s]')
hold off

%
%figure 6 - volume1 displaced over distance and inflating angle:

%length = rot90(length);

remap = fliplr(length);
remap_mat1 = [remap length];
remap_mat1 = flipud(remap_mat1);

% x_vecL = -delta_x;
% x_vecR = fliplr(delta_x);
% x_range = [x_vecL x_vecR];
% 
% contour(remap_mat1)


figure
surf(remap_mat1,'LineStyle','none');
% xticks([0 5 10])
% xticklabels('d = 0','d = 5','d = 10')

yticks([1 round(n/4) round(n/2) round(n - n/4) n])
y_tick5 = sprintf('%.1f',angle(1));
y_tick4 = sprintf('%.1f',angle(round(n/4)));
y_tick3 = sprintf('%.1f',angle(round(n/2)));
y_tick2 = sprintf('%.1f',angle(round(n - n/4)));
y_tick1 = sprintf('%.1f',angle(n));
yticklabels(char(y_tick1 , y_tick2 , y_tick3, y_tick4, y_tick5))

xticks([1 round(2*n/4) round(n) round(2*n - 2*n/4) 2*n])
x_tick1 = sprintf('%.1f',-inflating_angle(n));
x_tick2 = sprintf('%.1f',-inflating_angle(round(n/2)));
x_tick3 = sprintf('%.1f',inflating_angle(1));
x_tick4 = sprintf('%.1f',inflating_angle(round(n/2)));
x_tick5 = sprintf('%.1f',inflating_angle(n));
xticklabels(char(x_tick1 , x_tick2 , x_tick3, x_tick4, x_tick5))


zticks([0 .5*d_i d_i])
zticklabels({'0','5','10'})
% xlabel('Roller angle')
% ylabel('Angle of engagement')
% zlabel('Length value')

%zlab = sprintf('Length value (\x2113)','interpreter','latex');

xlabel('Roller angle ($\lambda$)','FontSize',16,'interpreter','latex');
ylabel('Angle of engagement ($\theta$)','FontSize',16,'interpreter','latex');
zlabel('Length value ($l$)','FontSize',16,'interpreter','latex');

hold on
xcx = zeros(1,2*n);
xcx(1:2*n) = n;
%Large x-axis limits
plot3(1:2*n, xcx, remap_mat1(n,:),'k--','LineWidth',.2)
plot3(1:2*n, xcx/2, remap_mat1(round(n/2),:),'k--','LineWidth',.2)
plot3(1:2*n, xcx/4, remap_mat1(round(n/4),:),'k--','LineWidth',.2)
plot3(1:2*n, (xcx - xcx/4), remap_mat1(n -round(n/4),:),'k--','LineWidth',.2)
plot3(1:2*n, (xcx - xcx), remap_mat1(1,:),'k--','LineWidth',.2)

%Large y-axis limits
ycy(1:n) = n;
plot3(ycy, 1:n , remap_mat1(:,n)','k--','LineWidth',.2)
plot3(ycy/2,1:n, remap_mat1(:,round(n/2))','k--','LineWidth',.2)
plot3(ycy/4,1:n, remap_mat1(:,round(n/4))','k--','LineWidth',.2)
plot3(ycy + ycy/2, 1:n, remap_mat1(:,round(n + n/2))','k--','LineWidth',.2)
plot3(2*ycy - ycy/4,1:n, remap_mat1(:,round(2*n - n/4))','k--','LineWidth',.2)

% figure
% for i = 1:21
%     plot_length(:,i) = length(:,i*20);
%     plot_angle(:,i) = angle(:,i*20);
% end
% 
% plot(plot_length)

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%   RESULT DISPLAY   %%%%%%%%%%%%%%%%%%%%%%%%

%Creating a result matrix for user use, interpolation is used to stretch or
%squeeze the vector to the desireable size
A(1,:) = interp1(1:numel(angle),angle,linspace(1,n,numel(Q_flow)));
A(2,:) = interp1(1:numel(inflating_angle),inflating_angle,linspace(1,n,numel(Q_flow)));
A(3,:) = interp1(1:numel(d1),d1,linspace(1,n,numel(Q_flow)));
A(4,:) =  Q_flow;
A(5,:) = interp1(1:numel(volume),volume,linspace(1,n,numel(Q_flow)));



%Display values obtained to user
elapsed_time = toc;

%Informative texts to indicate the user of the values obtained

txt = sprintf('Flow_solver complete at %.4f seconds',elapsed_time);
disp(txt)
% display = 'Recomended curve - %s, polynomial order %g, upper value = %.2f mL/s, lower value  %.2f mL/s \n';
% txt4 = sprintf(display, txt_recomend,poly,ymax,ymin);



end

