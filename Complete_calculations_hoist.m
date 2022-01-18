clc,clear
%Ranya project 2019
%Insert load data:

max_load=50000 %kg spreader+container
lenght_cables=90 %meter
weight_cables=5.63 %kg/m for 6 strands 19 wire cable 36mm thick
mass_cable=weight_cables*lenght_cables
m=0.5*(max_load+mass_cable) %because 4 cables and every cable has a pulley



%insert motor data:
n_motor=4;
rpm_m=1500;
omega_m=rpm_m*0.104719755;
j_ratio=2;
P_extra=1.1; %take extra power to account for losses.


%drum calc
layers_on_drum=1 %less layers is higher service life
drum_diameter=0.750 %meter
r_D=0.375; %m. Assumed to be the same size as a train wheel.
drum_radius=drum_diameter/2
drum_circumference=drum_diameter*pi %meter
lenght_cable=90%meter
diameter_cable=0.028 %meter
windings_on_drum=lenght_cable/drum_circumference %number
width_of_rolled_rope=windings_on_drum*diameter_cable
width_of_drum=width_of_rolled_rope/layers_on_drum

v=2*2; %m/s  120 m/min but *2 since because of the pulley 2 times the distance of rope

j_L=m*(r_D)^2;
omega_D=v/r_D;

i=omega_m/omega_D

j_ref=j_L/i^2;
j_m=j_ref/j_ratio;

%Calculation of forces for torque
F=9.8*m;



T_L=F*r_D;
P_L=T_L*omega_D;
P_m=P_L*P_extra;
T_m=(T_L*omega_D)/(omega_m*1.1);

P_mn=P_m/n_motor;
T_mn=T_m/n_motor;
%Display values:
disp(['Per motor (of which there are ',num2str(n_motor), ') the required values are:'])
disp(['Torque ', num2str(T_mn), ' Nm'])
disp(['Speed ', num2str(rpm_m), ' rpm'])
disp(['Moment of inertia ', num2str(j_m), ' kg/m'])
disp(['Power ',num2str(P_mn/1000), ' kW'])

%%
%Approximate the gear ratio as best as possible:
clf

step=.25;
init_i1_vector=(i/8):step:8;
[i_result,z1,z2,z3,z4,i1,i2]=i_approx(i,init_i1_vector);

Ipos=find(min(abs(i-i_result)));
i_closest=i_result(Ipos)
z1_closest=z1(Ipos) %gear at the engine side
z2_closest=z2(Ipos)
z3_closest=z3(Ipos)
z4_closest=z4(Ipos)  %gear at drum side
i1_closest=i1(Ipos);
i2_closest=i2(Ipos);
i_horizontal=ones(1,length(i_result))*i;
x_vec=1:length(i_result);
figure(1), hold on
plot(x_vec,i_result,'b')
plot(x_vec,i_horizontal,'k')
xlabel('Iteration')
ylabel('Transmission Ratio')
legend({'Calculated ratio', 'Goal ratio'}), hold off

%Display found values neatly:
disp(' ')
disp(['There are two gear pairs. The closest total transmission ratio is ', num2str(i_closest), ' where the goal was ', num2str(i)])
disp(['Gears one and two will both have 21 teeth, so both Z1 and Z2 are 21.'])
disp(['Wheras Z2 will have ', num2str(z2_closest) ' teeth and Z4 will have ', num2str(z4_closest), ' teeth.'])

%%
% Recalculate the MoI ratio:

j_ref_actual=j_L/i_closest^2;
j_m_actual=6.16;
t_m_actual=1795;
j_ratio_actual=j_ref_actual/j_m_actual

disp(['The actual moment of inertia ratio will be ', num2str(j_ratio_actual), ' where the goal was ',num2str(j_ratio)])
%%
%Spur gear design:

%Defining module options
mod_pref=[0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.25,1.5,2.0,2.5,3,4,5,6,8,10,12,16,20,25,32,40,50]; % series of standardized modules for gears. In mm.
%Pitch diameter vectors:
d1=z1_closest.*mod_pref; %all in mm
d2=z2_closest.*mod_pref;
d3=d1;
d4=z4_closest.*mod_pref;


figure(2), hold on
plot(mod_pref, d1,'.', 'Markersize', 10)
plot(mod_pref, d2,'.', 'Markersize', 10)
plot(mod_pref, d4,'.', 'Markersize', 10)
xlabel('Standardized modules')
ylabel('Pitch gear diameter in mm')
legend({'Gear 1 and 3','Gear 2','Gear 4'})
hold off

%Gear width calculation
width_dia_ratio=1.3; %Assumed to be a normalized material, with gears placed in asymmetrical positions.
b1=d1.*width_dia_ratio;      %Gear width in mm.
b2=d2.*width_dia_ratio;
b3=d3.*width_dia_ratio;
b4=d4.*width_dia_ratio;
%% Using a polynomial to find Yfa and Ysa for profile shift =0:
z_polyvector=[17,18,19,20,25,30,40,60,80,100,150,200,800];
Yfa_polyvector=[2.97,2.91,2.85,2.80,2.62,2.52,2.40,2.28,2.22,2.18,2.14,2.12,2.063];
Ysa_polyvector=[1.52,1.53,1.54,1.55,1.59,1.63,1.67,1.73,1.77,1.79,1.83,1.87,1.966];

%5th degree polynomial was chosen because it seems to have the most
%accurate results.
polyfit_Yfa=polyfit(z_polyvector,Yfa_polyvector,5);
polyfit_Ysa=polyfit(z_polyvector,Ysa_polyvector,5);

%Using the found polynomial coefficients to find out the Yfa and Ysa values
%for the given number of teeth:
Yfa_z1=polyval(polyfit_Yfa,z1_closest);
Yfa_z2=polyval(polyfit_Yfa,z2_closest);
Yfa_z3=polyval(polyfit_Yfa,z3_closest);
Yfa_z4=polyval(polyfit_Yfa,z4_closest);

Ysa_z1=polyval(polyfit_Ysa,z1_closest);
Ysa_z2=polyval(polyfit_Ysa,z2_closest);
Ysa_z3=polyval(polyfit_Ysa,z3_closest);
Ysa_z4=polyval(polyfit_Ysa,z4_closest);

%%
% Proof of tooth strength and contact stress
%insert limits and other material properties given by external sources:
bend_stress_limit=106;      %slides, for normalized carbon steel
contact_stress_limit=269;   %slides, for normalized carbon steel
Emod=205;                   %CES Edupack, for Carbon steel, AISI 1118, normalized

%Spur gear force calculation
%calculating torque per gear
T4=T_L   %Torque on drum shaft
T3=T4/i2_closest %torque on shaft with only 2 gears
T2=T3 
T1=T4/i_closest %Torque on motor shaft

%Calculating tangential force per gear
Ft_1=(2*T1)./(d1/1000); %/1000 because d is in mm
Ft_2=(2*T2)./(d2/1000);
Ft_3=(2*T3)./(d3/1000);
Ft_4=(2*T4)./(d4/1000);

%Calculating bending stresses per gear:
sig_bend_z1=(Ft_1./(b1.*mod_pref)).*(Yfa_z1*Ysa_z1);
sig_bend_z2=(Ft_2./(b2.*mod_pref)).*(Yfa_z2*Ysa_z2);
sig_bend_z3=(Ft_3./(b3.*mod_pref)).*(Yfa_z3*Ysa_z3);
sig_bend_z4=(Ft_4./(b4.*mod_pref)).*(Yfa_z4*Ysa_z4);


%Plotting the bending stresses
bend_stress_limit=ones(1,length(mod_pref)).*bend_stress_limit;

figure(3),hold on
plot(mod_pref,bend_stress_limit,'k')
plot(mod_pref,sig_bend_z1)
plot(mod_pref,sig_bend_z2)
plot(mod_pref,sig_bend_z3)
plot(mod_pref,sig_bend_z4)
xlabel('Module')
ylabel('Calculated Bending stress')
legend({'Given limit','Gear 1','Gear 2', 'Gear 3','Gear 4'})
hold off

%Calculating contact stresses per gear:
Z_curve=2.5;        %for standard spur gears, lecture
Z_elasticity=sqrt(0.175*Emod);
k_factor=1;                %No clue what K is supposed to be, so for now it is 1.
u1_2=i1_closest;            %Teeth ratio is the same as the gear ratio calculated before.
u3_4=i2_closest;

sig_cont_z1=Z_curve*Z_elasticity.*sqrt(((k_factor*Ft_1)./(b1.*d1)).*((u1_2+1)/(u1_2)));  %Formula taken from lecture slides
sig_cont_z2=Z_curve*Z_elasticity.*sqrt(((k_factor*Ft_2)./(b2.*d2)).*((u1_2+1)/(u1_2)));
sig_cont_z3=Z_curve*Z_elasticity.*sqrt(((k_factor*Ft_3)./(b3.*d3)).*((u3_4+1)/(u3_4)));
sig_cont_z4=Z_curve*Z_elasticity.*sqrt(((k_factor*Ft_4)./(b4.*d4)).*((u3_4+1)/(u3_4)));

%Plotting the contact stresses:
contact_stress_limit=ones(1,length(mod_pref)).*contact_stress_limit;

figure(4),hold on
plot(mod_pref,contact_stress_limit,'k')
plot(mod_pref,sig_cont_z1)
plot(mod_pref,sig_cont_z2)
plot(mod_pref,sig_cont_z3)
plot(mod_pref,sig_cont_z4)
xlabel('Module')
ylabel('Calculated Contact stress')
legend({'Given limit','Gear 1','Gear 2', 'Gear 3','Gear 4'})
hold off

%%
%Output of the minimum module and gear dimensions that come with it:
%for now just fill in a module yourself. Assumption is made that all gears
%in the transmission use the same module.

m_target=6 %old target was 2 so calc below this needs to be adjusted 

i_m_target=find(mod_pref==m_target);
%Gear 1:
D1=d1(i_m_target)
B1=b1(i_m_target)

%Gear 2:
D2=d2(i_m_target)
B2=b2(i_m_target)

%Gear 3:
D3=d3(i_m_target)
B3=b3(i_m_target)

%Gear 4:
D4=d4(i_m_target)
B4=b4(i_m_target)

% tangential force with chosen module
Ftg_1=Ft_1(i_m_target)
Ftg_2=Ft_2(i_m_target)
Ftg_3=Ft_3(i_m_target)
Ftg_4=Ft_4(i_m_target)

%% shaft calc
%pressure angle = 20 degrees = 0.3490658504 rad
%for bending moment Fr and Fn need to be known
a=0.3490658504
Fr1=Ftg_1*tan(a)
Fr2=Ftg_2*tan(a)
Fr3=Ftg_3*tan(a)
Fr4=Ftg_4*tan(a)

Fn1=Ftg_1/cos(a)
Fn2=Ftg_2/cos(a)
Fn3=Ftg_3/cos(a)
Fn4=Ftg_4/cos(a)

% only one shaft calc needed, so middle shaft chosenFr2
width_bearing=20 %mm
room_between_gears=10 %mm
L_shaft=B2+B3+2*width_bearing+room_between_gears
L_shaft_man=190 %for manual calculation and finding FBD FBD will be in report to see.
T_shaft=T3
B_Bearing=40 % width of bearing in mm
Room=10  % room between components in mm
L1=B2/2+B_Bearing/2+Room %mm distance between center of bearing and first gear
L2=B2/2+B3/2+Room %mm distance between center of both gears
L3=B3/2+B_Bearing/2+Room %mm distance between center of second gear and bearing
L12=L1+L2
L23=L2+L3
L123=L1+L2+L3
% from calculations on paper we know that for the middle shaft
% bearing force
fb2=((L12*Fr3)+(L1*Fr2))/L123
fb1=Fr2+Fr3-fb2            
%calc for bending moments: together with FBDs on paper
Bend1_min=(fb1*0)/1000 %Nm
Bend1_max=(fb1*L1)/1000 %Nm
Bend2_min=((fb1*L1)-(Fr2*(L1-L1)))/1000 %Nm
Bend2_max=((fb1*L12)-(Fr2*(L12-L1)))/1000 %Nm
Bend3_min=((fb1*L12)-(Fr2*(L12-L1))-(Fr3*(L12-L12)))/1000
Bend3_max=((fb1*L123)-(Fr2*(L123-L1))-(Fr3*(L123-L12)))/1000

%calc for shear forces
S1=fb1
S2=fb1-Fr2
S3=fb1-Fr2-Fr3
S4=fb1-Fr2-Fr3+fb2

Shear_line=[S1 S1 S2 S2 S3 S3 S4]

% calc for torsion line

Torsion_line=[0 0 T3 T3 0 0 0]
figure(5), hold on 
Bend_line=[Bend1_min Bend1_max Bend2_min Bend2_max Bend3_min Bend3_max Bend3_max]
Location=[0 L1 L1 L2+L1 L2+L1 L3+L2+L1 L3+L2+L1]
title('Shear Force')
plot(Location, Shear_line)
xlabel('Location (mm)')
ylabel('Shear Force (N)')

figure(6), hold on 
title('Bending Moment')
plot(Location, Bend_line)
xlabel('Location (mm)')
ylabel('Bending moment (N/m)')

figure(7), hold on 
title('Torsion Line')
plot(Location, Torsion_line)
xlabel('Location (mm)')
ylabel('Torque N/m)')

max_bend=(max(Bend_line))/1000 %Nmm to nm
max_tors=(max(Torsion_line)) %Nm

%minimal shaft diameter calc
% for Stainless steel, martensitic, AISI 410, intermediate temper. 
N=2.5 %design factor
Kt=2.5 %stress concentration factor for sharp inner corners
Ys= 850*10^6% Mpa to pa
Ys_1=425*10^6 %Mpa to pa

D_shaft=(((32*N)/pi)*sqrt(((Kt*max_bend)/Ys_1)^2+0.75*(max_tors/Ys)^2))^(1/3) %diameter in meter
 
%% bearing calculation

L10=10000 %hours
nShaft=(omega_D/0.104719755)*i2_closest %rpm shaft
L10m=(L10*60*nShaft)/(10^6) %revolutions

%bearing 1
X=0.56 %from table 11.1
R1=fb1
P1=X*R1
e=3 %for rolling contact bearings
C1=P1*((L10m)^(1/e))
% no suitable ball bearing for this load so cylindrical bearing needed unless
% making shaft far to thick
%SKF 319426 B-2LS with 130 mm diamter in stead om 100 mm


%bearing 2
X=0.56 %from table 11.1
R2=fb2
P2=X*R2
e=3 %for rolling contact bearings
C2=P2*((L10m)^(1/e))
% no suitable ball bearing for this load so cylindrical bearing needed unless
% making shaft far to thick
%SKF 319426 B-2LS with 130 mm diamter in stead om 100 mm

%% key calculation

%only one key calculation needed so chosen for gear 2 on the calculated
%shaft, since the shaft diameter needs to be increased to 160mm to
%accomodate for the right bearing diameter of the shaft can be taken to be 160mm
%but no values for this so for now we use that the shaft has a diameter of 80mm 
% we hope to revise this later when we know where our calculations went
% wrong

%table 4.11 from the book provides:
b_key=22 %mm
h_key=14 %mm
L_key=110 %mm one size bigger than the width of gear 2 = B2
l_key=L_key-b_key %for rounded key
allowable_bearing_stress= 150 %MPa
crushing_stress_real=((4*T3)*10^3)/(h_key*l_key*(D_shaft*1000))
 %much bigger key needed, no information given about this so own design
 %or using double keys
 
b_key=22 %mm %increase width of key by 4mm, should be possible since using big shaft
h_key=14 %mm %increase height of key by 2mm, should be possible since using big shaft
L_key=110 %mm one size bigger than the width of gear 2 = B2
l_key=L_key %for square key in stead of rounded key
allowable_bearing_stress= 150 %MPa
crushing_stress_real=(((4*T3)*10^3)/(h_key*l_key*(D_shaft*1000)))/2  %dived by 2 because of using 2 keys 180 degrees apart
%so to make it work we use 2 keys of the above sizes.
