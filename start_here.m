%% THIS is the MAIN-file %% DO NOT EDIT HERE!
clear all, clc

%% Variables
run variables.m

%% Coordinates

% Koordinaten Start
P.lat1 = 38.9613433; % Pjönjang, Nordkorea
P.lon1 = 125.8279959;
P.h1 = 29; % Abschusshöhe, zur Seehöhe gemessen [m]

%Koordinaten Ziel
P.lat2 = 48.203530; % St. Pölten, Österreich
P.lon2 = 15.638170;
P.h2 = 267; % Zielhöhe, zur Seehöhe gemessen [m]
%% state space vector / initial state
tspan = [0, 10000]; % row-vector; when to start and end Integration
r0    = [0; 0; 0]; % [x(0); y(0); z(0)];
v0    = [0; 0; 0]; % [xdot(0); ydot(0); zdot(0)]
s0    = [r0; v0]; % initial state-vector
%%% WICHTIG: WIE sieht STATE-SPACE-VECTOR aus?!
% s = [x; xdot;] aber mit x= [x(0);y(0);z(0)] & xdot =[xdot(0);ydot(0);zdot(0)]
% https://en.wikipedia.org/wiki/State_vector_(navigation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ODE Solver DJO

%options = odeset('MaxStep',1,'Events',@(t,s) bodenkontakt(s,P)); % hier kann man noch Stepgröße einstellen etc. 
%[t,s] = ode45(@(t,s) flightODE(t,s,P), tspan, s0, options); %die function ODE soll von state variable s zu sdot gehen

% s & s0 sind state Vektoren
% @ ruft func 'flightODE' auf; (t,s) sind immer noch Variablen; 
% der flightODE werden (t,s,P) übergeben (geht das mit struct?)
% sol= [t_out,s_out]; t_out... alle Zeitpunkte; s_out...state variables bei Zeitpunkten 
% (s_out->[x@t1,xdot@t1; x@t2, xdot@t2;...]; Spaltenanzahl=Statevar.anzahl)

fH_r = @(t,x) Bewegungsgleichung(t,x,P);
opt_contact = odeset('MaxStep', 0.1,"Events", @(t,x) event_bodenkontakt(t,x,P));

sol = ode45(fH_r, [0, P.t_max], x_0 ,opt_contact)

%% ODE Solver Phillip
tmax = 3000; % maximale Flugzeit
timespan = 0:tmax;
init_cond = [0, 0, 0, 0, 0, 0]; % Anfangsbedingungen im ENU
options = odeset('Events',@event_contact);
[t,s] = ode45(@SPS,timespan,init_cond,options);

Location_ENU = lla2enu([s(1),s(3),s(5)],[city1(1),city1(2),city1(3)]); %[xEast,yNorth,zUp] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots

hold on

xlabel('Flugweite (km)')
ylabel('Flughoehe (km)')
comet(sol.y(1,:), sol.y(2,:))
axis equal
xlim([0 sol.y(1,end)])
ylim([0 sol.y(2,end)])

hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function - SPS - MIKE
function sdot = flightODE(t,s,P) %in der Fkt wird beschrieben wie man von s zu sdot kommt!
% variablen: (kann man auch vor func definieren und als input einfügen
m = massconsumption(t,P); %Schwierigkeit: ist nun in einer externen Funktion "massconsumption.m" ; aufrufen?!
P.F_G = m*g_ECEF; % muss VEKTOR sein! -> von DJO 'Umrechnungen_ecef.m' hab ichs als function geschrieben
[g_ecef_x, g_ecef_y, g_ecef_z] = g_ecef;
P.thrust = thrust(t,P); 
v = [s(4);s(5);s(6)]; % Geschwindigkeit wird über state-vektor beschrieben
% v =10
P.drag = drag_rocket(v,P); % alle Kräfte auch als Funktion!
    % P.height wirft noch Fehler auf... wie implementieren? 
P.Fcor = 0;

sdot= zeros(6,1); % allgemein definieren als Spaltenvektor
sdot(1) = s(4); % xdot
sdot(2) = s(5); % ydot
sdot(3) = s(6); % zdot

% Bewegungsgleichungen ref. Newton 2nd law (SPS): -> a(xddot;yddot;zddot) = F(t)/m(t)

% xddot:
sdot(4) = P.FG(1) + P.thrust(1) + P.drag(1) + P.Fcor(1) + P.Fzentri(1); 
% das 'durch m(t)' fehlt noch

% yddot:
sdot(5) = P.FG(2) + P.thrust(2) + P.drag(2) + P.Fcor(2) + P.Fzentri(2);

% zddot:
sdot(6) = P.FG(3) + P.thrust(3) + P.drag(3) + P.Fcor(3) + P.Fzentri(3);

end

%% Functions - forces
%m(t)
function Mass_total = massconsumption(t, P)
    for n = 1:length(t)
        if t <= P.burn_time1
            Mass_total = P.m_stage1_dry + (P.m_stage2_fuel + P.m_stage2_dry) + P.m_sprengkopf + (P.m_stage1_fuel-P.m_stage1_fuel/P.burn_time1)*t;
        elseif t <= P.burn_time2 && t > P.burn_time1
            Mass_total = P.m_stage1_dry + P.m_sprengkopf + P.m_stage2_dry + (P.m_stage2_fuel-P.m_stage2_fuel/P.burn_time2)*t;
        else
            Mass_total = P.m_stage1_dry + P.m_sprengkopf + P.m_stage2_dry;
        end
    end
end

%thrust im ECEF:
function [thrust_ecef_x, thrust_ecef_y, thrust_ecef_z] = thrust(P)

    thrust = P.g_0 * P.I_sp_1 * P.m_dot;
    
    [thrust_x, thrust_y, thrust_z] = aer2enu(P.startwinkel_azimuth,P.startwinkel_elevation,thrust);
    
    thrust_enu = [thrust_x, thrust_y, thrust_z];
    
    [thrust_ecef_x, thrust_ecef_y, thrust_ecef_z] = enu2ecefv(thrust_x,thrust_y,thrust_z,P.lat,P.lon);
end

%gravity im ECEF:
function [g_ecef_x, g_ecef_y, g_ecef_z] = gravity_force(P)

    F_G = gravitywgs84(P.height,P.lat,'None');   %  Rechnet Erdbeschleunigung aus;
                                                 % 'None' Befehl stellt die eingebaute Warnung ab 20000m Hoehe ab
    
    [g_ecef_x, g_ecef_y, g_ecef_z] = enu2ecefv(0,0,F_G,P.lat,P.lon);

end

%drag im ECEF:
function [drag_ecef_x, drag_ecef_y, drag_ecef_z] = drag_force(velocity,P)

    if P.height > 25000 % upper stratosphere
        T     = -131.21 + 0.00299*P.height;
        p     = 2.488*((T+273.1)/216.6)^-(11.388);
        P.rho = p/(0.2869*(T+273.1));
    elseif P.height < 11000 %troposhere
        T     = 15.04-0.00649*P.height;
        p     = 101.29*((T+273.1)/288.08)^5.256;
        P.rho = p/(0.2869*(T+273.1));
    else  %lower stratosphere
        T     = -56.46;
        p     = 22.65*exp(1.73-0.000157*P.height);
        P.rho = p/(0.2869*(T+273.1));
    end
    
%     P.velocity = 10;        % hier brauchen wir noch eine Funktion die den Wert ...
%                             % nach der Zeit ändert
                            
    P.reference_area = pi * 3^2 / 4;
    P.drag_coefficient = 0.75;      % aus den NASA files, sollte gute Näherung sein
    
    drag = 0.5 * P.rho * velocity^2 * P.reference_area * P.drag_coefficient;
    
    % B = [cosd(startwinkel_elevation)*cosd(startwinkel_azimuth);...
    %     cosd(startwinkel_elevation)*sind(startwinkel_azimuth);...
    %     sind(startwinkel_elevation)];
    % 
    % P.drag_enu1 = drag * B;
    % P.drag_enu  = transpose(P.drag_enu1);
    
    % drag im ENU dargestellt, zeigt in die gleiche Richtung wie thrust, in der
    % Bewegungsgleichung einfach negativ nehmen
    [drag_x, drag_y, drag_z] = aer2enu(P.startwinkel_azimuth,P.startwinkel_elevation,drag);
    
    [drag_ecef_x, drag_ecef_y, drag_ecef_z] = enu2ecefv(drag_x,drag_y,drag_z,P.lat,P.lon);
%     y = enu2ecefv(drag_x,drag_y,drag_z,P.lat,P.lon);

end

%Coriolis im ECEF:
function y = corioliskraft(vel,P)       % hier muss man die Geschwindigkeit in x oder y eingfügen

    y = 2 * P.omega_erde * vel;
end

%Zentrifugalkraft im ECEF:
function y = zentrifugalkraft(abstand,P)        % Abstand ist in x oder in y Richtung

    y = P.omega_erde.^2 * abstand;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Bewegungsgleichung - Phillip
function f = SPS(P,rocket_mass,s)

    for t = 1:3000
       
        if t <= P.burn_time1
            f = [(P.m_dot1/rocket_mass)*P.I_spec1+(2*P.w_earth*s(4)/rocket_mass)+(P.w_earth^2*s(1)/rocket_mass),...
                (P.m_dot1/rocket_mass)*P.I_spec1+(2*P.w_earth*s(2)/rocket_mass)+(P.w_earth^2*s(1)/rocket_mass),...
                (P.m_dot1/rocket_mass)*P.I_spec1+gravity+(2*P.w_earth*s(4)/rocket_mass)+(P.w_earth^2*s(1)/rocket_mass),...
                ];
    
        elseif t > P.burn_time1 & t <= P.burn_time1+P.burn_time2
            f = [(P.m_dot2/rocket_mass)*P.I_spec2+(2*P.w_earth*s(4)/rocket_mass)+(P.w_earth^2*s(1)/rocket_mass),...
                (P.m_dot2/rocket_mass)*P.I_spec2+(2*P.w_earth*s(2)/rocket_mass)+(P.w_earth^2*s(1)/rocket_mass),...
                (P.m_dot2/rocket_mass)*P.I_spec2+gravity+(2*P.w_earth*s(4)/rocket_mass)+(P.w_earth^2*s(1)/rocket_mass),...
                ];
        else
            f = [(2*P.w_earth*ydot/rocket_mass)+(P.w_earth^2*s(1)/rocket_mass),...
                (2*P.w_earth*xdot/rocket_mass)+(P.w_earth^2*s(1)/rocket_mass),...
                gravity+(2*P.w_earth*s(4)/rocket_mass)+(P.w_earth^2*s(1)/rocket_mass),...
                ];
        end
    
    end

end
%% Function Rocket Mass - Phillip
function m = rocket_mass(t,P)

    if t <= P.burn_time1
        m = P.m_ges - P.m_dot1 * t;

    elseif t > P.burn_time1 & t <= P.burn_time1+P.burn_time2

        m = P.m_ges - P.m_stage1_full - P.m_dot2 * (t-P.burn_time1);

    else 
        m = P.m_warhead; 
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Events - MIKE
% damit man keine Zeitspanne vorgeben muss, sondern Simulation bei einem Event abbricht

function [value,isterminal,direction] = bodenkontakt(s,P)

P.earthradius_equator = 6378137; % [m] Erdradius am Äquator -> vl Zielstädte wie: "Quito" (20km südlich vom Äquator)
value      = s(3)-P.earthradius_equator; % 3.Zeile des State Vektors -> z soll 0 werden
isterminal = 1; % zwei Optionen: 0... flag; 1... simulation stops
direction  = 0; % drei Optionen: -1...negativ slope ; 0... dont care ; 1... positive slope
end


%% Events - Phillip
function f = event_contact(Location_ENU)

f = [value,isterminal,direction];

s = [Location_ENU(1) Location_ENU(2) Location_ENU(3)];
value = R2 - norm(s(1),s(2),s(3));
isterminal = 1;
direction = 0;

end