%% In dieser File werden ausschließlich konstante/bekannte Variablen definiert

%% Parameter Rakete - aus "ICBM_1_Version_Djordje.m"
%für drag:
P.reference_area = pi * 3^2 / 4;
P.drag_coefficient = 0.75; % beides aus den NASA file, sollte gute Näherung sein
% Daten von LGM - 25C Titan II
%%%%%%%%%%%%%%% Variablen durch Phipsis 'Test.m' ersetzt %%%%%%%%%%%%%%%%%
% P.burn_time_stage_1 = 156;  % in s
% P.burn_time_stage_2 = 180;
% P.I_sp_stage_1 = 258;       % LR 87
% P.I_sp_stage_2 = 316;       % LR 91
% P.phi_0 = 85;       % Startwinkel, Grad
% P.phi = P.phi_0 * pi/180;
% P.g = 9.81;                 % m/s^2
% P.rho_norm = 1.2250;        % kg / m^3
% P.C_d = 0.75;
% P.t_max = 1000;
% P.q_m_stage_1 = (P.m_stage_1_dry + P.m_stage_1_fuel)/P.burn_time_stage_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test.m
% Parameter Erde
P.m_earth = 5.972E24; % Erdmasse [kg]
P.r_earth_pole = 6356752; % Polradius [m]
P.r_earth_eq = 6378137; % Äquatorradius [m]
P.r_earth_m  =(P.r_earth_eq+P.r_earth_pole) / 2; % mittlerer Erdradius [m]
P.e_num = 0.08181919; % Numerische Exzentrizität der Erde [-] für das WGS84
P.w_earth = 7.292115e-5; % Rotationsgeschwindigkeit der Erde [rad/s]
P.G = 6.67430E-11; % Gravitationskonstante [m^3/(kg*s^2)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Rakete
% Werte der Titan II (https://de.wikipedia.org/wiki/Titan_(Interkontinentalrakete)
P.m_warhead = 4000; %Masse Sprengkopf [kg]

% Phase 1
P.m_stage1_dry = 4320; % [kg]
P.m_stage1_full = 116880; % [kg]
P.I_spec1 = 258; % Spezifischer Impuls [m/s]
P.burn_time1 = 156; % Burn time [s]
P.m_dot1 = (P.m_stage1_full-P.m_stage1_dry) / P.burn_time1; % Betrag Massestrom durch Verbrennung [kg/s]
P.thrust1 = 1893400; % Schub erste Stufe [N]

% Phase 2
P.m_stage2_dry = 2300; % [kg]
P.m_stage2_full = 26100; % [kg]
P.I_spec2 = 316; % Spezifischer Impuls [m/s]
P.burn_time2 = 180; % Burn time [s]
P.m_dot2 = (P.m_stage2_full-P.m_stage2_dry) / P.burn_time2; % Betrag Massestrom durch Verbrennung [kg/s]
P.thrust2 = 444819; % Schub zweite Stufe [N]

P.m_ges = P.m_stage1_full + P.m_stage2_full + P.m_warhead; %[kg]
%% Berechnung der Orthodrome zwischen 2 Punkten auf der Erdoberfläche, Genauigkeit 50m: 
%https://de.wikipedia.org/wiki/Orthodrome

P.f = 1/298.257223563; % Abplattung der Erde
P.F = (P.lat1+P.lat2)/2;
P.G = (P.lat1-P.lat2)/2;
P.l = (P.lon1-P.lon2)/2;

P.S = sind(P.G)^2*cosd(P.l)^2+cosd(P.F)^2*sind(P.l)^2;
P.C = cosd(P.G)^2*cosd(P.l)^2+sind(P.F)^2*sind(P.l)^2;
P.w = atan(sqrt(P.S/P.C));
P.D = 2*P.w*P.r_earth_eq;

P.T = sqrt(P.S*P.C)/P.w;
P.H1 = (3*P.T-1)/(2*P.C);
P.H2 = (3*P.T+1)/(2*P.S);

orthodrome = P.D*(1+P.f*P.H1*sind(P.F)^2*cosd(P.G)^2-P.f*P.H2*cosd(P.F)^2*sind(P.G)^2);

city1 = [P.lat1 ; P.lon1 ; P.h1]; % Abschussort
city2 = [P.lat2 ; P.lon2 ; P.h2]; % Zielort 

%Winkel zwischen zwei Punkten auf der Erdoberfläche, von Erdmittelpunkt aus
zeta = acosd(sind(P.lat1)*sind(P.lat2)+cosd(P.lat1)*cosd(P.lat2)*cosd(P.lon2-P.lon1));

%% Startrichtung 
[X_1, Y_1, Z_1] = geodetic2ecef(wgs84Ellipsoid,city1(1),city1(2),city1(3),"degrees");
[X_2, Y_2, Z_2] = geodetic2ecef(wgs84Ellipsoid,city2(1),city2(2),city2(3),"degrees");

Vektor1 = [X_1 Y_1 Z_1];
Vektor2 = [X_2 Y_2 Z_2];

R1 = norm(Vektor1); % Betrag des Vektors zum Startpunkt, entspricht Erdradius am Abschussort
R2 = norm(Vektor2); % Betrag des Vektors zum Zielpunkt, entspricht Erdradius am Zielort

r = (dot(Vektor1,Vektor1) - dot(Vektor1,Vektor2)) / dot(Vektor1, Vektor1); %Umrechnungsfaktor von Zielvektor auf Startvektor

alpha0 = 45; %Anfangswert für iterative Berechnung des Abschusswinkels
n = (norm(Vektor2)+norm(Vektor1)*(r-1)) / (r*cosd(alpha0)*norm(Vektor1)) - norm(Vektor2) / norm(Vektor1)*r + 1/r;
E = Vektor2 + r * (Vektor1 * n);

Start = E - Vektor1; % Richtung des Startvektors (Azimut, Elevation), tangential zur Erdoberfläche

lla = [P.lat1 P.lon1 P.h1];
%% Umrechnungen_ecef.m
%für Umrechnung_ecef:
P.lat = 0;
P.lon = 0;

% lat = 48.200714;        % Breitengrad
% lon = 16.363496;        % Laengengrad
P.height = 175;     % Hoehe von der Erdoberflaeche weg
%BA-Gebäude: 48.200714, 16.363496, 175 [m] Lösung stimmt überein mit https://www.oc.nps.edu/oc2902w/coord/llhxyz.htm

% Anfangswerte für den Abschuss;
% spaeter soll das Programm automatisch durch die Eingabe der Koordinaten
% vom Zielort, die Startwinkel ausgeben

P.startwinkel_azimuth = 0;        % von yNorth weg gemessen im ENU
P.startwinkel_elevation = 80;     % von der xEast und yNorth Ebene weg gemessen im ENU

P.omega_erde = 7.292115e-5;       % rad/s
% thrust im ecef:
P.g_0 = 9.81;           % m/s^2
P.I_sp_1 = 180;         % s, specific impulse 1 stage
P.m_dot = 515;          % kg/s , Massenstrom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%