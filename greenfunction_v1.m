% Solve total beam response using green's function for beam
% Euler-bernoulli beam element

% Single element dynamic model of steel beem
b = .03; % m, cross section base
h = .03; % m, cross section height
L = .15; % m, overall length

% Identify Modes
E = 200 * 10^9;  % GPa, Young's modulus of steel
p = 8000; % Kg/m^3, density of steel
I = 1/12 * b * h^3; % m^4, moment of inertia from rectangular face
A = b*h; % m^2, area of cross section

EI = E*I; % Fluxeral rigidity (Youngs Modulus * Cross section Inertia
MU = p*A; % Mass per unit length

ra = 0; % External damping
ri = 0; % Internal damping

FHZ = 0.9; % Hz Temporal frequncy of load
Om = FHZ * 2 * pi; % Temporal frequency (rad/s) of load

dt_plot = .01; % s, how often to plot time soln
dx_plot = .01;% m, how often to plot space soln

distplot = [0:dx_plot:L]; % Points of beam to plot
timeplot = [0:dt_plot:10]; % Points in time to plot

% Define beam parameter k
k4 = (MU*Om^2 - 1i*ra*Om) / (EI + 1i*ri*Om);
k3 = k4 ^ (3/4);
k2 = k4 ^ (1/2);
k  = k4 ^ (1/4);

syms x phi1 phi2 phi3 phi4 T fone ftwo ftre ffor Tlow mff chi

% Define basis functions
phi1 =  0.5 * (cosh(k*x) + cos(k*x));
phi2 =  0.5 * (sinh(k*x) + sin(k*x));
phi3 =  0.5 * (cosh(k*x) - cos(k*x));
phi4 =  0.5 * (sinh(k*x) - sin(k*x));
U = 1.0*(x>0);
%plot(phi1(0:.01:1))

% Calculate all IC's (from fixed-free conditions)
% Complete matrix is given, take subset that matters
T = [subs(phi1,L)  , subs(phi2,L)/k , subs(phi3,L)/k2, subs(phi4,L)/k3;
    k*subs(phi4,L) , subs(phi1,L)   , subs(phi2,L)/k , subs(phi3,L)/k2;
    k2*subs(phi3,L), k*subs(phi4,L) , subs(phi1,L)   , subs(phi2,L)/k ;
    k3*subs(phi2,L), k2*subs(phi3,L), k*subs(phi4,L) , subs(phi1,L)];
% Note that unknown IC's will be expressed in terms of a function 
% in the frequency domain!
fone =  subs(phi4,L-chi)/(k3*(EI+1i*ri*Om));
ftwo =  subs(phi3,L-chi)/(k2*(EI+1i*ri*Om));
ftre =  subs(phi2,L-chi)/(k *(EI+1i*ri*Om));
ffor =  subs(phi1,L-chi)/(1 *(EI+1i*ri*Om));
% Ignore rows corresponding to unkown BC's at x=L
% Take rows corresponding to known BC's at x=L
% for us, fixed-free
Tlow = T(3:4,3:4);
mff = -[ftre;ffor];

% Finally get IC's (changes based on boundary conditions)
syms x0vec xdd0 xddd0

x0vec = inv(Tlow) * mff;
xdd0 = x0vec(1); % initial condition captured as function of chi
xddd0 = x0vec(2);% initial condition captured as function of chi

x0 = 0;  % Fixed end
xd0 = 0; % Fixed end

%G =   subs(phi4,x-chi)*subs(U,x-chi) / (k3*(EI+1i*ri*Om)) ...
%    + x0*phi1        + (xd0/k)*phi2...
%    + (xdd0/k2)*phi3 + (xddd0/k3)*phi4;
Gxgtchi =   subs(phi4,x-chi)/ (k3*(EI+1i*ri*Om)) ...
          + x0*phi1        + (xd0/k)*phi2...
          + (xdd0/k2)*phi3 + (xddd0/k3)*phi4;
Gxltchi =   x0*phi1        + (xd0/k)*phi2...
          + (xdd0/k2)*phi3 + (xddd0/k3)*phi4;


% If our forcing function is a delta at the end, then
% integration from 0 to L of f(chi)G(x,chi) d chi
% will just be the value of G at chi=L times the magnitude of f

force_mag = 1000000;
SpatialResponseMode = force_mag * subs(Gxltchi,chi,L);
SpatialResponsePlot = double(subs(SpatialResponseMode,distplot));
TimeComponent = exp(1i*Om*timeplot);
TotalResponsePlot = real(TimeComponent.' * SpatialResponsePlot);

% Plot displacement
exaggeration = 1;
figure();
surf(distplot, timeplot, TotalResponsePlot)
xlim([distplot(1), distplot(end)]);
ylim([timeplot(1), timeplot(end)]);
zlim([-distplot(end)/exaggeration, distplot(end)/exaggeration]);
ylabel('Time')
xlabel('Space')

%ResponseIndef = int(G, chi);
%Response = subs(ResponseIndef,chi,L) - subs(ResponseIndef, chi,0);










