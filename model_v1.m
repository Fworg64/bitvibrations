% Single element dynamic model of steel beem
b = .03; % m, cross section base
h = .03; % m, cross section height
l = .15; % m, overall length

% Identify Modes
E = 200;  % GPa, Young's modulus of steel
p = 8000; % Kg/m^3, density of steel
I = 1/12 * b * h^3; % m^4, moment of inertia from rectangular face
A = b*h; % m^2, area of cross section

% Bn from fixed/free beam
% Zeros to cos(Bl) * cosh(Bl) = -1
% weighted natural freqs for fixed/free beam
Bl = linspace(0,15,15000);
determ = cos(Bl).*cosh(Bl) +1;
zeros_indices_set = abs(diff(sign(determ))) > 0.5;
zeros = Bl(zeros_indices_set);
n = (length(zeros)+1):100;
Bnl = [zeros, (2*n-1)*pi / 2]
Wn = (Bnl/l).^2 .*sqrt(E*I/(p*A))
freqs_hz = Wn/(2*pi)
