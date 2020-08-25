% Simulation illustrating effect of scanning an OT
%
% The trap is a tightly focussed Gaussian beam (OTT).  For this simulation
% the particle is constrained to move only in the XY plane.  The particle
% is moving at terminal velocity (i.e., v = F/gamma) and there is random
% noise simulating the effect of Brownian motion.
%
% Copyright Isaac Lenton (aka ilent2) 2020

addpath('../../Fedora/ott');

% Setup the beam (CP Gaussian)
beam = ott.BscPmGauss('NA', 1.2, 'index_medium', 1.33, ...
  'wavelength0', 1.0);

% Setup the particle
particle = ott.TmatrixMie(0.4, 'index_medium', 1.33, ...
  'index_particle', 1.59, 'wavelength0', 1.0);

% Calculate the forces in the XY plane
x = linspace(-2, 2, 80);
y = linspace(-2, 2, 80);
[xx, yy] = meshgrid(x, y);
xyz = [xx(:), yy(:), 0*xx(:)].';
fxyz = ott.forcetorque(beam, particle, 'position', xyz);

% Fields for interpolation (faster)
FX = reshape(fxyz(1, :), size(xx));
FY = reshape(fxyz(2, :), size(xx));

%% Generate visualisation of forces and beam

figure();
beam.visualise('range', [2, 2]);
colormap hot;
hold on;
qh = quiver(xyz(1, :), xyz(2, :), fxyz(1, :), fxyz(2, :));
qh.Color = 'white';
hold off;

%% Simulate the particle's trajectory
% Takes about 25 seconds

scan_radius = 1.0;   % In wavelengths
dtheta = 2*pi/2e4;   % Scane step size (radians per time step)
dt = 1.0e-5;         % Simulation time step (seconds)
max_time = 0.5;

t = 0:dt:max_time;   % Array of simulation times

% Array of particle positions
xy = zeros(2, numel(t));
% xy(:, 1) = [scan_radius; 0];   % Start initially in trap
xy(:, 1) = [0; 1.5];

% Calculate beam positions
bxy = scan_radius .* [cos(dtheta.*(0:numel(t))); sin(dtheta.*(0:numel(t)))];

invGamma = 2e8;   % Inverse drag [SI]
kbT = 300.0 .* 1.3806e-23;    % 300 kelvin [SI]
nPc = 0.01.*1.33./3e8;  % 200 mW in water [SI]
uConv = 1e6;    % Convert from SI to simulation distance (sans. n_medium?)

tic

for ii = 2:numel(t)
  
  % Calculate force acting on particle
  fx = interp2(xx, yy, FX, xy(1, ii-1)-bxy(1, ii-1), xy(2, ii-1)-bxy(2, ii-1));
  fy = interp2(xx, yy, FY, xy(1, ii-1)-bxy(1, ii-1), xy(2, ii-1)-bxy(2, ii-1));
  f = [fx; fy] .* nPc;
  
  % Calcualte change in position (wth BM)
  dx = invGamma .* f + sqrt(2*kbT.*invGamma./dt).*randn(2, 1);
  
  % Update position
  xy(:, ii) = xy(:, ii-1) + dx .* dt .* uConv;
  
end

toc

%% Generate plot of positions

figure();
plot(t, bxy(:, 1:end-1), t, xy)

%% Generate a visualisation

% Generate the beam data for visualisatino
imxrange = linspace(-2, 2, 80);
imyrange = linspace(-2, 2, 80);
im = beam.visualise('range', {imxrange, imyrange});

stride = 50;    % Output frame stride

% Setup the figure
figure('position', [0, 0, 480, 480]);
imh = imagesc(imxrange, imyrange, im);
colormap hot;
axh = gca();
set(axh, 'Color', [0,0,0]);
set(axh, 'Position', [0, 0, 1, 1]);
axis image;

hold on;
% Add a particle indicator
circtheta = linspace(0, 2*pi, 100);
circx = particle.radius.*cos(circtheta);
circy = particle.radius.*sin(circtheta);
circh = plot(circx, circy);
circh.LineWidth = 2;

% Add a particle trajectory marker
trajh = plot(0, 0);
trajh.Color = 'white';
hold off;

% Add text annotation for time
timean = annotation(gcf(),'textbox',...
  [0.602666666666666 0.895 0.4 0.075],'Color',[1 1 1],...
  'String',{'t = 0 ms'},...
  'FontWeight','bold',...
  'FontSize',18,...
  'FitBoxToText','off');

%% Generate video

trajh.LineWidth = 0.5;
circh.LineWidth = 2;

clear Frames
Frames = struct('cdata', {}, 'colormap', {});

for ii = 1:stride:numel(t)

  % Update beam position
  imh.XData = imxrange + bxy(1, ii);
  imh.YData = imyrange + bxy(2, ii);
  
  % Update trajectory
  set(trajh, 'XData', xy(1, 1:ii), 'YData', xy(2, 1:ii));
  
  % Update particle position
  circh.XData = circx + xy(1, ii);
  circh.YData = circy + xy(2, ii);
  timean.String = {['t = ' num2str(round(t(ii)*1000)) ' ms']};
  
  axis([-2, 2, -2, 2]);
  drawnow;
  
  Frames(end+1) = getframe();

end

% v = VideoWriter('ScannedBeam', 'Archival');
% open(v);
% writeVideo(v, Frames);
% close(v);

%% Select three nice stills to show in the paper

ii = 800;   % Particle initially untrapped
% ii = 1500;    % Particle pulled into trap
% ii = numel(t); % Particle orbiting with beam
% ii = 10000; % Particle orbiting with beam

circh.LineWidth = 4;
trajh.LineWidth = 1;
timean.String = {['t = ' num2str(round(t(ii)*1000)) ' ms']};

% Update beam position
imh.XData = imxrange + bxy(1, ii);
imh.YData = imyrange + bxy(2, ii);

% Update trajectory
set(trajh, 'XData', xy(1, 1:ii), 'YData', xy(2, 1:ii));

% Update particle position
circh.XData = circx + xy(1, ii);
circh.YData = circy + xy(2, ii);

axis([-2, 2, -2, 2]);
drawnow;
