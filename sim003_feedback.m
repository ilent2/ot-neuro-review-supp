% Simulate a trap with variable trap stiffness created using feedback
% and coutner propagating optical tweezers.
%
% Uses two weakly focussed beams for the count. prop. setup, opposite pol.
%
% Copyright Isaac Lenton (aka ilent2) 2020

addpath('../../Fedora/ott');

% Setup the beams
NA = 0.8;
beam1 = ott.BscPmGauss('NA', NA, 'index_medium', 1.33, ...
  'wavelength0', 1.0, 'polarisation', [1, 0]);
beam2 = ott.BscPmGauss('NA', NA, 'index_medium', 1.33, ...
  'wavelength0', 1.0, 'polarisation', [0, 1]);

% Setup the particle
particle = ott.TmatrixMie(0.4, 'index_medium', 1.33, ...
  'index_particle', 1.59, 'wavelength0', 1.0);

%% Calculate force along beam axis for different alphas
% Also calculate FZ(alpha, z) for axial interpolation

% Odd number of points so we have zero at the centre
alpha = linspace(0.3, 0.7, 41);
z = linspace(-1, 1, 81);

FZ = zeros(numel(alpha), numel(z));

for ii = 1:numel(alpha)

  beam = sqrt(alpha(ii))*beam1 + sqrt(1-alpha(ii))*beam2.rotateY(pi);
  
  fz = ott.forcetorque(beam, particle, 'position', [0;0;1].*z);
  FZ(ii, :) = fz(3, :);
  
% Check that the beam radial trap stiffness is mostly uniform
%   fz = ott.forcetorque(beam, particle, 'position', [1;0;0].*z + [0;0;1]);
%   FZ(ii, :) = fz(1, :);
  
end

figure();
imagesc(z, alpha, FZ);
xlabel('Z Position');
ylabel('\alpha');

% Look at the maximum force we can achieve for any position (cons. 0.07?)
figure();
plot(z, max(-z.*FZ./abs(z), [], 1));
xlabel('Z Position');

%% Calculate radial force
% This should be approximately uniform as long as the particle is
% near the equilibrium (z+/-1).  FX(x)

x = linspace(-0.3, 0.3, 50);
beam = sqrt(0.5)*beam1 + sqrt(0.5)*beam2.rotateY(pi);
fx = ott.forcetorque(beam, particle, 'position', [1;0;0].*x);
FX = fx(1, :);

%% Simulate particle moving around focus
% Simulate motion in just XZ, could also do XYZ
%
% We assume the beam ratio can be updated almost instantaniously and
% the driving range is between +/- 0.07 (from above graphs).
%
% We also assume the radial trap stiffness doesn't vary significantly
% along the beam axis.

stiff = 0.07;    % Moderate stiffness
% stiff = 0.02;  % Low stiffness

% Simulation time step of 1e-4 seems small enough for these two stiff.
dt = 1.0e-4;         % Simulation time step (seconds)
max_time = 2.5;

t = 0:dt:max_time;   % Array of simulation times

% Array of particle positions
xz = zeros(2, numel(t));
xz(:, 1) = [0; 0];

invGamma = 2e8;   % Inverse drag [SI]
kbT = 300.0 .* 1.3806e-23;    % 300 kelvin [SI]
nPc = 0.01.*1.33./3e8;  % 200 mW in water [SI]
uConv = 1e6;    % Convert from SI to simulation distance (sans. n_medium?)

sim_alphas = zeros(1, numel(t));   % Alphas needed for beams
sim_alphas(1) = 0.5;

tic

for ii = 2:numel(t)
  
  % Choose the force we want axially (linear potential)
%   target_fz = -xz(2, ii-1) .* stiff;

  % Double well potential
  if xz(2, ii-1) < -0.05
    target_fz = -(xz(2, ii-1)+0.1) .* stiff;
  elseif xz(2, ii-1) > 0.05
    target_fz = -(xz(2, ii-1)-0.1) .* stiff;
  else
    target_fz = xz(2, ii-1) .* stiff;
  end
  
  % Calculate alpha which produces this stiffness
  alpha_column = interp1(z.', FZ.', xz(2, ii-1)).';
  sim_alphas(ii) = interp1(alpha_column, alpha, target_fz);
  
  % Calculate force acting on particle
  fx = interp1(x, FX, xz(1, ii-1));
  f = [fx; target_fz] .* nPc;
  
  % Calcualte change in position (wth BM)
  dx = invGamma .* f + sqrt(2*kbT.*invGamma./dt).*randn(2, 1);
  
  % Update position
  xz(:, ii) = xz(:, ii-1) + dx .* dt .* uConv;
  
end

toc

%% Generate plot of trajectory and alphas

figure();
plot(t, xz, t, sim_alphas);
legend({'X Position', 'Z Position', 'Trap Ratio'});

figure();
hist(xz(2, :), 30);

%% Plot double well potential

dbz = linspace(-0.2, 0.2, 100);
target_fz = dbz .* stiff;
target_fz(dbz < -0.05) = -(dbz(dbz < -0.05)+0.1) .* stiff;
target_fz(dbz > 0.05) = -(dbz(dbz > 0.05)-0.1) .* stiff;

figure();
plot(dbz, target_fz);

%% Visualise the beams
% Generate a visualisation showing each beam with different colours.
% This isn't quite acurate (since both beams are coherent), but it makes
% it easier to visualise what is going on.
%
% This would be better but we can't see whats happening
% figure();
% beam = sqrt(our_alpha)*beam1 + sqrt(1-our_alpha)*beam2.rotateY(pi);
% beam.visualise('axis', 'y', 'range', [0.5, 0.5]);

rng = linspace(-0.5, 0.5, 60);

im1 = beam1.visualise('axis', 'y', 'range', {rng, rng});
im2 = beam2.rotateY(pi).visualise('axis', 'y', 'range', {rng, rng});

% our_alpha = min(sim_alphas);
our_alpha = 1;
imc_max = max(im1);

% Setup figure
figure('Position', [0, 0, 480, 480]);
coffset = 0.4;
cscale = 0.5 - coffset;
imc = cat(3, (our_alpha-coffset).*im1./cscale, 0*im1, (1-our_alpha+coffset).*im2./cscale)./imc_max.*2;
imh = image(rng, rng, imc);
axis image;
axh = gca();
set(axh, 'Position', [0, 0, 1, 1]);

hold on;
% Add a particle indicator
circtheta = linspace(0, 2*pi, 100);
circx = particle.radius.*cos(circtheta);
circy = particle.radius.*sin(circtheta);
circh = plot(circx, circy);
circh.LineWidth = 2;
circh.LineStyle = '--';
circh.Color = [0,1,0];

% Add a particle trajectory marker
trajh = plot(0, 0);
trajh.Color = [1,1,1,0.3];
hold off;

edges = linspace(-0.5, 0.5, 30);
ax2 = axes('position', [0, 0, 1, 0.2]);
histh = histogram(ax2, xz(2, :), edges);
histh.FaceColor = [1,1,1];
histh.FaceAlpha = 1;
ax2.Color = 'none';
ax2.Box = 'off';
ax2.YTick = [];

%% Generate video

trajh.LineWidth = 0.5;
circh.LineWidth = 2;

clear Frames
Frames = struct('cdata', {}, 'colormap', {});

for ii = 1:stride:numel(t)
  
  % Update beam colours (intensities)
  imc = cat(3, (sim_alphas(ii)-coffset).*im1./cscale, 0*im1, (1-sim_alphas(ii)-coffset).*im2./cscale)./imc_max;
%   imc = cat(3, sim_alphas(ii).*im1, 0*im1, (1-sim_alphas(ii)).*im2)./imc_max.*2;
  imh.CData = imc;
  
  % Update trajectory
  set(trajh, 'YData', xz(1, 1:ii), 'XData', xz(2, 1:ii));
  
  % Update particle position
  circh.YData = circx + xz(1, ii);
  circh.XData = circy + xz(2, ii);
  
  % Update histogram data
  histh.Data = xz(2, 1:ii);
  
  axis(axh, [-0.5, 0.5, -0.5, 0.5]);
  drawnow;
  
  Frames(end+1) = getframe(gcf());

end

% % v = VideoWriter('DoubleTrap', 'Archival');
% % v = VideoWriter('ModStiff', 'Archival');
% v = VideoWriter('LowStiff', 'Archival');
% open(v);
% writeVideo(v, Frames);
% close(v);

%% Create a colorbar for the figure

figure();
our_alphas = linspace(0, 1, 100);
colors = [our_alphas.', 0.*our_alphas.', flip(our_alphas).'].*2;
colormap(min(colors, ones(size(colors))));
imagesc();
h = colorbar()
% h.Location = 'SouthOutside';
h.Ticks = h.Limits;
h.TickLabels = {'Left', 'Right'};
h.FontSize = 12;
set(gcf(), 'Color', 'white');