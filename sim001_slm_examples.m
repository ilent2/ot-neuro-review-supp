% Figure illustrating different SLM patterns
%
% Uses OTSLM and OTT for simulating fields focussed by a high numerical
% aperture objective.
%
% Copyright Isaac Lenton (aka ilent2) 2020

addpath('../../Fedora/otslm');
addpath('../../Fedora/ott');

slm_sz = [256, 256];
NA = 1.0;
gauss = otslm.simple.gaussian(slm_sz, 50);

rng = [5, 5];

figure();
imagesc(gauss);
axis image;

%% Single beam, uniform phase

im = zeros(slm_sz);
Efield = gauss .* exp(1i*2*pi*im);
beam = ott.BscPmParaxial(NA, Efield, 'Nmax', 25);

figure();
imagesc(mod(im, 1));
caxis([-pi, pi]);
axis image;

figure();
beam.visualise('axis', 'y', 'range', rng);
colormap hot;

%% Radial Displaced beam

im = otslm.simple.linear(slm_sz, 50);
Efield = gauss .* exp(1i*2*pi*im);
beam = ott.BscPmParaxial(NA, Efield, 'Nmax', 25);

figure();
imagesc(mod(im, 1));
axis image;

figure();
beam.visualise('axis', 'y', 'range', rng);
colormap hot;

%% Axial Displaced beam

im = otslm.simple.spherical(slm_sz, 128)*2;
Efield = gauss .* exp(1i*2*pi*im);
beam = ott.BscPmParaxial(NA, Efield, 'Nmax', 25);

figure();
imagesc(mod(im, 1));
axis image;

figure();
beam.visualise('axis', 'y', 'range', rng);
colormap hot;

%% Axial and Radial Displaced beam

im = otslm.simple.spherical(slm_sz, 128)*2 + otslm.simple.linear(slm_sz, 50);
Efield = gauss .* exp(1i*2*pi*im);
beam = ott.BscPmParaxial(NA, Efield, 'Nmax', 25);

figure();
imagesc(mod(im, 1));
axis image;

figure();
beam.visualise('axis', 'y', 'range', rng);
colormap hot;

%% Multiple traps

im1 = otslm.simple.spherical(slm_sz, 128)*1.8 + otslm.simple.linear(slm_sz, 40);
im2 = otslm.simple.spherical(slm_sz, 128)*2 - otslm.simple.linear(slm_sz, 90);
im3 = otslm.simple.spherical(slm_sz, 128)*-2 + otslm.simple.linear(slm_sz, 90);
im = otslm.tools.combine({im1, im2, im3}, 'method', 'rsuper');
Efield = gauss .* exp(1i*2*pi*im);
beam = ott.BscPmParaxial(NA, Efield, 'Nmax', 25);

figure();
imagesc(mod(im, 1));
axis image;

figure();
beam.visualise('axis', 'y', 'range', rng);
colormap hot;

%% Annular beam

r1 = 70;

im = otslm.simple.checkerboard(slm_sz) .* ~otslm.simple.aperture(slm_sz, [r1, 110], 'shape', 'ring');
% Efield = gauss .* exp(1i*2*pi*im) .* otslm.simple.aperture(slm_sz, [r1, 110], 'shape', 'ring');
Efield = exp(1i*2*pi*im) .* otslm.simple.aperture(slm_sz, [r1, 110], 'shape', 'ring');
beam = ott.BscPmParaxial(NA, Efield, 'Nmax', 30);

figure();
imagesc(mod(im, 1));
axis image;

figure();
beam.visualise('axis', 'y', 'range', 1.5*rng);
colormap hot;

%% Single beam, LG phase

im = otslm.simple.lgmode(slm_sz, 3, 0);
Efield = gauss .* exp(1i*2*pi*im);
beam = ott.BscPmParaxial(NA, Efield, 'Nmax', 25);

figure();
imagesc(mod(im, 1));
axis image;

figure();
beam.visualise('axis', 'z', 'range', rng);
colormap hot;

%% TOW beam

lin = otslm.simple.linear([256, 128], 30);
prism = [fliplr(lin), lin];
lens = otslm.simple.spherical([256, 256], 512, 'type', '1d', 'angle_deg', 90)*20;
Eff = gauss .* exp(1i*2*pi.*(0.6.*prism -  lens));

figure();
imagesc(mod(prism +  lens, 1));
axis image;

beam = ott.BscPmParaxial(0.95, cat(3, Eff, zeros(size(Eff))), 'Nmax', 40);

figure();
beam.visualise('axis', 'y', 'range', 2*rng);
colormap hot;

%% Chiral optical field

lmodes = 1:3:10;
im = otslm.simple.chiralAsvpp(slm_sz, lmodes, 'eta', 0.01, 'radius', 128, ...
    'background', 'checkerboard');
Efield = gauss .* exp(1i*2*pi*im);
beam = ott.BscPmParaxial(NA, Efield, 'Nmax', 30);

figure();
imagesc(mod(im, 1));
axis image;

figure();
beam.visualise('axis', 'z', 'range', rng);
colormap hot;

%% Line shaped field

im = otslm.simple.sinc(slm_sz, 30, 'type', '1d');
im = otslm.tools.finalize(zeros(size(im)), 'amplitude', 0.8*im)./(2*pi);
Efield = gauss .* exp(1i*2*pi*im);
beam = ott.BscPmParaxial(NA, Efield, 'Nmax', 30);

figure();
imagesc(mod(im, 1));
axis image;

figure();
beam.visualise('axis', 'z', 'range', rng);
colormap hot;
