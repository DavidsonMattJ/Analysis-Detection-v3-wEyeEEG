% makeGabor_stimulus
%% Define parameters for the Gabor patch
lambda = 20;        % Wavelength of the sinusoidal carrier (in pixels)
theta = 45;         % Orientation of the sinusoidal carrier (in degrees)
sigma = 23;          % Standard deviation of the Gaussian envelope
phase = 0;          % Phase of the carrier wave
sizex = 128;         % Size of the Gabor patch (e.g., 128x128 pixels)

% Create a mesh grid of coordinates centered at the middle of the patch
[x, y] = meshgrid(-sizex/2:sizex/2-1, -sizex/2:sizex/2-1);

% Calculate the Gabor patch
gabor_patch = exp(-(x.^2 + y.^2) / (2 * sigma^2)) .*cos(2 * pi * x / lambda + deg2rad(theta) + phase);

% Define parameters for the circle
center = [0,0];  % Center of the circle
radius = sizex/2;           % Radius of the circle

% Calculate the distance of each point to the center
distance_to_center = sqrt((x - center(1)).^2 + (y - center(2)).^2);

% Create a logical mask to identify points outside the circle
outside_circle = distance_to_center > radius;
% make a noise patch
noise= rand(length(x), length(y));
% increase pixelation:


% make the stimuli.
stim1 = gabor_patch+noise;
stim2= noise;

% Set values outside the circle to NaN
stim1(outside_circle) = 1;
stim2(outside_circle)= 1; % 1 is white colour

% Display the Gabor patch using imshow
subplot(121)
imshow(stim1);
subplot(122);
imshow(stim2); 
% axis square
shg
%% course noise:
sizey = 48;         % Size of the Gabor patch (e.g., 128x128 pixels)
% Create a mesh grid of coordinates centered at the middle of the patch
[x, y] = meshgrid(-sizey/2:sizey/2-1, -sizey/2:sizey/2-1);
noisey= round(rand(length(x), length(y)));

% Calculate the distance of each point to the center
distance_to_center = sqrt((x - center(1)).^2 + (y - center(2)).^2);
radius= sizey/3;
% Create a logical mask to identify points outside the circle
outside_circle = distance_to_center > radius;
inside_circle = distance_to_center < radius*.5;
figure(2);
% Display the Gabor patch using imshow

% Set values outside the circle to NaN
noisey(outside_circle) = .5;
noisey(inside_circle)= .5; % 1 is white colour
subplot(121)
imshow(noisey);
shg
%% make auditory stimulus:

% sine wave:
% white noise
t=0:.1:1000;
Fs= 2000;
sw = sin(2*pi*t./Fs);
clf

plot(t,sw.^2);
shg
%
carrierFs = 40;
cw = sin(2*pi*t./carrierFs);
clf
hold on;
subplot(1,2,1);
% plot some white noise.
wnoise= rand([1,length(t)]);

plot(t, wnoise-.5);
hold on;
plot(t,(cw.*sw),'linew',2)
ylabel('Amplitude [a.u.]')
set(gca,'fontsize',24,'XColor','w' )
box off
subplot(1,2,2);
% plot some white noise.
wnoise= rand([1,length(t)]);
box off
plot(t, wnoise-.5);
% hold on;
% plot(t,abs(cw.*sw),'linew',2)
ylim([-1 1]);
% set(gca,'ytick',[-1 0 1]); 
set(gca,'fontsize',12)
box off
set(gca,'fontsize',24,'XColor','w')