close all;
clear all;
clc;

%Fixed Variables
re = 0.35; % Ω·kcm 
ri = 0.11; % Ω·kcm 
cm = 2.5; % uF/cm^2 
ar = 0.0007;% Axon Radius (cm) 
w = 0.0001; % Node of Ranvier width 
dx = 0.1; % internodal gap (cm) 
msc = 1445; % Max Sodium Conductance 
mlc = 128; % M Leak C 
snp = 115; % Sodium Nernst Potential 
lnp = -0.01 % Leak N P 

%Adjustable Variables
x = 90; %Length of Axon scaled
f = 0.001; 
t = 0:f:2; %0.2 ms
is_cat = -4000; %current density trigger
is_an = -3 * is_cat; % Comment out for #3
ist = zeros(length(t), x);
ist(100:350, 50) = is_an; % Comment out for #3
ist(100:350, 45) = is_cat;
%dw = 0.1; %Adjust for distances away
%dw = 0.2;
%dw = 0.4;
dw = 0.8;

%Pre-simulation prep
v = zeros(length(t), x);
dvdt = zeros(length(t),x);
m = zeros(length(t), x);
dmdt = zeros(length(t),x);
h = zeros(length(t), x);
dhdt = zeros(length(t),x);
act = zeros(length(t),x);
cable = zeros(length(t),x);
im = zeros(length(t),x);
ve = zeros(length(t),x);
Vm = v(1,1);

alpha_m = (97+0.363*Vm)/(1+exp((31-Vm)/5.3));
beta_m = alpha_m/exp((Vm-23.8)/4.17);
beta_h = 15.6/(1+exp((24-Vm)/10));
alpha_h = beta_h/exp((Vm-5.5)/5);
m(:,:) = alpha_m / (alpha_m + beta_m);
h(:,:) = alpha_h / (alpha_h +beta_h);

% Update Extracellular voltage respect to axon diameter 
for i = 1:x
    for j = 1:x
        ve(:,i) = ve(:, i) + (ist(:,j) * re / (4 * pi * sqrt(0.1 * (i-j)^2 + dw^2)));
    end
end

for i = 1:length(t)
    for j = 2 : (x-1) 
        Vm = v(i,j);

        alpha_m = (97+0.363*Vm)/(1+exp((31-Vm)/5.3));
        beta_m = alpha_m/exp((Vm-23.8)/4.17);
        beta_h = 15.6/(1+exp((24-Vm)/10));
        alpha_h = beta_h/exp((Vm-5.5)/5);

        act(i,j) = (2*ar) * (ve(i,j-1) - (2 * ve(i,j)) + ve(i,j+1)) / (4 * ri * w * dx);
        cable(i,j) = (2*ar) * (v(i,j-1) - (2 * v(i,j)) + v(i,j+1)) / (4 * ri * w * dx); 

        im(i,j) = msc * m(i,j)^2 * h(i,j) * (Vm - snp) + mlc * (Vm - lnp); 
        dmdt(i,j) = (-1 * (alpha_m + beta_m) * m(i,j) + alpha_m);
        dhdt(i,j) = (-1 * (alpha_h + beta_h) * h(i,j) + alpha_h);
       
        dvdt(i,j) = (-im(i,j) + cable(i,j) + act(i,j)) / cm;
        end
        
    m(i+1,:) = m(i,:) + dmdt(i,:)*f; 
    h(i+1,:) = h(i,:) + dhdt(i,:)*f; 
    v(i+1,:) = v(i,:) + dvdt(i,:)*f; 
end

%plot

figure
xaxis = (1:x) / 10;
y = 0:f:2+f;
surf(xaxis,y,v,'edgecolor','none','facecolor','interp')
xlabel('Axonial distance (cm)')
ylabel('Time (ms)')
zlabel('Voltage (mV)')
%title('Bi-directional AP')
title('Uni-directional AP')
set(gca,'fontweight','light','fontsize',15); 


plt = figure('units','normalized','outerposition',[0 0 1 1]);
video = VideoWriter('Directional Propagating AP', 'MPEG-4');
open(video);

xaxis=1:x;

for i = 1:10:length(t) 
    plot(xaxis/10 , v(i,xaxis) )
    ylim([-100 100])
    %title('Bi-Directional AP')
    title('Uni-directional AP')
    xlabel('Axonial Distance (cm)')
    ylabel('Voltage (mV)')
    set(gca,'fontweight','light','fontsize',15)
    mov = getframe(plt);
    writeVideo(video,mov)
end
close(video)











