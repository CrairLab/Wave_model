%% Compare activations of rhodopsin/ rods under different conditions

S_relative = Lamb(499, 'rhodopsin');
u = Lei;
D = 2.07; %mm mouse lense thickness
UV_incident = Incident_Intensity(1);
Blue_incident = Incident_Intensity(2);
IR_incident = Incident_Intensity(3);
TP_incident = Incident_Intensity(4);


UV_transmitted = UV_incident.* 10.^(-u * D);
plotRelativeInensity(UV_transmitted, 'UV transmitted');
Blue_transmitted = Blue_incident.* 10.^(-u * D);
plotRelativeInensity(Blue_transmitted, 'Blue transmitted');
IR_transmitted = IR_incident.* 10.^(-u * D);
plotRelativeInensity(IR_transmitted, 'IR transmitted');
TP_transmitted = TP_incident.* 10.^(-u * D);
plotRelativeInensity(TP_transmitted, '2p transmitted');

UV_activation = sumActivation(UV_transmitted, S_relative);
Blue_activation = sumActivation(Blue_transmitted, S_relative);
IR_activation = sumActivation(IR_transmitted, S_relative);
TP_activation = sumActivation(TP_transmitted, S_relative);


function activation = sumActivation(input_transmitted, S_relative)

inter_transmitted = (input_transmitted(1:end-1) + input_transmitted(2:end))/2;
inter_s = (S_relative(1:end-1) + S_relative(2:end))/2;
activation = sum(inter_transmitted .* inter_s);


end


function plotRelativeInensity(input_trace, input_title)

L = 350:1:1050;
figure; plot(L, input_trace); set(gca, 'YScale', 'log')
xlabel('Wavelength (nm)')
ylabel('Relative intensity')
title(input_title)

end

function S_relative = Lamb(Lmax, opsin_name)
%Lamb 1995 Opsin Template http://www.cvrl.org/database/text/pigments/lamb.htm

A = 0.880;
B = 0.924;
C = 1.104;
D = 0.655;
a = 70;
b = 28.5;
c = -14.1;
L = 350:1:1050;

S_relative = 1./(exp(a*(A - Lmax./L)) + exp(b*(B - Lmax./L)) + exp(c*(C - Lmax./L)) + D);
figure; plot(L, S_relative); set(gca, 'YScale', 'log')
xlabel('Wavelength (nm)')
ylabel('Relative sensitivity')
title(['Spectral sensitivity of ' opsin_name ' (max lambda = ' num2str(Lmax) 'nm)'])

end





function u = Lei
%Spectral attenuation https://www.sciencedirect.com/science/article/pii/S0014483506001825
% Omit the water term as the difference it might cause is negligible for
% the range considered in this study

% C57/B6 mice
a = 0.01939;
b = 2.426 / 1000;
c = 27277;
d = -2.27;
l0 = 353.6;
u0 = 0.00821;
L = 350:1:1050;

%{
% Albino mice%}
a = 0.06341;
b = 5.951 / 10000;
c = 961310;
d = -3.855;
l0 = 320.1;
u0 = 0.00534;
L = 350:1:1050;
%}

u = a./(1+b.*(L - l0).^2) + c.*L.^d + u0;
figure; plot(L, u); set(gca, 'YScale', 'log')
xlabel('Wavelength (nm)')
ylabel('Attenuation mm-1')
title('Attenuation coefficient (mouse lens)')


end



function I_incident = Incident_Intensity(flag, filter_set, power)
% 1: UVX, 2: BDX, 3: RLX, 4: 920nm Laser
%https://www.researchgate.net/figure/A-TiSapphire-laser-Mai-Tai-HP-from-Spectra-Physics-was-used-for-comparison-with-the_fig5_318560821

L = 350:1:1050;

if nargin < 3
    
    % Define bandpass range
    passed = zeros(1, length(L));
    ini = L(1)+1;
    blue = passed; blue(450-ini:490-ini) = 1;
    filter_set.blue = blue;% Chroma 49002 (470/40)
    uv = passed; uv(382-ini:408-ini) = 1;
    filter_set.uv = uv; %Chroma 49028 (395/25)
    %ir = passed; ir(610-ini:640-ini) = 1; %Chroma ET 625/30
    ir = passed; ir(580-ini:630-ini) = 1; %Chroma ET 605/50
    filter_set.ir = ir; 
    tp = passed; tp(900-ini:940-ini) = 1;
    filter_set.tp = tp; % Two photon
    
    % Define illumination power
    power.uv = 3; %3mW at 5%
    power.blue = 4; %4mW at 5%
    power.ir = 100; %100mW at 100%
    power.tp = 100; %100mW 
 
end




switch flag
    
    case 1
        sigma = 6.25;
        u = 387.5;
        I_max = 20;
        filter = filter_set.uv;
        power_percent = power.uv;
    case 2
        sigma = 11.25;
        u = 472.5;
        I_max = 22;
        filter = filter_set.blue;
        power_percent = power.blue;
    case 3
        sigma = 10;
        u = 635;
        I_max = 25;
        filter = filter_set.ir;
        power_percent = power.ir;
    case 4
        sigma = 5;
        u = 920;
        Intensity_factor = 16; %Compared to mesoscope (2mm/ 500um)^2
        I_max = 100 * Intensity_factor; %100 mW * intensity power
        filter = filter_set.tp;
        power_percent = power.tp;
end



I_incident = I_max .* exp(-0.5*((L - u)./(2.*sigma)).^2) .* filter .* power_percent;

end