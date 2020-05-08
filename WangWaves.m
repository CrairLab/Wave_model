%% Model based on Keane&Gong paper: https://www.jneurosci.org/content/jneuro/35/4/1591.full.pdf
%Date:05/07/2020; Author: yixiang.wang@yale.edu

%% Set parameters.
dt = 0.5; %time step 0.5ms
c = 1; %capacitance 1uF 
gl = 50; %Leak conductance 50uS 
vl = -70; %Leak voltage -70mV
vr = -70; %Resting poteintial -70mV
ve = 0; %Excitatory reversal potential 0mV
vi = -80; %Inhibitory reversal potential -80mV
vt = -55; %spiking threshold -55mV
tau_rf = 5; %refractory period after spiking 5ms
fe = 15; %constrant external input (excitatory) 15uS 
fi = 2; %constant external input (inhibitory) 2uS 
tau_re = 0.5; %rise time (excitatory) 0.5ms
tau_ri = 0.5; %rise time (inhibitory) 0.5ms
tau_de = 2.0; %decay time (excitatory) 2.0ms
tau_di = 7.0; %decay time (inhibitory) 7.0ms
re = 5; %radius of excitatory connection
ri = 15; %radius of inhibitory connection
de = 12; %spatial scale for Gaussian decay of excitatory couplings
we = 0.23; %excitatory coupling strength
wi = 0.23; %inhibitory coupling strength (balanced btw 0.23-0.35)
width = 200; %width of the square network
max_memory = 30; %memory of recent spikes, most recent 25ms;
max_steps = max_memory/dt; %corresponding maximum memory in steps
runTime = 7.5*10^3; %Run 7.5s



%% Construct networks
dis_i = sqrt(5); %distance between inhibitory neurons, evenly distributed. E:I = 4:1;
ID = zeros(width,width); %Creat a width*width grid identity matrix
ID(1:dis_i:width,1:dis_i:width) = 1; %1: inhibitory neurons; 0: excitatory neurons
ID = logical(ID);


%% Construct K filters (coupling)

%Excitatory coupling filter
x_idx = [-re:re]';
y_idx = [-re:re];
dis_square = x_idx.^2 + y_idx.^2;
dis_square(dis_square > re^2) = nan;
KEf = we.*exp(-dis_square./de);
KEf(isnan(KEf)) = 0;
KEf(re+1,re+1) = 0;

%Inhibitory coupling filter
x_idx = [-ri:ri]';
y_idx = [-ri:ri];
dis_square = x_idx.^2 + y_idx.^2;
KIf = wi.*(dis_square<=ri^2);
%KIf(:,1:ri) = KIf(:,1:ri) * 0.95;
%KIf(:,ri+2:end) = KIf(:,ri+2:end) * 1.1;
%KIf(1:ri,1:ri) = KIf(1:ri,1:ri) * 0.5;
KIf(ri+1,ri+1) = 0;


%{
KE_all = zeros(2*re+1,2*re+1,width,width,'single');
KI_all = zeros(2*ri+1,2*ri+1,width,width,'single');
for i = 1:width
    for j = 1:width
        %current neighboorhood
        curNb_e = ID_p(i+rp-re:i+rp+re, j+rp-re:j+rp+re);   
        curNb_i = ID_p(i+rp-ri:i+rp+ri, j+rp-ri:j+rp+ri); 
        %construct all KE matrices
        KE_all(:,:,i,j) = KEf.*~curNb_e;
        %construct all KI matrices
        KI_all(:,:,i,j) = KIf.*curNb_i;
    end
end

K_p = conv2(ID_p, KEf, 'same') + conv2(ID_p, KIf, 'same');
K = K_p(rp+1:rp+width, rp+1:rp+width);
%}

%% Calculate time course of post-synaptic conductance
t = dt:dt:max_memory;
ge_t = (exp(-t/tau_de)-exp(-t/tau_re))./(tau_de-tau_re);
ge_t = ge_t';
gi_t = (exp(-t/tau_di)-exp(-t/tau_ri))./(tau_di-tau_ri);
gi_t = gi_t';


%% Run iteration
V = rand(width, 'single').*(vt-vr) + vr; %Initialize Voltage matrix;
GE = fe.*ones(width, 'single'); %initialize excitatory Conductance matrix
GI = fi.*ones(width, 'single'); %initialize inhibitory Conductance matrix
T = false(width*width,max_steps); %initialize spike matrix, only memorize the recent 50 steps
Refractory = zeros(width, 'single'); %initialize refractory matrix
s = 0; %step 0

%Calculate extended/padded identity matrix
ID_p = repmat(ID,[3,3]); %extend 2D identity matrix
rp = max(re,ri); %padding radius
ID_p = ID_p(1+width-rp:2*width+rp, 1+width-rp:2*width+rp); %padded identity matrix


for t = 0:dt:runTime
    
    %update refractory matrix
    Refractory = Refractory - dt; 
    
    %select active neurons at this time
    Active = Refractory<= 0; 
    
    %Current excitatory/inhibitory postsynaptic conductance matrices
    curGe = reshape(T*ge_t, [width,width]);
    curGi = reshape(T*gi_t, [width,width]);
      
    %Padding postsynaptic conductance matrices
    curGe_p = repmat(curGe,[3,3]); %extend 2D identity matrix
    curGe_p = curGe_p(1+width-rp:2*width+rp, 1+width-rp:2*width+rp); %padded identity matrix  
    curGi_p = repmat(curGi,[3,3]); %extend 2D identity matrix  
    curGi_p = curGi_p(1+width-rp:2*width+rp, 1+width-rp:2*width+rp); %padded identity matrix
    
    %Convolve postsynaptic conductance with couplting strength matrices
    %KG_p = conv2(~ID_p.*curGe_p, KEf, 'same') + conv2(ID_p.*curGi_p, KIf, 'same');
    KGe_p = conv2(~ID_p.*curGe_p, KEf, 'same');
    KGe = KGe_p(rp+1:rp+width, rp+1:rp+width);
    KGi_p = conv2(ID_p.*curGi_p, KIf, 'same');
    KGi = KGi_p(rp+1:rp+width, rp+1:rp+width);
    
    %The final conductance matrices at this moment
    GE = fe + KGe*1e3;
    GI = fi + KGi*1e3;
    
    %calculate change of voltages for active neruons
    dV = Active.*(1e-3.*dt*(-gl*(V-vl) - GE.*(V-ve) - GI.*(V-vi))/c); 
    
    %update V
    V = V+dV; 
    
    %Display setteing
    imshow(mat2gray(V))
    %hold on
    colormap jet
    pause(0.001)
    
    %Neurons that spike at this time
    curT = (V>=-55);
    %curT = Active&(V>=-55); %unecessary
       
    %Neurons fired have V back to vr
    V(curT) = vr;
    
    %Mark neurons enter refractory periods
    Refractory(curT) = tau_rf;
    
    %Update spike matrix
    T(:,2:max_steps) = T(:,1:max_steps-1); 
    T(:,1) = curT(:);
    
    %Show progress
    disp([num2str(t) ' ms']);
end
