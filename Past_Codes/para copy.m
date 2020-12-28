clear
clc

% Notes #################################################################
%
% 1. Particle Reynolds number Re_par = rho*gamma*r^2/eta = 0.1
% 2. The effective RBC diameter, Deffm = sqrt(area / pi), is calculated.
% 3. Foppl-von Karman number is chosen as FvK = G*Deff/kb = 680
% 4. The bending constant kb is determined by FvK.
% 5. kv_m  = 25*Gm 
% 6. kag_m = 50*Gm
% 7. kal_m = 0
% 8. wallConst = 10*kb_m
% 9. wallForceDis = 0.7
%
%########################################################################


% Physical parameters
%###############################################################
Dp           = 7.82e-6;               % meter
Yp           = 18.9e-6;               % N/m
%Gp           = 6.3e-6;                % N/m
%kb_p         = 2e-19;                 % J 
eta_p        = 1.2e-3;                % Pa s
%kv_p         = 248.1832;              % J/m^3
%kag_p        = 2.3598e-4;             % J/m^2
%kal_p        = 4.8159e-6;             % J/m^2
kT300        = (1.38064852e-23)*300;  % thermal energy at T=300K (J)
%shear_p_base = 15.5;                  % s-1 this value is chosen such that Re(D)~0.5 
%shear_p      = 250;
%factor       = shear_p_base / shear_p; % scaling factor 
area = 134e-12;                        % m^2
Deff = sqrt(area/(pi));                % effective RBC diameter (m)

% Model parameters
%###############################################################

max_x    = 144;  % box size
max_y    = 24;
max_z    = 144;
cellSize = 8.0;   % for setting the neighbor list
level    = -3; % set particle shape and resolution
volFrac  = 3;  % particle volume fraction
config   = 1;  % set the intial configuration of particles

% particle diameter
%        n3     n2    n4
%RBC     15.64  7.82  31.28
%sphere  13.4

% particle volume
%        n4    n3        n2  
%RBC     6016  746.3163  90.91
%sphere        1248.991

if level == -3
  Dm     = 15.64;
  Rm     = 0.5*Dm;
  particle_vol = 746.3163;
  particle_area = 533.0397;
end
if level == 3
    Dm = 13.4;
    Rm = 0.5*Dm;
    particle_vol = 1248.991;
end
if level ==2
    Rm=3.35;
    Dm=2*Rm;
    particle_vol=152.1502;
end


Deffm = sqrt(particle_area/pi);
numPar  = floor(max_x*max_y*max_z*(volFrac/100.0)/particle_vol);


rho=36.0;        % model fluid density
nu=1.0/6;        % model kinematic viscosity
eta_m = rho*nu;  % model dynamic viscosity


% set wall velocity different
Re_par          = 0.1;
shear_m         = Re_par * eta_m / (rho*Rm*Rm);
d_u             = max_y*shear_m;                 % wall velocity diff.
u_max = 1.5*d_u;     % for setting Poiseuille flow
% set particle forces
Ca              = 0.2;
Gm              = eta_m*shear_m*Rm / Ca;  % model shear modulus 
FvK             = 680;
kb_m            = (Gm*Deffm*Deffm)/FvK;   % model bending modulus
kv_m            = 25*Gm;                  % vol. force strength
kag_m           = 50*Gm;                   % global area force strength
kal_m           = 0;                      % local area force strength 
%epsWCA          = 0.01*kb_m;              % strength of the WCA potential
wallConst       = 10*kb_m;                % wall force constant
wallForceDis    = 0.7;                    % parameter for wall force


springType  = 1;      % set the elastic potential for the spring network
x0          = 2.2;    % a parameter in WLC potential
engScale    = 1.0;    % not used now


attractionType = 5; % set the type of interparticle interaction
if attractionType == 5
  Dmorse = 0;
else
  Dmorse = 2.1; % set the depth of the Morse potential
end

eqWCA = 0.8;
eqWCA = 2^(1/6);

%if attractionType == 5
%  eqWCA  = 2^(1/6);     % zero force distance for the WCA potential
%else
%  eqWCA = 0.8;
%end

poiseuille = 0;    % set the flow profile
wallFlag   = 1;    % set the box boundary
if poiseuille == 1
  ux = u_max;
else
  ux = d_u;
end
uy = 0;
uz = 0;


%Gm     = 1.0;   % model shear modulus 
%kv_m   = 25.8;  % vol. force strength
%kag_m  = 50.0;  % global area force strength
%kal_m  = 1.;    % local area force strength 


% derivative
springConst = 2/(3^(0.5))*Gm;     % for harmonic springs (not used now) 
K  = 2*Gm + kag_m + kal_m;          % compression modulus
Ym = 4*K*Gm/(K+Gm);                 % Young's modulus
Poisson_ratio = (K-Gm)/(K+Gm);      % Poisson ratio
Ym_deviation  = (4*Gm-Ym)/(4*Gm);


t_lbe          = 1.0;    % LB time step
growthStep     = 20000;  % for generating a dense suspension
frictionalMass = 1000;   % for generating a dense suspension

    
nlistRenewal = 0.8;   % set the criterion for renewing the neighbor list
 

% Model scales
%##########################################################################
Ls = Dp/Dm;
Es = (Yp/Ym)*(Dp/Dm)^2;
Ns = (Yp/Ym)*(Dp/Dm);
Ts = (Dp/Dm)*(Ym/Yp)*(eta_p/eta_m);


% Use the model energy scale to estimate or set some quantities
%##############################################################################
Gp      = Gm*(Es/(Ls*Ls));      % physical shear modulus
kb_p    = kb_m * Es;            % model bending modulus
kv_p    = kv_m*(Es/(Ls*Ls*Ls)); % physical vol. force strength
kag_p   = kag_m*(Es/(Ls*Ls));   % physical global area force strength
kal_p   = kal_m*(Es/(Ls*Ls));   % physical local area force strength
kT300_m = kT300 / Es;           % model thermal energy at T=300k
shear_p = shear_m /Ts;
%shear_m = shear_p*Ts;           % model shear rate
%FvK     = (Gp*Deff*Deff)/kb_p;  %


% Set the flow velocity according to the channel height and 
% the desired shear rate 
%##############################################################################
epsLJ    = 0;
eqLJ     = 0;
cutoffLJ = 0;

epsWCA = kT300_m; % set the strength of the WCA potential

epsMorse    = Dmorse*kT300_m; % set the depth of the Morse potential
widthMorse  = 0.75;           % set the width of the Morse potential
eqMorse     = 0.8;            % set the zero force distance
cutoffMorse = 2.2;            % set the interaction cutoff distance

if attractionType == 5        % set the cutoff distance for the neighbor list
  nlistCutoff = eqWCA+1;
elseif attractionType == 4
  nlistCutoff = cutoffMorse +1;
end


strain = 1.0;           
step_1strain = strain / shear_m; % the # of steps for a single strain
numCycle = 1;  % the frequency of writing checkpoint files
numStep  = round(100*step_1strain);  % the total # of steps
writeProps = floor(1*step_1strain);  % the frequency of writing data
writeConfig = floor(2*step_1strain); % the frequency of writring particle configuration in vtk format
writeFluid = floor(100*step_1strain); % the frequency of writing flow field in vtk format 


% Estimate Ca and Re
%########################################
%Ca      = eta_m * shear_m * Dm / Gm;
%Re_par  = rho * shear_m * Dm*Dm / eta_m;
%Re_chan = rho * u_avg * h_eff /eta_m;


fprintf('Ls = %e\n', Ls);
fprintf('Es = %e\n', Es);
fprintf('Ns = %e\n', Ns);
fprintf('Ts = %e\n\n', Ts);

fprintf('Gp      = %e\n', Gp);
fprintf('kb_p    = %e\n', kb_p);
fprintf('kv_p    = %e\n', kv_p);
fprintf('kag_p   = %e\n', kag_p);
fprintf('kal_p   = %e\n', kal_p);
fprintf('shear_p = %f\n\n', shear_p);

fprintf('FvK           = %f\n', FvK);
fprintf('Poisson ratio = %f\n', Poisson_ratio);
fprintf('Ym deviation  = %f\n\n', Ym_deviation);

fprintf('Gm                      = %f\n', Gm);
fprintf('spirng const.           = %e\n', springConst);
fprintf('kb_m                    = %e\n', kb_m);
fprintf('kv_m                    = %e\n', kv_m);
fprintf('kag_m                   = %e\n', kag_m);
fprintf('kal_m                   = %e\n', kal_m);
fprintf('kwall_m                 = %e\n', wallConst);
fprintf('wall vel difference_m   = %e\n', d_u);
fprintf('max_vel_m (Poiseuille)  = %e\n', 1.5*d_u);
fprintf('shear_rate_m            = %e\n', shear_m);
fprintf('kT300_m                 = %e\n\n', kT300_m);

fprintf('Re(R) = %f\n', Re_par);
fprintf('Ca(R) = %f\n', Ca);
fprintf('Re(D) = %f\n', 4*Re_par);
fprintf('Ca(D) = %f\n\n', 2*Ca);

fprintf('1 stain steps             = %f\n', step_1strain);
fprintf('# of steps for each cycle = %f\n', numStep);
fprintf('total steps               = %d\n', numCycle*numStep);
fprintf('writeProps                = %d\n', writeProps);
fprintf('writeConfig               = %d\n', writeConfig);
fprintf('writeFluid                = %d\n', writeFluid);
fprintf('# of particles            = %d\n\n', numPar);

% Output file name
%#####################
%trial=1;  % for the naming the data folder
%name=['phi' num2str(volFrac) '-D' num2str(Dmorse) '-gamma' num2str(shear_p) '-' num2str(trial)]; % set the name of the data folder

name=['h' num2str(max_y) 'phi' num2str(volFrac) 'Re' num2str(Re_par) 'Ca' num2str(Ca) 'D' num2str(Dmorse) 'eqWCA' num2str(eqWCA)]; 
filename  = [name '.dat'];


numCycle    = 5;  % the frequency of writing checkpoint files
numStep     = 1000000;  % the total # of steps
writeProps  = 2000;  % the frequency of writing data
writeConfig = 4000; % the frequency of writring particle configuration in vtk format
writeFluid  = 100000; % the frequency of writing flow field in vtk format 


fw = fopen(filename,'w');
fprintf(fw,'NumParticle    level    configuration\n');
fprintf(fw,'%d             %d       %d\n', numPar, level, config);

fprintf(fw,'Lx    Ly    Lz\n');
fprintf(fw,'%d    %d    %d\n', max_x, max_y, max_z);

fprintf(fw,'NumCycle    NumStep    t_lbe\n');
fprintf(fw,'%d          %d         %f\n', numCycle, numStep, t_lbe);

fprintf(fw,'Poiseuille    WallFlag    ux    uy    uz\n');
fprintf(fw,'%d            %d          %e    %e    %e\n',poiseuille, wallFlag, ux, uy, uz);

fprintf(fw,'GrowthStep    FrictionalMass\n');
fprintf(fw,'%d            %f\n', growthStep, frictionalMass);

fprintf(fw,'WriteProps    WriteConfig    WriteFluid\n');
fprintf(fw,'%d            %d             %d\n', writeProps, writeConfig, writeFluid);

fprintf(fw,'Gm    x0\n');
fprintf(fw,'%f    %f\n', Gm, x0);

fprintf(fw,'EngScale    SpringType    AttractionType\n');
fprintf(fw,'%f          %d            %d\n', engScale, springType, attractionType);

fprintf(fw,'SpringConst    kb_m    kv_m    kag_m    kal_m\n');
fprintf(fw,'%e             %e      %e      %e       %e\n', springConst, kb_m, kv_m, kag_m, kal_m);

fprintf(fw,'epsLJ    eqLJ    cutoffLJ\n');
fprintf(fw,'%e       %f      %f\n', epsLJ, eqLJ, cutoffLJ);

fprintf(fw,'epsWCA    eqWCA\n');
fprintf(fw,'%e        %f\n', epsWCA, eqWCA);

fprintf(fw,'Dmorse    widthMorse    eqMorse    cutoffMorse\n');
fprintf(fw,'%e        %f            %f         %f\n', epsMorse, widthMorse, eqMorse, cutoffMorse);

fprintf(fw,'nlistCutoff    CellSize    nlistRenewal\n');
fprintf(fw,'%f             %f          %f\n', nlistCutoff, cellSize, nlistRenewal);

fprintf(fw,'WallConst    WallForceDis\n');
fprintf(fw,'%f           %f\n', wallConst, wallForceDis);

fprintf(fw, '\n');
fprintf(fw, '============================================================================\n');
fprintf(fw, 'some information\n\n');

fprintf(fw,'Ls = %e\n',   Ls);
fprintf(fw,'Es = %e\n',   Es);
fprintf(fw,'Ns = %e\n',   Ns);
fprintf(fw,'Ts = %e\n\n', Ts);

fprintf(fw,'Compression modulus (K) = %e\n', K);
fprintf(fw,'Ym                      = %e\n', Ym);
fprintf(fw,'Poisson ratio           = %f\n', Poisson_ratio);
fprintf(fw,'Ym deviation            = %f\n', Ym_deviation);
fprintf(fw,'FvK                     = %f\n\n', FvK);

fprintf(fw,'Ca(R) = %f\n', Ca);
fprintf(fw,'Re(R) = %f\n', Re_par);
fprintf(fw,'Ca(D)   = %f\n', 2*Ca);
fprintf(fw,'Re(D)   = %f\n\n', 4*Re_par);

fprintf(fw,'Dp             = %e\n', Dp);
fprintf(fw,'Yp             = %e\n', Yp);
fprintf(fw,'Gp             = %e\n', Gp);
fprintf(fw,'kb_p           = %e\n', kb_p);
fprintf(fw,'eta_p          = %e\n', eta_p);
fprintf(fw,'kv_p           = %e\n', kv_p);
fprintf(fw,'kag_p          = %e\n', kag_p);
fprintf(fw,'kal_p          = %e\n', kal_p);
%fprintf(fw,'shear_p_base   = %e\n', shear_p_base);
fprintf(fw,'shear_p        = %e\n', shear_p);
%fprintf(fw,'Scaling factor = %f\n', factor);
fprintf(fw,'Vol fraction   = %f\n\n', volFrac);

fprintf(fw,'Gm                      = %f\n', Gm);
fprintf(fw,'kb_m                    = %e\n', kb_m);
fprintf(fw,'kv_m                    = %e\n', kv_m);
fprintf(fw,'kag_m                   = %e\n', kag_m);
fprintf(fw,'kal_m                   = %e\n', kal_m);
fprintf(fw,'eff_wall_vel_m          = %e\n', d_u);
fprintf(fw,'max_vel_m (Poiseuille)  = %e\n', 1.5*d_u);
fprintf(fw,'shear_rate_m            = %e\n', shear_m);
fprintf(fw,'kT300_m                 = %e\n', kT300_m);
fprintf(fw,'# of steps for 1 strain = %f\n\n', step_1strain);

