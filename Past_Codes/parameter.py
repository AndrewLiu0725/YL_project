import math
import sys

# The format of the argument will be aggregation_?KT_flow_?_ncycle_?
cycle_number, aggregation_intensity, shear_rate = 50, 1, 0.0
components = sys.argv[1].split('_')
for setting in ['aggregation', 'flow', 'ncycle']:
    try:
        if setting == 'aggregation':
            aggregation_intensity = float(components[components.index(setting)+1][:-2])
        elif setting == 'ncycle':
            cycle_number = int(components[components.index(setting)+1])
        elif setting == 'flow':
            shear_rate = float(components[components.index(setting)+1])
    except:
        pass

# Physical parameters
################################################################
Dp           = 7.82e-6                # meter
Yp           = 18.9e-6                # N/m
Gp           = 6.3e-6                 # N/m
kb_p         = 2e-19                  # J 
eta_p        = 1.2e-3                 # Pa s
kv_p         = 248.1832               # J/m^3
kag_p        = 2.3598e-4              # J/m^2
kal_p        = 4.8159e-6              # J/m^2
kT300        = (1.38064852e-23)*300   # thermal energy at T=300K (J)
shear_p_base = 37.6                   # s-1 this value is chosen such that Re(D)~0.5 
shear_p      = shear_rate             # desired shear rate, need assigned
try: factor       = shear_p_base / shear_p  # scaling factor 
except: factor = 1
area = 134e-12                         # m^2
Deff = math.sqrt(area/(math.pi))                 # effective RBC diameter (m)

# Model parameters
################################################################
max_x = 96   # box size
max_y = 64 
max_z = 96  
level = -3  # set particle shape and resolution
numPar = 2 
config  = 2   # set the intial configuration of particles

if level == -3:
    Dm     = 15.64 
    Rm     = 0.5*Dm 
    particle_vol = 746.3163 

elif level == 3:
    Dm = 13.4 
    Rm = 0.5*Dm 
    particle_vol = 1248.991 

volFrac = numPar*particle_vol/(max_x*max_y*max_z)

Gm     = 1.0 * factor     # model shear modulus 
kv_m   = 25.8  * factor   # vol. force strength
kag_m  = 50.0 * factor   # global area force strength
kal_m  = 1. * factor     # local area force strength 

# derivative
K  = 2*Gm + kag_m + kal_m           # compression modulus
Ym = 4*K*Gm/(K+Gm)                  # Youngs modulus
Poisson_ratio = (K-Gm)/(K+Gm)       # Poisson ratio
Ym_deviation  = (4*Gm-Ym)/(4*Gm) 
rho=36.0         # model fluid density
nu=1.0/6         # model kinematic viscosity
eta_m = rho*nu   # model dynamic viscosity
t_lbe          = 1.0     # LB time step
growthStep     = 20000   # for generating a dense suspension
frictionalMass = 1000    # for generating a dense suspension
x0          = 2.2     # a parameter in WLC potential
engScale    = 1.0     # not used now
springType  = 1       # set the elastic potential for the spring network
springConst = 0.      # for harmonic springs (not used now)  
cellSize     = 8.0    # for setting the neighbor list
nlistRenewal = 0.8    # set the criterion for renewing the neighbor list
wallConst    = kv_m   # wall force constant
wallForceDis = 1.0    # parameter for wall force

# Model scales
###########################################################################
Ls = Dp/Dm 
Es = (Yp/Ym)*(Dp/Dm)**2 
Ns = (Yp/Ym)*(Dp/Dm) 
Ts = (Dp/Dm)*(Ym/Yp)*(eta_p/eta_m) 

# Use the model energy scale to estimate or set some quantities
###############################################################################
Gp      = Gm*(Es/(Ls*Ls))       # physical shear modulus
kv_p    = kv_m*(Es/(Ls*Ls*Ls))  # physical vol. force strength
kag_p   = kag_m*(Es/(Ls*Ls))    # physical global area force strength
kal_p   = kal_m*(Es/(Ls*Ls))    # physical local area force strength
kb_m    = kb_p / Es             # model bending modulus
shear_m = shear_p*Ts            # model shear rate
kT300_m = kT300 / Es            # model thermal energy at T=300k
FvK     = (Gp*Deff*Deff)/kb_p   #

# Set the flow velocity according to the channel height and the desired shear rate 
###############################################################################
ux = max_y*shear_m  # for setting the wall velocity
uy = 0 
uz = 0 
poiseuille = 0

wallFlag   = 1     # set the box boundary

numCycle = cycle_number   # the frequency of writing checkpoint files, need assigned, recommended value = 50
numStep  = 25600   # the total # of steps

try:
    strain = 1.0            
    step_1strain = strain / shear_m  # the # of steps for a single strain
    writeProps = math.floor(1*step_1strain)   # the frequency of writing data
    writeConfig = math.floor(1*step_1strain)  # the frequency of writring particle configuration in vtk format
    writeFluid = 100*math.floor(1*step_1strain)  # the frequency of writing flow field in vtk format 
except:
    step_1strain = 0
    writeProps = 1280
    writeConfig = 1280
    writeFluid = 128000  

attractionType = 4  # set the type of interparticle interaction

epsLJ    = 0 
eqLJ     = 0 
cutoffLJ = 0 

Dmorse = aggregation_intensity  # set the depth of the Morse potential, need assigned

epsWCA = Dmorse*kT300_m  # set the strength of the WCA potential
eqWCA  = 0.8      # zero force distance for the WCA potential

epsMorse    = Dmorse*kT300_m  # set the depth of the Morse potential
widthMorse  = 0.75            # set the width of the Morse potential
eqMorse     = 0.8             # set the zero force distance
cutoffMorse = 2.5             # set the interaction cutoff distance

if attractionType == 5:        # set the cutoff distance for the neighbor list
    nlistCutoff = eqWCA+1 
elif attractionType == 4:
    nlistCutoff = cutoffMorse +1 

# Estimate Ca and Re
#########################################
Ca      = eta_m * shear_m * Dm / Gm 
Re_par  = rho * shear_m * Dm*Dm / eta_m 


# Write the file
#########################################
#filename = 'parameter/' +  sys.argv[1] + '.dat'
filename =  sys.argv[1] + '.dat'
f = open(filename, 'w+')

f.write( 'NumParticle    level    configuration\n') 
f.write( '{}             {}       {}\n'.format(numPar, level, config)) 
f.write( 'Lx    Ly    Lz\n') 
f.write( '{}    {}    {}\n'.format(max_x, max_y, max_z)) 
f.write( 'NumCycle    NumStep    t_lbe\n') 
f.write( '{}          {}         {}\n'.format(numCycle, numStep, t_lbe)) 
f.write( 'Poiseuille    WallFlag    ux    uy    uz\n') 
f.write( '{}            {}          {}    {}    {}\n'.format(poiseuille, wallFlag, ux, uy, uz)) 
f.write( 'GrowthStep    FrictionalMass\n') 
f.write( '{}            {}\n'.format(growthStep, frictionalMass)) 
f.write( 'WriteProps    WriteConfig    WriteFluid\n') 
f.write( '{}            {}             {}\n'.format(writeProps, writeConfig, writeFluid)) 
f.write( 'Gm    x0\n') 
f.write( '{}    {}\n'.format(Gm, x0)) 
f.write( 'EngScale    SpringType    AttractionType\n') 
f.write( '{}          {}            {}\n'.format(engScale, springType, attractionType)) 
f.write( 'SpringConst    kb_m    kv_m    kag_m    kal_m\n') 
f.write( '{}             {}      {}      {}       {}\n'.format(springConst, kb_m, kv_m, kag_m, kal_m)) 
f.write( 'epsLJ    eqLJ    cutoffLJ\n') 
f.write( '{}       {}      {}\n'.format(epsLJ, eqLJ, cutoffLJ)) 
f.write( 'epsWCA    eqWCA\n') 
f.write( '{}        {}\n'.format(epsWCA, eqWCA)) 
f.write( 'Dmorse    widthMorse    eqMorse    cutoffMorse\n') 
f.write( '{}        {}            {}         {}\n'.format(epsMorse, widthMorse, eqMorse, cutoffMorse)) 
f.write( 'nlistCutoff    CellSize    nlistRenewal\n') 
f.write( '{}             {}          {}\n'.format(nlistCutoff, cellSize, nlistRenewal)) 
f.write( 'WallConst    WallForceDis\n') 
f.write( '{}           {}\n'.format(wallConst, wallForceDis)) 

f.write('\n###########################################################################\n')
f.write( 'Ls = %e\n' % Ls) 
f.write( 'Es = %e\n' % Es) 
f.write( 'Ns = %e\n' % Ns) 
f.write( 'Ts = %e\n\n' % Ts) 
f.write( 'Compression modulus (K) = %e\n' % K) 
f.write( 'Ym                      = %e\n' % Ym) 
f.write( 'Poisson ratio           = %f\n' % Poisson_ratio) 
f.write( 'Ym deviation            = %f\n' % Ym_deviation) 
f.write( 'FvK                     = %f\n\n' % FvK) 
f.write( 'Ca(diameter) = %f\n' % Ca) 
f.write( 'Re(diameter) = %f\n' % Re_par) 
f.write( 'Ca(radius)   = %f\n' % (0.5*Ca)) 
f.write( 'Re(radius)   = %f\n\n' % (0.25*Re_par)) 
f.write( 'Dp             = %e\n' % Dp) 
f.write( 'Yp             = %e\n' % Yp) 
f.write( 'Gp             = %e\n' % Gp) 
f.write( 'kb_p           = %e\n' % kb_p) 
f.write( 'eta_p          = %e\n' % eta_p) 
f.write( 'kv_p           = %e\n' % kv_p) 
f.write( 'kag_p          = %e\n' % kag_p) 
f.write( 'kal_p          = %e\n' % kal_p) 
f.write( 'shear_p_base   = %e\n' % shear_p_base) 
f.write( 'shear_p        = %e\n' % shear_p) 
f.write( 'Scaling factor = %f\n' % factor) 
f.write( 'Vol fraction   = %f\n\n' % volFrac) 
f.write( 'Gm                      = %f\n' % (Gm/factor)) 
f.write( 'kb_m                    = %e\n' % kb_m) 
f.write( 'kv_m                    = %e\n' % kv_m) 
f.write( 'kag_m                   = %e\n' % kag_m) 
f.write( 'kal_m                   = %e\n' % kal_m) 
f.write( 'eff_wall_vel_m          = %e\n' % ux) 
f.write( 'max_vel_m (Poiseuille)  = %e\n' % (1.5*ux)) 
f.write( 'shear_rate_m            = %e\n' % shear_m) 
f.write( 'kT300_m                 = %e\n' % kT300_m) 
f.write('# of steps for 1 strain = %f\n' % step_1strain)
f.close()