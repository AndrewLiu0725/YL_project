# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 02/28/2021
# ===============================================================================
import math
import sys
import subprocess
import os

"""
This code is to make parameter.dat based on the given variables and then submit it.
It will also check whether this job is already done or if there is a shorter one.
"""

# The format of the argument will be (test_)phi_?_Re_?_Ca_?_aggregation_?KT_ncycle_?_np_?_angle_?
# Initilaization
cycle_number, aggregation_intensity = 10, 1.0
volume_fraction, particle_number = 5.0, 2
components = sys.argv[1].split('_')
config  = 2   # set the intial configuration of particles
test = 0

aggregation_intensity = float(components[components.index('aggregation')+1][:-2])
cycle_number = int(components[components.index('ncycle')+1]) # number of strains
Re_input = float(components[components.index('Re')+1])
Ca_input = float(components[components.index('Ca')+1])
volume_fraction = 0.01*float(components[components.index('phi')+1])
particle_number = int(components[components.index('np')+1])
if 'angle' in components: config  = 3 # use the rotated init.config
if 'test' in components: test = 1

# Physical parameters
###############################################################################
Dp           = 7.82e-6                # meter
Yp           = 18.9e-6                # N/m
eta_p        = 1.2e-3                 # Pa s
kT300        = (1.38064852e-23)*300   # thermal energy at T=300K (J)
area = 134e-12                         # m**2
Deff = math.sqrt(area/(math.pi))                 # effective RBC diameter (m)

# Model parameters
###############################################################################
max_y = 24   
cellSize = 8.0
level = -3  # set particle shape and resolution
numPar = particle_number

if level == -3:
    Dm     = 15.64 
    Rm     = 0.5*Dm 
    particle_vol = 746.3163 
    particle_area = 533.0397

elif level == 3:
    Dm = 13.4 
    Rm = 0.5*Dm 
    particle_vol = 1248.991 
    particle_vol = 1248.991

Deffm = math.sqrt(particle_area/math.pi)

z = math.sqrt(numPar*particle_vol/(volume_fraction*max_y))
max_zs, max_xs, volFracs = [], [], []
max_z = math.ceil(z)
max_x = math.ceil(z) # use the same x-to-z ratio
volFrac = numPar*particle_vol/(max_x*max_y*max_z)
error = abs(volFrac - volume_fraction )
if (abs(numPar*particle_vol/((max_x-1)*max_y*(max_z-1)) - volume_fraction ) < error):
	volFrac = numPar*particle_vol/((max_x-1)*max_y*(max_z-1))
	max_z = max_z - 1
	max_x = max_x - 1
rho = 36.0        # model fluid density
nu = 1.0/6        # model kinematic viscosity
eta_m = rho*nu  # model dynamic viscosity

# Set wall velocity different
###############################################################################
Re_par          = Re_input
shear_m         = Re_par * eta_m / (rho*Rm*Rm)
d_u             = max_y*shear_m                 # wall velocity diff.
u_max = 1.5*d_u     # for setting Poiseuille flow
poiseuille = 0    # set the flow profile
wallFlag   = 1    # set the box boundary
if poiseuille == 1: ux = u_max
else: ux = d_u
uy = 0
uz = 0

# Set particle forces
###############################################################################
Ca              = Ca_input
Gm              = eta_m*shear_m*Rm / Ca  # model shear modulus 
FvK             = 680
kb_m            = (Gm*Deffm*Deffm)/FvK   # model bending modulus
kv_m            = 25*Gm                  # vol. force strength
kag_m           = 50*Gm                  # global area force strength
kal_m           = 0                      # local area force strength 
wallConst       = 10*kb_m                # wall force constant
wallForceDis    = 0.7                    # parameter for wall force
springType  = 1      # set the elastic potential for the spring network
x0          = 2.2    # a parameter in WLC potential
engScale    = 1.0    # not used now


# Derivative
###########################################################################
springConst = 2/(3**(0.5))*Gm     # for harmonic springs (not used now) 
K  = 2*Gm + kag_m + kal_m          # compression modulus
Ym = 4*K*Gm/(K+Gm)                 # Young's modulus
Poisson_ratio = (K-Gm)/(K+Gm)      # Poisson ratio
Ym_deviation  = (4*Gm-Ym)/(4*Gm)


t_lbe          = 1.0    # LB time step
growthStep     = 20000  # for generating a dense suspension
frictionalMass = 1000   # for generating a dense suspension
nlistRenewal = 0.8   # set the criterion for renewing the neighbor list
 

# Model scales
###########################################################################
Ls = Dp/Dm
Es = (Yp/Ym)*(Dp/Dm)**2
Ns = (Yp/Ym)*(Dp/Dm)
Ts = (Dp/Dm)*(Ym/Yp)*(eta_p/eta_m)


# Use the model energy scale to estimate or set some quantities
###############################################################################
Gp      = Gm*(Es/(Ls*Ls))      # physical shear modulus
kb_p    = kb_m * Es            # model bending modulus
kv_p    = kv_m*(Es/(Ls*Ls*Ls)) # physical vol. force strength
kag_p   = kag_m*(Es/(Ls*Ls))   # physical global area force strength
kal_p   = kal_m*(Es/(Ls*Ls))   # physical local area force strength
kT300_m = kT300 / Es           # model thermal energy at T=300k
shear_p = shear_m /Ts


# Set interparticle interaction
###############################################################################
attractionType = 5 # 5 is for WCA, and 4 is for Morse
if attractionType == 5: Dmorse = 0
else: Dmorse = 2.1 # set the depth of the Morse potential

epsLJ    = 0
eqLJ     = 0
cutoffLJ = 0

eqWCA = 2**(1/6)
epsWCA = aggregation_intensity*kT300_m # set the strength of the WCA potential

epsMorse    = Dmorse*kT300_m # set the depth of the Morse potential
widthMorse  = 0.75           # set the width of the Morse potential
eqMorse     = 0.8            # set the zero force distance
cutoffMorse = 2.2            # set the interaction cutoff distance

if attractionType == 5:        # set the cutoff distance for the neighbor list
	nlistCutoff = eqWCA+1
elif attractionType == 4:
	nlistCutoff = cutoffMorse +1


# Set the frequency of writing file
###############################################################################
strain = 1.0           
step_1strain = strain / shear_m # the # of steps for a single strain
writeProps = 1*math.floor(step_1strain)  # the frequency of writing data
writeConfig = 2*math.floor(step_1strain) # the frequency of writring particle configuration in vtk format
writeFluid_frequency = 100
writeFluid = writeFluid_frequency*math.floor(step_1strain) # the frequency of writing flow field in vtk format 

if test == 1:
    numStep  = cycle_number 
    numCycle =  1 
else:
    numStep  = 200*writeProps   # write 200 strains in each cycle
    numCycle =  int(cycle_number/200)  # the total number of timesteps is ncycle strains


# Write the file
###############################################################################
root_folder = "/userdata4/ajliu/RBC_doublet/"
components[components.index('phi')+1] = '.'.join([str(volFrac*100).split('.')[0], str(volFrac*100).split('.')[1][0]]) # Correct the volume fraction here
filename = root_folder + 'parameter/' + '_'.join(components) + '.dat'
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
f.write( 'Effective RBC diameter = %f\n' % Deff)
f.write( 'Modal RBC diameter = %f\n' % Deffm)
f.write( 'Ca(diameter) = %f\n' % (2*Ca) ) 
f.write( 'Re(diameter) = %f\n' % (4*Re_par) ) 
f.write( 'Ca(radius)   = %f\n' % Ca) 
f.write( 'Re(radius)   = %f\n\n' % Re_par) 
f.write( 'Dp             = %e\n' % Dp) 
f.write( 'Yp             = %e\n' % Yp) 
f.write( 'Gp             = %e\n' % Gp) 
f.write( 'kb_p           = %e\n' % kb_p) 
f.write( 'eta_p          = %e\n' % eta_p) 
f.write( 'kv_p           = %e\n' % kv_p) 
f.write( 'kag_p          = %e\n' % kag_p) 
f.write( 'kal_p          = %e\n' % kal_p)
f.write( 'shear_p        = %e\n' % shear_p) 
f.write( 'Vol fraction   = %f\n\n' % volFrac) 
f.write( 'Gm                      = %f\n' % Gm) 
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

# Run the simulation
f_log = open("submit.log", "a+")
#print('\n-------------------------------------------------------------------------------------------------', flush = True)
job_name = '_'.join(components)
#print("Job name:", job_name, flush = True)
f_log.write('-------------------------------------------------------------------------------------------------\n')
f_log.write("Job name: {}\n".format(job_name))

# Check whether the job is already done
if os.path.isfile("/userdata4/ajliu/RBC_doublet/Data/" + job_name + "/data/u.{}.vtk".format(int(writeFluid*(cycle_number/writeFluid_frequency - 1)))):
	#print("This job is already done!", flush = True)
	f_log.write("This job is already done!\n")
	f_log.flush()
	f_log.close()
	exit() 

# Check whether there is a run with shorter ncycle
have_shorter = False
for fn in os.listdir(root_folder+'Data'):
	try:
		if (fn.split('_')[0:7] == components[0:7]) and (fn.split('_')[-2:] == components[-2:])and (int(fn.split('_')[9]) < cycle_number):
			folder_to_rename = fn
			have_shorter = True
			#print("Rename the folder {} to {}".format(folder_to_rename, job_name), flush = True)
			f_log.write("Rename the folder {} to {}\n".format(folder_to_rename, job_name))
			subprocess.call(['mv',root_folder+'Data/'+folder_to_rename, root_folder+'Data/'+job_name])
	except:
		pass

# Roate the init.config if it is necessary
if ('angle' in components) and (not have_shorter):
	try:
		for fn in os.listdir(root_folder+'Data'):
			if fn.split('_')[0] == 'test' and fn.split('_')[1:3] == components[0:2]: # as long as the volume fraction is the same, the init config would be the same
				folder_particle = fn
		subprocess.call(['cp',root_folder+'Data/'+folder_particle+'/init/init_config.dat', root_folder+'src/init/init_config.dat'])
		#print("Complete copying the config of particles from {} to src".format(folder_particle), flush = True)
		f_log.write("Complete copying the config of particles from {} to src\n".format(folder_particle))
		subprocess.call(['python3', root_folder+'Script/Rotate.py', '_'.join(components)])
		#print("Rotate the init_config.dat in src for {} degrees".format(components[components.index('angle')+1]), flush = True)
		f_log.write("Rotate the init_config.dat in src for {} degrees\n".format(components[components.index('angle')+1]))
	except:
		#print("Couldn't find the folder to copy!", flush = True)
		f_log.write("Couldn't find the folder to copy!\n")
		f_log.flush()
		f_log.close()
		quit()

# send the job
#print('Start running the simulation "{}"'.format(job_name), flush = True)
f_log.write('Start running the simulation "{}"\n'.format(job_name))
f_log.flush()
f_log.close()
subprocess.call(['sh', root_folder+'job.sh', job_name])