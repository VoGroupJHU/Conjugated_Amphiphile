import hoomd
from hoomd import md
import datetime
import coxeter as cx
import os
import numpy as np
import random
import gsd.hoomd

# Hoomd v3.11 simulation script for conjugated polymer systems.
# particle type: A - upright pentagon (conjugated core); B - inverted pentagon (conjugated core); C - side chain segment; D - bonding particle; E - bump sphere.

path = os.path.dirname(os.path.abspath(__file__))
path = path + "/"
index_f = 0
# Read file for bond configuration
file_name = path + 'init_bond.txt'
f = open(file_name,'r')
message = f.readlines()
bond_config = []
bond_type_config = []
for line in message:

  tmp_line = line
  bond_type,bond_i,bond_j = tmp_line.split()
  bond_tmp = [int(bond_i),int(bond_j)]  

  bond_config.append(bond_tmp)
  bond_type_config.append(bond_type)

f.close()

# read file for angle configuration
file_name = path + 'angle_init.txt'
f = open(file_name, 'r')
message = f.readlines()
angle_config = []
angle_type_config = []
for line in message:

  tmp_line = line
  angle_type,angle_i,angle_j,angle_k, Theta_value = tmp_line.split()
  angle_tmp = [int(angle_i),int(angle_j),int(angle_k)]  

  angle_config.append(angle_tmp)
  angle_type_config.append(angle_type)

f.close()
Theta_value = float(Theta_value)

# read file for dihedral configuration
file_name = path + 'dihedral_init.txt'
f = open(file_name, 'r')
message = f.readlines()
dihedral_config = []
dihedral_type_config = []
for line in message:

  tmp_line = line
  dihedral_type,dihedral_i,dihedral_j,dihedral_k, dihedral_l = tmp_line.split()
  dihedral_tmp = [int(dihedral_i),int(dihedral_j),int(dihedral_k),int(dihedral_l)]  

  dihedral_config.append(dihedral_tmp)
  dihedral_type_config.append(dihedral_type)

f.close()

# read file for side chain particle diameter, volume and moment of inertia
file_name = path + 'monomer_sphere.txt'
f = open(file_name, 'r')
message = f.readlines()
for line in message:
    tmp_line = line
    peptide_number, r_sphereC, vC_sphere, Ixx_sphere, Iyy_sphere, Izz_sphere = tmp_line.split()
f.close()
peptide_number = int(peptide_number)
r_sphereC = float(r_sphereC)
sphere_C_d = r_sphereC * 2.0

vC_sphere = float(vC_sphere)
Ixx_sphere = float(Ixx_sphere)
Iyy_sphere = float(Iyy_sphere)
Izz_sphere = float(Izz_sphere)

# read file for bonding particle diameter, volume and moment of inertia
file_name = path + 'bonding_sphere.txt'
f = open(file_name, 'r')
message = f.readlines()
for line in message:
    tmp_line = line
    r_sphereD, vD_sphere, IxxD_sphere, IyyD_sphere, IzzD_sphere = tmp_line.split()
f.close()
r_sphereD = float(r_sphereD)
sphere_D_d = r_sphereD * 2.0

vD_sphere = float(vD_sphere)
IxxD_sphere = float(IxxD_sphere)
IyyD_sphere = float(IyyD_sphere)
IzzD_sphere = float(IzzD_sphere)

# read file for bump sphere diameter, volume and moment of inertia
file_name = path + 'bump_sphere.txt'
f = open(file_name, 'r')
message = f.readlines()
for line in message:
    tmp_line = line
    n_5gon, r_sphereE, vE_sphere, IxxE_sphere, IyyE_sphere, IzzE_sphere = tmp_line.split()
f.close()
n_5gon = int(n_5gon)
r_sphereE = float(r_sphereE)
sphere_E_d = r_sphereE * 2.0

vE_sphere = float(vE_sphere)
IxxE_sphere = float(IxxE_sphere)
IyyE_sphere = float(IyyE_sphere)
IzzE_sphere = float(IzzE_sphere)

# Read file for A pentagon information, including vertices positions, insphere diameter, circumsphere diameter, moment of inertia
file_name = path + '5_gon.txt'
f = open(file_name,'r')
message = f.readlines()
vertices_Aup = []
for line in message:

  tmp_line = line
  x_tmp,y_tmp,z_tmp, insphere_A,sigma_prism_A,v_tmp_A,Ixx_A,Iyy_A,Izz_A,Rratio_A = tmp_line.split()
  pts_tmp = [float(x_tmp),float(y_tmp),float(z_tmp)]  

  vertices_Aup.append(pts_tmp)  

f.close()
insphere_A_d = float(insphere_A) * 2.0
v_A = float(v_tmp_A)

insphere_A = float(insphere_A)
Rratio_A = float(Rratio_A)

# Moment of inertia
Ixx_A = float(Ixx_A)
Iyy_A = float(Iyy_A)
Izz_A = float(Izz_A)

# Read file for B pentagon information, including vertices positions, insphere diameter, circumsphere diameter, moment of inertia
file_name = path + '5_gon_inverse.txt'
f = open(file_name,'r')
message = f.readlines()
vertices_Bdown = []
for line in message:

  tmp_line = line
  x_tmp,y_tmp,z_tmp, insphere_B,sigma_prism_B,v_tmp_B,Ixx_B,Iyy_B,Izz_B,Rratio_B = tmp_line.split()
  pts_tmp = [float(x_tmp),float(y_tmp),float(z_tmp)]  

  vertices_Bdown.append(pts_tmp)  

f.close()

v_B = float(v_tmp_B)

Rratio_B = float(Rratio_B)
insphere_B = float(insphere_B)

# Moment of inertia
Ixx_B = float(Ixx_B)
Iyy_B = float(Iyy_B)
Izz_B = float(Izz_B)

##Read A pentagon bonding particle positions relative to the A pentagon
file_name = path + 'rigid_ghost_up.txt'
f = open(file_name,'r')
message = f.readlines()
rigid_up = []
rigid_up_type = []
for line in message:

  tmp_line = line
  x_tmp,y_tmp,z_tmp = tmp_line.split()
  pts_tmp = [float(x_tmp),float(y_tmp),float(z_tmp)]  

  rigid_up.append(pts_tmp)  
  rigid_up_type.append('D')
  
f.close()

##Read B pentagon bonding particle positions relative to the A pentagonwn
file_name = path + 'rigid_ghost_down.txt'
f = open(file_name,'r')
message = f.readlines()
rigid_down = []
for line in message:

  tmp_line = line
  x_tmp,y_tmp,z_tmp = tmp_line.split()
  pts_tmp = [float(x_tmp),float(y_tmp),float(z_tmp)]  

  rigid_down.append(pts_tmp)  

f.close()

# Read the file for all particles positions, rigid body information, initial quaternion, volume, and box size. 
file_name = path + 'init.txt'
f = open(file_name,'r')
message = f.readlines()
pts_config = []
q_config = []
type_config = []
body_config = []
moment_config = []
mass_config = []
diameter_config = []
spherical_P_MI = 1.0
for line in message:

  tmp_line = line
  type_tmp,body_tmp,x_tmp,y_tmp,z_tmp,qs,qx,qy,qz,v_monomer,Lx,Ly,Lz = tmp_line.split()
  pts_tmp = [float(x_tmp),float(y_tmp),float(z_tmp)]    
  q_tmp = [float(qs),float(qx),float(qy),float(qz)]

  pts_config.append(pts_tmp)
  type_config.append(int(type_tmp))
  body_config.append(int(body_tmp))
  q_config.append(q_tmp)
  moment_factor = 5E-6
  if int(type_tmp) == 1:  
    monomer_add = [moment_factor*Ixx_B,moment_factor*Iyy_B,moment_factor*Izz_B]
    moment_config.append(monomer_add)
    mass_config.append(1.0)  
    diameter_config.append(0.0)
  elif int(type_tmp) == 0:
    monomer_add = [moment_factor*Ixx_A,moment_factor*Iyy_A,moment_factor*Izz_A]
    moment_config.append(monomer_add)
    mass_config.append(1.0)
    diameter_config.append(0.0)
  elif int(type_tmp) == 2:
    moment_add = [spherical_P_MI,spherical_P_MI,spherical_P_MI]
    moment_config.append(moment_add)
    mass_config.append(1.0)  
    diameter_config.append(float(sphere_C_d))
  elif int(type_tmp) == 3:
    moment_add = [0.0,0.0,0.0]
    moment_config.append(moment_add)
    mass_config.append(0.0)      
    diameter_config.append(float(sphere_D_d))
  elif int(type_tmp) == 4:
    moment_add = [0.0,0.0,0.0]
    moment_config.append(moment_add)
    mass_config.append(0.0)      
    diameter_config.append(float(sphere_E_d))  
    if int(body_tmp) == peptide_number:
        rigid_up.append([0.0,0.0,0.0])
        rigid_up_type.append('E')
    elif int(body_tmp) == peptide_number+1:
        rigid_down.append([0.0,0.0,0.0])       
f.close()

# Monomer total volume
v_monomer = float(v_monomer)

# Box size
Lx = float(Lx) * 2.0
Ly = float(Ly) * 2.0
Lz = float(Lz) * 0.2
 
#####################################
#####################################

###################
### Set up sims ###
###################

# Box factor
box_factor = 1

# Define simulation device
cpu = hoomd.device.CPU()

# Create simulation object
sim = hoomd.Simulation(device=cpu, seed=87)

######### make snapshot ################
snapshot = gsd.hoomd.Snapshot()

# box configuration
snapshot.configuration.box = [box_factor*Lx, box_factor*Ly, box_factor*Lz,0,0,0]

# partticles information
snapshot.particles.N = len(pts_config)
snapshot.particles.position = pts_config
snapshot.particles.orientation = q_config
snapshot.particles.moment_inertia = moment_config
snapshot.particles.typeid = type_config
snapshot.particles.diameter = diameter_config
snapshot.particles.body = body_config      
snapshot.particles.types = ['A', 'B', 'C', 'D', 'E']
snapshot.particles.mass = mass_config

# define rigid for A and B pentagon. Each rigid body includes one pentagon, two bonding particles and one bump sphere.
rigid = md.constrain.Rigid()
rigid.body['A'] = {
"positions": rigid_up,
"constituent_types": rigid_up_type,   
"orientations": [(1.0, 0.0, 0.0, 0.0), (1.0, 0.0, 0.0, 0.0), (1.0, 0.0, 0.0, 0.0)],
"charges": [0.0, 0.0, 0.0],
"diameters": [sphere_D_d, sphere_D_d, sphere_E_d]
}

rigid.body['B'] = {
"positions": rigid_down,
"constituent_types": rigid_up_type,   
"orientations": [(1.0, 0.0, 0.0, 0.0), (1.0, 0.0, 0.0, 0.0), (1.0, 0.0, 0.0, 0.0)],
"charges": [0.0, 0.0, 0.0],
"diameters": [sphere_D_d, sphere_D_d, sphere_E_d]
}


##Bonds information 
snapshot.bonds.N = len(bond_config)
snapshot.bonds.group = bond_config
snapshot.bonds.typeid = bond_type_config
snapshot.bonds.types = ['sphere-sphere', 'sphere-ghost', 'ghost-ghost']

## angle information 
snapshot.angles.N = len(angle_config)
snapshot.angles.group = angle_config
snapshot.angles.typeid = angle_type_config
snapshot.angles.types = ['5gon']

## dihedral information 
snapshot.dihedrals.N = len(dihedral_config)
snapshot.dihedrals.group = dihedral_config
snapshot.dihedrals.typeid = dihedral_type_config
snapshot.dihedrals.types = ['3-9-10-4']

# Replicate the system in x, y, z dimension
n_repx = 2
n_repy = 2
n_repz = 40

communicator = hoomd.communicator.Communicator()
snapshot_1 = hoomd.Snapshot.from_gsd_snapshot(snapshot, communicator)
snapshot_1 = snapshot_1.replicate(nx=n_repx,ny=n_repy,nz=n_repz)


sim.create_state_from_snapshot(snapshot_1, domain_decomposition=(2, 2, 3))

Lz_current = sim.state.box.Lz
Lx_current = sim.state.box.Lx
Ly_current = sim.state.box.Ly

# Compute volumes
V_total = (n_repx*n_repy*n_repz)*(v_A*n_5gon + vC_sphere*peptide_number*2 + vD_sphere*n_5gon*2)

# Define final volume density
rho_target = 0.02

# Box size
Lfinal = (V_total/rho_target)**(1.0/3.0)

# Compute final packing fraction
phi_final = V_total/(Lfinal**3.0)

###############################
### Define potential params ###
###############################

# create neighbor list
nl = md.nlist.Cell(buffer=0.0, exclusions=['body'])

# ALJ parameters
kbt = 1.0
eps_rep = 1.0
eps_att = 1.0

# Define cutoffs
cutoff_rep = 2.0**(1.0/6.0)
cutoff_att = 2.5

# Cutoff buffer
buffer_cutoff = 1.25 * 3.0
buffer_cutoff_sphere = 2.0
Rratio_A = 2.5
# Cutoffs
rcut_5gon = buffer_cutoff*cutoff_rep*insphere_A*Rratio_A
rcut_sphereC = buffer_cutoff_sphere*cutoff_rep*r_sphereC
rcut_sphereD = buffer_cutoff_sphere*cutoff_rep*r_sphereD
rcut_sphereE = buffer_cutoff_sphere*cutoff_rep*r_sphereE
rcut_sphereCD = (rcut_sphereC+rcut_sphereD)
rcut_sphereCD_E = rcut_sphereC + rcut_sphereE
rcut_5gon_sphereC = (rcut_sphereC+rcut_5gon) * 0.8
rcut_5gon_sphereD = (rcut_sphereD+rcut_5gon) * 0.8
rcut_5gon_sphereE = (rcut_sphereE+rcut_5gon) * 0.8

# Shape definition in ALJ
shapeA = cx.shapes.ConvexPolyhedron(np.array(vertices_Aup))
shapeB = cx.shapes.ConvexPolyhedron(np.array(vertices_Bdown))

test_E = eps_rep

alj = md.pair.aniso.ALJ(nl, default_r_cut=0)
alj.params[(['A','B'], ['A','B'])] = dict(epsilon=eps_rep,sigma_i=insphere_A_d,sigma_j=insphere_A_d, alpha=0)     
alj.params[(['A','B'], 'C')] = dict(epsilon=eps_rep,sigma_i=insphere_A_d,sigma_j=sphere_C_d, alpha=0)  
alj.params[(['A','B'], 'D')] = dict(epsilon=eps_rep,sigma_i=insphere_A_d,sigma_j=sphere_D_d, alpha=0) 
alj.params[(['A','B'], 'E')] = dict(epsilon=test_E,sigma_i=insphere_A_d,sigma_j=sphere_E_d, alpha=0)    
alj.params[('C', 'C')] = dict(epsilon=eps_rep,sigma_i=sphere_C_d,sigma_j=sphere_C_d, alpha=0)  
alj.params[('C', 'D')] = dict(epsilon=eps_rep,sigma_i=sphere_C_d,sigma_j=sphere_D_d, alpha=0) 
alj.params[('D', 'D')] = dict(epsilon=eps_rep,sigma_i=sphere_D_d,sigma_j=sphere_D_d, alpha=0) 
alj.params[('C', 'E')] = dict(epsilon=test_E,sigma_i=sphere_C_d,sigma_j=sphere_E_d, alpha=0)   
alj.params[('D', 'E')] = dict(epsilon=test_E,sigma_i=sphere_D_d,sigma_j=sphere_E_d, alpha=0)
alj.params[('E', 'E')] = dict(epsilon=test_E,sigma_i=sphere_E_d,sigma_j=sphere_E_d, alpha=0)   
alj.shape["A"] = dict(vertices = shapeA.vertices, rounding_radii=(0.0,0.0,0.0), faces = shapeA.faces)
alj.shape["B"] = dict(vertices = shapeB.vertices, rounding_radii=(0.0,0.0,0.0), faces = shapeB.faces)
alj.shape["C"] = dict(vertices = [], rounding_radii = r_sphereC, faces = [])
alj.shape["D"] = dict(vertices = [], rounding_radii = r_sphereD, faces = []) 
alj.shape["E"] = dict(vertices = [], rounding_radii = r_sphereE, faces = []) 

alj.r_cut[(['A','B'], ['A','B'])] = rcut_5gon
alj.r_cut[(['A','B'], 'C')] =  rcut_5gon_sphereC
alj.r_cut[(['A','B'], 'D')] =   rcut_5gon_sphereD
alj.r_cut[(['A','B'], 'E')] =   rcut_5gon_sphereE
alj.r_cut[('C', 'C')] = rcut_sphereC
alj.r_cut[('C', 'D')] = rcut_sphereCD
alj.r_cut[(['C','D'], 'E')] = rcut_sphereCD_E
alj.r_cut[('D', 'D')] = rcut_sphereD
alj.r_cut[('E', 'E')] = rcut_sphereE

# dipole dipole parameter setup parameter
cut_coefficient = 8.05
dd_cutoff = insphere_A_d*10.0
dipole = md.pair.aniso.Dipole(nl, default_r_cut=Lfinal/cut_coefficient)
dipole.params[(['A','B','C','D','E'], ['A','B','C','D','E'])] = dict(A = 0.0, kappa = 0.0) 
dipole.params[(['A','B'], ['A','B'])] = dict(A=10000.0, kappa=0.0)
dipole.mu['A'] = (0.0, 0.0, 8.0)
dipole.mu['B'] = (0.0, 0.0, 8.0)
dipole.mu[('C', 'D', 'E')] = (0.0, 0.0, 0.0)
dipole.r_cut[(['A','B','C','D','E'], ['C','D','E'])] = 0.0
dipole.r_cut[(['A','B'], ['A','B'])] = Lfinal/cut_coefficient

# Define FENEWCA bond potential
kspring_gg = 30.0
kspring_nn = 30.0
kspring_mg = 30.0
ro = 1.5
ro_gg = 0.8*ro
bond_CC_sigma = 1.15*(r_sphereC + r_sphereC)
bond_CD_sigma = 1.15*(r_sphereC + r_sphereD)
bond_DD_sigma = 1.15*(r_sphereD + r_sphereD)
fenewca = md.bond.FENEWCA()
ro_init = 50
fenewca.params['sphere-sphere'] = dict(k=kspring_gg, r0=ro_init*bond_CC_sigma, sigma=bond_CC_sigma, epsilon= 0, delta = 0.0)
fenewca.params['sphere-ghost'] = dict(k=kspring_gg, r0=ro_init*bond_CD_sigma, sigma=bond_CD_sigma, epsilon= 0, delta = 0.0)
fenewca.params['ghost-ghost'] = dict(k=kspring_gg, r0=ro_init*bond_DD_sigma, sigma=bond_DD_sigma, epsilon= 0, delta = 0.0) 

# angle potential setup
angleharmonic = md.angle.Harmonic()
angleharmonic.params['5gon'] = dict(k=5000, t0=Theta_value)
# dihedral potential setup
dihedralOPLS = md.dihedral.OPLS()
dihedralOPLS.params['3-9-10-4'] = dict(k1=5000.0, k2=0.0, k3=0.0, k4=0.0)

# set up integrator
integrator = hoomd.md.Integrator(dt=7.5E-4, integrate_rotational_dof=True)
sim.operations.integrator = integrator
integrator.rigid = rigid


# Define groups
all = hoomd.filter.All()
rigid_centers_and_free = hoomd.filter.Rigid(("center", "free"))
null = hoomd.filter.Null()

# set up the integrator
nvt = hoomd.md.methods.NVT(filter=rigid_centers_and_free, kT=kbt, tau=0.1)
integrator.methods.append(nvt)
integrator.forces.append(fenewca)       
integrator.forces.append(alj)
integrator.forces.append(angleharmonic) 
integrator.forces.append(dihedralOPLS)                            
# Assign to simulation
sim.operations.integrator = integrator

# Thermalize - assigns random values of momentum, velocity, and angular momentum consistent with ensemble properties
sim.state.thermalize_particle_momenta(filter=rigid_centers_and_free, kT=kbt)

ndump = 10000
ndump_gsd = ndump

# Define the writer to GSD output
logger_shape = hoomd.logging.Logger()
logger_shape.add(alj,quantities=['type_shapes'])
gsd_writer = hoomd.write.GSD(filename=path + 'test' + str(index_f) + '.gsd',
                             trigger=hoomd.trigger.Periodic(int(ndump_gsd)),
                            mode='wb',filter = all,
                            log=logger_shape)
sim.operations.writers.append(gsd_writer)

#output initial cfg
hoomd.write.GSD.write(state=sim.state, mode='wb', filename=path + 'initial.gsd')
                 
ro_space_array_SS = np.linspace(50,1.5,int(40))            
ro_space_array_SG = np.linspace(50,1.5,int(40))             
ro_space_array_GG = np.linspace(50,1.5,int(40))            
kspring_gg = 30.0 

print(Lfinal)
# initilize FENE bond  
for i in range(0,len(ro_space_array_GG)):
    if (cpu.communicator.rank == 0):
        print('Run: ',i)

    # Update bonds
    fenewca.params['sphere-sphere'] = dict(k=kspring_gg, r0=ro_space_array_SS[i]*bond_CC_sigma, sigma=bond_CC_sigma, epsilon= 0, delta = 0.0)
    fenewca.params['sphere-ghost'] = dict(k=kspring_gg, r0=ro_space_array_SG[i]*bond_CD_sigma, sigma=bond_CD_sigma, epsilon= 0, delta = 0.0)
    fenewca.params['ghost-ghost'] = dict(k=kspring_gg, r0=ro_space_array_GG[i]*bond_DD_sigma, sigma=bond_DD_sigma, epsilon= 0, delta = 0.0) 

    sim.run(10000)            

# zero momentum, and renormalize momentum
zero_trigger = hoomd.trigger.On(sim.timestep)
hoomd.md.update.ZeroMomentum(zero_trigger)
sim.state.thermalize_particle_momenta(filter=rigid_centers_and_free, kT=kbt)

##output thermodynamic properties to logger file.
thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(
    filter=rigid_centers_and_free)
sim.operations.computes.append(thermodynamic_properties)
## Define logger for tps printing ###
logger = hoomd.logging.Logger(categories=['scalar', 'string'])
logger.add(sim, quantities=['timestep', 'tps'])
logger.add(thermodynamic_properties, quantities = ['kinetic_temperature'])
class Status():
    def __init__(self, sim):
        self.sim = sim
    @property
    def seconds_remaining(self):
        try:
            return (self.sim.final_timestep - self.sim.timestep) / self.sim.tps
        except ZeroDivisionError:
            return 0
    @property
    def etr(self):
        return str(datetime.timedelta(seconds=self.seconds_remaining))
status = Status(sim)
logger[('Status', 'etr')] = (status, 'etr', 'string')
table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(period=int(ndump)),
                          logger=logger)
sim.operations.writers.append(table)

# compress the simulation box
n_compress = 1E6                               
initial_box = sim.state.box
final_box = hoomd.Box.cube(Lfinal)
box_resize_trigger = hoomd.trigger.Periodic(5)                     
ramp = hoomd.variant.Ramp(A=0, B=1, t_start=sim.timestep, t_ramp=int(n_compress))
box_resize = hoomd.update.BoxResize(box1=initial_box,
                                    box2=final_box,
                                    variant=ramp,
                                    trigger=box_resize_trigger)
sim.operations.updaters.append(box_resize)                                    
sim.run(n_compress)

integrator.forces.append(dipole)
##### zero momentum and simulation running at the targeted volume density
zero_trigger = hoomd.trigger.On(sim.timestep)
hoomd.md.update.ZeroMomentum(zero_trigger)
sim.state.thermalize_particle_momenta(filter=rigid_centers_and_free, kT=kbt)
sim.run(2E8)

