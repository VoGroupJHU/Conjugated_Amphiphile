The simulation script contains the following txt files:
Python script:
1. conjugated_polymer_new_run.py: the simulation script to initialize and start a new simulation.
2. continue_run_write_new_file.py: the script to continue the simulation.


3. 5_gon.txt: Contains information for type A pentagons, including vertex positions, insphere diameter, circumsphere diameter, and moments of inertia.
4. 5_gon_inverse.txt: Stores similar information for type B pentagons, including vertex positions, insphere diameter, circumsphere diameter, and moments of inertia.
5. angle_init.txt: Specifies the initial angular configuration for the simulation.
6. bonding_sphere.txt: Provides details for bonding particle parameters, including diameter, volume, and moments of inertia.
7. bump_sphere.txt: Contains information on bump spheres, including their diameter, volume, and moments of inertia.
8. dihedral_init.txt: Includes data for the dihedral configuration in the system.
9. init.txt: File that defines the initial conditions for the simulation, including: Positions of particles, Rigid body configurations, Initial quaternion orientations,
	Volume and Box size.
10. init_bond.txt: Defines the bond types and the specific particles that are bonded together.
11. monomer_sphere.txt: Contains parameters for side block monomers, such as particle diameter, volume, and moments of inertia.
12. rigid_ghost_down.txt: Includes positions of bonding particles for type B pentagons.
13. rigid_ghost_up.txt: Includes positions of bonding particles for type A pentagons.
