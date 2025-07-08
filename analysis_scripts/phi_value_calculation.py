# calculate phi value

import gsd.hoomd
import numpy as np
import math

# read pentagon information
file_name = '5_gon.txt'
f = open(file_name,'r')
message = f.readlines()
vertices_Aup = []
for line in message:

  tmp_line = line
  x_tmp,y_tmp,z_tmp, insphere_A,sigma_prism_A,v_tmp_A,Ixx_A,Iyy_A,Izz_A,Rratio_A = tmp_line.split()
  pts_tmp = [float(x_tmp),float(y_tmp),float(z_tmp)]  

  vertices_Aup.append(pts_tmp)  

f.close()
insphere_A = float(insphere_A)
thickness_5gon = insphere_A*2.0
# define the cut-off distance for the cluster defination
cut_off = 21.24 * 2.0 * 1.5   # 21.24: pentaon thickness         

######################################

# file name, number of pentagons per chain, number of side chain segments at each side per chain
filename = "test0_one_two_three_four"
N_5gon_perchain = 4        
N_peptide_eachside = 1
# path for the targeted gsd file
path_temp = '/scratch4/tvo12/xzhan357/project_7_conjugated_polymer/4mer_pep_len' + '/n=' + str(N_peptide_eachside) + '/'

# read the last snapshot from the gsd file
file_temp = path_temp + filename + '.gsd'
f = gsd.hoomd.open(name = file_temp, mode='rb')
frame_temp = f[len(f) - 1]
N_total = frame_temp.particles.N

# total number of particles in each chain
N_length = N_5gon_perchain * 4 + N_peptide_eachside * 2      
# total number of chains 
N_chains = int(N_total/N_length)

total_5gon = N_5gon_perchain*N_chains

N_frame = len(f)        # number of frames in the gsd file
N_file = 50             # totally average over "N_files"
N_skip = 5              # skip every "N_skip" files
N_start = N_frame - 1 - (N_file - 1)*  N_skip      #  the start file number
N_end = N_start + N_file * N_skip               # the end file number

# define the function to find whether an element is in the list
def is_element_in_nd_list(nd_list, target):
    if isinstance(nd_list, list):
        for element in nd_list:
            if is_element_in_nd_list(element, target):  # if "target" in the "element" list
                return True
    else:
        return nd_list == target
    return False
# define the function to find all pentagon clusters in the system
def find_association(ref_5gon, last_run, posi_x, posi_y, posi_z, cut_off, Lx, Ly, Lz):  
    # input a ref_5gon list
    last_temp = ref_5gon
    for ref in ref_5gon: 
        # if pentagon is already in the cluster, skip to the next pentagon
        if is_element_in_nd_list(last_run, ref):
            continue
        x_ref = posi_x[ref]
        y_ref = posi_y[ref]
        z_ref = posi_z[ref]      
        for i in range(len(posi_x)):
            if is_element_in_nd_list(ref_5gon, i):
                continue
            x_i = posi_x[i]
            y_i = posi_y[i]
            z_i = posi_z[i]
 # actual distance calculation for PBC            
            drx = x_i-x_ref
            dry = y_i-y_ref
            drz = z_i-z_ref
            drx = drx - Lx*round(drx/Lx, 0)
            dry = dry - Ly*round(dry/Ly, 0)
            drz = drz - Lz*round(drz/Lz, 0)
            dr = np.sqrt(drx**2 + dry**2 + drz**2)
            if (dr < cut_off):
                ref_5gon.append(int(i))
    last_run = last_temp
# if there are no updates for the cluster list, exit    
    if len(last_run) == len(ref_5gon):
        return ref_5gon
    else:
        return find_association(ref_5gon, last_run, posi_x, posi_y, posi_z, cut_off, Lx, Ly, Lz)

counting = 0
distribution = np.zeros(total_5gon + 1)
overall = np.zeros(N_file)
for i in range (N_start,N_end, N_skip):
    counting = counting + 1
    frame_temp = f[i]
    print(frame_temp.configuration.step)
    Lx = frame_temp.configuration.box[0]
    Ly = frame_temp.configuration.box[1]
    Lz = frame_temp.configuration.box[2]   
    posi = frame_temp.particles.position
    typeid = frame_temp.particles.typeid
    x_5gon = []
    y_5gon = []
    z_5gon = []
# read pentagon coordinates    
    count = 0
    for j in range(N_chains):
        for k in range(N_length):
            if N_5gon_perchain == 1:
                if (typeid[count] == 0):
                    x_5gon.append(float(posi[count][0]))
                    y_5gon.append(float(posi[count][1]))
                    z_5gon.append(float(posi[count][2]))                    
            else:    
                if (typeid[count] == 0 or typeid[count] == 1):
                    x_5gon.append(float(posi[count][0]))
                    y_5gon.append(float(posi[count][1]))
                    z_5gon.append(float(posi[count][2]))   
            count = count + 1
    
    fivegon_cluster = []
    exit_flag = -1
    # running until all of pentagons are scanned.
    while exit_flag < 0:
        for i in range(len(x_5gon)):
            if is_element_in_nd_list(fivegon_cluster, i):
                continue
            temp = []
            temp.append(int(i))
            last_run = []
            temp1 = find_association(temp, last_run, x_5gon, y_5gon, z_5gon, cut_off, Lx, Ly, Lz)
            fivegon_cluster.append(temp1)
        if (sum(len(sublist) for sublist in fivegon_cluster) == len(x_5gon)):
            exit_flag = 1
            
    temp = 0  
# calculate the distribution based on the phi equation.    
    for sublist in fivegon_cluster:
        distribution[len(sublist)] = distribution[len(sublist)] + 1
        temp = temp + (len(sublist)/N_chains)**2
    temp = temp/N_5gon_perchain    
    overall[counting-1] = temp

f.close()
print(counting)
distribution = distribution/counting
total = sum(distribution)
file_temp = path_temp + filename + '_2nd_5gon_dis_a.txt'
output_file = open(file_temp, 'w') 
# output final alpha value and the std.
output_file.write(f"{np.mean(overall)} {np.std(overall)}\n")
# output the distribution
for i in range(len(distribution)):
    output_file.write(f"{i} {distribution[i]} {distribution[i]/total} {1.0/N_5gon_perchain*(i/N_chains)**2*distribution[i]}\n")
output_file.close()    


