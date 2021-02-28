# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 02/28/2021
# ===============================================================================
import math
import sys

root_folder = "/userdata4/ajliu/RBC_doublet/"
components = sys.argv[1].split('_')
angle = float(components[components.index('angle')+1])*(math.pi)/180.0 # convert to radian
path_parameter = root_folder+'parameter/'+'_'.join(sys.argv[1].split('_'))+'.dat'
fp = open(path_parameter)
dim = fp.readlines()[3].split()
fp.close()
center = [float(i)/2.0 for i in dim]

rotation_matrix = [[math.cos(angle), math.sin(angle), 0], [-math.sin(angle), math.cos(angle), 0], [0, 0, 1]]
def rotate(r):
	r_rotated = [0, 0, 0]
	for j in range(3):
		for k in range(3):
			r_rotated[j] += ((r[k]-center[k]) * rotation_matrix[j][k])
		r_rotated[j] += center[j]
	return r_rotated

with open(root_folder+'src/init/init_config.dat') as f:
	lines = f.readlines()

for line_index, line in enumerate(lines):
	if '.' in line:
		pos = [float(coord) for coord in (line.split())]
		pos_rotated = rotate(pos)
		lines[line_index] = ' '.join([str(i) for i in pos_rotated])+'\n'

fn = open(root_folder+'src/init/init_config.dat', 'w')
fn.writelines(lines)
fn.close()

