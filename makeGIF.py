import imageio
import os

r = 1.0
files = []
for i in range(1, 21):
    files.append("./Pictures/IntrinsicRelativeViscosity_vs_phi_unseperated/TwoCellSystem_Viscosity_vs_phi_Ca_{}_fixed_ylim_r_{}_unseperated.png".format(i*0.01, r))

images = [imageio.imread(file) for file in files]
imageio.mimwrite('./Pictures/IntrinsicRelativeViscosity_vs_phi_unseperated/TwoCellSystem_Viscosity_vs_phi_r_{}_unseperated.gif'.format(r), images, fps = 2)