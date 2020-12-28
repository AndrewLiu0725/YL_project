import imageio
import os

files = []
#phi = 3.8991
phis = ['0.89978', '1.4996', '1.9495', '2.3994', '2.9993', '3', '3.4492', '3.8991', '4.4989', '4.5', '4.9488', '5', '5.3987', '5.9986', '6.4484']
Ca = 0.12
path = "./Pictures/Suspension/DoubletFraction_Histogram/" 
for phi in phis:
    filename = "h24phi{}Re0.1Ca{}WCA1zero0.8_histogram_wholetimeseries.png".format(phi, Ca)
    if os.path.isfile(path+filename):
        files.append(path+filename)

images = [imageio.imread(file) for file in files]
#imageio.mimwrite('./Pictures/phi_{}.gif'.format(phi), images, fps = 2)
imageio.mimwrite('./Pictures/Ca_{}.gif'.format(Ca), images, fps = 2)