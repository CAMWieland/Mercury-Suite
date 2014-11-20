import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
import warnings

warnings.simplefilter("ignore", UserWarning) #please shut up about the Agg renderer.

# mu0 = 22032
# AU3perkm3 = 3.348E24
# year2persec2 = 9.945E14

# mu = mu0 * AU3perkm3 / year2persec2
mu = 2.959122082855911E-4 # = k2

with open("MERCURY.aei") as f:
	data = []
	for line in f: #reads each line in the file into Python as a string, and splits that string at whitespaces
		singleline = []
		foo = line.split()
		for element in foo: 
			try: singleline.append(np.float64(element))
			except ValueError: singleline.append(element) #where it can, turns values into floats, passes where not (actually writes to a new array bc the original output is read-only)
		data.append(singleline)

data = data[4:]

formatted_data = []
time = []
for line in data:
	time.append(line[0])
	formatted_data.append([line[0], [line[1],line[2],line[3]], [line[4],line[5],line[6]], line[7], line[8], line[9] ]) #[Time , r vec., v vec., e, M, a]

ecc_vec = []
ecc_vec2 = []

for line in formatted_data:
	crossed = (np.cross(line[2], np.cross(line[1],line[2]))/mu) - (line[1]/np.linalg.norm(line[1]))
	ecc_vec.append(crossed)#/np.linalg.norm(crossed))

	# longform = ((np.linalg.norm(line[2])**2)*np.asarray(line[1])/mu) - (np.dot(line[1], line[2])*np.asarray(line[2])/mu) - np.asarray(line[1])/np.linalg.norm(line[1])
	# ecc_vec2.append(longform/np.linalg.norm(longform))


ecc_vec = np.transpose(ecc_vec)
# ecc_vec2 = np.transpose(ecc_vec2)

phi = np.arctan(ecc_vec[1]/ecc_vec[0])
theta = np.arccos(ecc_vec[2]/np.sqrt(ecc_vec[0]**2 + ecc_vec[1]**2 + ecc_vec[2]**2))

dphi = phi-phi[0]
dtheta = theta-theta[0]

moving_phi = []
moving_theta = []

for j in np.arange(len(time)-1000)+500:
	moving_theta.append(np.mean(dtheta[j-500:j+500]))
	moving_phi.append(np.mean(dphi[j-500:j+500]))


print len(time[500:-500])
print len(moving_phi)
# print formatted_data[:5]
# print '\n'
# print np.transpose(ecc_vec)[:5]
##########################################################################################

fig = plt.figure()
ax1 = fig.add_subplot(311)

ax1.set_title('Eccentricity Vector')

# ax1.plot(time, np.arccos(ecc_vec[2]/np.sqrt(ecc_vec[0]**2 + ecc_vec[1]**2 + ecc_vec[2]**2)), '-', color='r', label='e_theta', alpha=0.3)
ax1.plot(time, dphi, '-', color='b', label=r'$\Delta \left|{e}_{\phi}\right|$', alpha=0.3)
ax1.plot(time[500:-500], moving_phi, '-', color='r', label=r'moving average of $\phi$', alpha=0.3)
ax1.plot(time, dtheta, '-', color='g', label=r'$\Delta \left|{e}_{\theta}\right|$', alpha=0.3)
ax1.plot(time[500:-500], moving_phi, '-', color='k', label=r'moving average of $\theta$', alpha=0.3)
# ax1.plot(time, np.linalg.norm(ecc_vec, axis=0), '-', color='green', label='|e|', alpha=0.3)
# ax1.plot(time, np.arccos(ecc_vec2[2]/np.sqrt(ecc_vec2[0]**2 + ecc_vec2[1]**2 + ecc_vec2[2]**2)), '-', color='b', label='e_theta 2', alpha=0.3)
# ax1.plot(time, np.arctan(ecc_vec2[1]/ecc_vec2[0]), '-', color='y', label='e_phi 2', alpha=0.3)

ax1.set_xlabel('Time [years]', fontsize='small')
ax1.set_ylabel('Eccentricity components', fontsize='small')
# ax1.set_xscale('log')
# ax1.set_xlim([0.002, 2095])

ax1.grid()

# ax1.legend()

ax3 = fig.add_subplot(312)
ax3.set_xlim([0, 5])
ax3.set_title('Eccentricity Vector (Zoomed In)')

# ax3.plot(time, np.arccos(ecc_vec[2]/np.sqrt(ecc_vec[0]**2 + ecc_vec[1]**2 + ecc_vec[2]**2)), '-', color='r', label='e_theta', alpha=0.3)
ax3.plot(time, dphi, '-', color='b', label=r'$\Delta \left|{e}_{\phi}\right|$', alpha=0.3)
ax3.plot(time[500:-500], moving_phi, '-', color='r', label=r'moving average of $\phi$', alpha=0.3)
ax3.plot(time, dtheta, '-', color='g', label=r'$\Delta \left|{e}_{\theta}\right|$', alpha=0.3)
ax3.plot(time[500:-500], moving_phi, '-', color='k', label=r'moving average of $\theta$', alpha=0.3)
# ax3.plot(time, np.linalg.norm(ecc_vec, axis=0), '-', color='green', label='|e|', alpha=0.3)

# ax3.plot(time, np.arccos(ecc_vec2[2]/np.sqrt(ecc_vec2[0]**2 + ecc_vec2[1]**2 + ecc_vec2[2]**2)), '-', color='b', label='e_theta 2', alpha=0.3)
# ax3.plot(time, np.arctan(ecc_vec2[1]/ecc_vec2[0]), '-', color='y', label='e_phi 2', alpha=0.3)

ax3.set_xlabel('Time [years]', fontsize='small')
ax3.set_ylabel('Eccentricity components', fontsize='small')

ax3.grid()

ax3.legend()
##########################################################################################

ax2 = fig.add_subplot(313, projection='3d')
ax2.set_aspect("equal")
# ax2.dist = 5.

testu=np.linspace(0.0, 2.0*np.pi, 100)
testv=np.linspace(0.0, np.pi, 100)
# x, y, and z are the coordinates of the points for plotting
# each is arranged in a 100x100 array
testx=np.outer(np.cos(testu),np.sin(testv))
testy=np.outer(np.sin(testu),np.sin(testv))
testz=np.outer(np.ones(np.size(testu)),np.cos(testv))

ax2.plot_wireframe(testx,testy,testz,alpha=0.1)

ax2.plot(ecc_vec[0]/np.linalg.norm(ecc_vec, axis=0), ecc_vec[1]/np.linalg.norm(ecc_vec, axis=0), ecc_vec[2]/np.linalg.norm(ecc_vec, axis=0), '-', c='red')

ax2.axis([-1., 1., -1., 1.])
ax2.set_zlim3d(-1., 1.)

ax2.set_xlabel('X/|e|')
ax2.set_ylabel('Y/|e|')
ax2.set_zlabel('Z/|e|')

##########################################################################################

plt.tight_layout()

plt.savefig('images/ecc_vec.png')
























