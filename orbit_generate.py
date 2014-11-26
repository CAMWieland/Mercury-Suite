# -*- coding: utf-8 -*-
"""Bahcall–Wolf distribution
n(r)~r^-7/4
	in the center of the galaxy, maybe closer to ^-3/2 or even -1/2 for red giants

	code something like r^-a, pick a

	but that's 3d distribution; integrate twice for 1D density profile

	dn/dr ~ r^2-a

	e~0.97"""

#code checks: 
# mean anomaly goes [0,2pi]         !!!!!
# j vectors evenly distributed      !!!!!
# a histogram matches bacall wolf   !!!!!
# plot log rho vs log r             
# visualize orbits


import numpy as np
import scipy
import scipy.optimize
from scipy.optimize import curve_fit
from scipy import misc
import matplotlib as mpl
import matplotlib.pyplot as plt
import codecs, argparse, warnings, os, commands, sys
import pylab as py
import mpl_toolkits.mplot3d.axes3d as p3

from mayavi import mlab

warnings.simplefilter("ignore", RuntimeWarning) #this is not an efficient program. shush.

parser = argparse.ArgumentParser(description='Generates stars in a cluster.') #adds command line arguments - c.f. https://docs.python.org/2/howto/argparse.html
parser.add_argument("-t", "--test", action="store_true", help="turns on testing mode")
parser.add_argument("-v", "--verbose", action="store_true", help="turns on verbosity")
parser.add_argument("-m", "--mayavi", action="store_true", help="does the mayavi render")
args = parser.parse_args()

numstars = int(raw_input('Number of disrupted binary stars to generate: '))
numdisk = int(raw_input('Number of disk stars to generate: '))

# while True:
# 	foo = raw_input('Include intermediate-mass BH? [y/n]')
# 	if (foo.lower())[0] == 'n': 
# 		imbh = False
# 		break
# 	elif (foo.lower())[0] == 'y':
# 		imbh = True
# 		break
# 	else:
# 		print "Not a valid answer." 
# 		continue


if os.path.isfile('bigS.in'): os.remove('bigS.in')


Msun                     = 1.98892e33
Yr                       = 3.156926e7
Day                      = 86400
G_Newton                 = 6.67384e-8
pc                       = 3.08568025e18
AU                       = 1.49597871e13
AU_per_day               = AU/86400
km                       = 1e5
km_per_sec_to_AU_per_day = 1e5/AU_per_day
Mbh                      = 4.3e6*Msun
days_per_year            = 365.
r_influence_MW           = 3*pc #http://en.wikipedia.org/wiki/Sphere_of_influence_(astronomy)
R_earth                  = 6.3674447e8
M_earth                  = 5.9721986e27
M_sun_AUsq_per_day       = 5.150489e54 #erg seconds [angular momentum]

class mbody:
	def __init__(self, name="Unnamed Object", mass=0.0, rvec=[0.0, 0.0, 0.0], vvec=[0.0, 0.0, 0.0], euler=[0.0, 0.0, 0.0], a=0.0,  b=0.0, e=0.0, mean_anomaly=0.0, eccentric_anomaly=0.0, true_anomaly=0.0): 
		self.name              = name
		self.mass              = mass
		self.rvec              = rvec
		self.vvec              = vvec
		self.euler             = euler #a vector in angle space containing the three euler angles	
		self.a                 = a
		self.b                 = b
		self.e                 = e
		self.mean_anomaly      = mean_anomaly
		self.eccentric_anomaly = eccentric_anomaly
		self.true_anomaly      = true_anomaly

def write_star_data(star_name, r_vec, v_vec, mass_in_M_Sun=0., close_encounter_radius=5., density_in_d_h20=1., spin_J_vec=[0.,0.,0.]):
	starname = star_name.upper().replace(" ", "_").replace("-", "") #makes the name all caps, replaces spaces with underscores, and removes dashes.
	mstring  = 'm={0:.18E}'.format(mass_in_M_Sun)

	if 0.1 <= close_encounter_radius < 100: rstring  = 'r={0:.3f}d3'.format(close_encounter_radius)
	else: rstring  = 'r={0:.3E}'.format(close_encounter_radius)

	if 0.1 <= density_in_d_h20 < 100: dstring  = 'd={0:.3f}'.format(density_in_d_h20)
	else: dstring  = 'd={0:.3E}'.format(density_in_d_h20)
	
	data1 = [starname, mstring, rstring, dstring]

	infile.write(' ' + '{0[0]:<15}{0[1]:<28}{0[2]:<12}{0[3]:<12}\n'.format(data1)) #cf http://ebeab.com/2012/10/10/python-string-format/ and https://docs.python.org/2/library/string.html#formatstrings
	infile.write(' ' + '{0[0]: .18E} {0[1]: .18E} {0[2]: .18E}\n'.format(r_vec))
	infile.write(' ' + '{0[0]: .18E} {0[1]: .18E} {0[2]: .18E}\n'.format(v_vec))
	if spin_J_vec == [0.,0.,0.]:
		infile.write('  0. 0. 0.\n')
	else:
		infile.write(' ' + '{0[0]: .18E} {0[1]: .18E} {0[2]: .18E}\n'.format(spin_J_vec))

number_of_orbits = numstars #for historical reasons, this refers to the number of disrupted binaries created. numdisk is the number of stars in the disk.
number_of_points = 10*number_of_orbits/0.0627553 #10x the fraction of points that would ideally fall in the r**-7/4 distribution (to be safe)
v = []

while len(v) != number_of_orbits:
		x   = np.linspace(0.01/r_influence_MW,1/r_influence_MW,1000000)
		f   = lambda x: x**(-7/4)
		fx  = f(x)
		u1  = np.random.rand(number_of_points)*(1-0.05)+0.05 # uniform random samples scaled to the interval
		u2  = np.random.rand(number_of_points)*(0.05**(-7/4))
		idx = np.where(u2<=f(u1))[0]
		v   = u1[idx][:number_of_orbits] #takes only the first [# of orbits] of accepted a's; if it doesn't get enough points, it'll try again, refilling the distribution with 10x as many tries
		number_of_points *= 10

list_of_as = v*r_influence_MW
list_of_es = 1-np.random.beta(72/25,2328/25,number_of_orbits) #gives a mean of .97 and a mode of .98 #tested this with 'for x in 1-np.random.beta(72/25,2328/25,100000000):if x>=1: print x' none were printed, yet mercury is giving me e's>1 and I don't know why
# list_of_es = np.random.beta(72/25,2328/25,number_of_orbits)

disk_as = np.random.rand(numdisk)*5.0*pc
disk_es = np.random.rand(numdisk)

if args.test == True:
	list_of_stars = list() #used later
	listofthetas = list()  #    "
	fig = plt.figure()
	n, bins, patches = plt.hist(list_of_as/r_influence_MW, 10*int(np.log(numstars)), histtype='stepfilled')
	plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
	popt, pcov = curve_fit(lambda x, a: a*(x**(-7/4)), bins[:-1], n)
	y = popt[0]*(bins**(-7/4))
	l = plt.plot(bins, y, 'k--', linewidth=1.5)
	plt.savefig('bacal-wolf_test.png')
	plt.clf()

	fige = plt.figure()
	n, bins, patches = plt.hist(list_of_es, 10*int(np.log(numstars)), histtype='stepfilled')
	plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
	plt.savefig('eccentricity_test.png')
	plt.clf()

def eulerrotation(a, b, c, vec): 
	rotated_vec = np.dot(np.asarray([[np.cos(c), -np.sin(c), 0],[np.sin(c), np.cos(c), 0],[0, 0, 1]]), np.dot(np.asarray([[np.cos(b), 0, np.sin(b)],\
		[0, 1, 0],[-np.sin(b), 0, np.cos(0)]]), np.dot(np.asarray([[1, 0, 0],[0, np.cos(a), -np.sin(a)],[0, np.sin(a), np.cos(a)]]),vec)))
	return rotated_vec

def rotation(theta, phi, vec):
	vec = np.asarray(vec)
	aroundx = np.asarray([[1.0, 0.0, 0.0],[0.0, np.cos(theta), -np.sin(theta)],[0.0, np.sin(theta), np.cos(theta)]])
	aroundy = np.asarray([[np.cos(phi), 0.0, np.sin(phi)],[0.0, 1.0, 0.0],[-np.sin(phi), 0.0, np.cos(phi)]])
	aroundz = np.asarray([[np.cos(theta), -np.sin(theta), 0.0],[np.sin(theta), np.cos(theta), 0.0],[0.0, 0.0, 1.0]])
	rotated_vec = np.dot(aroundz, np.dot(aroundy,vec))
	return rotated_vec

def f(ecc, mean_anomaly):
	return ecc - e*np.sin(ecc) - mean_anomaly

def inv_kepler_eq(mean_anomaly, e, a, b):
	ecc_anomaly = scipy.optimize.root(f, np.pi, mean_anomaly).x
	if args.verbose == True: print '\nb = ', b, '\na = ', a, '\ne = ', e, '\necc anomaly = ', ecc_anomaly[0], '\nmean anomaly = ', mean_anomaly, ' =?= ', (ecc_anomaly - e*np.sin(ecc_anomaly))[0]
	new_x = a*(np.cos(ecc_anomaly)-e) 
	new_y = b*np.sin(ecc_anomaly)
	new_r = np.sqrt(new_x**2 + new_y**2)

	return new_x, new_y, new_r, ecc_anomaly[0]

with codecs.open('bigS.in', "a", encoding="utf-8") as infile:
	infile.write(""")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)
) Lines beginning with `)' are ignored.
)---------------------------------------------------------------------
 style (Cartesian, Asteroidal, Cometary) = Cartesian
 epoch (in days) = 0.0
)---------------------------------------------------------------------\n""")

	def rotate_and_save(a, e, orbit, starprefix='S', phi_range=(0.0, 2.*np.pi), theta_range=(0.0, 2*np.pi), retrograde_possible=True): #'orbit' == list index in list of as and es, for historical reasons
		# t = np.random.rand()*1000000*Yr
		phi_min, phi_max = phi_range[0], phi_range[1] 
		theta_min, theta_max = theta_range[0], theta_range[1] 

		m_obj = np.random.normal(20, 3)*Msun #unless you have a better way to generate them
		b = a*np.sqrt(1.-e**2)
		
		mean_anomaly = 2.*np.pi*np.random.rand()

		mu = G_Newton*(Mbh+m_obj) 
		p = a*(1. - e**2)
		specific_rel_ang_mom = np.sqrt(mu*p)

		x, y, r, ecc = inv_kepler_eq(mean_anomaly, e, a, b) #computes the eccentric anomaly from the mean anomaly
		
		theta = 2.*np.arctan( np.sqrt((1.+e)/(1.-e))*np.tan(ecc/2.) ) #radians:   http://en.wikipedia.org/wiki/Eccentric_anomaly#From_the_true_anomaly
		
		while True:
			if theta < 0.: 
				theta = 2.*np.pi - theta
				continue
			if theta > 2.*np.pi: 
				theta = theta - 2.*np.pi
				continue
			else: break

		if retrograde_possible == True:
			random_bool = bool(np.random.randint(2)) #orbits can be retrograde! (does this work?)
			if random_bool == True: theta = np.pi-theta #randomly switches the orbit to the other side of the sphere, but with the cartesian velocity vectors pointing the same way
			else: pass
		else: pass

		if args.verbose == True: 
			print 'generated eccentricity = ', e
			print 'theta = ', theta
			print 'cos theta = ', np.cos(theta)
			print 'x = ', x
			print 'y = ', y
			
		if args.test == True: globals()["listofthetas"].append(theta)

		planar_vec = np.asarray([x,y])
		r_vec = np.append(planar_vec, 0)

		ecc_dot = specific_rel_ang_mom/(a*b*(1.-e*np.cos(ecc)))
		vx = -a*np.sin(ecc)*ecc_dot
		vy = b*np.cos(ecc)*ecc_dot
		v_vec = np.asarray([vx,vy,0.0])

		rottheta = (np.random.rand()*(theta_max-theta_min))+theta_min

		# rotphi = 1.0*np.arccos(2.0*np.random.rand()-1.0)+np.pi/2 #phi ranges from 0 to pi
		rotphi = 9999999
		while True:
			if phi_min <= rotphi <= phi_max: break
			else:
				rotphi = 1.0*np.arccos(2.0*np.random.rand()-1.0)+np.pi/2.
				continue

		print "rotphi = ", rotphi
		print "rottheta = ", rottheta

		r_vec = rotation(rottheta, rotphi, r_vec)/AU
		v_vec = rotation(rottheta, rotphi, v_vec)/AU_per_day #and rotate orbit plane to new orientation

		print "r_vec = ", r_vec
		print "v_vec = ", v_vec

		star_name = starprefix + str(orbit+1)
		#disrupted binaries

		if args.test == True:
			globals()[star_name]                               = mbody(name=star_name)
			globals()[star_name].__dict__['mass']              = m_obj
			globals()[star_name].__dict__['rvec']              = r_vec
			globals()[star_name].__dict__['vvec']              = v_vec
			globals()[star_name].__dict__['rot_angles']        = [rottheta, rotphi]
			globals()[star_name].__dict__['a']                 = a
			globals()[star_name].__dict__['b']                 = b
			globals()[star_name].__dict__['e']                 = e
			globals()[star_name].__dict__['mean_anomaly']      = mean_anomaly
			globals()[star_name].__dict__['eccentric_anomaly'] = ecc
			globals()[star_name].__dict__['true_anomaly']      = theta

			globals()["list_of_stars"].append(star_name)

		write_star_data(star_name, r_vec, v_vec, mass_in_M_Sun=m_obj/Msun)

	if int(number_of_orbits) != 0:
		print "number of binaries = ", numstars
		for orbit in np.arange(number_of_orbits):
			a = list_of_as[orbit]
			e = list_of_es[orbit]
			rotate_and_save(a, e, orbit, starprefix='DB', phi_range=(0.0, 2*np.pi), theta_range=(0.0, 2*np.pi), retrograde_possible=True)
	else: pass

	if int(numdisk) != 0:
		print "number in disk = ", numdisk
		for orbit in np.arange(numdisk):
			a = disk_as[orbit]
			e = disk_es[orbit]
			rotate_and_save(a, e, orbit, starprefix='BG', phi_range=(np.pi*(1.-.005), np.pi*(1.+.005)), theta_range=(0.0, 2*np.pi), retrograde_possible=False)
	else: pass

	infile.write('\n')


if args.test == True:

	from matplotlib.patches import FancyArrowPatch
	from mpl_toolkits.mplot3d import proj3d

	print '"e"s = ', list_of_es
	print '"a"s = ', list_of_as
	print 'rvec = '
	for astyr in list_of_stars: print globals()[astyr].rvec*AU

	fig = plt.figure()
	n, bins, patches = plt.hist(listofthetas, 10*int(np.log(numstars)), histtype='stepfilled')
	plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
	plt.savefig('cos_dist_test.png')
	plt.clf()

	# u and v are parametric variables.
	testu=np.linspace(0.0, 2.0*np.pi, 100)
	testv=np.linspace(0.0, np.pi, 100)
	# x, y, and z are the coordinates of the points for plotting
	# each is arranged in a 100x100 array
	testx=np.outer(np.cos(testu),np.sin(testv))
	testy=np.outer(np.sin(testu),np.sin(testv))
	testz=np.outer(np.ones(np.size(testu)),np.cos(testv))

	fig2=plt.figure(figsize=(20.48, 15.36), dpi=100)
	ax2 = p3.Axes3D(fig2)
	ax2.plot_wireframe(testx,testy,testz,alpha=0.1)
	ax2.set_xlabel('X')
	ax2.set_ylabel('Y')
	ax2.set_zlabel('Z')

	test_vecs = []
	for astyr in list_of_stars:
		[rottheta, rotphi] = globals()[astyr].rot_angles
		rotated_vec = rotation(rottheta, rotphi, np.asarray([1,0,0]))
		test_vecs.append(rotated_vec)
	
	test_vecs = np.transpose(test_vecs)
	ax2.scatter(test_vecs[0][:numstars-1], test_vecs[1][:numstars-1], test_vecs[2][:numstars-1], c='OrangeRed', s=6, alpha=0.5)
	ax2.scatter(test_vecs[0][-numdisk:], test_vecs[1][-numdisk:], test_vecs[2][-numdisk:], c='black', s=6, alpha=0.5)
	ax2.set_title('ArcCos Covering Test')

	for angle in np.arange(360):
		ax2.azim = 0
		ax2.elev = angle
		savename = 'arccos-covering_test'+str(angle).zfill(3)+'.png'
		plt.savefig(savename)
		print '\r[{0: <10}] {1: <7.2%}: Rendering frame {2:} of 720'.format('#'*int(float(angle)/(720./10)), float(angle)/(720.), angle),
		sys.stdout.flush()
	for angle in np.arange(360):
		ax2.azim = angle
		ax2.elem = 0
		savename = 'arccos-covering_test'+str(angle+360).zfill(3)+'.png'
		plt.savefig(savename)
		print '\r[{0: <10}] {1: <7.2%}: Rendering frame {2:} of 720'.format('#'*int(float(angle+360)/(720./10)), float(angle+360)/(720.), angle+360),
		sys.stdout.flush()


	if commands.getstatusoutput("man ffmpeg")[0] == 0:
		syscommand = 'ffmpeg -f image2 -r 24 -i arccos-covering_test%03d.png -vcodec mpeg4 -y arccos-covering_test.mp4'	#this condition will be true, and you'll have the option to turn them into a movie with it. Otherwise, you're on your own.
		os.system(syscommand)
		os.system('rm arccos-covering_test*.png')

	plt.clf()

	fig2=plt.figure(figsize=(20.48, 15.36), dpi=100)
	ax2 = p3.Axes3D(fig2)
	# ax2.plot_wireframe(testx,testy,testz,alpha=0.1)
	ax2.set_xlabel('X')
	ax2.set_ylabel('Y')
	ax2.set_zlabel('Z')

	for astyr in list_of_stars:
		tarange = np.linspace(0.0,2.0*np.pi,200)
		[rottheta, rotphi] = globals()[astyr].rot_angles
		orba = globals()[astyr].a
		orbe = globals()[astyr].e
		orbx = (orba*(1.0-orbe**2)/(1.0+orbe*np.cos(tarange)))*np.cos(tarange)
		orby = (orba*(1.0-orbe**2)/(1.0+orbe*np.cos(tarange)))*np.sin(tarange)
		vecs = np.transpose([orbx,orby,np.zeros(200)])
		rotated_vecs  = [] 
		for orbvec in vecs: rotated_vecs.append(rotation(rottheta, rotphi, orbvec))
		rotated_vecs = np.transpose(rotated_vecs)
		ax2.plot(rotated_vecs[0], rotated_vecs[1], rotated_vecs[2], alpha=0.5)
		ax2.scatter(globals()[astyr].rvec[0]*AU, globals()[astyr].rvec[1]*AU, globals()[astyr].rvec[2]*AU, c='red', s=5)

		class Arrow3D(FancyArrowPatch):
		    def __init__(self, xs, ys, zs, *args, **kwargs):
		        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
		        self._verts3d = xs, ys, zs

		    def draw(self, renderer):
		        xs3d, ys3d, zs3d = self._verts3d
		        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
		        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
		        FancyArrowPatch.draw(self, renderer)

		print globals()[astyr].rvec[0]*AU, globals()[astyr].rvec[0]*AU+globals()[astyr].vvec[0]*AU_per_day*200000000*days_per_year, globals()[astyr].vvec[0]*AU_per_day*100*days_per_year
		print globals()[astyr].rvec[1]*AU, globals()[astyr].rvec[1]*AU+globals()[astyr].vvec[1]*AU_per_day*200000000*days_per_year, globals()[astyr].vvec[1]*AU_per_day*100*days_per_year
		print globals()[astyr].rvec[2]*AU, globals()[astyr].rvec[2]*AU+globals()[astyr].vvec[2]*AU_per_day*200000000*days_per_year, globals()[astyr].vvec[2]*AU_per_day*100*days_per_year
		print '\n'

		a = Arrow3D([globals()[astyr].rvec[0]*AU ,globals()[astyr].rvec[0]*AU+globals()[astyr].vvec[0]*AU_per_day*100000000*days_per_year],[globals()[astyr].rvec[1]*AU ,globals()[astyr].rvec[1]*AU+globals()[astyr].vvec[1]*AU_per_day*100000000*days_per_year],[globals()[astyr].rvec[2]*AU ,globals()[astyr].rvec[2]*AU+globals()[astyr].vvec[2]*AU_per_day*100000000*days_per_year],\
		 mutation_scale=15, lw=1, arrowstyle="-|>", color="k")
		ax2.add_artist(a)
	
	ax2.set_title('Unperturbed Orbits')

	i = -1
	for angle in np.arange(360):
		ax2.azim = 0
		ax2.elev = angle
		i += 1
		savename = 'unperturbed_orbits'+str(i).zfill(3)+'.png'
		plt.savefig(savename)
		print '\r[{0: <10}] {1: <7.2%}: Rendering frame {2:} of 720'.format('#'*int(float(i)/(720./10)), float(i)/(720.), i),
		sys.stdout.flush()
	for angle in np.arange(360):
		ax2.azim = angle
		ax2.elem = 0
		i += 1
		savename = 'unperturbed_orbits'+str(i).zfill(3)+'.png'
		plt.savefig(savename)
		print '\r[{0: <10}] {1: <7.2%}: Rendering frame {2:} of 720'.format('#'*int(float(i)/(720./10)), float(i)/(720.), i),
		sys.stdout.flush()

	if commands.getstatusoutput("man ffmpeg")[0] == 0:
		syscommand = 'ffmpeg -f image2 -r 24 -i unperturbed_orbits%03d.png -vcodec mpeg4 -y unperturbed_orbits.mp4'	#this condition will be true, and you'll have the option to turn them into a movie with it. Otherwise, you're on your own.
		os.system(syscommand)
		os.system('rm unperturbed_orbits*.png')

	plt.clf()








