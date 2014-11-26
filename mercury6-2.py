import matplotlib as mpl
mpl.use('Agg') #hopefully fixes a problem with ssh-ing 
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import matplotlib.pyplot as plt
import matplotlib.animation as animation
try:
   import cPickle as pickle
except:
   import pickle
import random, os, sys, glob, traceback, argparse, zlib, warnings, commands
warnings.simplefilter("ignore", RuntimeWarning) #this is not an efficient program either.

# from mayavi import mlab

termy = int(os.popen('stty size', 'r').read().split()[0]) #prints enough newlines to clear the terminal
print termy*'\n'


parser = argparse.ArgumentParser(description='Extracts and plots data from a Mercury6 run. If the user supplies no commandline arguments, all possible subroutines will be output.') #adds command line arguments - c.f. https://docs.python.org/2/howto/argparse.html
parser.add_argument("-v", "--verbose", action="store_true", help="increases verbosity")
parser.add_argument("-c", "--clo", help="plot eccentricities and (de/dt)/E_init's and mark close encounters for each object individually", action="store_true")
parser.add_argument("-e", "--ecc", help="plot all objects' eccentricities on the same plot", action="store_true")
parser.add_argument("-D", "--dedt", help="plot all objects' (de/dt)/E_init's on the same plot", action="store_true")
parser.add_argument("-m", "--moviemake", nargs='*', help="renders frames for compilation into a movie (with ffmpeg or equivalent) of the orbits of all objects together")
parser.add_argument("-o", "--orbitmovies", nargs='*', help="renders frames for compilation into a movie (with ffmpeg or equivalent) of the orbits of all objects, individually. Not specifying an object will render all objects.")
parser.add_argument("-a", "--anonymous", action="store_true", help="Will not label plotted orbits, allowing for cleaner graphing of large clusters")
parser.add_argument("-rl", "--reload", action="store_true", help="uses previously extracted data rather than recalculating")

#should add name options to --orbitmovies to allow specific orbits to be requested

if not len(sys.argv) > 1: args = parser.parse_args(['-c', '-e', '-D', '-m', '1', '-o', '1']) #if you don't supply any commands, it'll run all of them (without loading preexisting data)
elif (len(sys.argv) == 2) and ('-v' in sys.argv): args = parser.parse_args(['-v', '-c', '-e', '-D', '-m', '1', '-o', '1']) #if you get a '-v' and nothing else, run everything, verbosely
elif (len(sys.argv) == 2) and ('-rl' in sys.argv): args = parser.parse_args(['-rl', '-c', '-e', '-D', '-m', '1', '-o', '1']) #if you get a '-rl' and nothing else, run everything, from pickle
elif (len(sys.argv) == 3) and ('-v' and '-rl' in sys.argv): args = parser.parse_args(['-rl', '-v', '-c', '-e', '-D', '-m', '1', '-o', '1']) #if you get '-rl -v' yeah
elif (len(sys.argv) == 3) and ('-v' and '-a' in sys.argv): args = parser.parse_args(['-v', '-c', '-e', '-D', '-m', '1', '-o', '1', '-a']) #if you get '-v -a' yeah
else: args = parser.parse_args()
# print 'reached line 23'
############ conversion tables ############
Msun = 1.98892e33
Yr   = 3.156926e7
Day = 86400
G_Newton = 6.67384e-8
pc   = 3.08568025e18
AU = 1.49597871e13
AU_per_day = AU/86400
km = 1e5
km_per_sec_to_AU_per_day = 1e5/AU_per_day
Mbh = 4.3e6*Msun
days_per_year = 365.
###########################################


#### create an 'mbody' class to store ####
########## info about an object ##########
class mbody:
	def __init__(self, name="Unnamed Object", x=[0.0], y=[0.0], z=[0.0], vx=[0.0], vy=[0.0], vz=[0.0], a1=[0.0], a2=[0.0], a3=[0.0], a=[],  b=[], d=[1.0], e=[], f=[], g=[], i=[], l=[],  m=[0.0], n=[], o=[], p=[], q=[], r=[], s=[], t=[0.0], En=[], J=[], Jx=[], Jy=[], Jz=[], djdt=[], Jae=[], clo=False, clo_t=[], clo_obj=[], clo_a1=[], clo_a2=[], clo_e1=[], clo_e2=[], clo_i1=[], clo_i2=[], clo_dmin=[], dedt=[], timeinc='Days', clo_timeinc='Days', randcol=None): 
		if randcol is None: 
			randcol = (np.random.rand(), np.random.rand(), np.random.rand())
		self.name = name
		self.x    = x
		self.y    = y
		self.z    = z	
		self.vx   = vx
		self.vy   = vy
		self.vz   = vz	
		self.a1   = a1
		self.a2   = a2
		self.a3   = a3		
		self.a    = a
		self.b    = b
		self.d    = d
		self.e    = e
		self.f    = f
		self.g    = g
		self.i    = i
		self.l    = l
		self.mass = m
		self.n    = n
		self.o    = o
		self.p    = p
		self.q    = q
		self.r    = r
		self.spin = s
		self.Time = t
		self.energy = En
		self.j = J
		self.jx = Jx
		self.jy = Jy
		self.jz = Jz
		self.dJdt_Jinit = djdt
		self.J_ae = Jae
		self.vtot = []
		self.timeincrement = timeinc
		#The next few are relating to close encounters rather than constantly evaluated parameters#
		self.clo_Occur = clo
		self.clo_Time = clo_t
		self.clo_Object = clo_obj
		self.clo_dmin = clo_dmin
		self.clo_a1 = clo_a1
		self.clo_a2 = clo_a2
		self.clo_e1 = clo_e1
		self.clo_e2 = clo_e2
		self.clo_i1 = clo_i1
		self.clo_i2 = clo_i2
		self.dEdt_Einit = dedt
		self.clo_timeincrement = clo_timeinc
		self.colortriplet = randcol

"""
  a = semi-major axis (in AU)
  b = apocentre distance (in AU, b is short for Big q)
  d = density (g per cm^3)
  e = eccentricity
  f = true anomaly (degrees)
  g = argument of perihelion (degrees)
  i = inclination (degrees)
  l = mean anomaly (degrees)
  m = mass (solar masses)
  n = longitude of ascending node
  o = obliquity (degrees)
  p = longitude of perihelion (degrees)
  q = pericentre distance (AU)
  r = radial distance (AU)
  s = spin period (days)
  t = time (days)
  x, y or z = Cartesian coordinates x, y or z
  u, v or w = Cartesian velocities vx, vy or vz
"""
##########################################
##########################################

########### get and process all ###########
######### outputted files from M6 #########

def orbitload(filename, output_files, body_names):
	filenumb = filename
	outrigger = output_files
	andnobody = body_names
	# print 'In OrbitLoad'
	# print 'output_files: ', output_files
	# print 'filename: ', filenumb
	# print outrigger[filenumb]
	# print outrigger[filenumb][:-4]
	bodyname = outrigger[filenumb][:-4]
	globals()['body_names'].append(bodyname)

	globals()[bodyname] = mbody(name=bodyname) 	#initializes a bunch of body classes to the names of the objects

	if globals()[bodyname].name[:2] == 'DB':
		globals()[bodyname].colortriplet = "DarkBlue"
	elif globals()[bodyname].name[:2] == 'BG':
		globals()[bodyname].colortriplet = "Red"

	with open(outrigger[filenumb]) as f:
		data = []
		for line in f: #reads each line in the file into Python as a string, and splits that string at whitespaces
			singleline = []
			foo = line.split()
			for element in foo: 
				try: singleline.append(np.float64(element))
				except ValueError: singleline.append(element) #where it can, turns values into floats, passes where not (actually writes to a new array bc the original output is read-only)
			data.append(singleline)

	data = data[3:] #removes the top three lines because they're boring
	
	# try:
	# 	print '\r[{0}] {1}%'.format('#'*int(float(filename)/(len(output_files)/100)), float(filename)/(len(output_files)/100.)),  #writes a progress bar of the render
	# 	sys.stdout.flush()
	# except ZeroDivisionError: pass

	if data[0][1] == '(days)':
		time_increment = 'Days'
	elif data[0][1] == '(years)':
		time_increment = 'Years'
		
	globals()[bodyname].__dict__['timeincrement'] = time_increment

	try:
		data[0][1] = data[0][0]
		data[0] = data[0][1:] #turns 'Time [Days]' into just 'Time'
		data = np.transpose(data)
		# goes from 	"Time 	x 	y 	z 	..."	to 	"Time 	0 	1 	2 	..."
		# 					0 	1 	2 	3 					x 	1 	2 	4
		# 					1 	2 	4 	6 					y 	2 	4 	8
		# 					2 	4 	8 	9 					z 	3 	6 	9
		# 						(etc.) 							(etc.)
	except MemoryError:
		print bodyname, ': You have too much data either this computer or Numpy to handle. Decimating results...'
		decimations = 1
		while True: 
			try:
				array_len = len(data) #should just get the number of rows while using less memory than np.shape
				indexes = np.append(np.arange(int(array_len/10))*10, array_len-1)
				foodata = []
				index = 0
				while index <= array_len:
					if index == 0: 
						foodata.append(data[0])
					else: 
						foodata.append(data[10])
						del data[:10]
					index += 10
				foodata.append(data[-1])
				data = foodata
				data[0][1] = data[0][0]
				data[0] = data[0][1:] #turns 'Time [Days]' into just 'Time'
				data = np.transpose(data)
			except MemoryError: 
				print traceback.format_exc()
				print '\nTry ' + str(decimations) + ' failed. Re-decimating... (current length = '+str(array_len)+' rows.)'
				decimations += 1
				continue
			else: break
		print decimations, ' decimations were performed on ', globals()[bodyname]


	for line in data:
		datakey = line[0] #...which was useful because now I can interpret the first element in each line/list as the key to a property of the object and set that property to the remaining values
		values = []
		for val in line[1:]:
			try: values.append(np.float64(val))
			except ValueError: values.append(np.float64(0.0))
		# print datakey, values
		
		values = np.asarray(values)
		globals()[bodyname].__dict__[datakey] = values

	del data

def cloload(clo_filename, close_encounter_files):
	bodyname = clo_filename[:-4]
	# print bodyname
	
	with open(clo_filename) as f:
		data2 = []
		for line in f:
			singleline = []
			foo = line.split()
			for element in foo: 
				try: singleline.append(np.float64(element))
				except ValueError: singleline.append(element)
			data2.append(singleline)

	data2 = data2[3:]

	# try:
	# 	print '\r[{0}] {1}%'.format('#'*int(float(clo_filename)/(len(close_encounter_files)/100)), float(clo_filename)/(len(close_encounter_files)/100.)), #writes a progress bar of the render
	# 	sys.stdout.flush()
	# except ZeroDivisionError: pass

	if data2[0][1] == '(days)':
		clo_time_increment = 'Days'
	elif data2[0][1] == '(years)':
		clo_time_increment = 'Years'
		
	globals()[bodyname].__dict__['clo_timeincrement'] = clo_time_increment

	data2[0][1] = data2[0][0] 
	data2[0] = data2[0][1:] 	#takes 'Time (days)' to just 'Time'
	data2[0].pop(3)			#takes 'dmin (AU)' to just 'dmin'
	foo2 = []
	for element in data2[0]: foo2.append('clo_'+element)
	data2[0] = foo2
	data2 = np.transpose(data2)

	if len(data2[0][1:])>0: 
		globals()[bodyname].__dict__['clo_Occur'] = True
		
		for line in data2:
			datakey = line[0]
			values = []
			for val in line[1:]:
				try: values.append(np.float64(val))
				except ValueError: values.append(val)

			values = np.asarray(values)
			globals()[bodyname].__dict__[datakey] = values
	else: pass 
	
	del data2

def convertandj(bodyname):
	for step in np.arange(len(globals()[bodyname].Time)): #not sure this step works, agh why is python being dumb (I think it works though...)
		globals()[bodyname].mass[step] *= Msun #anyway, M6 outputs in stupid units so convert them to CGS
		globals()[bodyname].x[step] *=  AU
		globals()[bodyname].y[step] *=  AU
		globals()[bodyname].z[step] *=  AU
		globals()[bodyname].a[step] *=  AU
		globals()[bodyname].vx[step] *=  AU/Day
		globals()[bodyname].vy[step] *=  AU/Day
		globals()[bodyname].vz[step] *=  AU/Day
	
	if globals()[bodyname].clo_Occur == True: #same conversions, but only if a close encounter actually occurs
		if globals()[bodyname].timeincrement == 'Years' and globals()[bodyname].clo_timeincrement == 'Days':
			# print globals()[bodyname].clo_Time[:10]
			globals()[bodyname].clo_Time *= 1/days_per_year
			globals()[bodyname].__dict__['clo_timeincrement'] = 'Years'
		elif globals()[bodyname].timeincrement == 'Days' and globals()[bodyname].clo_timeincrement == 'Years':
			globals()[bodyname].clo_Time *= days_per_year
			globals()[bodyname].__dict__['clo_timeincrement'] = 'Days'

		for occurance in np.arange(len(globals()[bodyname].clo_Time)):
			try:
				globals()[bodyname].clo_a1[occurance] *=  AU
				globals()[bodyname].clo_a2[occurance] *=  AU
				for dmon in np.arange(len(globals()[bodyname].clo_dmin)):
					try: globals()[bodyname].clo_dmin[dmon] *=  AU
					except TypeError, ValueError: pass
			except TypeError: pass
			except Exception, err:
				print traceback.format_exc()
	else: pass

	if args.verbose == True: print globals()[bodyname].name + ': Close encounter == ', globals()[bodyname].clo_Occur

	currentr = np.sqrt(np.asarray(globals()[bodyname].x)**2 + np.asarray(globals()[bodyname].y)**2 + np.asarray(globals()[bodyname].z)**2) #get the radial coordinate
	curvtot = np.sqrt(np.asarray(globals()[bodyname].vx)**2 + np.asarray(globals()[bodyname].vy)**2 + np.asarray(globals()[bodyname].vz)**2) #and |v|

	globals()[bodyname].__dict__['r'] = np.ndarray.tolist(currentr)
	globals()[bodyname].__dict__['vtot'] = np.ndarray.tolist(curvtot)
	globals()[bodyname].__dict__['energy'] = np.ndarray.tolist(-(G_Newton*Mbh/np.abs(currentr)) + 0.5*(curvtot**2)) 

	# semimajor_E = G_Newton*(Mbh + np.asarray(globals()[bodyname].mass))/(2*np.asarray(globals()[bodyname].a)) #turns out it's the same as the Newtonian E
	
	dEdT_E = []
	for n in np.arange(len(globals()[bodyname].Time)-1): 
		dEdT_E.append(((globals()[bodyname].energy[n+1]-globals()[bodyname].energy[n])/(globals()[bodyname].Time[n+1]-globals()[bodyname].Time[n]))/globals()[bodyname].energy[0])
	globals()[bodyname].__dict__['dEdt_Einit'] = dEdT_E
	# try: 
		# ax.plot(globals()[bodyname].Time, globals()[bodyname].energy, '-', color='0.5', label="From Gravitational + Kinetic")
		# ax.plot(globals()[bodyname].Time, globals()[bodyname].energy, '--r', label="From Kepler's Laws")

	J_ae = []

	# print len(globals()[bodyname].a), len(globals()[bodyname].e)

	for n in np.arange(len(globals()[bodyname].Time)): #technically this is J/J_circ
		try: J_ae.append(np.sqrt(G_Newton * Mbh * np.abs(globals()[bodyname].a[n]) * (1 - globals()[bodyname].e[n]**2))/np.sqrt(G_Newton * Mbh * np.abs(globals()[bodyname].a[n])))
		except: 
			print 'a = ', globals()[bodyname].a[n], '\ne = ', globals()[bodyname].e[n] 
			J_ae.append(0.0)
	globals()[bodyname].__dict__['J_ae'] = J_ae

	J_rxv = []
	J_rxv_vec = []
	for n in np.arange(len(globals()[bodyname].Time)): 
		J_rxv_vec.append(np.cross([globals()[bodyname].x[n], globals()[bodyname].y[n], globals()[bodyname].z[n]], [globals()[bodyname].vx[n], globals()[bodyname].vy[n], globals()[bodyname].vz[n]])/np.sqrt(G_Newton * Mbh * np.abs(globals()[bodyname].a[n])))
		J_rxv.append(np.linalg.norm(np.cross([globals()[bodyname].x[n], globals()[bodyname].y[n], globals()[bodyname].z[n]], [globals()[bodyname].vx[n], globals()[bodyname].vy[n], globals()[bodyname].vz[n]]))/np.sqrt(G_Newton * Mbh * np.abs(globals()[bodyname].a[n])))
	globals()[bodyname].__dict__['j'] = J_rxv
	globals()[bodyname].__dict__['jx'] = np.transpose(J_rxv_vec)[0] 
	globals()[bodyname].__dict__['jy'] = np.transpose(J_rxv_vec)[1]
	globals()[bodyname].__dict__['jz'] = np.transpose(J_rxv_vec)[2]


	dJdT_J = []
	for n in np.arange(len(globals()[bodyname].Time)-1): 
		dJdT_J.append(((globals()[bodyname].j[n+1]-globals()[bodyname].j[n])/(globals()[bodyname].Time[n+1]-globals()[bodyname].Time[n]))/globals()[bodyname].j[0])
	globals()[bodyname].__dict__['dJdt_Jinit'] = dJdT_J

# def loaddata(body_names, output_files, close_encounter_files):
# 	bodayyy = body_names
# 	outputty = output_files
# 	thethirdkind = close_encounter_files
def loaddata(bodayyy, outputty, thethirdkind):
	# print 'in LoadData: '
	# print 'body_names: ', bodayyy
	# print 'output_files: ', outputty
	# print 'close_encounter_files: ', thethirdkind

	for filename in np.arange(len(outputty)): #note that 'filename' is actually a number
		#print '1. output_files: ', outputty
		namae = outputty[filename][:-4]
		#print '2. output_files: ', outputty
		try:
			print '\r[{0: <10}] {1: <7.2%} {2}: {3: <37}'.format('#'*int(float(filename)/(len(outputty)/10)), float(filename)/(len(outputty)), namae, 'Importing orbit data...'),
			sys.stdout.flush()
		except ZeroDivisionError: pass
		orbitload(filename, outputty, bodayyy)
		# print output_files[filename][:-4]
		# print 'Importing orbit data...'
		#orbitload(filename, bodayyy, outputty)
		clo_filename = outputty[filename][:-4]+'.clo'
		if clo_filename in thethirdkind: 
			try:
				print '\r[{0: <10}] {1: <7.2%} {2}: {3: <37}'.format('#'*int(float(filename)/(len(outputty)/10)), float(filename)/(len(outputty)), namae, 'Importing close encounter data...'),
				sys.stdout.flush()
			except ZeroDivisionError: pass
			# print 'Importing close encounter data...'
			cloload(clo_filename, thethirdkind)
		try:
			print '\r[{0: <10}] {1: <7.2%} {2}: {3: <37}'.format('#'*int(float(filename)/(len(outputty)/10)), float(filename)/(len(outputty)), namae, 'Calculating unit conversions and J...'),
			sys.stdout.flush()
		except ZeroDivisionError: pass
		convertandj(namae)
		try:
			print '\r[{0: <10}] {1: <7.2%} {2}: {3: <37}'.format('#'*int(float(filename)/(len(outputty)/10)), float(filename)/(len(outputty)), namae, 'Saving...'),
			sys.stdout.flush()
		except ZeroDivisionError: pass
		if os.path.isfile('body_names.pickle'):
			pickle.dump(bodayyy, file('body_names_foo.pickle','w'))
			os.remove('body_names.pickle')
			os.rename('body_names_foo.pickle', 'body_names.pickle') #to be absolutely sure we get everything, we're going to write a whole new version before deleting the old; this will be our method from hereon
		else: pickle.dump(bodayyy, file('body_names.pickle','w'))
		
		picklename = namae + '.gz'
		with open(picklename, 'wb') as fp:
			fp.write(zlib.compress(pickle.dumps(globals()[namae], pickle.HIGHEST_PROTOCOL),9))

		with open('loaded.txt', 'wb') as out:
			out.write(namae + '\n') #This writes after the save is finished, so anything on that list has already been saved; anything not on it has not

	print '\r[{0: <10}] {1: <7.2%} {2: <37}'.format('#'*10, 1., 'Import Complete.'),
	sys.stdout.flush()
	print '\n'
	
	with open('loaded.txt', 'wb') as out: out.write('!!##IMPORT_FINISHED##!!')


output_files          = glob.glob('*.aei')
close_encounter_files = glob.glob('*.clo')
# print output_files

if args.reload == False:
	
	if os.path.isfile('loaded.txt'): os.remove('loaded.txt')
	else: pass

	body_names       = []
	close_encounters = []

	print 'Loading data...\n'
	# loaddata(globals()['body_names'], globals()['output_files'], globals()['close_encounter_files'])
	loaddata(body_names, output_files, close_encounter_files)
	globals()['time_increment'] = globals()[body_names[0]].timeincrement
	globals()['integration_length'] = max(globals()[body_names[0]].Time)


else: #i.e. if args.reload == True:
	while True: #just allows me to use a continue statement to go directly to the 'else' after dealing with loaded.txt's existence
		if not os.path.isfile('loaded.txt'): 
			assumecomplete = raw_input('Error: there is no completion list; assume all zipped pickles are complete? [y/n]')
			while True:
				if (assumecomplete.lower())[0] == 'y':
					prevloaded = [x[:-3] for x in glob.glob('*.gz')]
					with open('loaded.txt', 'wb') as out:
						for nam in prevloaded: out.write(nam, '\n')
					break
				elif (assumecomplete.lower())[0] == 'n':
					startover = raw_input('Reload failed. Start from scratch? [y/n]')
					while True:
						if (assumecomplete.lower())[0] == 'y':
							print 'Loading data...\n'
							loaddata(body_names, output_files, close_encounter_files)
							break
						elif (assumecomplete.lower())[0] == 'n' or (assumecomplete.lower())[0] == 'q':
							print 'Quitting...'
							sys.exit()
						else: 
							assumecomplete = raw_input('Error: not a valid response. Please answer [y/n]')
							continue
					loaddata(body_names, output_files, close_encounter_files)
				else: 
					assumecomplete = raw_input('Error: not a valid response. Please answer [y/n]')
					continue
			continue #after it deals with all this, it'll go on to the next part
		else:
		    with open('loaded.txt') as f: 
		    	prevloaded = [x.rstrip('\n') for x in f.readlines()]
			
			if os.path.isfile('body_names.pickle'): 
				globals()['body_names'] = pickle.load(file('body_names.pickle'))
			else: 
				print 'Warning: body_names.pickle not found.'
				globals()['body_names'] = prevloaded

			if prevloaded[-1] == '!!##IMPORT_FINISHED##!!':
				print 'Importing files from pickles...'
				loadnum = 0
				prevloaded = globals()['body_names'] = prevloaded[:-1]
				for bodyname in prevloaded:
					picklename = bodyname + '.gz'
					with open(picklename, 'rb') as fp:
						dcdata = zlib.decompress(fp.read())
						globals()[bodyname] = pickle.loads(dcdata)
					# print '\r[{0}] {1}%'.format('#'*int(loadnum/(len(body_names)/100.)), float(loadnum)/(len(body_names)/100.)),  #progress bar
					print '\r[{0: <10}] {1: <7.2%}: Frame {2:} of {3}'.format('#'*int(float(loadnum)/(len(body_names)/10)), float(loadnum)/(len(body_names)), float(loadnum), len(body_names)),
					sys.stdout.flush()
					loadnum += 1
				print '\nImport complete.'
				time_increment = globals()[body_names[0]].timeincrement
				integration_length = max(globals()[body_names[0]].Time)	

			else:
				print 'Warning: Not all objects completed saving: Attempting to complete reading of unfinished objects before continuing.'
				prevfiles = [obj+'.aei' for obj in prevloaded]
				needs_loading = []
				for obj in output_files:
					if obj not in prevfiles: needs_loading.append(obj)
				loaddata(body_names, needs_loading, close_encounter_files) #load the new ones from the aei files
				
				print 'Loading previously saved objects...'
				for bodyname in prevloaded: #then load the previously saved files (hopefully at least this will get them all saved as class objects)
					picklename = bodyname + '.gz'
					with open(picklename, 'rb') as fp:
						dcdata = zlib.decompress(fp.read())
						globals()[bodyname] = pickle.loads(dcdata)
					# print '\r[{0}] {1}%'.format('#'*int(loadnum/(len(body_names)/100.)), float(loadnum)/(len(body_names)/100.)),  #progress bar
					print '\r[{0: <10}] {1: <7.2%}: Frame {2:} of {3}'.format('#'*int(float(loadnum)/(len(body_names)/10)), float(loadnum)/(len(body_names)), float(loadnum), len(body_names)),
					sys.stdout.flush()
					loadnum += 1
				print '\nImport complete.'
				time_increment = globals()[body_names[0]].timeincrement
				integration_length = max(globals()[body_names[0]].Time)	

			break

	# print body_names
	globals()['time_increment'] = globals()[body_names[0]].timeincrement
	globals()['integration_length'] = max(globals()[body_names[0]].Time)


if not os.path.isdir('images'): os.makedirs('images') #creates the directory if it doesn't already exist.
if not os.path.isdir('movies'): os.makedirs('movies')

for bodyname in body_names:
	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111)
	fig2.suptitle(globals()[bodyname].name + 'J/J_circ comparison')
	ax2.set_xlabel('Time [' + time_increment + ']')
	ax2.set_ylabel('J')
	if args.verbose == True: print bodyname, len(globals()[bodyname].Time), len(globals()[bodyname].J_ae), len(globals()[bodyname].j)
	limiting_dim_jae = min([len(globals()[bodyname].J_ae), len(globals()[bodyname].Time)])
	limiting_dim_j = min([len(globals()[bodyname].j), len(globals()[bodyname].Time)])
	ax2.plot(globals()[bodyname].Time[:limiting_dim_jae], globals()[bodyname].J_ae[:limiting_dim_jae], '-', alpha=0.5, label='from a and e')
	ax2.plot(globals()[bodyname].Time[:limiting_dim_j], globals()[bodyname].j[:limiting_dim_j], '-', alpha=0.5, label='from r x v')
	leg = plt.legend(loc='best', prop={'size':6})
	leg.get_frame().set_alpha(0.5)
	plt.savefig('images/J_comp-'+globals()[bodyname].name)
	plt.clf()

if args.clo:
	print 'Plotting eccentricities, energies, and angular momenta...'
	for bodyname in body_names:
		##############################################################################
		################ PLOT THE ECCENTRICITIES AND CLOSE ENCOUNTERS ################
		##############################################################################
		fig = plt.figure(figsize=(6,12), dpi=100)
		ax = fig.add_subplot(311)
		# ax0 = ax.twinx()
		axe = fig.add_subplot(312, sharex=ax)
		axj = axe.twinx()
		ax1 = fig.add_subplot(313, sharex=ax)

		ax.set_title(globals()[bodyname].name + ' Energy and Angular Momentum')

		ax.plot(globals()[bodyname].Time[:-1], globals()[bodyname].dEdt_Einit, '-', alpha=0.5, color='b', label='(dE/dt)/E_0')
		ax.plot(globals()[bodyname].Time[:-1], globals()[bodyname].dJdt_Jinit, '-', alpha=0.5, color='g', label='(dJ/dt)/J_0')

		ax1.plot(globals()[bodyname].Time, globals()[bodyname].e, '-', color='0.5')

		lns1 = axe.plot(globals()[bodyname].Time, globals()[bodyname].energy, '-', alpha=0.5, color='b', label='Energy')
		lns2 = axj.plot(globals()[bodyname].Time, globals()[bodyname].j, '-', alpha=0.5, color='g', label='J/J_circ')


		if globals()[bodyname].clo_Occur == True: #all this ugliness with annotation position is just trying to get the labels to not crash into each other

			# annotation_spacing_top = (max(dEdT_E)-min(dEdT_E))/len(globals()[bodyname].clo_Time) #(np.max(globals()[bodyname].energy)-np.min(globals()[bodyname].energy))/(len(globals()[bodyname].clo_Time))
			# annotation_spacing_bottom = (np.max(globals()[bodyname].e)-np.min(globals()[bodyname].e))/(len(globals()[bodyname].clo_Time))
			# annotation_number = 0
			# toporbottom = 1

			for occurance in np.arange(len(globals()[bodyname].clo_Time)):
				ax.axvline(globals()[bodyname].clo_Time[occurance], color='r', alpha=0.25)
				ax1.axvline(globals()[bodyname].clo_Time[occurance], color='r', alpha=0.25)
				axe.axvline(globals()[bodyname].clo_Time[occurance], color='r', alpha=0.25)

				# pickaYtop = (annotation_number*annotation_spacing_top)+np.min(dEdT_E)
				# pickaYbottom = (annotation_number*annotation_spacing_bottom)+np.min(globals()[bodyname].e)
				# annotation_number += 1
				
				# clo_data = str(globals()[bodyname].clo_Object[occurance]) + '\n' + str(globals()[bodyname].clo_Time[occurance]) + '\n' + 'minimum distance: ' + str(globals()[bodyname].clo_dmin[occurance]) + '\n' + 'eccentricity: ' + str(globals()[bodyname].clo_e1[occurance]) + ' --> ' + str(globals()[bodyname].clo_e2[occurance]) + '\n' + 'semi-major axis: ' + str(globals()[bodyname].clo_a1[occurance]) + ' --> ' + str(globals()[bodyname].clo_a2[occurance])
				# if toporbottom == 1:
				# 	ax.annotate(clo_data, xy=(globals()[bodyname].clo_Time[occurance], pickaYtop), xytext=(5, 5), color='b', size='xx-small', ha='left', textcoords='offset points')
				# elif toporbottom == -1:
				# 	ax1.annotate(clo_data, xy=(globals()[bodyname].clo_Time[occurance], pickaYbottom), xytext=(5, 5), color='b', size='xx-small', ha='left', textcoords='offset points')
				# toporbottom *= -1

		ax1.set_xlabel('Time [' + str(time_increment) + ']', fontsize='small')
		ax.set_ylabel('[unitless]', fontsize='small')
		# ax0.set_ylabel('[unitless]', fontsize='small')
		ax1.set_ylabel('Eccentricity', fontsize='small')
		axe.set_ylabel('E [Joules/unit mass]', fontsize='small')
		axj.set_ylabel('J [Joule seconds/unit mass]', fontsize='small')

		ax.grid()
		axe.grid()
		ax1.grid()

		ax.legend(loc='best', prop={'size':6})
		# leg.get_frame().set_alpha(0.5)

		lns = lns1+lns2
		labs = [l.get_label() for l in lns]
		axe.legend(lns, labs, loc='best', prop={'size':6})

		# lns2 = lns3+lns4
		# labs2 = [l.get_label() for l in lns2]
		ax.legend(loc='best', prop={'size':6})

		plt.tight_layout()

		plt.savefig('images/' + globals()[bodyname].name + '_fractionalEnergy')
		# plt.clf()
		plt.close('all')
		# fig = plt.figure()
		# ax = fig.add_subplot(211)
		# ax1 = fig.add_subplot(212, sharex=ax)
		# ax1.title(globals()[bodyname].name + ' Energy and Angular Momentum')
		# except ValueError: pass
		##############################################################################
		##############################################################################
		##############################################################################
# print 'reached line 248'
if args.ecc:
	##################################################################################
	############################# PLOT ALL ECCENTRICITIES ############################
	##################################################################################
	print 'Plotting eccentricities...'
	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111)
	fig2.suptitle('S Stars: All Eccentricities')
	ax2.set_xlabel('Time [' + str(time_increment) + ']')
	ax2.set_ylabel('e')
	for body in body_names:
			ax2.plot(globals()[body].Time, globals()[body].e, '-', label=globals()[body].name)
			if globals()[body].clo_Occur == True:
				ax2.plot(globals()[body].clo_Time, globals()[body].clo_e1, 'o', markersize=2, alpha=0.2) #putting a point on the line each time a close encounter occurs
				# for close in np.arange(len(globals()[body].clo_Time)): #and label it with the colliding body
					# ax2.annotate(globals()[body].clo_Object[close], xy=(globals()[body].clo_Time[close], globals()[body].clo_e1[close]), xytext=(-5, 5), ha='right', textcoords='offset points')
	leg = plt.legend(loc='best', prop={'size':6})
	leg.get_frame().set_alpha(0.5)
	ax2.grid()
	plt.savefig('images/All_Eccentricities')
	plt.clf()
	##################################################################################
	##################################################################################
	##################################################################################
# print 'reached line 307'
if args.dedt:
	##################################################################################
	############################ PLOT ALL (dE/dt)/E_init #############################
	##################################################################################
	print 'Plotting differential energy evolution...'
	fig3 = plt.figure()
	ax3 = fig3.add_subplot(111)
	fig3.suptitle('S Stars: All (dE/dt)/E_init')
	ax3.set_xlabel('Time [' + str(time_increment) + ']')
	ax3.set_ylabel('(dE/dt)/E_init')

	i = 1
	for body in body_names:
			# print np.shape(globals()[body].Time), np.shape(globals()[body].dEdt_Einit)
			ax3.plot(globals()[body].Time[:-1], globals()[body].dEdt_Einit, '-', label=globals()[body].name, alpha=0.5)

			if globals()[body].clo_Occur == True:
				for close in np.arange(len(globals()[body].clo_Time)):
					ax3.axvline(globals()[body].clo_Time[close], color='black', alpha=0.1) #the annotation gets a bit messy but I didn't want to put in another ugly and questionably functional location picker like above.
					# ax3.annotate(str(globals()[body].name + '\n&\n' + globals()[body].clo_Object[close]), xy=(globals()[body].clo_Time[close], i*max(globals()[body].dEdt_Einit)*0.75), xytext=(-5, 5), ha='right', textcoords='offset points', size='xx-small')
					# i *= -1 #so the most I did was have the annotation flip sides of the axis with each occurance.
	leg = plt.legend(loc='best', prop={'size':6})
	leg.get_frame().set_alpha(0.5)
	ax3.grid()
	for zoom in [0.05, 0.25, 2]:
		ax3.axis([0, np.max(globals()[body_names[0]].Time[:-1]), -zoom, zoom])
		name = 'images/All_dedtfracs-yrange'+str(zoom)+'.png'
		plt.savefig(name)
	plt.clf()
	##################################################################################
	##################################################################################
	##################################################################################
# print 'reached line 335'

needed_frames = int(np.log10(len(globals()[body_names[0]].Time)))+1

if args.moviemake is not None:
	##################################################################################
	############################# ALL ORBIT MOVIE-MAKER ##############################
	##################################################################################
	print '\n \n Now rendering movie of all orbits...'
	
	frame = 0

	# print args.moviemake
	# if args.moviemake in globals(): print 'yes'
	# else: print 'no'

	if len(args.moviemake) == 0: mov_framerate = 1
	else: 
		if len(args.moviemake) == 1:
			try: mov_framerate = int(args.moviemake[0])
			except ValueError: print 'Error: The framerate, if given, must be an integer.'
		else: 
			print 'Error: It seems you have specified multiple frame rates. Only the first integer value will be used.'
			for argument in np.arange(len(args.moviemake)):
				try: 
					mov_framerate = int(args.moviemake[argument])
					break
				except ValueError: pass
			if mov_framerate not in globals():
				print 'Error: none of the supplied arguments were viable framerates. Setting frame rate to 1 (Default).'
				mov_framerate = 1

	movfig = plt.figure(figsize=(36, 36), dpi=100)
	movfig.suptitle('Generated Cluster: ' + str(integration_length) + '-' + str(time_increment) + ' integration')
	movax1 = movfig.add_subplot(3, 3, 1)
	movax2 = movfig.add_subplot(3, 3, 2)
	movax3 = movfig.add_subplot(3, 3, 3)
	movax4 = movfig.add_subplot(3, 3, 4, projection='3d') #pyplot is not super great at 3d, and it becomes apparent in the 3d plot here; should probably switch to mayavi at some point
	movax5 = movfig.add_subplot(3, 3, 5)
	movax6 = movfig.add_subplot(3, 3, 6)
	movax4.set_aspect("equal")
	boxlen = 0.75E20 #0.5E17
	movax1.axis([-boxlen, boxlen, -boxlen, boxlen])
	movax2.axis([-boxlen, boxlen, -boxlen, boxlen])
	movax3.axis([-boxlen, boxlen, -boxlen, boxlen])
	movax4.axis([-boxlen, boxlen, -boxlen, boxlen])
	movax4.set_zlim3d(-boxlen, boxlen)

	movax5.axis([-boxlen, boxlen, 0, 1])
	movax6.axis([-boxlen, boxlen, 0, 1])
	# mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0,0,0))#, size=(1800,950))
	# mlab.disable_render = True

	def circles(x, y, s, c='b', ax=None, vmin=None, vmax=None, **kwargs):
	    """
	    Make a scatter of circles plot of x vs y, where x and y are sequence 
	    like objects of the same lengths. The size of circles are in data scale.

	    Parameters
	    ----------
	    x,y : array_like, shape (n, )
	        Input data
	    s : scalar or array_like, shape (n, ) 
	        Radius of circle in data scale (ie. in data unit)
	    c : color or sequence of color, optional, default : 'b'
	        `c` can be a single color format string, or a sequence of color
	        specifications of length `N`, or a sequence of `N` numbers to be
	        mapped to colors using the `cmap` and `norm` specified via kwargs.
	        Note that `c` should not be a single numeric RGB or
	        RGBA sequence because that is indistinguishable from an array of
	        values to be colormapped.  `c` can be a 2-D array in which the
	        rows are RGB or RGBA, however.
	    ax : Axes object, optional, default: None
	        Parent axes of the plot. It uses gca() if not specified.
	    vmin, vmax : scalar, optional, default: None
	        `vmin` and `vmax` are used in conjunction with `norm` to normalize
	        luminance data.  If either are `None`, the min and max of the
	        color array is used.  (Note if you pass a `norm` instance, your
	        settings for `vmin` and `vmax` will be ignored.)

	    Returns
	    -------
	    paths : `~matplotlib.collections.PathCollection`

	    Other parameters
	    ----------------
	    kwargs : `~matplotlib.collections.Collection` properties
	        eg. alpha, edgecolors, facecolors, linewidths, linestyles, norm, cmap

	    Examples
	    --------
	    a = np.arange(11)
	    circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')
	    """
	    from matplotlib.patches import Circle
	    from matplotlib.collections import PatchCollection
	    import pylab as plt
	    #import matplotlib.colors as colors

	    if ax is None:
	        ax = plt.gca()    

	    if isinstance(c,basestring):
	        color = c     # ie. use colors.colorConverter.to_rgba_array(c)
	    else:
	        color = None  # use cmap, norm after collection is created
	    kwargs.update(color=color)

	    if isinstance(s, (int, long, float)):
	        patches = [Circle((x_,y_), s) for x_,y_ in zip(x,y)]
	    else:
	        patches = [Circle((x_,y_), s_) for x_,y_,s_ in zip(x,y,s)]
	    collection = PatchCollection(patches, **kwargs)

	    if color is None:
	        collection.set_array(np.asarray(c))
	        if vmin is not None or vmax is not None:
	            collection.set_clim(vmin, vmax)

	    ax.add_collection(collection)
	    return collection #draws circles. Will be used later.

	def find_nearest(array,value):
	    idx = (np.abs(array-value)).argmin()
	    return array[idx] #finds the nearest value of an array to a given value

	for time in np.arange(len(globals()[body_names[0]].Time), step=mov_framerate)+1:
		# print len(globals()[body_names[0]].Time)
		# print time, len(globals()[body_names[0]].Time)
		# print '\r[{0}] {1}%'.format('#'*(time/(len(globals()[body_names[0]].Time)/100)), float(time)/(len(globals()[body_names[0]].Time)/100)), #writes a progress bar of the render
		try:
			print '\r[{0: <10}] {1: <7.2%}: Frame {2:} of {3}'.format('#'*int(float(time)/(len(globals()[body_names[0]].Time)/10)), float(time)/(len(globals()[body_names[0]].Time)), int(time), len(globals()[body_names[0]].Time)),
			sys.stdout.flush()
		except ZeroDivisionError: pass

		timelabel = 'Time [' + str(time_increment) + ']: ' + str(globals()[body_names[0]].Time[time-1])
		
		movfig.suptitle('Generated Cluster: ' + str(integration_length) + '-' + str(time_increment) + ' integration\n'+timelabel)

		movax1.set_xlabel('Y [cm]')
		movax1.set_ylabel('Z [cm]')
		# ax1.text(0.05, 0.1, timelabel, horizontalalignment='left', verticalalignment='center', transform = ax2.transAxes)
	 
		movax2.set_xlabel('X [cm]')
		movax2.set_ylabel('Z [cm]')
		# ax2.text(0.05, 0.1, timelabel, horizontalalignment='left', verticalalignment='center', transform = ax2.transAxes)

		movax3.set_xlabel('X [cm]')
		movax3.set_ylabel('Y [cm]')
		# ax3.text(0.05, 0.1, timelabel, horizontalalignment='left', verticalalignment='center', transform = ax2.transAxes)

		movax4.set_xlabel('X [cm]')
		movax4.set_ylabel('Y [cm]')
		movax4.set_zlabel('Z [cm]')
		# ax4.text(0.05, 0.1, 0.1, timelabel, horizontalalignment='left', verticalalignment='center', transform = ax2.transAxes)


		#MIGHT PUT THIS IN A DIFFERENT PLOT LATER BUT FOR NOW THIS IS THE QUICKEST WAY TO INCORPORATE THIS FUNCTIONALITY
		movax5.set_xlabel('a [cm]')
		movax5.set_ylabel('e')

		movax6.set_xlabel('projected radial position [cm]')
		movax6.set_ylabel('h statistic')
		# mlab.clf()
		# view = mlab.view()

		for body in body_names:
			# print body
			try:
				movax1.scatter(0,0,c='grey',s=40,marker='+')
				movax2.scatter(0,0,c='grey',s=40,marker='+')
				movax3.scatter(0,0,c='grey',s=40,marker='+')
				movax4.scatter(0,0,0,c='grey',s=40,marker='+') #puts a + at (0,0,0)

				movax1.plot(globals()[body].y[:time], globals()[body].z[:time], '-', c=globals()[body].colortriplet, label=globals()[body].name)
				movax2.plot(globals()[body].x[:time], globals()[body].z[:time], '-', c=globals()[body].colortriplet, label=globals()[body].name)
				movax3.plot(globals()[body].x[:time], globals()[body].y[:time], '-', c=globals()[body].colortriplet, label=globals()[body].name)
				if -boxlen<globals()[body].z[time-1]<boxlen:
					movax4.plot(globals()[body].x[:time], globals()[body].y[:time], globals()[body].z[:time], '-', c=globals()[body].colortriplet, label=globals()[body].name)
				else: pass #my workaround for matplotlib being shitty and not letting me control the z limits

				movax1.plot(globals()[body].y[time-1], globals()[body].z[time-1], 'o', c=globals()[body].colortriplet, label=globals()[body].name)
				movax2.plot(globals()[body].x[time-1], globals()[body].z[time-1], 'o', c=globals()[body].colortriplet, label=globals()[body].name)
				movax3.plot(globals()[body].x[time-1], globals()[body].y[time-1], 'o', c=globals()[body].colortriplet, label=globals()[body].name)
				if -boxlen<globals()[body].z[time-1]<boxlen:
					movax4.scatter(globals()[body].x[time-1],globals()[body].y[time-1],globals()[body].z[time-1], c=globals()[body].colortriplet)
				else: pass #my workaround for matplotlib being shitty and not letting me control the z limits
				
				# if globals()[body].clo_Occur == True:
				# 	for clotime in np.arange(len(globals()[body].clo_Time)):
				# 		besty = find_nearest(globals()[body].Time, globals()[body].clo_Time[clotime])
				# 		besttime = [besty, np.where(globals()[body].Time==besty)[0][0]]
				# 		#not 100% sure that clotimes are subsets of general time (and not interpolated), so this takes care of that
				# 		#time should be 1d and monotonic so I think this will work
				# 		try:
				# 			initialscale = float((globals()[body].clo_dmin)[clotime])
				# 		except ValueError: initialscale = 0.1
				# 		u = np.linspace(0, 2 * np.pi, 20)
				# 		v = np.linspace(0, np.pi, 20)
				# 		shiftmatrix = np.ndarray(shape=(20,20), dtype=float)
				# 		shiftmatrix.fill(1.0)
				# 		if (np.abs(besttime[0]-globals()[body].Time[time]) <= mov_framerate) and (np.abs(besttime[0]-globals()[body].Time[time]) <= np.abs(besttime[0]-globals()[body].Time[time+mov_framerate])):
				# 			#i.e., if this is the frame closest to your 'besttime,' taking framerate into account (it's less than one frame away from the best time, and closer than the next frame)
				# 			xc = initialscale * np.outer(np.cos(u), np.sin(v)) + (globals()[body].x[time-1] * shiftmatrix)
				# 			yc = initialscale * np.outer(np.sin(u), np.sin(v)) + (globals()[body].y[time-1] * shiftmatrix)
				# 			zc = initialscale * np.outer(np.ones(np.size(u)), np.cos(v)) + (globals()[body].z[time-1] * shiftmatrix)
				# 			movax4.plot_surface(xc, yc, zc,  rstride=4, cstride=4, color='r', linewidth=0, alpha=1.0, antialiased=False)
				# 			# movax3.scatter(globals()[body].x[time-1],globals()[body].y[time-1], c=(1.,0.,0.), s=initialscale, alpha=1., marker='o')
				# 			# movax1.scatter(globals()[body].y[time-1],globals()[body].z[time-1], c=(1.,0.,0.), s=initialscale, alpha=1., marker='o')
				# 			# movax2.scatter(globals()[body].x[time-1],globals()[body].z[time-1], c=(1.,0.,0.), s=initialscale, alpha=1., marker='o')
				# 		elif globals()[body].Time[time]-besttime[0] <= 10*mov_framerate: #if it's less than 10 frames away
				# 			framesaway = float(((globals()[body]).Time[time]-besttime[0])/mov_framerate)
				# 			xc = initialscale * (framesaway**2) * np.outer(np.cos(u), np.sin(v)) + (globals()[body].x[time-1] * shiftmatrix)
				# 			yc = initialscale * (framesaway**2) * np.outer(np.sin(u), np.sin(v)) + (globals()[body].y[time-1] * shiftmatrix)
				# 			zc = initialscale * (framesaway**2) * np.outer(np.ones(np.size(u)), np.cos(v)) + (globals()[body].z[time-1] * shiftmatrix)
				# 			clo_alpha = (1.0-(framesaway/20.))
				# 			if 1.0 >= clo_alpha >= 0.0: pass
				# 			else: clo_alpha = 0.0
				# 			movax4.plot_surface(xc, yc, zc,  rstride=4, cstride=4, color='r', linewidth=0, alpha=clo_alpha, antialiased=False)
				# 			# movax3.scatter(globals()[body].x[time-1],globals()[body].y[time-1], c=(1.,0.,0.), s=(initialscale*(framesaway**2)), alpha=(1.0-(framesaway/10.)), marker='o')
				# 			# movax1.scatter(globals()[body].y[time-1],globals()[body].z[time-1], c=(1.,0.,0.), s=(initialscale*(framesaway**2)), alpha=(1.0-(framesaway/10.)), marker='o')
				# 			# movax2.scatter(globals()[body].x[time-1],globals()[body].z[time-1], c=(1.,0.,0.), s=(initialscale*(framesaway**2)), alpha=(1.0-(framesaway/10.)), marker='o')
				# 		else: pass #this should work but I haven't tested it yet

				if args.anonymous == False:
					movax1.annotate(globals()[body].name, xy=(globals()[body].y[time-1],globals()[body].z[time-1]), xytext=(-5, 5), ha='right', textcoords='offset points')
					movax2.annotate(globals()[body].name, xy=(globals()[body].x[time-1],globals()[body].z[time-1]), xytext=(-5, 5), ha='right', textcoords='offset points')
					movax3.annotate(globals()[body].name, xy=(globals()[body].x[time-1],globals()[body].y[time-1]), xytext=(-5, 5), ha='right', textcoords='offset points')
				else: pass
				

				movax5.scatter(globals()[body].a[time-1],globals()[body].e[time-1],c=globals()[body].colortriplet,marker='o')

				#calculate projected radius in the y-z plane (assuming edge-on to disk).
				proj_r = np.sqrt((globals()[body].y[time-1])**2 + (globals()[body].z[time-1])**2)
				#calculate the h statistic in the y-z plane. Might put this in the class attributes later.
				h_stat = ((globals()[body].y[time-1]*globals()[body].vz[time-1])-(globals()[body].z[time-1]*globals()[body].vy[time-1]))/np.sqrt(G_Newton*Mbh*p)

				movax6.scatter(proj_r,h_stat,c=globals()[body].colortriplet,marker='o')


				# x2, y2, _ = proj3d.proj_transform(globals()[body].x[time-1],globals()[body].y[time-1],globals()[body].z[time-1], movax4.get_proj())

				# if (-boxlen<globals()[body].x[time-1]<boxlen) and (-boxlen<globals()[body].y[time-1]<boxlen) and (-boxlen<globals()[body].z[time-1]<boxlen):
				# 	label = movax4.annotate(globals()[body].name, xy = (x2, y2), xytext = (-5, 5), textcoords = 'offset points', ha = 'right', va = 'bottom')#,
				# 	    # bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
				# 	    # arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
				# else: pass
					# ax4.annotate(globals()[body].name, xy=(x2, y2), xytext=(-5, 5), ha='right', textcoords='offset points')
			except IndexError: pass #should stop it from complaining if an object gets ejected

		movax4.azim = frame
		# axes = mlab.axes()

		# arr = mlab.screenshot()

		# movax4.imshow(arr)
		# movax4.axis('off')

		fileindex = 'movies/'+str(frame).zfill(needed_frames)+'.png'
		plt.savefig(fileindex)
		plt.clf()

		# mlab.view(azimuth=mlab.view()[0]+1)
		frame += 1

		movfig.suptitle('Generated Cluster: ' + str(integration_length) + '-' + str(time_increment) + ' integration')
		movax1 = movfig.add_subplot(3, 3, 1)
		movax2 = movfig.add_subplot(3, 3, 2)
		movax3 = movfig.add_subplot(3, 3, 3)
		movax4 = movfig.add_subplot(3, 3, 4, projection='3d') #pyplot is not super great at 3d, and it becomes apparent in the 3d plot here; should probably switch to mayavi at some point
		movax5 = movfig.add_subplot(3, 3, 5)
		movax6 = movfig.add_subplot(3, 3, 6)
		movax4.set_aspect("equal")
		movax1.axis([-boxlen, boxlen, -boxlen, boxlen])
		movax2.axis([-boxlen, boxlen, -boxlen, boxlen])
		movax3.axis([-boxlen, boxlen, -boxlen, boxlen])
		movax4.axis([-boxlen, boxlen, -boxlen, boxlen])
		movax4.set_zlim3d(-boxlen, boxlen)

		movax5.axis([-boxlen, boxlen, 0, 1])
		movax6.axis([-boxlen, boxlen, 0, 1])

	if commands.getstatusoutput("man ffmpeg")[0] == 0: #there are probably classier ways to do this but it's late and I'm not thinking of them. Basically, if ffmpeg is installed on your computer, it's not on the astro servers and I don't think I can install it for them;) 
		syscommand = 'ffmpeg -f image2 -r 24 -i movies/%0'+str(needed_frames)+'d.png -vcodec mpeg4 -y movies/movie.mp4'	#this condition will be true, and you'll have the option to turn them into a movie with it. Otherwise, you're on your own.
		delframes = raw_input("Delete individual frames after creating movie? [y/n] ")
		os.system(syscommand)
		while True:
			if (delframes.lower())[0] == 'y':
				for movframe in list(np.arange(int(globals()['frame']))): 
					try: os.remove('movies/'+str(movframe).zfill(int(globals()['needed_frames']))+'.png')
					except IOError: pass
				break
			elif (delframes.lower())[0] == 'n': break
			else: 
				print 'That response is not an option.'
				delframes = raw_input("Delete individual frames after creating movie? [y/n] ")
				continue
	else: print 'FFMPEG does not seem to be installed on this computer; the user is left to their own devices to combine the frames into a video.'
		

	##################################################################################
	##################################################################################
	##################################################################################

def singleorbit(star, framerate, needed_frames=6):
	print '\n \n Now rendering a movie of the orbit of ' + star + '...'
	frame = 0
	movfig = plt.figure(figsize=(36, 19), dpi=100)
	# movfig.suptitle(star + ': 1000-year integration \nTime [days]: 0.00')
	movax1 = movfig.add_subplot(2, 2, 1)
	movax2 = movfig.add_subplot(2, 2, 2)
	movax3 = movfig.add_subplot(2, 2, 3)
	movax4 = movfig.add_subplot(2, 2, 4, projection='3d') #pyplot is not super great at 3d, and it becomes apparent in the 3d plot here; should probably switch to mayavi at some point
	movax1.axis([1.5*np.min(globals()[star].y), 1.5*np.max(globals()[star].y), 1.5*np.min(globals()[star].z), 1.5*np.max(globals()[star].z)])
	movax2.axis([1.5*np.min(globals()[star].x), 1.5*np.max(globals()[star].x), 1.5*np.min(globals()[star].z), 1.5*np.max(globals()[star].z)])
	movax3.axis([1.5*np.min(globals()[star].x), 1.5*np.max(globals()[star].x), 1.5*np.min(globals()[star].y), 1.5*np.max(globals()[star].y)])
	movax4.axis([1.5*np.min(globals()[star].x), 1.5*np.max(globals()[star].x), 1.5*np.min(globals()[star].y), 1.5*np.max(globals()[star].y)])

	for time in np.arange(len(globals()[body_names[0]].Time), step=framerate)+1:
		# print len(globals()[body_names[0]].Time)
		# print time, len(globals()[body_names[0]].Time)
		# print '\r[{0}] {1}%'.format('#'*(time/(len(globals()[body_names[0]].Time)/100)), float(time)/(len(globals()[body_names[0]].Time)/100)), #writes a progress bar of the render
		print '\r[{0: <10}] {1: <7.2%}: Frame {2:} of {3}'.format('#'*int(float(time)/(len(globals()[body_names[0]].Time)/10)), float(time)/(len(globals()[body_names[0]].Time)), int(time), len(globals()[body_names[0]].Time)),
		sys.stdout.flush()

		timelabel = 'Time [days]: ' + str(globals()[body_names[0]].Time[time-1])
			# movfig.suptitle('Generated Cluster: ' + str(integration_length) + '-' + str(time_increment) + ' integration')

		movfig.suptitle(star + ': '+ str(integration_length) + '-' + str(time_increment) + ' integration \n'+timelabel)

		movax1.set_xlabel('Y [cm]')
		movax1.set_ylabel('Z [cm]')
		# ax1.text(0.05, 0.1, timelabel, horizontalalignment='left', verticalalignment='center', transform = ax2.transAxes)
	 
		movax2.set_xlabel('X [cm]')
		movax2.set_ylabel('Z [cm]')
		# ax2.text(0.05, 0.1, timelabel, horizontalalignment='left', verticalalignment='center', transform = ax2.transAxes)

		movax3.set_xlabel('X [cm]')
		movax3.set_ylabel('Y [cm]')
		# ax3.text(0.05, 0.1, timelabel, horizontalalignment='left', verticalalignment='center', transform = ax2.transAxes)

		movax4.set_xlabel('X [cm]')
		movax4.set_ylabel('Y [cm]')
		movax4.set_zlabel('Z [cm]')
		# ax4.text(0.05, 0.1, 0.1, timelabel, horizontalalignment='left', verticalalignment='center', transform = ax2.transAxes)
		
			# print body
		movax1.plot(globals()[star].y[:time], globals()[star].z[:time], '-', c=globals()[star].colortriplet, label=globals()[star].name)
		movax2.plot(globals()[star].x[:time], globals()[star].z[:time], '-', c=globals()[star].colortriplet, label=globals()[star].name)
		movax3.plot(globals()[star].x[:time], globals()[star].y[:time], '-', c=globals()[star].colortriplet, label=globals()[star].name)
		movax4.plot(globals()[star].x[:time], globals()[star].y[:time], globals()[star].z[:time], '-', c=globals()[star].colortriplet, label=globals()[star].name)



		movax1.plot(globals()[star].y[time-1], globals()[star].z[time-1], 'o', c=globals()[star].colortriplet, label=globals()[star].name)
		movax2.plot(globals()[star].x[time-1], globals()[star].z[time-1], 'o', c=globals()[star].colortriplet, label=globals()[star].name)
		movax3.plot(globals()[star].x[time-1], globals()[star].y[time-1], 'o', c=globals()[star].colortriplet, label=globals()[star].name)
		movax4.scatter(globals()[star].x[time-1],globals()[star].y[time-1],globals()[star].z[time-1], c=globals()[star].colortriplet)

		# movax1.annotate(globals()[body].name, xy=(globals()[body].y[time-1],globals()[body].z[time-1]), xytext=(-5, 5), ha='right', textcoords='offset points')
		# movax2.annotate(globals()[body].name, xy=(globals()[body].x[time-1],globals()[body].z[time-1]), xytext=(-5, 5), ha='right', textcoords='offset points')
		# movax3.annotate(globals()[body].name, xy=(globals()[body].x[time-1],globals()[body].y[time-1]), xytext=(-5, 5), ha='right', textcoords='offset points')
		# ax4.annotate(globals()[body].name, xy=(x2, y2), xytext=(-5, 5), ha='right', textcoords='offset points')

		fileindex = 'movies/'+ star + '_' + str(frame).zfill(needed_frames)+'.png'
		plt.savefig(fileindex)
		plt.clf()
		frame += 1

		# movfig.suptitle(star + ': 1000-year integration \n'+timelabel)
		movax1 = movfig.add_subplot(2, 2, 1)
		movax2 = movfig.add_subplot(2, 2, 2)
		movax3 = movfig.add_subplot(2, 2, 3)
		movax4 = movfig.add_subplot(2, 2, 4, projection='3d')
		movax1.axis([1.125*np.min(globals()[star].y), 1.125*np.max(globals()[star].y), 1.125*np.min(globals()[star].z), 1.125*np.max(globals()[star].z)])
		movax2.axis([1.125*np.min(globals()[star].x), 1.125*np.max(globals()[star].x), 1.125*np.min(globals()[star].z), 1.125*np.max(globals()[star].z)])
		movax3.axis([1.125*np.min(globals()[star].x), 1.125*np.max(globals()[star].x), 1.125*np.min(globals()[star].y), 1.125*np.max(globals()[star].y)])
		movax4.axis([1.125*np.min(globals()[star].x), 1.125*np.max(globals()[star].x), 1.125*np.min(globals()[star].y), 1.125*np.max(globals()[star].y)]) #the code here is almost identical to that in --moviemake; see there for documentation. !!!!!should add framerate option.
# print args.orbitmovies
# if args.orbitmovies: print 'yes'
# else: print 'no'
if args.orbitmovies is not None:
	##################################################################################
	############################ SINGLE ORBIT MOVIE-MAKER ############################
	##################################################################################
	# print args
	if args.orbitmovies == None:
		for star in body_names:
			singleorbit(star, 1, needed_frames) #if you don't give it any commands, it'll do all of them.

		if commands.getstatusoutput("man ffmpeg")[0] == 0: 
			for star in body_names:
				syscommand = 'ffmpeg -f image2 -r 24 -i movies/'+str(star)+'_%0'+str(needed_frames)+'d.png -vcodec mpeg4 -y movies/'+str(star)+'_movie.mp4'
				delframes = raw_input("Delete individual frames after creating movies? [y/n] ")
				os.system(syscommand)
				while True:
					if (delframes.lower())[0] == 'y':
						for movframe in list(np.arange(int(globals()['frame']))): 
							try: os.remove('movies/'+star+'_'+str(movframe).zfill(int(globals()['needed_frames']))+'.png')
							except IOError: pass
						break
					elif (delframes.lower())[0] == 'n': break
					else: 
						print 'That response is not an option.'
						delframes = raw_input("Delete individual frames after creating movies? [y/n] ")
						continue
		else: print 'FFMPEG does not seem to be installed on this computer; the user is left to their own devices to combine the frames into a video.'
	
	else: 
		orbarg_integer = []
		for orbarg in np.arange(len(args.orbitmovies)):
			try: 
				frame_fraction = int((args.orbitmovies)[orbarg])
				orbarg_integer.append(frame_fraction) #if some of the arguments you give it are integers, remove them from the list of objects to plot and assume they are intended as frame rates.
				(args.orbitmovies).pop(orbarg)
			except Exception: 
				pass

		if len(orbarg_integer) <= 1: #if there's only one or zero assumed frame rate(s), we're good; otherwise it'll display the error message below.
			if len(orbarg_integer) == 1: orbarg_integer = orbarg_integer[0]
			else: orbarg_integer = 1
			
			if len(args.orbitmovies) == 0: #once again, if we only have the frame rate, it'll plot all of them at that frame rate (eventually)
				for star in body_names:
					singleorbit(star, orbarg_integer, needed_frames)
				if commands.getstatusoutput("man ffmpeg")[0] == 0: 
					for star in body_names:
						syscommand = 'ffmpeg -f image2 -r 24 -i movies/'+str(star)+'_%0'+str(needed_frames)+'d.png -vcodec mpeg4 -y movies/'+str(star)+'_movie.mp4'
						delframes = raw_input("Delete individual frames after creating movies? [y/n] ")
						os.system(syscommand)
						while True:
							if (delframes.lower())[0] == 'y':
								for movframe in list(np.arange(int(globals()['frame']))): 
									try: os.remove('movies/'+star+'_'+str(movframe).zfill(int(globals()['needed_frames']))+'.png')
									except IOError: pass
								break
							elif (delframes.lower())[0] == 'n': break
							else: 
								print 'That response is not an option.'
								delframes = raw_input("Delete individual frames after creating movies? [y/n] ")
								continue
				else: print 'FFMPEG does not seem to be installed on this computer; the user is left to their own devices to combine the frames into a video.'
			
			else: 
				plotted_stars = []
				for starname in args.orbitmovies:	#if we have a frame rate as well as specified stars, it'll plot those at that frame rate
					starname = starname.upper()
					if starname in body_names:
						plotted_stars.append(starname)
						singleorbit(starname, orbarg_integer, needed_frames)
					else: print 'Error: A star by the name of ' + starname + ' was not included in the integration.'
				if commands.getstatusoutput("man ffmpeg")[0] == 0: 
					for star in plotted_stars:
						syscommand = 'ffmpeg -f image2 -r 24 -i movies/'+str(star)+'_%0'+str(needed_frames)+'d.png -vcodec mpeg4 -y movies/'+str(star)+'_movie.mp4'
						delframes = raw_input("Delete individual frames after creating movies? [y/n] ")
						os.system(syscommand)
						while True:
							if (delframes.lower())[0] == 'y':
								for movframe in list(np.arange(int(globals()['frame']))): 
									try: os.remove('movies/'+star+'_'+str(movframe).zfill(int(globals()['needed_frames']))+'.png')
									except OSError: pass
								break
							elif (delframes.lower())[0] == 'n': break
							else: 
								print 'That response is not an option.'
								delframes = raw_input("Delete individual frames after creating movies? [y/n] ")
								continue
				else: print 'FFMPEG does not seem to be installed on this computer; the user is left to their own devices to combine the frames into a video.'
			

		else: print 'Error: You have either specified the frame rate multiple times or refered to stars by name only, confusing the program. Please rectify this.'

	##################################################################################
	##################################################################################
	##################################################################################



















