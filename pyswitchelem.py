import glob, os

elements = glob.glob('element*.in')
if ('element-edited.in' in elements) and ('element.in' in elements):
	os.rename('element.in', 'element-original.in')
	os.rename('element-edited.in', 'element.in')
elif ('element-original.in' in elements) and ('element.in' in elements):
	os.rename('element.in', 'element-edited.in')
	os.rename('element-original.in', 'element.in')
else: print "This program only switches element-edited.in and/or element-original.in"