#  Open binary file to plot Cartesian co-ords 
#  Note: current *.aei files output in Kepler elements

import numpy as np
import matplotlib.pyplot as plt
import struct

file = open("xv.out", "rb")
a = np.fromfile(file, dtype=np.dtype('f8') )


print a
print a[0], a[1], a[2]
#plt.plot(a)
#plt.show()