# CAUTIONS:
# compile the code before running
# DON'T touch/compile the code once you submit this script

import os
from subprocess import Popen, PIPE
import time

gtCases = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5]
#gtCases = [1.1, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5]
waitTime = 60

for gtCase in gtCases:
	#os.system("./Tunegt %s" %gtCase)
	#os.system("mv readPDF/INPUT_auto/Tunegt.root readPDF/INPUT_auto/DiffCroSec%s.root" %gtCase)
	#time.sleep(waitTime)

	os.system("python readPDF/2d.py readPDF/INPUT_auto/DiffCroSec%s.root readPDF/OUTPUT_auto/reweighing%s.root 2000 3000" %(gtCase, gtCase))
	time.sleep(waitTime)


