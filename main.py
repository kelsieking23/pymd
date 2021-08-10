import os
import sys
sys.path.append(os.getcwd())

from pymd.mdanalysis.system import System

test = System(root='/work/cascades/aksharp/OR2T7_MD/WT/MDrun/', reps=4, name='WT')
test.run.dssp(dssp_path='/home/aksharp7/software/bin/dssp', email='aksharp7@vt.edu', start=0, stop=1000)