import os
import subprocess
from pymd.mdanalysis.analysis import Analysis
import pymd
DIRECTORY = pymd.dir_path

class Gmx(Analysis):

    def __init__(self, inp, top, ndx='index.ndx', cluster='infer', version='2020.4', 
                 partition='p100', email=None, **kwargs):
        self._inp = inp
        self._topfile = top
        self.ndx = ndx
        self.cluster = cluster
        self.version = version
        self.email = email
        self.partition = partition
        self.__dict__.update(kwargs)
    
    def getHeader(self, job_name, walltime, nodes=1, ntasks_per_node=24):
        header = []
        base = os.path.join(DIRECTORY, 'mdrun', 'sh', self.cluster, self.partition, '{}.sh'.format(self.version))
        f = open(base, 'r')
        for line in f:
            line_parts = line.split()
            try:
                if ('-t' in line_parts[-1]) and ('all' not in line_parts[-1]):
                    line = line.strip() + ' ' + str(walltime) + '\n'
                    header.append(line)
                elif '--job-name=' in line_parts[-1]:
                    line = line.strip() + job_name + '\n'
                    header.append(line)
                elif '--mail-user=' in line_parts[-1]:
                    if self.email is not None:
                        line = line.strip() + self.email + '\n'
                        header.append(line)
                elif '--nodes=' in line_parts[-1]:
                    line = line.strip() + str(nodes) + '\n'
                    header.append(line)
                elif '--ntasks-per-node=' in line_parts[-1]:
                    if ntasks_per_node is None:
                        continue
                    else:
                        line = line.strip() + str(ntasks_per_node) + '\n'
                        header.append(line)
                else:
                    header.append(line)
            except:
                header.append(line)
        f.close()
        return header
    
    def mindist(self, job_name='mindist', walltime='50:00:00', nodes=1, ntasks_per_node=24, )
gmx = Gmx('string', 'string', email='kelsieking23@vt.edu')
header = gmx.getHeader('mindist', '72:00:00')
for line in header:
    print(line)