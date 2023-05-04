import os

class Infer:

    def __init__(self):
        pass
    

    def getHeader(self):
        header = []
        f = open('infer_header.sh', 'r')
        for line in f:
            line_parts = line.split()
            try:
                if ('-t' in line_parts[-1]) and ('all' not in line_parts[-1]):
                    line = line.strip() + ' ' + str(walltime) + '\n'
                    header.append(line)
                elif '--job-name=' in line_parts[-1]:
                    line = line.strip() + job_name + '\n'
                    header.append(line)
                else:
                    header.append(line)
            except:
                header.append(line)
        f.close()
        return header
    