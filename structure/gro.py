class GroHandler:

    def __init__(self, gro):
        self.gro = gro
        self.data = None

    
    def fixData(self):
        f = open(self.gro, 'r')
        contents = f.readlines()
        f.close()
        data = []
        i = 0
        m = 0
        for line in contents:
            if i < 2:
                i += 1
                continue
            elif line == contents[-1]:
                break
            else:
                line_parts = line.split()
                res_name = line_parts[0][-3:]
                res_num = line_parts[0][-3]
                if res_name == 'SOL':
                    break
                atom_type = line_parts[1]
                x, y, z = line_parts[3:6]
                
                m += 1
