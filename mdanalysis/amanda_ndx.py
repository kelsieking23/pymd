from pymd.mdanalysis.ndx import Ndx


f = open('pd1_mut1.ndx','r')
indeces = {}
keys = []
for line in f:
    if line.startswith('['):
        key = line[1:-1].strip()
        indeces[key] = []
        if key not in keys:
            keys.append(key)
    else:
        nums = line.split()
        for num in nums:
            indeces[key].append(num.strip())
    
f.close()
for key, value in indeces.items():
    print(key,value)

n = Ndx('pd1_mut1.gro')
sele = n.select('backbone')
print(sele)
for index in keys:
    indeces['{}_backbone'.format(index)] = [i for i in indeces[index] if i in sele]

n.selections = indeces
n.makeNDX('pd1_mut1_backbone.ndx')