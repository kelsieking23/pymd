from ndx import Ndx


# # string = '(ri1 or ri2 or ri3) and backbone'
# string = 'ri1-ri5'
n = Ndx('mn11.gro')
for key, value in n.types.items():
    print(key)
n.select('p1')
n.select('p2-p3')
n.makeNDX('index.ndx')

# print(parsed_selection)

