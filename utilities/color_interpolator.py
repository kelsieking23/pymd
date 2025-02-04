import warnings
warnings.filterwarnings('ignore')
import pandas as pd
from pymd.structure.library import get_scale
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))



def plot_interpolation(c1, c2, n, mix=0):
    fig, ax = plt.subplots(figsize=(8, 5))
    for x in range(n+1):
        ax.axvline(x, color=colorFader(c1,c2,x/n), linewidth=4) 
        print(colorFader(c1,c2,x/n))
    plt.show()

def interploate_hydrophobicity(scale, c1, c2, mix=0, print_pymol_cmd = False):
    scale_list = get_scale(scale, get_as='list')
    sorted_list = sorted(scale_list, key = lambda x:x[1])
    interpolation_data = []
    for i, item in enumerate(sorted_list):
        code = item[0]
        val = item[1]
        interpolation_data.append([code, val, colorFader(c1, c2, i/19)])
    df = pd.DataFrame(interpolation_data, columns=['aa', 'score', 'hex'])
    if print_pymol_cmd:
        for index in df.index:
            aa = df.loc[index, 'aa']
            hex = df.loc[index, 'hex']
            rgb = hex_to_rgb(hex)
            print(f"cmd.set_color('color_{aa.lower()}X', {rgb})")
            print(f"cmd.color('color_{aa.lower()}X', '('+s+' and resn {aa})')")
    return df


def create_interpolations_rainbow():
    color_list = [
        ['#e62615', '#fadcd9'], # fire truck red -> pale pink
        ['#f2921d', '#f7cf9e'], # dark orange -> peach
        ['#f5ee1d', '#f7f5a8'], # bright yellow -> pale yellow
        ['#3d911d', '#7c9174'], # med green -> seafoam
        ['#0a6ed1', '#bfd9f2'], # primary blue -> baby blue
        ['#571270', '#b695c2']  # dark blue -> periwinkle

    ]
    for item in color_list:
        interploate_hydrophobicity('kyte-doolittle', item[0], item[1], print_pymol_cmd=True)
        print('****************************************')
        print('\n\n')

# create_interpolations_rainbow()