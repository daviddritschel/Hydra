#!/usr/bin/env python
from matplotlib.colors import LinearSegmentedColormap
from numpy import *
import numpy as np
from scipy.optimize import ridder
from matplotlib.pyplot import cm

# Inspired by http://stackoverflow.com/questions/4296249/how-do-i-convert-a-hex-triplet-to-an-rgb-tuple-and-back
# Available under Creative Commons Attribution License
_NUMERALS = '0123456789abcdefABCDEF'
_HEXDEC = {v: int(v, 16) for v in (x+y for x in _NUMERALS for y in _NUMERALS)}
LOWERCASE, UPPERCASE = 'x', 'X'

# See https://github.com/hetland/octant/blob/master/octant/sandbox/plotting.py
# Available under BSD
def cmap_brightened(cmap,factor=0.5):
    """
    Brightens colormap cmap with using a saturation factor 'factor'
    (0.5 by default).
    """
    return cmap_map(lambda x: (1.-factor) + factor*x, cmap)

# See https://github.com/hetland/octant/blob/master/octant/sandbox/plotting.py
# And scipy cookbook. Scipy is under BSD
def cmap_map(function,cmap):
    """ Applies function (which should operate on vectors of shape 3:
    [r, g, b], on colormap cmap. This routine will break any discontinuous
    points in a colormap.
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red','green','blue'):         step_dict[key] = map(lambda x: x[0], cdict[key])
    step_list = reduce(lambda x, y: x+y, step_dict.values())
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : np.array(cmap(step)[0:3])
    old_LUT = np.array(map( reduced_cmap, step_list))
    new_LUT = np.array(map( function, old_LUT))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i,key in enumerate(('red','green','blue')):
        this_cdict = {}
        for j,step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j,i]
            elif new_LUT[j,i]!=old_LUT[j,i]:
                this_cdict[step] = new_LUT[j,i]
        colorvector=  map(lambda x: x + (x[1], ), this_cdict.items())
        colorvector.sort()
        cdict[key] = colorvector
    return LinearSegmentedColormap('colormap',cdict,1024)

def rgb(triplet):
    return _HEXDEC[triplet[0:2]]/255., _HEXDEC[triplet[2:4]]/255., _HEXDEC[triplet[4:6]]/255.
      
def circ_ave(wght,a0,a1):
    # determine how to obtain minimal distance between a0 and a1
    # (which gives averaging trajectory)
    # by moving a0 around the circle an integer number of times
    # inspired by http://stackoverflow.com/questions/1416560/hsl-interpolation
    # Available under Creative Commons Attribution License
    d0=abs(a1-a0)
    a0min=a0-1.0
    d1=abs(a0min-a1)
    a0plus=a0+1.0
    d2=abs(a0plus-a1)
    if(d0<d1 and d0<d2):
       return (wght*a0+(1.0-wght)*a1)
    elif(d1<d2):
       return (wght*a0min+(1.0-wght)*a1)%1.0     
    else:
       return (wght*a0plus+(1.0-wght)*a1)%1.0              

def calcluminosity(r,g,b):
    return 0.2126*r+0.7152*g+0.0722*b

def lumfunc(hue,lum,sat,targetlumin):
    import colorsys   
    r,g,b=colorsys.hls_to_rgb(hue,lum,sat)
    lumin=calcluminosity(r,g,b)
    return lumin-targetlumin
    
# Root finding: Ridder's algorithms (stable when root is bracketed)
# Used to calculate the distibution of levels
def ridderlum(hue,sat,targetlumin):   
    def lf(lum):
        return lumfunc(hue,lum,sat,targetlumin)
    return ridder(lf,0.,1.)
                                 
def hlsinterp(listin):
    import colorsys
    rgbs=[rgb(i) for i in listin]
    hsls=[colorsys.rgb_to_hls(i[0],i[1],i[2]) for i in rgbs]
    lumins=[calcluminosity(i[0],i[1],i[2]) for i in rgbs]
    hues=[i[0] for i in hsls]
    lums=[i[1] for i in hsls]
    sats=[i[2] for i in hsls]
    lenn=len(hues)
    hueindex=linspace(0,lenn-1.0-1e-12,255)
    hueinterp=zeros(255)
    for i in range(255):
        ilower=int(hueindex[i])
        iupper=ilower+1
        wght=iupper-hueindex[i]
        # gives divergent colormaps a hue discontinuity in the middle
        if(lumins[ilower]>0.97 or sats[ilower]<0.02):
           hueinterp[i]=hues[iupper]
        if(lumins[iupper]>0.97 or sats[iupper]<0.02):
           hueinterp[i]=hues[ilower]
        else:         
           hueinterp[i]=circ_ave(wght,hues[ilower],hues[iupper])
    luminterpguess=interp(linspace(0,lenn-1,255),range(lenn),lums)   
    satinterp=interp(linspace(0,lenn-1,255),range(lenn),sats)
    luminindex=linspace(0,lenn-1.0-1e-12,255)
    lumininterp=zeros(255)
    for i in range(255):
        ilower=int(luminindex[i])
        iupper=ilower+1
        wght=iupper-luminindex[i]
        # reduce weights close to zero
        if(lumins[ilower]>0.9 and lumins[iupper]<0.9):
           wght=wght**2
        if(lumins[iupper]>0.9 and lumins[ilower]<0.9):
           wght=1.0-(1.0-wght)**2
        lumininterp[i]=wght*lumins[ilower]+(1.0-wght)*lumins[iupper]   
    # adjust luminosity using Ridder's method.
    # idea: luminance matching
    lumout=zeros(255)
    for i in range(255):     
        lumout[i]=ridderlum(hueinterp[i],satinterp[i],lumininterp[i])
    return [colorsys.hls_to_rgb(hueinterp[i],lumout[i],satinterp[i]) for i in range(255)]

# Inspried by http://stackoverflow.com/questions/3279560/invert-colormap-in-matplotlib
# Available under Creative Commons Attribution License
def reverse_colourmap(cmap, name = 'my_cmap_r'):
    """
    In: 
    cmap, name 
    Out:
    my_cmap_r

    Explanation:
    t[0] goes from 0 to 1
    row i:   x  y0  y1 -> t[0] t[1] t[2]
                   /
                  /
    row i+1: x  y0  y1 -> t[n] t[1] t[2]

    so the inverse should do the same:
    row i+1: x  y1  y0 -> 1-t[0] t[2] t[1]
                   /
                  /
    row i:   x  y1  y0 -> 1-t[n] t[2] t[1]
    """        
    reverse = []
    k = []   

    for key in cmap._segmentdata:    
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:                    
            data.append((1-t[0],t[2],t[1]))            
        reverse.append(sorted(data))    

    LinearL = dict(zip(k,reverse))
    my_cmap_r = LinearSegmentedColormap(name, LinearL) 
    return my_cmap_r
    
sb_cmaps={}
sb_cmaps['yrb1']=LinearSegmentedColormap.from_list('sj',hlsinterp(['F3E300','FF8C44','E43B19','7D3D46','143134']))
sb_cmaps['yrb2']=LinearSegmentedColormap.from_list('sj',hlsinterp(['ECF65A','F4BA34','EA4B43','8D4454','3A5178','0E2D37']))
sb_cmaps['silver']=LinearSegmentedColormap.from_list('sj',hlsinterp(['D9D7DF','95A2E4','0089B3','7D651F','6E211C','390E2D']))
sb_cmaps['fullcircle']=LinearSegmentedColormap.from_list('sj',hlsinterp(['FFC637','AFC600','00B3A1','3162E8','77349C','7E2A30','443011']))
sb_cmaps['spec1']=LinearSegmentedColormap.from_list('sj',hlsinterp(['FF7F16','9EA700','00A97C','196AE3','863AAF','933F63','873E3A']))
sb_cmaps['spec2']=LinearSegmentedColormap.from_list('sj',hlsinterp(['FFF9B3','CDE000','00D11C','12909E','1340AA','251F5C']))
sb_cmaps['spec3']=LinearSegmentedColormap.from_list('sj',hlsinterp(['FFBAA3','FF964B','B09C3B','699511','007D53','0B5665','121A68']))
sb_cmaps['spec4']=LinearSegmentedColormap.from_list('sj',hlsinterp(['F0C775','BEBE3F','56B656','3E9089','5E46A4','762742','630002']))
sb_cmaps['golden']=LinearSegmentedColormap.from_list('sj',hlsinterp(['050505','BD0A0D','E23512','EC7018','EAA545','E8CD7B','EEEABB']))
sb_cmaps['br1']=LinearSegmentedColormap.from_list('sj',hlsinterp(['9EDFFF','5DB6F6','4D96E0','808080','C00000','862800','542300']))
sb_cmaps['cloud']=LinearSegmentedColormap.from_list('sj',hlsinterp(['FFFFFF','6a91a0','1a1e8e','2E0F20']))
sb_cmaps['lcloud']=LinearSegmentedColormap.from_list('sj',hlsinterp(['FFFFFF','BFCFD6','8AA3D8','8885AD']))
sb_cmaps['br2']=LinearSegmentedColormap.from_list('sj',hlsinterp(['1916C4','0C44EA','0677F9','4FA2F6','FFFFFF','F77C63','EB4747','C43131','9A2626']))
sb_cmaps['ljet']=cmap_brightened(cm.jet)
sb_cmaps['lspectral']=cmap_brightened(cm.Spectral,0.8)

for key in sb_cmaps.keys():
    sb_cmaps[key+'_r']=reverse_colourmap(sb_cmaps[key])

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import savefig
    import os
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))
    
    # inspired by http://matplotlib.org/examples/color/colormaps_reference.html
    fig,axes = plt.subplots(nrows=len(sb_cmaps))
    fig.subplots_adjust(top=0.95, bottom=0.01, left=0.2, right=0.99)
    axes[0].set_title('Experimental colormaps', fontsize=14)

    for ax,key in zip(axes, sb_cmaps.keys()):
        ax.imshow(gradient, aspect='auto', cmap=sb_cmaps[key])
        pos = list(ax.get_position().bounds)
        x_text = pos[0] - 0.01
        y_text = pos[1] + pos[3]/2.
        fig.text(x_text, y_text, key, va='center', ha='right', fontsize=10)

    # Turn off *all* ticks & spines, not just the ones with colormaps.
    for ax in axes:
        ax.set_axis_off()
    
    savefig('ccmaps.png')
    try:
        os.system('convert ccmaps.png -colorspace gray ccmaps_g.png')
    except:
        pass
