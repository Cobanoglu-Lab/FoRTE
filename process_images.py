import sys, os
import re
import multiprocessing as mp
from collections import namedtuple, defaultdict

import pandas as pd
import numpy as np

import scipy
from scipy.ndimage import find_objects

import skimage
from skimage import measure, io, color, data, img_as_float
from skimage.filters import try_all_threshold
from skimage.filters import gaussian, threshold_yen

NUM_CPU = 32
PLOT_IMG = True

def illum_filter(img,gauss_sigma=100):
    if img is None:
        return None
    fltr = gaussian(img,gauss_sigma)
    img = img-fltr
    img[img<0] = 0
    return img

def read_image(path,xy,z,c):
    try:
        if xy == '00':
            img = img_as_float(
                io.imread(path+"z{0:s}c{1:d}.tif".format(z,int(c))))
        else:
            img = img_as_float(
                io.imread(path+"xy{0:s}z{1:s}c{2:d}.tif".format(xy,z,int(c))))
    except (FileNotFoundError, ValueError) as e:
        return None
    try:
        return img.mean(2)
    except IndexError:
        print("ERROR in image: shape:{0:s}, xy:{1:s}, z:{2:s}, c:{3:s}".format( 
				str(img.shape), str(xy), str(z), str(c)))
        return None

def nuclei_area(thr_img, min_area=10, max_area=None):
    if max_area is None:
        max_area = np.prod(thr_img.shape)
    cell_labels = measure.label(thr_img, background=0, connectivity=1)
    tot_area = 0
    img_mask = np.zeros(thr_img.shape)
    for n in find_objects(cell_labels):
        block = thr_img[n]
        area = (block>0).sum()
        if area >= min_area and area <= max_area:
            tot_area += area
            img_mask[n][block>0]=1
    return tot_area, img_mask

def get_img_thr(img, ratio=10): 
    thr = threshold_yen(img)
    background = img[img<thr]
    if thr < ratio*np.mean(background):
        print("thr:{0:.5f}, ratio*mean:{1:.5f}".format(thr,ratio*img.mean()))
        return 0
    return thr

def get_file(input_tuple):
    path,xy,z,c = input_tuple
    #z = str(z).rjust(3,'0') TODO: MAKE THIS DYNAMIC
    img = illum_filter(read_image(path, xy, z, c))
    return img

def get_thr(path, xy, c, z_min=84, z_range=32, z_len=3):
    with mp.Pool(NUM_CPU) as p:
        images = p.map(get_file,[(path, xy, str(z).rjust(z_len,'0'), c) for z in range(z_min,z_min+z_range)])
    print([(path, xy, str(z).rjust(z_len,'0'), c) for z in range(z_min,z_min+z_range)])
    images = [img for img in images if img is not None]
    print("xy:{0:s}, c:{1:s}, #images={2:d}".format(xy,str(c),len(images)))
    if len(images)==0:
        return 0
    thresholds = []
    with mp.Pool(NUM_CPU) as p:
        thresholds = p.map(get_img_thr, images)
    thresholds = [thr for thr in thresholds if thr>0]
    if len(thresholds)==0:
        return 0
    sel_ind = np.argmin(thresholds)
    sel_thr = thresholds[sel_ind]
    pix_area, _ = nuclei_area(images[sel_ind]>=sel_thr)
    if pix_area > 0:
        return sel_thr
    else:
        return 0

def get_area(input_tuple):
    "thr_d, path, xy, z = input_tuple"
    thr_d, path, xy, z = input_tuple
    images = {}
    for c in name.keys():
        img = read_image(path,xy,z,c)
        if img is None: return None
        else:
            images[c] = illum_filter(img)
    results = {}
    for ind, chn in name.items():
        tot_area, nuc_mask = nuclei_area(images[ind]>=thr_d[ind])
        if PLOT_IMG:
            plt.figure(figsize=(9,9))
            plt.imshow(np.maximum(nuc_mask*np.max(images[ind]),images[ind]))
            plt.savefig(output+os.sep+"xy{0:s}z{1:s}c{2:s}.png".format(xy, z, chn))
        results[chn] = PixelRes(tot_area, thr_d[ind], nuc_mask, images[ind])
    pixel_counts = {'xy':xy,'z':z}
    for i,chn1 in enumerate(names):
        pixel_counts[chn1] = results[chn1].tot_pixels
        for j,chn2 in enumerate(names[i+1:]):
            pixel_counts[chn1+'_and_'+chn2] = (results[chn1].nuc_mask+\
                                               results[chn2].nuc_mask==2).sum() 
    return pixel_counts

def get_available_z(path,xy):
    z_dict = defaultdict(int)
    z_set = set()
    for im_f in os.listdir(path):
        xy_curr, z_curr, c_curr = check_pattern(im_f)
        if xy == xy_curr:
            z_dict[z_curr]+=1
            if z_dict[z_curr] == len(names):
                z_set.add(z_curr)
    return sorted(list(z_set))

def get_available_files(path):
    files = {}
    for im_f in os.listdir(path):    
        xy_curr, z_curr, c_curr = check_pattern(im_f)
        if xy_curr not in files:
            files[xy_curr] = get_available_z(path, xy_curr)
    return files

def check_pattern(im_f):
    if im_f.startswith('xy'):
        pattern = re.compile('xy([0-9]+)z([0-9]+)c([0-9]+)\.tif')
        return pattern.match(im_f).groups()
    else:
        pattern = re.compile('z([0-9]+)c([0-9]+)\.tif')
        z_curr, c_curr = pattern.match(im_f).groups()
        xy_curr = '00'
        return xy_curr, z_curr, c_curr

path = sys.argv[1]
if not(os.path.isdir(path)):
    print("Input folder {0:s} does not exist. Terminating.".format(path))
    exit()
output = sys.argv[2]
dirs = output.split(os.sep)
if len(dirs) > 0:
    for d in dirs:
        is not(os.path.isdir(d)):
            os.mkdir(d)
print(path)

#name = {1:'Nuc', 2:'Apop', 3:'EtHd'}; names = list(name.values())
name = {}
for ind, chn_name in enumerate(sys.argv[3:]):
    name[ind+1] = chn_name
names = list(name.values())
print(name,names)

PixelRes = namedtuple('PixelResults', ['tot_pixels','threshold','nuc_mask','image'])

final_res = defaultdict(list)
files = get_available_files(path)
for xy, z_avail in files.items():
    print('xy:{:s}. Thresholding...'.format(xy))
    thr_d = {}
    z_range = NUM_CPU
    z_min = int(len(z_avail)/2-z_range/2)
    if z_min < 1:
        z_min = 1
    z_len = len(z_avail[0])
    for ind in name.keys():
        thr_d[ind] = get_thr(path,xy,ind,z_min,z_range,z_len)
        #thr_d[ind] = get_thr(path,xy,ind,1,20,2)
    print(thr_d)
    print('Threshold complete, running through slices')
    try:
        with mp.Pool(NUM_CPU) as p:
            for r in p.map(get_area, [(thr_d, path, xy, z) for z in z_avail]):
                if r is not None:
                    for k,v in r.items():
                        final_res[k].append(v)
    except IndexError:
        print("IndexError:", thr_d, path, xy, z_avail)
pd.DataFrame(final_res).to_excel(output+'.xlsx',index=False)

