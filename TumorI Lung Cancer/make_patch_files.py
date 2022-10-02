'''
Make small image files containing patches of large image files.
'''

from PIL import Image
Image.MAX_IMAGE_PIXELS = None
import numpy as np
import json

infile = '../raw/I1.DF1.09.tif'
outprefix = 'I1.DF1.09.'
coordsfile = '../labeled/I1.DF1.09.json'

im = Image.open(infile)
ima = np.array(im)   # convert to numpy

coordshandle = open(coordsfile)
coords = json.load(coordshandle)

num_patches = len ( coords[0]['annotations'] )
print (num_patches, "patches")

patch_size=863  # read this from file?

# Assume the coords file annotates only one image.
patch_num=0
for patch in coords[0]['annotations']:
    patch_num += 1
    outfile = '%s%03d.tif'%(outprefix,patch_num)
    x=int(patch['coordinates']['x'])
    y=int(patch['coordinates']['y'])
    print ("x y=",x,y,outfile)
    patch1 = ima[x:x+patch_size,y:y+patch_size,]
    im1 = Image.fromarray(patch1)
    im1.save(outfile)
