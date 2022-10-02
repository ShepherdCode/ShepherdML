'''
Make small image files containing patches of large image files.
'''

from PIL import Image
Image.MAX_IMAGE_PIXELS = None
import numpy as np

infile = 'I1.DF1.09.tif'
im = Image.open(infile)
ima = np.array(im)   # convert to numpy
patch1 = ima[3662:3662+863,3969:3969+863,]
im3 = Image.fromarray(patch1)
im3.save('im3.tif')
# Opened with gwenview, I saw the zoomed in region
