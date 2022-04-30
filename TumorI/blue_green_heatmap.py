import time
import sys
import os
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
import numpy as np

def load_pixel_array(filename,verbose=False):
    im = Image.open(filename)
    ima = np.array(im)   # convert to numpy
    if verbose:
        print(filename, ima.size, ima.shape)
    return ima

def pixel_to_heatmap(green,blue):   # TO DO: this is slow, need a hash function
    bins=[10,20,30,40,50,60,70,80,90,256]
    gbin=None
    bbin=None
    for bin in range(0,10):
        if gbin is None and green<=bins[bin]:
            gbin=bin
        if bbin is None and blue<=bins[bin]:
            bbin=bin
        if gbin is not None and bbin is not None:
            return gbin,bbin
    return gbin,bbin

def accumulate_pixels(imary,verbose=False):  # TO DO: nested for loop is very slow
    heatmap=np.zeros( (10,10), dtype=np.int32)
    if verbose:
        print("accumulate",imary.shape)
    nrows,ncols,nchannel=imary.shape
    for row in range(0,nrows):
        for col in range(0,ncols):
            pixel = imary[row,col]
            #red = pixel[0]
            green = pixel[1]
            blue = pixel[2]
            gbin,bbin = pixel_to_heatmap(green,blue)
            heatmap[gbin,bbin] += 1
    return heatmap

def main(argv):
    filename = argv[0]
    pixels = load_pixel_array(filename)
    hm = accumulate_pixels(pixels)
    print(hm)

if __name__ == "__main__":
    main(sys.argv[1:])
