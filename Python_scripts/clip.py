from os import listdir
from os.path import isdir, isfile, join
from PIL import Image
import glob
import os.path

for infile in glob.glob("*.png"):
    file, ext = os.path.splitext(infile)
    if file.endswith('_clipped'): #don't clip a clippa
        continue
    im = Image.open(infile)
    im2 = im.crop((46, 88, 953, 991)) #left, upper, right, lower
    im2.save(file+'_clipped.png')
