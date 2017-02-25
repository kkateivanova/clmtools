from os.path import basename, isfile, join

def valid_file(d, f):
  fn = join(data_folder, d, f)
  return isfile(fn) and not (f.startswith(".") or f.lower().endswith(".md"))

# mirror scheme
mirror_map = [6,5,4,3,2,1,0,11,10,9,8,7,14,13,12,17,16,15]

# path for drawing face
path = {\
   'normal' : [\
      [0,1,2,3,4,5,6,7,8,9,10,11,0,12,13,14,6,15,16,17,0],
   ], \
   'vertices' : [\
      [0,1,17,0],\
      [1,2,17,1],\
      [2,3,17,2],\
      [3,17,16,3],\
      [3,15,16,3],\
      [3,4,15,3],\
      [4,5,15,4],\
      [5,6,15,5],\
      [6,7,14,6],\
      [7,8,14,7],\
      [8,13,14,8],\
      [8,9,13,8],\
      [9,10,13,9],\
      [10,12,13,10],\
      [10,11,12,10],\
      [0,11,12,0],\
   ]\
}

# list of new positions of array 1
num_patches = 18

# wanted width of pdm
# a width of 100 will give ocular distance of approximately ?
#modelwidth = 400
modelwidth = 65 # default for creating pdm
#modelwidth = 40 # default for creating detection filters

# wanted patchsize, must be odd
patch_size = 11 # default for creating SVM filters
#patch_size = 16 # default for creating MOSSE filters

# raw image folder
data_folder = "./data/"
images = "./data/images/"
annotations = "./data/annotations.csv"
