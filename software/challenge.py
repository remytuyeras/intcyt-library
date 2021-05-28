from intcyt import *
#------------------
#Libraries to open datasets
#------------------
import gzip
import numpy as np
import scipy.io as sio
#-------------------------------------------------------------------------------------------
#Loading the dataset for learning
#-------------------------------------------------------------------------------------------
#Open the SVHN dataset
if len(sys.argv) > 1 and sys.argv[1] == "SVHN":
  dic_svhn = sio.loadmat("data/train_32x32.mat")
  #Access to the dictionary
  dicim = dic_svhn["X"] #shape: 32, 32, 3, 73257
  dictr = dic_svhn["y"] #shape: 73257, 1
  image_size = [32,32,3] #height, width, depth
#------------------
#Open the MNIST dataset
elif len(sys.argv) > 1 and sys.argv[1] == "MNIST":
  fim = gzip.open("data/train-images-idx3-ubyte.gz","r")
  ftr = gzip.open("data/train-labels-idx1-ubyte.gz","r")
  ftr.read(8)
  fim.read(16)
  image_size = [28,28,1] #height, width, depth
#------------------
#Open the fashion-MNIST dataset
elif len(sys.argv) > 1 and sys.argv[1] == "fashion-MNIST":
  fim = gzip.open("data/t10k-images-idx3-ubyte.gz","r")
  ftr = gzip.open("data/t10k-labels-idx1-ubyte.gz","r")
  ftr.read(8)
  fim.read(16)
  image_size = [28,28,1] #height, width, depth
  categories = ["t-shirt","trousers","pullover", \
                "dress","jacket","sandal","shirt", \
                "sneaker","bag","ankle-boot"]
else:
  exit(0)
#-------------------------------------------------------------------------------------------
#Parameters  
#-------------------------------------------------------------------------------------------
dim = image_size[0] * image_size[1] * image_size[2]
flatten_width = image_size[1] * image_size[2]
try:
  ary = int(sys.argv[2])
except:
  print "Error:challenge.py:parameter [ary] is missing its value"
  exit() 
#-------------------------------------------------------------------------------------------
#Randomly select data in the dataset.
#-------------------------------------------------------------------------------------------
init = list()
while len(init) < ary:
  n = random.randint(0,999)
  if not(n in init): 
    init.append(n)
init = sorted(init)
print init
#-------------------------------------------------------------------------------------------
#Save the randomly picked data in a file
#-------------------------------------------------------------------------------------------
fself = gzip.open("result-load/load_initial.gz","w")
#-------------------------------------------------------------------------------------------
random_init = list()
if len(sys.argv) > 3 and sys.argv[3][:6] == "-right":
  random_init.extend([0, flatten_width/2-1])
elif len(sys.argv) > 3 and sys.argv[3][:5] == "-left":
  random_init.extend([flatten_width/2, flatten_width-1])
#-------------------------------------------------------------------------------------------
for i in range(1000):
  #------------------
  #Get labels and inputs
  #------------------
  if sys.argv[1] == "SVHN":
    label = dictr[i]
    inputs = dicim[:,:,:,i].reshape(dim)
  #------------------
  if sys.argv[1] == "MNIST":
    buf_lab = ftr.read(1)
    buf_inp = fim.read(dim)
    label = np.frombuffer(buf_lab, dtype=np.uint8).astype(np.int64)
    inputs = np.frombuffer(buf_inp, dtype=np.uint8).astype(np.float32)
  #------------------
  if sys.argv[1] == "fashion-MNIST":
    buf_lab = ftr.read(1)
    buf_inp = fim.read(dim)
    label = map(lambda x: categories[x], np.frombuffer(buf_lab, dtype=np.uint8).astype(np.int64))
    inputs = np.frombuffer(buf_inp, dtype=np.uint8).astype(np.float32) 
  #------------------  
  if i in init:
    #------------------  
    debug_time.set("Data picked: " + ", ".join(map(str,label)))
    #------------------
    if len(random_init) == 2:
      for k in range(len(inputs)):
        if random_init[0] <= k % flatten_width <= random_init[1]:
          inputs[k] = random.randint(0,250)
        elif len(sys.argv) > 3 and sys.argv[3] in ["-right-noisy","-left-noisy"] and inputs[k] == 0:
          if len(sys.argv) == 4:
            inputs[k] = random.randint(0,250)
          elif len(sys.argv) > 5:
            if 1 <= random.randint(0,int(sys.argv[5])) <= int(sys.argv[4]):
              inputs[k] = random.randint(0,250)
    #------------------
    linputs = inputs.tolist()
    usf.print_data(linputs,image_size,sys.stdout,option = "display")
    usf.print_data(linputs,image_size,fself,option = "save")
    fself.write("\n")
#-------------------------------------------------------------------------------------------
fself.close()
if sys.argv[1] in ["MNIST","fashion-MNIST"]:
  fim.close()
  ftr.close()
