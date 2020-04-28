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
  #------------------
  dic_svhn = sio.loadmat("data/train_32x32.mat")
  #Access to the dictionary
  dicim = dic_svhn["X"] #shape: 32, 32, 3, 73257
  dictr = dic_svhn["y"] #shape: 73257, 1
  #------------------
  image_size = [32,32,3] #height, width, depth
#------------------
#Open the MNIST dataset
elif len(sys.argv) > 1 and sys.argv[1] == "MNIST":
  #------------------
  fim = gzip.open("data/train-images-idx3-ubyte.gz","r")
  ftr = gzip.open("data/train-labels-idx1-ubyte.gz","r")
  ftr.read(8)
  fim.read(16)
  #------------------
  image_size = [28,28,1] #height, width, depth
#------------------
#Open the fashion-MNIST dataset
elif len(sys.argv) > 1 and sys.argv[1] == "fashion-MNIST":
  #------------------
  fim = gzip.open("data/t10k-images-idx3-ubyte.gz","r")
  ftr = gzip.open("data/t10k-labels-idx1-ubyte.gz","r")
  ftr.read(8)
  fim.read(16)
  #------------------
  image_size = [28,28,1] #height, width, depth
  categories = ["t-shirt","trousers","pullover", \
                "dress","jacket","sandal","shirt", \
                "sneaker","bag","ankle-boot"]
#------------------
elif len(sys.argv) > 1 and sys.argv[1] == "BINARY":
  #------------------
  fim = gzip.open("data/binary_training.gz","r")
  listim = np.array(usf.get_memory(fim,1,[0,1]))
  #------------------
  image_size = [28,28,1] #height, width, depth
  categories = ["image A","image B"]
#------------------
else:
  exit(0)
#-------------------------------------------------------------------------------------------
#Initializing the super cell
#-------------------------------------------------------------------------------------------
#Main parameters defining the initial super cell
#------------------
dim = image_size[0] * image_size[1] * image_size[2]
operad = Operad(dim)
ary = 20
#------------------
#Initializing the organelles and the cytosolic content (either new ones or saved ones)
#------------------
if len(sys.argv) > 2 and sys.argv[2][:5] == "-load":
  #Load a previously generated tree (height 0)
  load_index = 0
  fload_name = "result-load/load_initial.gz"
  #------------------
  #To load the index of the last super cell saved in fload
  #------------------
  if len(sys.argv) > 3 and sys.argv[3] == "-final":
    fload = gzip.open(fload_name,"r")
    load_index = usf.get_last_cycle(fload,ary)
    fload.close()
  #------------------
  fload = gzip.open(fload_name,"r")
  initial_organelles = usf.get_memory(fload,ary,[load_index])
  fload.close()
  c = Cell(dim,0,[0]*dim,initial_organelles[0])
  sc = SuperCell(c)
else:
  #Generate a random tree of height 0
  sc = operad.generate(levels = 0, arity = ary, key = lambda : random.randint(0,30))
#------------------
#The initial sate of the super cell is saved in the memory
#------------------
fsave_roo = gzip.open("save_roo.gz","w")
fsave_org = gzip.open("save_org.gz","w")
fsave_tre = gzip.open("save_tre.gz","w")
usf.print_root(sc,image_size,fsave_roo,option="save")
usf.print_organelles(sc,image_size,fsave_org,option="save")
usf.print_tree(sc,image_size,fsave_tre,option="save")
fsave_roo.flush() 
fsave_org.flush()
fsave_tre.flush()
#-------------------------------------------------------------------------------------------
#Parameters for self-supervision challenges
#-------------------------------------------------------------------------------------------
selfsup = list()
if len(sys.argv) > 2 and sys.argv[2] == "-load-selfsup-right":
  selfsup.append([image_size[1]/2, #middle of image
                  image_size[1]-1, #right of image
                  image_size[1]])  #width
elif len(sys.argv) > 2 and sys.argv[2] == "-load-selfsup-left":
  selfsup.append([0,                  #left of image
                  image_size[1]/2-1, #middle of image
                  image_size[1]])    #width
#-------------------------------------------------------------------------------------------
#Running the Learning algorithms
#-------------------------------------------------------------------------------------------
for i in range(10000):
  #------------------
  #Get labels and inputs
  #------------------
  if sys.argv[1] == "SVHN":
    label = ", ".join(map(str,dictr[i]))
    inputs = dicim[:,:,:,i].reshape(dim)
  #------------------
  if sys.argv[1] == "MNIST":
    buf_lab = ftr.read(1)
    buf_inp = fim.read(dim)
    label = ", ".join(map(str,np.frombuffer(buf_lab, dtype=np.uint8).astype(np.int64)))
    inputs = np.frombuffer(buf_inp, dtype=np.uint8).astype(np.float32)
  #------------------
  if sys.argv[1] == "fashion-MNIST":
    buf_lab = ftr.read(1)
    buf_inp = fim.read(dim)
    label = map(lambda x: categories[x], np.frombuffer(buf_lab, dtype=np.uint8).astype(np.int64))
    inputs = np.frombuffer(buf_inp, dtype=np.uint8).astype(np.float32)
  #------------------
  if sys.argv[1] == "BINARY":
    label = categories[i%2]
    inputs = listim[i%2][0]
  #------------------  
  #Normalization of the input vector
  #------------------
  inputs = 0.0001*inputs/max(inputs)
  vector = inputs.tolist()
  #------------------
  #The algorithm
  #------------------  
  debug_time.set("Learning data labeled as " + label)
  #------------------
  if sys.argv[1] in ["SVHN","MNIST","fashion-MNIST"]:
    start = 4
    epoch = 4
    fission_events = [0]
    fusion_events = [2]
    compose_events = [1,3]
    #------------------
    filtering = [1.5, #fission
                 1.5]  #fusion
  #------------------ 
  if sys.argv[1] in ["BINARY"]:
    start = 16
    epoch = 4
    fission_events = (i % 16 in [0]) * range(epoch)
    fusion_events = (i % 16 in [i for i in range(1,16) if i%2 == 0]) * range(epoch)
    compose_events = (i % 16 in [i for i in range(1,16) if i%2 == 1]) * range(epoch)
    #------------------
    filtering = [2, #fission
                 .7]  #fusion
  #------------------
  events = [start,epoch,fission_events,fusion_events,compose_events]
  #------------------
  for k in range(epoch):
    debug_time.set("TREE")
    sc.stdout(vector)
    #------------------
    #brightness: intensity levels (length m)
    #profiles: intervals for intensity counts (length n*m)
    #scores: agreement scores for a given profile (length n)
    if sys.argv[1] == "SVHN":
      #------------------
      brightness = [.1,.25,.5,.75,.9]
      profiles = [[(0,.6),(0,.375),(0,.15),(0,.05),(0,.005)]] #expert
      scores = [.82] #expert
      E = 13
      F = 25
    #------------------
    if sys.argv[1] == "MNIST":
      #------------------
      brightness = [.1,.25,.5,.75,.9]
      #the average profile for MNIST is [(0,.175),(0,.156),(0,.133),(0,.106),(0,.088)]
      profiles = [[(0,.6),(0,.25),(0,1),(0,1),(0,.1)], #student
                  [(0,.2),(0,.15),(0,.1),(0,.1),(0,.1)]]  #expert
      scores = [.7, #student
                .8] #expert
      E = 12.5
      F = 20
    #------------------
    '''
    #for self-supervision
    if sys.argv[1] == "MNIST":
      #------------------
      brightness = [.1,.25,.5,.75,.9]
      #the average profile for MNIST is [(0,.175),(0,.156),(0,.133),(0,.106),(0,.088)]
      profiles = [[(0,.45),(0,.25),(0,.1),(0,.05),(0,.01)], #student
                  [(0,.3),(0,.3),(0,.1),(0,.05),(0,.01)]]  #expert
      scores = [.7, #student
                .8] #expert
      E = 14
      F = 25
    '''
    #------------------
    if sys.argv[1] == "fashion-MNIST":
      #------------------
      brightness = [.1,.25,.5,.75,.9]
      profiles = [[(0,.7),(0,.5),(0,.3),(0,.1),(0,.01)],  #student
                  [(0,.5),(0,.3),(0,.2),(0,.15),(0,.1)]]  #expert
      scores = [.5, #student
                .8] #expert
      E = 13.5
      F = 20
    #------------------
    '''
    #for self-supervision
    if sys.argv[1] == "fashion-MNIST":
      #------------------
      brightness = [.1,.25,.5,.75,.9]
      profiles = [[(0,.7),(0,.5),(0,.3),(0,.1),(0,.01)],  #student
                  [(0,.5),(0,.3),(0,.2),(0,.15),(0,.1)]]  #expert
      scores = [.5, #student
                .8] #expert
      E = 13.5
      F = 22.5
    '''
    #------------------
    if sys.argv[1] == "BINARY":
      #------------------
      brightness = [.1,.25,.5,.75,.9]
      #the average profile for MNIST is [(0,.175),(0,.156),(0,.133),(0,.106),(0,.088)]
      profiles = [[(0,.1),(0,.1),(0,.1),(0,.1),(0,.1)], #student
                  [(0,.1),(0,.1),(0,.1),(0,.1),(0,.1)],
                  [(0,.1),(0,.1),(0,.1),(0,.1),(0,.1)],
                  [(0,.1),(0,.1),(0,.1),(0,.1),(0,.1)],
                  [(0,.1),(0,.1),(0,.1),(0,.1),(0,.1)],
                  [(0,.1),(0,.1),(0,.1),(0,.1),(0,.1)],
                  [(0,.1),(0,.1),(0,.1),(0,.1),(0,.1)],
                  [(0,.1),(0,.1),(0,.1),(0,.1),(0,.1)],
                  [(0,.1),(0,.1),(0,.1),(0,.1),(0,.1)],
                  [(0,.1),(0,.1),(0,.1),(0,.1),(0,.1)]]  #expert
      scores = [.1,
                .2,
                .3,
                .4, #student
                .5,
                .6,
                .7,
                .8,
                .9,
                .95] #expert
      E = 13
      F = 35
    #------------------
    gamma_parameter = usf.gamma(E,F,brightness,profiles,scores,*selfsup)
    intcyt(operad,sc,epoch*i+k,events,vector,gamma_parameter,filtering)
    #------------------
    #The sate of the super cell is saved in the memory
    #------------------
    #usf.print_root(sc,image_size,sys.stdout,option="display")
    #usf.print_organelles(sc,image_size,sys.stdout,option="display")
    #usf.print_tree(sc,image_size,sys.stdout,option="display")
    usf.print_root(sc,image_size,fsave_roo,option="save")
    usf.print_organelles(sc,image_size,fsave_org,option="save")
    usf.print_tree(sc,image_size,fsave_tre,option="save")
    fsave_roo.flush() 
    fsave_org.flush()
    fsave_tre.flush()
#------------------
fsave_roo.close() 
fsave_org.close()
fsave_tre.close()
if sys.argv[1] in ["MNIST","fashion-MNIST"]:
  fim.close()
  ftr.close()
if sys.argv[1] in ["BINARY"]:
  fim.close()
