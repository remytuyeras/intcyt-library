from intcyt import *
import gzip
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
#-------------------------------------------------------------------------------------------
def set_legends(fig,ax,image_size,query,ysubtitle,title):
  #remove frame
  bottom_side = ax.spines["bottom"]
  bottom_side.set_visible(False) 
  left_side = ax.spines["left"]
  left_side.set_visible(False)
  right_side = ax.spines["right"]
  right_side.set_visible(False)
  top_side = ax.spines["top"]
  top_side.set_visible(False)
  #display tickers
  ax.xaxis.set_major_locator(ticker.NullLocator())
  numlabel = map(lambda x :image_size[0]/2+(x)*(image_size[0]+1),range(len(query)))
  ax.yaxis.set_major_locator(ticker.FixedLocator(numlabel))
  labels = map(lambda x :x,list(sorted(query)))
  ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
  #display legends
  plt.suptitle(title,fontsize=15)
  plt.title("Memory",loc='right',fontsize=13)
  if ysubtitle == "":
    plt.ylabel("Inputs",fontsize=15)
    #left_side.set_visible(True)
  else:
    plt.title("Inputs",loc='left',fontsize=13)
    plt.ylabel(ysubtitle,fontsize=15)
#-------------------------------------------------------------------------------------------
def visualize(fload,inputs,batch,query,image_size,title):
  memory = usf.get_memory(fload,batch,query)
  image = usf.make_rgb_grid(memory,image_size)
  extension = usf.make_rgb_grid(map(lambda x: [x],inputs),image_size,grayscale = True)
  image = usf.join_images(extension,image,image_size[1])
  fig, ax = plt.subplots(1)
  #display labels and legends
  set_legends(fig,ax,image_size,query,"Learning cycles",title)
  #display the image
  ax.imshow(image)
  plt.show()
#-------------------------------------------------------------------------------------------
def visualize_tree(fload,inputs,query,image_size,title):
  trees = usf.get_trees(fload,query)
  memory = list()
  tree = list()
  for i in range(len(trees)):
    if i != len(trees)-1:
      memory.append(trees[i][0])
    else:
      tree = trees[i]
  image1 = usf.make_rgb_grid(memory,image_size)
  image2 = usf.make_rgb_tree(tree,image_size)
  image = image1 + image2
  extension = usf.make_rgb_grid(map(lambda x: [x],inputs),image_size,grayscale = True)
  image = usf.join_images(extension,image,image_size[1])
  fig, ax = plt.subplots(1)
  #display labels and legends
  set_legends(fig,ax,image_size,query,"Learning cycles",title)
  #display the image
  ax.imshow(image)
  plt.show()
#-------------------------------------------------------------------------------------------
def visualize_forest(fload,inputs,query,image_size,title):
  trees = usf.get_trees(fload,query)
  image = list()
  panel_width = 0
  for i in range(len(trees)):
    tree_image = usf.make_rgb_tree(trees[i],image_size)
    if i == 0:
      panel_width = len(tree_image[0])
    extension = usf.make_rgb_panel([inputs[i]],image_size,grayscale = True)
    tree_image = usf.join_images(extension,tree_image,image_size[1],explicit = 2,center = True)
    image = image + tree_image
    if i != len(trees)-1:
      arrow_image = usf.draw_arrow(panel_width,length=6)
      extension = [[[500,500,500]]*image_size[1]]
      arrow_image = usf.join_images(extension,arrow_image,image_size[1],explicit = 0)
      image = image + arrow_image
  fig, ax = plt.subplots(1)
  #display legends
  set_legends(fig,ax,image_size,[],"",title)
  #display the image
  ax.imshow(image)
  plt.show()
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
if sys.argv[2][0] != "-":
  r = 2
else:
  r = 3
#-------------------------------------------------------------------------------------------
query = list()
for i in range((len(sys.argv)-2)/2):
  query = query + range(int(sys.argv[2*i+r]),int(sys.argv[2*i+r+1])) 
#-------------------------------------------------------------------------------------------
#To visualize learning from SVHN dataset
if sys.argv[1] == "SVHN":
  image_size = [32,32,3] #rgb
  batch_roo = 1
  batch_org = 20
  epoch = 4
  
  dim = image_size[0] * image_size[1] * image_size[2]
  dic_svhn = sio.loadmat("data/train_32x32.mat")
  dicim = dic_svhn["X"]
    
  query_by_epoch = map(lambda x : x/epoch,query)
  collect_inputs = list()
  for i in range(max(query)+1):
    inputs = dicim[:,:,:,i].reshape(dim)
    if i in query_by_epoch:
      collect_inputs.append([i,inputs.tolist()])
  
  inputs = list()
  for i in range(len(query_by_epoch)):
    for j in range(len(collect_inputs)):
      if query_by_epoch[i] == collect_inputs[j][0]:
        inputs.append(collect_inputs[j][1])
#-------------------------------------------------------------------------------------------
#To visualize learning from MNIST dataset
if sys.argv[1] == "MNIST":
  image_size = [28,28,1] #bicolor
  batch_roo = 1
  batch_org = 20
  epoch = 4
  
  dim = image_size[0] * image_size[1] * image_size[2]
  fim = gzip.open("data/train-images-idx3-ubyte.gz","r")
  fim.read(16)
  
  query_by_epoch = map(lambda x : x/epoch,query)
  collect_inputs = list()
  for i in range(max(query)+1):
    buf_inp = fim.read(dim)
    inputs = np.frombuffer(buf_inp, dtype=np.uint8).astype(np.float32)
    if i in query_by_epoch:
      collect_inputs.append([i,inputs.tolist()])
  
  inputs = list()
  for i in range(len(query_by_epoch)):
    for j in range(len(collect_inputs)):
      if query_by_epoch[i] == collect_inputs[j][0]:
        inputs.append(collect_inputs[j][1])
#-------------------------------------------------------------------------------------------
if sys.argv[1] == "fashion-MNIST":
  image_size = [28,28,1] #bicolor
  batch_roo = 1
  batch_org = 20
  epoch = 4
  
  dim = image_size[0] * image_size[1] * image_size[2]
  fim = gzip.open("data/t10k-images-idx3-ubyte.gz","r")
  fim.read(16)
  
  query_by_epoch = map(lambda x : x/epoch,query)
  collect_inputs = list()
  for i in range(max(query)+1):
    buf_inp = fim.read(dim)
    inputs = np.frombuffer(buf_inp, dtype=np.uint8).astype(np.float32)
    if i in query_by_epoch:
      collect_inputs.append([i,inputs.tolist()])
  
  inputs = list()
  for i in range(len(query_by_epoch)):
    for j in range(len(collect_inputs)):
      if query_by_epoch[i] == collect_inputs[j][0]:
        inputs.append(collect_inputs[j][1])
#-------------------------------------------------------------------------------------------
#To visualize learning from MNIST dataset
if sys.argv[1] == "BINARY":
  image_size = [28,28,1] #bicolor
  batch_roo = 1
  batch_org = 20
  epoch = 4
  
  dim = image_size[0] * image_size[1] * image_size[2]
  fim = gzip.open("data/binary_training.gz","r")
  listim = usf.get_memory(fim,1,[0,1])
  
  query_by_epoch = map(lambda x : x/epoch,query)
  collect_inputs = list()
  for i in range(max(query)+1):
    inputs = listim[i%2][0]
    if i in query_by_epoch:
      collect_inputs.append([i,inputs])
  
  inputs = list()
  for i in range(len(query_by_epoch)):
    for j in range(len(collect_inputs)):
      if query_by_epoch[i] == collect_inputs[j][0]:
        inputs.append(collect_inputs[j][1])
#-------------------------------------------------------------------------------------------
if sys.argv[2][0] != "-":
  fload_roo = gzip.open("save_roo.gz","r")
  fload_org = gzip.open("save_org.gz","r")
  #--------------
  visualize(fload_roo,inputs,batch_roo,query,image_size,"Evolution of the root")
  visualize(fload_org,inputs,batch_org,query,image_size,"Evolution of the "+str(batch_org)+" organelles")
  #--------------
  fload_roo.close()
  fload_org.close()
#-------------------------------------------------------------------------------------------
elif sys.argv[2] == "-tree":
  fload_tre = gzip.open("save_tre.gz","r")
  #--------------
  visualize_tree(fload_tre,inputs,query,image_size,"Evolution of the "+str(batch_org)+" organelles + last super cell")
  #--------------
  fload_tre.close()
#-------------------------------------------------------------------------------------------
elif sys.argv[2] == "-forest":
  fload_tre = gzip.open("save_tre.gz","r")
  #--------------
  visualize_forest(fload_tre,inputs,query,image_size,"Evolution of the super cell")
  #--------------
  fload_tre.close()
#-------------------------------------------------------------------------------------------
