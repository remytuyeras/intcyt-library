from intcyt import *
import numpy as np
import matplotlib.pyplot as plt
import hdbscan
import gzip
import matplotlib.ticker as ticker
#------------------------------------------------------------------------------
ary = 20
fload = gzip.open("result-load/load_initial.gz","r")
initial_organelles = usf.get_memory(fload,ary,[0])
fload.close()
image_size = [28,28,1]
#------------------------------------------------------------------------------
data = initial_organelles[0]
#------------------------------------------------------------------------------
for i in range(len(data)):
  usf.print_data(data[i],image_size,sys.stdout)
#------------------------------------------------------------------------------
clusterer = hdbscan.HDBSCAN(min_cluster_size=2,gen_min_span_tree=True).fit(data)
print clusterer.labels_
#------------------------------------------------------------------------------
l = set(clusterer.labels_)
table = list()
for i in range(len(l)):
  table.append([])
max_length = [0]*len(l)

for i in range(len(clusterer.labels_)):
  table[clusterer.labels_[i]+1].append(data[i])
  max_length[clusterer.labels_[i]+1] += 1

m = max(max_length)
final_image = list()
for i in range(len(table)):
  image = usf.make_rgb_panel(table[i],image_size)
  for j in range(len(image)):
    image[j] = image[j] + [[1,1,1]] + [[500,500,500]]*(28+1)*(m-max_length[i])
  final_image = final_image + image
  if i != len(table)-1:
    final_image = final_image + [[[1,1,1]]*len(image[0])]
 
fig, ax = plt.subplots(1)
#remove axis
#ax.axis("off")
ax.xaxis.set_major_locator(ticker.NullLocator())
numlabel = map(lambda x :14+(x+1)*29,list(sorted(l)))
ax.yaxis.set_major_locator(ticker.FixedLocator(numlabel))
labels = map(lambda x :x,list(sorted(l)))
ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
plt.title("HDBSCAN.labels_")
#plt.xlabel("images")
plt.ylabel("clusters")
#display the image
ax.imshow(final_image)
plt.show()
#------------------------------------------------------------------------------
#clusterer.single_linkage_tree_.plot(cmap="viridis", colorbar=True)
#plt.show()
#------------------------------------------------------------------------------
#clusterer.condensed_tree_.plot()
#plt.show()
