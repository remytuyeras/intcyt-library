from intcyt import *
import gzip
import math
import scipy.stats as scstat
import matplotlib.pyplot as plt
#---------------------------------------------------------------------------------------
#This function get various types of information from the DREAM3 files
#title= contains the titles indidicating the fields of the columns (ID, gene name, times)
#genes = contains the gene names
#expressions = contains the expression values from a certain index start_expr 
#(start_expr is not the same in every file: see the class Data below)
def get_file_structure(fdoc,start_expr):
  title = list()
  genes = list()
  expressions = list()
  for i,s in enumerate(fdoc):
    l = s.rstrip().split("\t")
    if i == 0:
      title = l
    else:
      genes.append(l[0])
      expr = l[start_expr:]
      new_expr = list()
      for j in range(len(expr)):
        try:
          new_expr.append(float(expr[j]))
        except:
          new_expr.append(expr[j])
      expressions.append(new_expr)
  return title, genes, expressions
#---------------------------------------------------------------------------------------
#Get the gene expression of the genes that are not part of the targets/golden standard 
#for the DREAM3 challenge
def get_complementary_genes(genes,target):
  antitarget = list()
  for i in range(len(genes)):
    if not(genes[i] in target):
      antitarget.append(genes[i])
  return antitarget
#---------------------------------------------------------------------------------------
#Organize the data of the DREAM3 challenge
class Data(object):
  def __init__(self,fdata,fchal,fansw):
    #Get the different fields of the DREAM3 files
    #self.expressions contain the training data
    self.title,   self.genes,   self.expressions  = get_file_structure(fdata,3)
    #self.targets contain the name of the genes to impute
    _,            self.targets,  _                = get_file_structure(fchal,1)
    #self.answer contains the ranked imputation (golden standard values for the challenge)
    _,            _,            self.answer       = get_file_structure(fansw,2)
    #Gene expressions that are not part of the tagerts
    self.antitargets = get_complementary_genes(self.genes,self.targets)
    #Turn characters into integers
    for i in range(len(self.answer)):
      self.answer[i] = map(int,self.answer[i])
#---------------------------------------------------------------------------------------
#Allows to get a subset of a matrix of gene expressions
#data = full matrix of gene expressions
#tagerts = restrict the genes to be selected (can be a list of lists for parallel computation)
#times = restrict the time points to be selected
#strains = restrict the yeast strains to be selected
def get_targets_expressions(targets,data,times,strains,function = lambda x : x):
  expr_list = list()
  for k in range(len(targets)):
    expr = list()
    for j in range(len(targets[k])):
      for i in range(len(data.genes)):
        #print data.genes[i], targets[k][j], data.genes[i] in targets[k][j]
        if data.genes[i] in targets[k][j]:
          e = list()
          for s in range(len(strains)):
            e.extend(map(lambda x : function(data.expressions[i][x+strains[s]*8]),times)) 
          expr.append(e)
    expr_list.append(expr)
  return expr_list
#---------------------------------------------------------------------------------------
#Turn gene expressions into positive values by replacing the expressions with a pair of values
#where the left component contains the absolute value of a negative expression
#and the right compoenent contains the value of a positive expression
def separate_up_down(a_list,function = lambda x: x):
  output = list()
  for i in range(len(a_list)):
    try:
      v = float(a_list[i])
      if v < 0:
        output.extend([function(-v),0])
      else:
        output.extend([0,function(v)])
    except:
      print "[Warning] separate_up_down: "+str(a_list[i])+" is not a float."
      output.extend([function(1),function(1)])
  return output
#---------------------------------------------------------------------------------------
#Turn a tuple (x,y) of positive values back to a real value by using x-y
def reunite_up_down(a_list,function = lambda x: x):
  output = list()
  for i in range(len(a_list)/2):
    output.append(function(a_list[2*i+1]-a_list[2*i]))
  return output
#---------------------------------------------------------------------------------------
#Euclidean distance
def distance(list1,list2):
  d = 0
  for i in range(len(list1)):
    d += (list1[i] - list2[i])**2
  return math.sqrt(d)
#---------------------------------------------------------------------------------------
#Chose the closest vector to a specific vector from a lists of vector using the Euclidean distance
def closest_vector(a_list,list_of_lists):
  distances = list()
  for i in range(len(list_of_lists)):
    d = distance(a_list,list_of_lists[i])
    distances.append((i,d))
  distances = sorted(distances, key = lambda x :x[1])
  return distances[0]
#---------------------------------------------------------------------------------------
#Agreement: see definition in the mathematical development
def agreement(vector1,vector2):
  agree = 0
  norm1 = 0
  norm2 = 0
  for i in range(len(vector1)):
    agree += vector1[i] * vector2[i]
    norm1 += vector1[i]**2
    norm2 += vector2[i]**2
  return agree/math.sqrt(norm1*norm2)
#---------------------------------------------------------------------------------------  
#Merge two images of the same height with a black vertical black line in the middle
def merge(im1,im2):
  im = list()
  for i in range(min(len(im1),len(im2))):
    im.append(im1[i]+[-1]+im2[i])
  return im
#---------------------------------------------------------------------------------------
#Merge two images of different height with a vertical white strip in the middle  
def merge2(im1,im2):
  im = list()
  for i in range(min(len(im1),len(im2))):
    im.append(im1[i]+[[300,300,300]]*10+im2[i])
  if len(im1) < len(im2):
    for i in range(len(im1),len(im2)):
      im.append([[300,300,300]]*len(im1[0])+[[300,300,300]]*10+im2[i])
  if len(im2) < len(im1):
    for i in range(len(im2),len(im1)):
      im.append(im1[i]+[[300,300,300]]*10+[[300,300,300]]*len(im2[0]))
  return im
#--------------------------------------------------------------------------------------- 
def return_answer(matrix,n,m):
  answer = list()
  for i in range(n):
    row = list()
    for k in range(m):
      row.append(0)  
    answer.append(row)
    
  for k in range(m):
    order = list()
    for i in range(n):
      order.append([i,matrix[i][k]])
    order = sorted(order, key = lambda x:-x[1])
    for j in range(n):
      answer[order[j][0]][k] = j+1
      
  return answer
#---------------------------------------------------------------------------------------
def gene_accuracy(gstandard,predicted,n,m):
  correlations = list()
  pvalues = list() 
  for k in range(m):
    column1 = list()
    column2 = list()
    
    for i in range(n):
      column1.append(gstandard[i][k])
      column2.append(predicted[i][k])
    
    #rho, pval = scstat.pearsonr(column1,column2)
    rho, pval = scstat.spearmanr(column1,column2)
    correlations.append(rho)
    pvalues.append(pval/2) #two tails
      
  return correlations, pvalues
#---------------------------------------------------------------------------------------
def time_accuracy(gstandard,predicted,n,m):
  correlations = list()
  pvalues = list()   
  for k in range(n):
    row1 = list()
    row2 = list()
    
    for i in range(m):
      row1.append(gstandard[k][i])
      row2.append(predicted[k][i])

    #rho, pval = scstat.pearsonr(row1,row2)
    rho, pval = scstat.spearmanr(row1,row2)
    correlations.append(rho)
    pvalues.append(pval/2) #two tails
      
  return correlations, pvalues
#--------------------------------------------------------------------------------------- 
def contrast_profiles(brightness,vector):
  maximum = 0
  for u in range(len(vector)):
    if maximum < vector[u]:
      maximum = vector[u]

  counts = [0] * len(brightness)
  total = 0
  for u in range(len(vector)):
    total += 1.
    pixel_u = vector[u]/float(maximum)
    for j in range(len(counts)):
      if pixel_u > brightness[j]:
        counts[j] += 1

  percents = map(lambda x: x/total,counts)
  
  percents = map(lambda x: str(int(10000.*x)/10000.),percents)
  print "[brightness]  " + " ".join(percents)
#---------------------------------------------------------------------------------------
def geometric_mean(pvalues):
  gmean = 1
  for i in range(len(pvalues)):
    gmean = gmean * pvalues[i]
  return gmean **(1/float(len(pvalues)))
#---------------------------------------------------------------------------------------
def arithmetic_mean(coorelations):
  amean = 0
  for i in range(len(coorelations)):
    amean = amean + coorelations[i]
  return  amean/float(len(coorelations))
#---------------------------------------------------------------------------------------
def zero_list(length, val = 0):
  support = list()
  for j in range(length):
    support.append(val)
  return support
#---------------------------------------------------------------------------------------
def cumulate(a_list):
  cumul = list()
  c = 0
  for i in range(len(a_list)):
    c += a_list[i]
    cumul.append(c)
  return cumul
#---------------------------------------------------------------------------------------
def discrete_pdf(vector,subdivisions):
  pdf     = zero_list(subdivisions)
  values  = zero_list(subdivisions)
  step    = 1/float(subdivisions)
  M     = max(vector)
  m     = min(vector)
  line  = lambda s: m*(1-s)+M*s
  for i in range(len(vector)):
    for j in range(subdivisions):
      values[j] = line(1-(j+1)*step)
      if line(1-(j+1)*step) <= vector[i] < line(1-j*step):
        #print "[line]", line(1-(j+1)*step), vector[i], line(1-j*step), M, m , step, [1-(j+1)*step,1-j*step]
        pdf[j] += 1
        break
  return pdf, values
#---------------------------------------------------------------------------------------
def flexible_discrete_pdf(vector,subdivisions):
  pdf     = zero_list(len(subdivisions))
  values  = zero_list(len(subdivisions))
  line  = lambda s: m*(1-s)+M*s
  for i in range(len(vector)):
    for j in range(len(subdivisions)):
      print "-------["+str(i)+","+str(j)+"]--------"
      if subdivisions[j][0] <= vector[i] < subdivisions[j][1]:
        values[j] = subdivisions[j][0]
        pdf[j] += 1
  return pdf, values
#---------------------------------------------------------------------------------------
def data_binning(vector,bins):
  values  = zero_list(len(bins),1)
  for j in range(len(bins)):
    for i in range(len(vector)):
      if bins[j][0] <= i <= bins[j][1]:
        values[j] = min(values[j],vector[i])
        #values[j] += vector[i]
    #values[j] = values[j] / float(bins[j][1]-bins[j][0]+1)
  return values
#---------------------------------------------------------------------------------------
def distribution(vector,subdivisions,option = "pdf"):
  pdf, values = discrete_pdf(vector,subdivisions)  
  
  if option == "cumul":
    cumul = cumulate(pdf)
  else:
    cumul = pdf
    
  colors = list()
  count = 0
  threshold = 0
  for i in range(len(cumul)):
    if cumul[i] == 0:
      colors.append("red")
    else:
      colors.append("blue")
      count += 1
      if count == 1:
        threshold = values[i]
  return (range(subdivisions),cumul,colors,threshold)
#---------------------------------------------------------------------------------------
def plot_bar(values,color,tittle,xlabel,xticks,ylabel,precision = lambda t : int(100*t)/100.0):
  fig, ax = plt.subplots(1)
  barplot = ax.bar(range(len(values)),values, color = color)
  #plt.bar(x, y)
  for bar in barplot:
    yval = bar.get_height()
    #print yval
    plt.text(bar.get_x() + bar.get_width()/2.0-.02*len(str(precision(yval))), yval,precision(yval), va='bottom') #va: vertical alignment y positional argument
  
  plt.xticks(range(len(xticks)),xticks)
  plt.title(tittle)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
#---------------------------------------------------------------------------------------  
def smoothen(a_list,left=1,right=1):
  average = list()
  for u in range(len(a_list)):
    s = 0
    n = 0
    for t in range(u-left,u+right):
      if t in range(len(a_list)):
        s += a_list[t]
        n += 1
    s = s/float(n)
    average.append(s)
  return average
#--------------------------------------------------------------------------------------- 
def average_curves(list1,list2):
  average = list()
  for i in range(min(len(list1),len(list2))):
    average.append((list1[i]+ list2[i])/2.0)
  return average
#---------------------------------------------------------------------------------------  
def approximate(a_list):
  approx = a_list
  for i in range(len(a_list)):
    r = len(a_list)-i
    approx = average_curves(smoothen(a_list,r,r),approx)
  return approx
#---------------------------------------------------------------------------------------  
def index_interval(a_list,e):
  for i in range(len(a_list)-1):
    #print a_list[i] , e , a_list[i+1]
    if a_list[i] <= e <= a_list[i+1] or a_list[i] >= e >= a_list[i+1]:
      return i
  return len(a_list)  
#---------------------------------------------------------------------------------------
def plot_bar2(new_integration,vertical,horizontal,interval,color,option = "raw",precision = lambda t : int(100*t)/100.0):
  plt.rc('xtick',labelsize=8)
  plt.rc('ytick',labelsize=8)
  fig, ax = plt.subplots(vertical,horizontal)
  fig.tight_layout(pad=.3)
  #------------------------------------
  pdflist = list()
  vallist = list()
  #------------------------------------
  for i in interval:
    pdf, val = discrete_pdf(new_integration,i)
    pdflist.append(pdf)
    vallist.append(val)
  #------------------------------------
  for j in range(vertical):
    for k in range(horizontal):
      pdf = pdflist[horizontal*j+k]
      val = vallist[horizontal*j+k]
      #------------------------------------
      ax[j][k].bar(range(len(pdf)),pdf, color = color, alpha = 0.8)
      xticks = map(lambda t: str(int(100*t)/100.0),val)
      plt.sca(ax[j, k])
      plt.xticks(range(len(xticks)), xticks)
      plt.setp(plt.xticks()[1], rotation=90)
      #------------------------------------
      if option == "approximate" or option == "model":
        smooth = approximate(approximate(pdf))
        ax[j][k].plot(range(len(pdf)),smooth,c = "purple")
      if option == "model":
        h     = max(smooth)
        mu    = index_interval(val,val[smooth.index(h)])
        sigma = abs((mu - index_interval(val,.9))/3.0)
        c = scstat.norm.pdf(range(len(smooth)),mu,sigma)
        c = [h/float(max(c)) *i for i in c]
        plt.plot(range(len(smooth)),c)
  #------------------------------------    
  plt.rc('xtick',labelsize=12)
  plt.rc('ytick',labelsize=12)
#---------------------------------------------------------------------------------------
#Display the expression values of a table using green to blue colors
#Takes care of coloring the missing data "PREDICT" in black
def pixel_format(image):
  new_image = list()
  for i in range(len(image)):
    row = list()
    for j in range(len(image[i])):
      if image[i][j] == "PREDICT":
        row.append([30,30,30])
      else:
        row.append([50,int(130*float(image[i][j])),100])
    new_image.append(row)
  return new_image
#---------------------------------------------------------------------------------------
#Display a hashed window using the colors [40,0,80] and [35,0,72].
def space(height,width):
  sp = list()
  for i in range(height):
    s = list()
    for k in range(width):
      if (k+i)%4 != 0 and (k-i)%4 != 0 and not(i in [0,height-1]):
        s.append([40,0,80])
      else:
        s.append([35,0,72])
    sp.append(s)
  return sp
#---------------------------------------------------------------------------------------
fdata = open("data/expression_challenge/ExpressionData_UPDATED.txt","r")
fchal = open("data/expression_challenge/TargetList.txt","r")
fansw = open("data/expression_challenge/ExpressionChallenge.txt","r")
data = Data(fdata,fchal,fansw)
fdata.close()
fchal.close()
fansw.close()
#---------------------------------------------------------------------------------------
#This is a test for the use of scstat.norm.pdf, which is used in the function plot_bar2 (see above).
#The function plot_bar2 is used to define moving averages (see tutorial, section 2.8).
#we plot the cumulative gaussian distribution with mean = .9 and standard deviation = .3.
if "gauss" in sys.argv[1:]:
  plt.plot(range(10),scstat.norm.pdf(map(lambda t:t/10.0, range(10)),.9,.3))
  plt.show()
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
if "figure1" in sys.argv[1:]:
  #The DREAM3 challenge data is organized within the class Data.
  #data.antitargets is list of complete gene expressions (no missing values)
  #data.targets is the list of gene name with missing values
  #get_targets_expressions will ouput the first 8 time points for all strains and
  #for each of the genes in data.antitargets and in data.targets (the output is a lists of lists)
  queries = [data.antitargets,data.targets]
  training, load_initial = get_targets_expressions(queries,data,range(0,8),[0,1,2,3])
 
  #This is the first image shown in section 2.8 of the documentation
  #The variable "images" is a matrix of rgb vectors
  images = pixel_format(training[:40]) + space(10,8*4) + pixel_format(training[200:220])+ pixel_format(load_initial)  + pixel_format(training[500:520]) + space(10,8*4) + pixel_format(training[len(training)-40:])
  
  #Dipslay the matrix of pixels
  fig, ax = plt.subplots(1)
  ax.imshow(images)
  plt.xlabel("expressions")
  plt.ylabel("genes")
  plt.axis("off")
  plt.show()
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
if "figure2" in sys.argv[1:]:
  #The DREAM3 challenge data is organized within the class Data.
  #data.antitargets is list of complete gene expressions (no missing values)
  #data.targets is the list of gene name with missing values
  #get_targets_expressions will ouput the first 8 time points of strains 0, 2 and 3,
  #for each of the genes in data.antitargets and in data.targets (the output is a lists of lists)
  queries = [data.antitargets,data.targets]
  training, load_initial = get_targets_expressions(queries,data,range(0,8),[0,2,3])

  #This image is use din the second figure shown in section 2.8 of the documentation
  #The variable "images" is a matrix of rgb vectors
  images = pixel_format(training[:40]) + space(10,8*3) + pixel_format(training[200:220]) + pixel_format(training[500:520]) + space(10,8*3) + pixel_format(training[len(training)-40:])
  
  #Display an image representing the database without the incomplete gene expressions
  fig, ax = plt.subplots(1)
  ax.imshow(images)
  plt.xlabel("expressions")
  plt.ylabel("genes")
  plt.axis("off")
  plt.show()
  
  #get_targets_expressions will ouput the first 8 time points of all strains (in the order 1,0,2,3)
  #for each of the genes in data.targets (the output is a lists of lists)
  queries = [data.targets]
  load_initial = get_targets_expressions(queries,data,range(0,8),[1,0,2,3])[0]
  
  #This image is use din the second figure shown in section 2.8 of the documentation
  #The variable "images" is a matrix of rgb vectors
  images = pixel_format(load_initial)
  
  #Display an image representing the gene expressions with missing values
  fig, ax = plt.subplots(1)
  ax.imshow(images)
  plt.xlabel("expressions")
  plt.ylabel("genes")
  plt.axis("off")
  plt.show()
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#Same as figure3, but where the two images are merged as one
if "figure3" in sys.argv[1:]:
  queries = [data.antitargets,data.targets]
  training, load_initial = get_targets_expressions(queries,data,range(0,8),[1,0,2,3])
    
  image1 = pixel_format(training[:40]) + space(10,8*4) + pixel_format(training[200:220])+ space(10,8*4)  + pixel_format(training[500:520]) + space(10,8*4) + pixel_format(training[len(training)-40:])    
    
  image2 = pixel_format(load_initial)
  
  images = merge2(image1,image2)
  
  fig, ax = plt.subplots(1)
  ax.imshow(images)
  plt.xlabel("expressions")
  plt.ylabel("genes")
  plt.axis("off")
  plt.show()
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#This is a ad-hoc machine learning method that uses agreements/cos-similarities to impute 
#the missing gene expression. This method works well, but cannot be generalized to other problems.
if "easy_answer" in sys.argv[1:]:
  queries = [data.antitargets,data.targets]
  training, load_initial = get_targets_expressions(queries,data,range(0,8),[0,2,3])
  
  model = get_targets_expressions([data.antitargets],data,range(0,8),[1])[0]
  
  print "load"
  
  train = list()
  for i in range(len(training)):
    train_i = separate_up_down(training[i], lambda x : 100*x)
    train.append(train_i)
  
  intensities = list()
  memor = list()
  for i in range(len(load_initial)):
    memor_i = separate_up_down(load_initial[i], lambda x : 100*x)
    memor.append(memor_i)

  nearest = list()
  for i in range(len(memor)):
    print "----------"+str(i)+"---------"
    search_max = 0
    row = list()
    for j in range(len(train)):
      agree = agreement(memor[i],train[j])
      if search_max < agree:
        search_max = agree
        row = model[j]
      print j, agree  
    nearest.append(row)

  answer = return_answer(nearest,50,8)
  
  #--------------------------------------------------------------#
  
  fig, ax = plt.subplots(1)
  ax.imshow(merge(answer,data.answer))
  plt.show()
  
  corr, pval = gene_accuracy(answer,data.answer,50,8)
  average_corr_gene = arithmetic_mean(corr)
  average_pval_gene = geometric_mean(pval)
  print "\033[93maverage correlation for gene-profile accuracy", average_corr_gene, "\033[0m"
  print "\033[93maverage p-value for gene-profile accuracy", average_pval_gene, "\033[0m"
  
  fig, ax = plt.subplots(1)
  ax.plot(range(len(corr)), corr, '-o', color='black')
  ax.plot(range(len(corr)), [average_corr_gene]*len(corr), color='red')
  plt.show()
  
  fig, ax = plt.subplots(1)
  ax.plot(range(len(pval)), pval, '-o', color='gray')
  
  print "-log(pval) = ", map(lambda t: int(-100*math.log10(t))/100.0,pval)
  
  plt.show()
  
  corr, pval = time_accuracy(answer,data.answer,50,8)
  average_corr_time = arithmetic_mean(corr)
  average_pval_time = geometric_mean(pval)
  #for i in range(len(corr)):
  #  print i, corr[i]
  print "\033[93maverage correlation for time-profile accuracy", average_corr_time, "\033[0m"
  print "\033[93maverage p-value for time-profile accuracy", average_pval_time, "\033[0m"
  
  fig, ax = plt.subplots(1)
  ax.plot(range(len(corr)), corr, '-o', color='black')
  ax.plot(range(len(corr)), [average_corr_time]*len(corr), color='red')
  plt.show()
  
  fig, ax = plt.subplots(1)
  ax.plot(range(len(pval)), pval, '-o', color='gray')
  
  print "-log(pval) = ", map(lambda t: int(-100*math.log10(t))/100.0,pval)
  
  plt.show()
  
  #---------------------------  

  sc = -0.5*math.log10(average_pval_gene*average_pval_time)
  print "\033[91mscore = ", sc, "\033[0m"
#---------------------------------------------------------------------------------------
#---------------------------------------METHOD------------------------------------------
#--------------------------------------------------------------------------------------- 
# Analyses shown in section 2.8.4 of the documentation
if "method" in sys.argv[1:]:
  queries = [data.antitargets,data.targets]
  training, load_initial = get_targets_expressions(queries,data,range(0,8),[0,2,3])

  train = list()
  for i in range(len(training)):
    train_i = separate_up_down(training[i], lambda x : 100*x)
    train.append(train_i)
  
  whole = list()
  intensities = list()
  memor = list()
  for i in range(len(load_initial)):
    memor_i = separate_up_down(load_initial[i], lambda x : 100*x)
    memor.append(memor_i)
    row = list()
    for j in range(len(train)):
      row.append(agreement(memor_i,train[j]))
    print "[row "+str(i)+"] done"
    whole = whole + row
    intensities.append(row)
  
  #-----------------------------------------------
  
  #fig, ax = plt.subplots(10,5)
  #for i in range(10):
  #  for j in range(5):
  #    ax[i][j].imshow([intensities[5*i+j]]*600)
  #plt.show()
  
  #-----------------------------------------------
  
  thresholds = list()
  fig, ax = plt.subplots(10,5)
  for i in range(10):
    for j in range(5):
      print "[distribution at "+str(i)+","+str(j)+"] done"
      x, y, colors, threshold = distribution(intensities[5*i+j],100)
      thresholds.append(threshold)
      ax[i][j].scatter(x, y,s=1,color = colors)
  fig.tight_layout(pad=.3)
  
  plt.show()

  #-----------------------------------------------
  
  integration = list()
  agree_graph = list()
  x = []
  for i in range(len(memor)):
    print "----------"+str(i)+"---------", thresholds[i]
    agree_list = list()
    for j in range(len(train)):
      agree = agreement(memor[i],train[j])
      if agree >= thresholds[i]:
        agree_list.append(agree)
        print j, agree
    agree_graph.extend(agree_list)
    integration.append(arithmetic_mean(agree_list))
    x =  x + [i] * len(agree_list)  
  
    #print len(x), len(agree_graph)
  
  
  fig, ax = plt.subplots(1)
  ax.scatter(x,agree_graph,s=2, color = "gray")
  ax.scatter(range(len(integration)),integration,s=5, color = "red")
  
  fig, ax = plt.subplots(1)
  plt.axhspan(.9, max(agree_graph), xmin=0, xmax=50, facecolor='red', alpha=0.4)
  ax.scatter(x,agree_graph,s=2, color = "black")
  
  
  fig, ax = plt.subplots(1)
  plt.axhspan(.75, max(agree_graph), xmin=0, xmax=50, facecolor='blue', alpha=0.4)
  ax.scatter(x,agree_graph,s=2, color = "black")
  
  plt.show()
  
  #-----------------------------------------------
  
  new_integration = integration[:]
  new_integration.sort()
  M = max(new_integration)
  m = min(new_integration)
  I = (M-m)/5.0
  
  fig, ax = plt.subplots(1)
  ax.scatter(x,agree_graph,s=2, color = "gray")
  ax.scatter(range(len(integration)),integration,s=5, color = "red")
  for k in range(6):
    plt.hlines(m+k*I,0,len(integration)-1, color = "#58cd4c")
    
  plt.show()
    
  #-----------------------------------------------
  
  interval = [5*(11-i) for i in range(0,11)]+[3]
  plot_bar2(new_integration,4,3,interval,color = "#7fcdac")
  plot_bar2(new_integration,4,3,interval,color = "#7fcdac",option = "approximate")
  #plot_bar2(new_integration,4,3,interval,color = "#7fcdac",option = "model")
  
  plt.show()

  #pdf, val = flexible_discrete_pdf(whole,[(.9,1),(.85,1),(.8,1)])
  #plot_bar(pdf,"#7fcdac","Proportions","",map(str,val),"")
  #plt.show()
  
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
if "training" in sys.argv[1:]:
  queries = [data.antitargets,data.targets]
  training, load_initial = get_targets_expressions(queries,data,range(0,8),[1,0,2,3])
  
  ftrain = gzip.open("data/dream3_training.gz","w")
  for i in range(len(training)):
    train_i = separate_up_down(training[i], lambda x : 100*x)
    ftrain.write("\t".join(map(str,train_i))+"\n")
    ftrain.flush()
  ftrain.close()
    
  
  fmemor = gzip.open("result-load/load_initial.gz","w")
  image = list()
  for i in range(len(load_initial)):
    memor_i = separate_up_down(load_initial[i], lambda x : 100*x)
    image.append(memor_i)
    fmemor.write("\t".join(map(str,memor_i))+"\n")
    fmemor.flush()
  fmemor.close()
  
  print "training = ", len(training)
  print "ary = ", len(load_initial)
  print "dim = ", len(image[0])
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
if "result" in sys.argv[1:]:
  ary = 50
  load_index = 0
  fload_name = "save_org.gz"
  #---------------------------
  if len(sys.argv) == 2:
    fload = gzip.open(fload_name,"r")
    load_index = usf.get_last_cycle(fload,ary)
    print "loading organelles at index "+str(load_index)+" from save_org.gz"
    fload.close()
  elif len(sys.argv) > 2:
    load_index = int(sys.argv[2])
  #---------------------------
  try:
    interval = int(sys.argv[3])
  except:
    interval = 1000
  index_history = list()
  for i in range(int(load_index/float(interval))):
    index_history.append(interval*i)
  index_history.append(load_index)
  #---------------------------
  fload = gzip.open(fload_name,"r")
  list_of_organelles = usf.get_memory(fload,ary,index_history)
  fload.close()
  #---------------------------
  recovery_of_memory = get_targets_expressions([data.targets],data,range(0,8),[0,2,3],float)[0]
  list_of_scores = list()
  #---------------------------
  for k in range(len(list_of_organelles)):
    organelles = list_of_organelles[k]
    #---------------------------
    expression = list()
    for i in range(len(organelles)):
      r = reunite_up_down(organelles[i][:16], lambda t: t/100.0)
      x = reunite_up_down(organelles[i][16:], lambda t: t/100.0)
      j, d = closest_vector(x,recovery_of_memory)
      expression.append((j,r,i))
    #---------------------------
    ordered_expression = sorted(expression, key = lambda t: t[0])
    matrix = list()
    for i in range(len(organelles)):
      matrix.append(ordered_expression[i][1])
    #---------------------------
    answer = return_answer(matrix,50,8)
  
    corr_gene, pval_gene = gene_accuracy(answer,data.answer,50,8)
    average_corr_gene = arithmetic_mean(corr_gene)
    average_pval_gene = geometric_mean(pval_gene)
    print "\033[93maverage correlation for gene-profile accuracy", average_corr_gene, "\033[0m"
    print "\033[93maverage p-value for gene-profile accuracy", average_pval_gene, "\033[0m"
   
    corr_time, pval_time = time_accuracy(answer,data.answer,50,8)
    average_corr_time = arithmetic_mean(corr_time)
    average_pval_time = geometric_mean(pval_time)
    print "\033[93maverage correlation for time-profile accuracy", average_corr_time, "\033[0m"
    print "\033[93maverage p-value for time-profile accuracy", average_pval_time, "\033[0m"
    
    #---------------------------
    
    dim = 64
    selfsup = [[16,dim-1,dim]]
    brightness = [.98]           
    profiles = [ [(.85,1)], [(0,.85)] ]
    contrast = usf.contrast(brightness,profiles)
    
    
    #---------------------------  
    if average_pval_gene*average_pval_time != 0:
      sc = -0.5*math.log10(average_pval_gene*average_pval_time)
      list_of_scores.append(sc)
    else:
      list_of_scores.append(float("nan"))
    print "\033[91mscore = ", sc, "\033[0m"
    #---------------------------
  
    #See https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0008944
    #to compare with the results of Table 1
    try:
      if k == len(list_of_organelles)-1:
        print "----------"
        for i in range(len(organelles)):
          print i, organelles[ordered_expression[i][2]][:16], "...", int(-100*math.log10(pval_time[i]))/100.0
          #contrast(organelles[ordered_expression[i][2]],*selfsup)
        print "\n----------\n"
        for i in range(len(organelles)):
          print i, map(lambda x : int(1000000*x),matrix[i]), "...", int(-100*math.log10(pval_time[i]))/100.0
        print "\n----------\n"
        for i in range(len(organelles)):
          print i, answer[i], "...", int(-100*math.log10(pval_time[i]))/100.0, "...", data.answer[i]
        print "\n----------\n"
  
        fig, ax = plt.subplots(1)
        ax.imshow(merge(answer,data.answer))
        plt.show()

        fig, ax = plt.subplots(1)
        ax.plot(range(len(corr_gene)), corr_gene, '-o', color='black')
        ax.plot(range(len(corr_gene)), [average_corr_gene]*len(corr_gene), color='red')
        plt.show()
        
        fig, ax = plt.subplots(1)
        ax.plot(range(len(pval_gene)), pval_gene, '-o', color='gray')
        print "-log(pval_gene) = ", map(lambda t: int(-100*math.log10(t))/100.0,pval_gene)
        plt.show()

        fig, ax = plt.subplots(1)
        ax.plot(range(len(corr_time)), corr_time, '-o', color='black')
        ax.plot(range(len(corr_time)), [average_corr_time]*len(corr_time), color='red')
        plt.show()
        
        fig, ax = plt.subplots(1)
        ax.plot(range(len(pval_time)), pval_time, '-o', color='gray')
        print "-log(pval_time) = ", map(lambda t: int(-100*math.log10(t))/100.0,pval_time)
        plt.show()
      
    except:
        
      if k == len(list_of_organelles)-1:
        print "----------"
        for i in range(len(organelles)):
          print i, organelles[ordered_expression[i][2]][:16]
          #contrast(organelles[ordered_expression[i][2]],*selfsup)
        print "\n----------\n"
        for i in range(len(organelles)):
          print i, map(lambda x : int(1000000*x),matrix[i])
        print "\n----------\n"
        for i in range(len(organelles)):
          print i, answer[i], "...", data.answer[i]
        print "\n----------\n"
        
        fig, ax = plt.subplots(1)
        ax.imshow(merge(answer,data.answer))
        plt.show()

        fig, ax = plt.subplots(1)
        ax.plot(range(len(corr_gene)), corr_gene, '-o', color='black')
        ax.plot(range(len(corr_gene)), [average_corr_gene]*len(corr_gene), color='red')
        plt.show()
        

        fig, ax = plt.subplots(1)
        ax.plot(range(len(corr_time)), corr_time, '-o', color='black')
        ax.plot(range(len(corr_time)), [average_corr_time]*len(corr_time), color='red')
        plt.show()

    #---------------------------
    
  colors = list()
  for i in range(int(load_index/float(interval))+1):
    if i in map(lambda x: int(1000*x/float(interval)),[36,37,38,73,74,75,119,120,121]):
      colors.append("red")
    else:
      colors.append("black")
      
  fig, ax = plt.subplots(1)
  x = map(lambda s : interval * s, range(int(load_index/float(interval)))) + [load_index]
  ax.scatter(x, list_of_scores, s = 5, color=colors)
  plt.show()  
#---------------------------------------------------------------------------------------
