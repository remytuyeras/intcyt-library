import sys
sys.path.insert(0, '../intcyt/useful/')
from useful import *
import math
#-----------------------------------------------------------------
class Cell(object):
  def __init__(self,dimension,residual,cytosol,organelles):
    #~~~~~~~~~~~~~~~~~~~
    valid = (dimension == len(cytosol))
    for x in organelles:
      valid = valid and (dimension == len(x))
    if not(valid):
      print "Error: in Cell.__init__: list dimensions do not equal "+ str(dimension)
      exit(0)
    #~~~~~~~~~~~~~~~~~~~  
    self.dimension = dimension #int
    self.residual = residual #float
    self.cytosol = cytosol #list(int)
    self.organelles = organelles #list(list(int))
    #~~~~~~~~~~~~~~~~~~~
    self.K = list() #fast access to content
    self.SK = 0 #fast access to sum(content)
    self.Sorg = [0]*len(self.organelles) #fast access to the values sum(organelles[i])
    self.well_defined = (self.residual >= 0.0) and (self.residual + sum(self.cytosol) >= 0.0)
    #~~~~~~~~~~~~~~~~~~~
    for u in range(self.dimension):
      K_u = self.cytosol[u]
      for i in range(len(self.organelles)):
        #~~~~~~~~~~~~~~~~~~~
        self.well_defined = self.well_defined and (self.organelles[i][u] >= 0.0)
        K_u = K_u + self.organelles[i][u]
        self.Sorg[i] += self.organelles[i][u]
        #~~~~~~~~~~~~~~~~~~~
      self.K.append(K_u)
      self.SK += K_u
    #~~~~~~~~~~~~~~~~~~~
    self.well_defined_cytosol = list()
    for u in range(self.dimension):
      self.well_defined_cytosol.append(self.cytosol[u] >= 0.0)
#-----------------------------------------------------------------
  def check(self):
    print "residual:", self.residual >= 0.0
    print "residual + sum(cytosol):", self.residual + sum(self.cytosol) >= 0.0
    for u in range(self.dimension):
      print "-> cytosol["+str(u)+"]:", self.cytosol[u] >= 0.0
    for i in range(len(self.organelles)):
      for u in range(self.dimension):
        print "organelles["+str(i)+"]["+str(u)+"]:", self.organelles[i][u] >= 0.0
    return self
#-----------------------------------------------------------------  
  def copy(self):
    organelles = map(lambda x: x[:],self.organelles)
    return Cell(self.dimension,self.residual,self.cytosol[:],organelles)
#-----------------------------------------------------------------  
  def stdout(self):
    print "residual: " + str(self.residual)
    for u in range(self.dimension):
      print "cytosol[" + str(u+1)+"]: "+ str(self.cytosol[u]) 
    for i in range(len(self.organelles)):
      print str(i+1)+"-th organelle: " + str(self.organelles[i])
    return self
#-----------------------------------------------------------------  
  def content(self):
    the_content = list()
    for u in range(self.dimension):
      K = self.cytosol[u]
      for i in range(len(self.organelles)):
        if self.organelles[i] != []:
          K = K + self.organelles[i][u]
      the_content.append(K)
    return the_content
#-----------------------------------------------------------------  
  def compose(self,index,a_cell):
    #~~~~~~~~~~~~~~~~~~~
    if self.dimension != a_cell.dimension:
      print "Error: in Cell.factorize: dimensions are not equal"
      exit(0)
    #~~~~~~~~~~~~~~~~~~~ 
    new_residual = self.residual + a_cell.residual
    new_cytosol = list()
    for u in range(self.dimension):
      new_cytosol.append(self.cytosol[u]+ a_cell.cytosol[u]) 
    new_organelles = list()
    for i in range(len(self.organelles)):
      if i == index:
        new_organelles.extend(a_cell.organelles)
      else:
        new_organelles.append(self.organelles[i])
    return Cell(self.dimension,new_residual,new_cytosol,new_organelles)
#-----------------------------------------------------------------  
  def compose_all(self,cells):
    #~~~~~~~~~~~~~~~~~~~
    if len(cells) != len(self.organelles):
      print "Error: in simultaneously_compose: input needs", len(self.organelles), "cells"
      exit(0)
    #~~~~~~~~~~~~~~~~~~~
    the_cell = self
    composition_index = 0
    for k in range(len(self.organelles)):
      the_cell = the_cell.compose(composition_index,cells[k])
      composition_index += len(cells[k].organelles)
    return the_cell 
#-----------------------------------------------------------------
  def left(self,alpha_var,kappa_var,lambda_var):
    #~~~~~~~~~~~~~~~~~~~
    if not(len(alpha_var) == len(kappa_var) == len(lambda_var)):
      print "Error: in Cell.left: dimensions of variables are not equal"
      exit(0)
    #~~~~~~~~~~~~~~~~~~~
    new_cytosol = self.cytosol[:]
    new_residual = self.residual
    for i in range(len(lambda_var)):
      new_residual -= alpha_var[i]
      for u in range(self.dimension):
        new_cytosol[u] -= kappa_var[i][u]
    return Cell(self.dimension,new_residual,new_cytosol,lambda_var)
#-----------------------------------------------------------------  
  def right(self,mu_var,alpha_var,kappa_var):
    #~~~~~~~~~~~~~~~~~~~
    if not(sum(mu_var) == len(self.organelles) 
    and len(mu_var) == len(alpha_var) == len(kappa_var)):
      print "Error: in Cell.right: tilling exceed the number of organelles",sum(mu_var), len(self.organelles)
      exit(0)
    #~~~~~~~~~~~~~~~~~~~
    output = list()
    L = 0
    for i in range(len(mu_var)):
      c = Cell(self.dimension,alpha_var[i],kappa_var[i],self.organelles[L:L+mu_var[i]])
      L += mu_var[i]
      output.append(c)
    return output
#-----------------------------------------------------------------
  def spontaneous_reaction(self):
    new_residual = sum(self.cytosol)
    self.residual += new_residual
    #~~~~~~~~~~~~~~~~~~~
    if self.residual < 0:
      print "Error: in Cell.spontaneous_reaction: self.residual is not positive"
      self.stdout()
      exit(0)
    #~~~~~~~~~~~~~~~~~~~ 
    for u in range(self.dimension):
      self.K[u] -= self.cytosol[u] #fast update
    self.SK = self.SK - new_residual #fast update
    self.cytosol = [0] * self.dimension
    return self
#-----------------------------------------------------------------
  def action(self,matrix):
    #~~~~~~~~~~~~~~~~~~~
    if len(matrix) != len(self.organelles):
      print "Error: in Cell.action: the input list should contain "+str(len(self.organelles))+" list(s)"
      print matrix
      self.stdout()
      exit(0)
    #~~~~~~~~~~~~~~~~~~~ 
    the_action = [0]*self.dimension
    if self.SK != 0:
      for i in range(len(self.organelles)):
        #~~~~~~~~~~~~~~~~~~~
        if len(matrix[i]) != self.dimension:
          print "Error: in Cell.action: the lists of the input should contain "+str(self.dimension)+" float items"
          self.stdout()
          exit(0)
        #~~~~~~~~~~~~~~~~~~~ 
        for u in range(self.dimension):
          the_action[u] += (self.organelles[i][u]/float(self.SK)) * matrix[i][u]
    return the_action
#------------------------------------------------------------------------------
  def algebra_operator(self,matrix):
    #~~~~~~~~~~~~~~~~~~~
    if len(matrix) != len(self.organelles):
      print "Error: in Cell.algebra_operator: the input list should contain "+str(len(self.organelles))+" list(s)"
      self.stdout()
      exit(0)
    #~~~~~~~~~~~~~~~~~~~ 
    algebra_operator = [0] * self.dimension
    if self.SK != 0:
      for u in range(self.dimension):
        for k in range(len(self.organelles)):
          #~~~~~~~~~~~~~~~~~~~
          if len(matrix[k]) != self.dimension:
            print "Error: in Cell.algebra_operator: the lists of the input should contain "+str(self.dimension)+" float items"
            self.stdout()
            exit(0)
          #~~~~~~~~~~~~~~~~~~~ 
          dividend = self.Sorg[k]-self.organelles[k][u]  
          algebra_operator[u] += (dividend/float(self.SK)) * matrix[k][u]
    return algebra_operator
#-----------------------------------------------------------------   
  def allostasis(self,matrix,weight,org_index,dim_index):
    #~~~~~~~~~~~~~~~~~~~
    if len(matrix) != len(self.organelles) or len(matrix[org_index]) != self.dimension:
      print "Error: in Cell.allostasis: the size of the first input should be "+str(len(self.organelles))+"x"+str(self.dimension)+"-list"
      self.stdout()
      exit(0)
    #~~~~~~~~~~~~~~~~~~~ 
    diff = 0
    if self.SK != 0:
      for u in range(self.dimension):
        if u != dim_index:
          ratio = self.organelles[org_index][u]/float(self.Sorg[org_index])
          term1 = matrix[org_index][dim_index] * weight[dim_index]
          term2 = matrix[org_index][u] * weight[u]
          diff += ratio * (term1 - term2)
      diff = -diff/float(self.SK)
    return diff
#-----------------------------------------------------------------
  def agreement(self,index,vector,*challenge):
    orga_i = self.organelles[index]
    agreement = 0
    norm_o = 0
    norm_v = 0
    for u in range(self.dimension):
      if len(challenge) == 1 and not(challenge[0][0] <= u % challenge[0][2] <= challenge[0][1]):
        continue
      agreement += orga_i[u] * vector[u]
      norm_o += orga_i[u]**2
      norm_v += vector[u]**2
    norm_o = math.sqrt(norm_o)
    norm_v = math.sqrt(norm_v)
    if norm_o != 0 and norm_v != 0: 
      return agreement/float(norm_o*norm_v)
    else:
      return 0
#------------------------------------------------------------------------------
  def merge(self,list_of_organelles,order="non-sorted"):
    if order !="sorted":
      list_of_organelles.sort()
    if len(list_of_organelles) != 0:
      #From here, j should be the minimun element of list_of_organelles
      j = list_of_organelles[0]
      new_organelles = list()
      merged_organelles = [0] * self.dimension
      for i in range(len(self.organelles)):
        if i in list_of_organelles:
          for u in range(self.dimension):
            merged_organelles[u] += self.organelles[i][u]
          if i == j:
            new_organelles.append(merged_organelles)
        else:
          new_organelles.append(self.organelles[i]) 
      self.organelles = new_organelles
    return self
#------------------------------------------------------------------------------
  def divide(self,list_of_organelles,order="non-sorted"):
    #~~~~~~~~~~~~~~~~~~~
    if self.cytosol != [0] * self.dimension:
      print "Error: in Cell.divide: cytosol is not zero -- cannot divide"
      exit(0)
    #~~~~~~~~~~~~~~~~~~~
    new_organelles_1 = list()
    new_organelles_2 = list()
    if order !="sorted":
      list_of_organelles.sort()
    j = 0
    for i in range(len(self.organelles)):
      if j < len(list_of_organelles) and i == list_of_organelles[j]: 
        new_organelles_1.append(self.organelles[i])
        j += 1
      else:
        new_organelles_2.append(self.organelles[i])
    c_1 = Cell(self.dimension,0.5*self.residual,[0] * self.dimension,new_organelles_1) 
    c_2 = Cell(self.dimension,0.5*self.residual,[0] * self.dimension,new_organelles_2)   
    return [c_1,c_2]
#------------------------------------------------------------------------------
  def organelle_proportions(self):
    prop = list()
    for i in range(len(self.organelles)):
      prop.append(map(lambda x : float(x)/self.Sorg[i] if self.Sorg[i]!=0 else 0, self.organelles[i]))
    return prop
#------------------------------------------------------------------------------
  def content_proportions(self):
    s = sum(self.Sorg)
    prop = map(lambda x : float(x)/s if s!=0 else 0, self.Sorg)
    return prop
#------------------------------------------------------------------------------
  def best_compartment(self,cliques):
    result = list()
    for k in range(len(cliques)):
      result_k = 0
      for u in range(self.dimension):
        dividend1 = 0
        divisor1 = 0
        dividend2 = 0
        divisor2 = 0
        for i in range(len(self.organelles)):
          if i in cliques[k]:
            dividend1 += self.organelles[i][u]
            divisor1 += self.Sorg[i]
          else:
            dividend2 += self.organelles[i][u]
            divisor2 += self.Sorg[i]
        if divisor1 * dividend2 != 0 and dividend1/float(divisor1) > dividend2/float(divisor2):   
          result_k += (dividend1*divisor2)/float(dividend2*divisor1)
          #result_k += 1
      result.append([k,result_k])
      result.sort(key = lambda x :-x[1])
    if result[0][1] != 0:
      return cliques[result[0][0]]
    else:
      return list() 
#------------------------------------------------------------------------------      
  def proposed_clustering(self,matrix,option,filtering = 1.5):
    #~~~~~~~~~~~~~~~~~~~
    if len(matrix) != len(self.organelles):
      print "Error: in Cell.proposed_clustering: the input list should contain " + \
             str(len(self.organelles))+" list(s)"
      self.stdout()
      exit(0)
    #~~~~~~~~~~~~~~~~~~~
    #Computing prop increases the running time; instead we use self.Sorg[i]/float(self.SK)
    #prop = self.content_proportions()
    graph = usf.zero_matrix(len(self.organelles))
    #We use "nonempty" to know if there exists at least one relationship
    nonempty = False
    for u in range(self.dimension):
      #~~~~~~~~~~~~~~~~~~~
      barycenter = 0
      if self.SK != 0:
        for i in range(len(self.organelles)):
          #~~~~~~~~~~~~~~~~~~~
          if len(matrix[i]) != self.dimension:
            print "Error: in Cell.proposed_clustering: the "+str(i)+ \
                  "-th list of the input list should be of length "+str(self.dimension)
            self.stdout()
            exit(0)
          #~~~~~~~~~~~~~~~~~~~
          #barycenter += matrix[i][u] * prop[i]
          barycenter += matrix[i][u] * (self.Sorg[i]/float(self.SK))
      #~~~~~~~~~~~~~~~~~~~
      graph_u = list()
      for i in range(len(self.organelles)):
        if option == "division":
          if matrix[i][u] > filtering * barycenter:
            graph_u.append(i)
        elif option == "merging":
          if matrix[i][u] < filtering * barycenter:
            graph_u.append(i)
        else:
          print "Warning: Error: in Cell.proposed_clustering: no valid option given (option = " + \
                str(option)+")"
      #~~~~~~~~~~~~~~~~~~~
      for i in range(len(graph_u)):
        for j in range(i+1,len(graph_u)):
          graph[graph_u[i]][graph_u[j]] += 1
          nonempty = True
      #~~~~~~~~~~~~~~~~~~~
    if nonempty:
      cliques =  usf.cliques(graph)
      return self.best_compartment(cliques)
      #~~~~~~~~~~~~~~~~~~~
    else:
      return list() 
#------------------------------------------------------------------------------
