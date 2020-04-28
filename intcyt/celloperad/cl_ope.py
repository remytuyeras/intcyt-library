import random

from cl_sup import SuperCell
from cl_cel import Cell
#-----------------------------------------------------------------
class Operad(object):
  def __init__(self,dimension):
    self.dimension = dimension
#-----------------------------------------------------------------  
  def identity(self,vector):
    if len(vector) != self.dimension:
      print "Error: in Operad.identity: dimensions do not match"
      exit(0) 
    return Cell(self.dimension,0,[0]*self.dimension,[vector])
#-----------------------------------------------------------------  
  def tensor(self,list_of_cells):
    for i in range(len(list_of_cells)):
      if self.dimension != list_of_cells[i].dimension:
        print "Error: in Operad.tensor: dimensions do not match"
        exit(0)
    tensor_cell = Cell(self.dimension,0,[0]*self.dimension,[[0]*self.dimension]*len(list_of_cells))
    return tensor_cell.compose_all(list_of_cells)
#-----------------------------------------------------------------  
  def pseudo_tensor(self,list_of_super_cells):
    for i in range(len(list_of_super_cells)):
      if self.dimension != list_of_super_cells[i].cell.dimension:
        print "Error: in Operad.merging_tensor: dimensions do not match"
        exit(0)
        
    organelles = list()
    for i in range(len(list_of_super_cells)):
      organelles.append(list_of_super_cells[i].cell.content())

    tensor_cell = Cell(self.dimension,0,[0]*self.dimension,organelles)
    return SuperCell(tensor_cell,list_of_super_cells)
#-----------------------------------------------------------------  
  def merging_tensor(self,list_of_super_cells):
    tensor =  self.pseudo_tensor(list_of_super_cells)
    #inner composition
    if not(tensor.is_leaf):
      for i in range(len(tensor.innercells)):
        tensor.innercells[i].compose_state = True
    return tensor   
#-----------------------------------------------------------------  
  def dividing_tensor(self,list_of_super_cells):
    tensor =  self.pseudo_tensor(list_of_super_cells)
    #outer composition
    tensor.compose_state = True
    return tensor 
#-----------------------------------------------------------------
  def generate(self,levels,arity,key):
    cells = list()
    sup_cells = list()
    for i in range(arity**levels):
      first_cells = list()
      for j in range(arity):
        random_cell = list()
        for u in range(self.dimension):
          random_cell.append(key())
        first_cells.append(random_cell)
        
      c = Cell(self.dimension,0,[0]*self.dimension,first_cells)
      cells.append(c) 
      sup_cells.append(SuperCell(c))
      
    while len(sup_cells) > 1:
    
      tmp_cells = list()
      tmp_sup_cells = list()
      i = 0
      while i < len(cells):
        next_cells = list()
        next_sup_cells = list()
        
        for j in range(arity):
          next_cells.append(cells[i+j].content()[:])
          next_sup_cells.append(sup_cells[i+j])
          
        c = Cell(self.dimension,0,[0]*self.dimension,next_cells)
        sc = SuperCell(c,next_sup_cells)
        
        tmp_cells.append(c)
        tmp_sup_cells.append(sc)
        
        i += arity
        
      cells = tmp_cells
      sup_cells = tmp_sup_cells
      
    return sup_cells[0]
#-----------------------------------------------------------------
