import sys
#------------------------------------------------------------------------------
class SuperCell(object):
#------------------------------------------------------------------------------
  def __init__(self,cell,*innercells):
    self.depth = 0
    self.level = 0
    self.cell = cell
    self.is_leaf = (len(innercells) == 0) or (innercells[0] == [])
    self.compose_state = False
    self.pre_action = list()
    if not(self.is_leaf):
      self.innercells = innercells[0]
      for i in range(len(self.innercells)):
        self.innercells[i].set_level(1)
      for i in range(len(self.innercells)):
        self.depth = max(self.depth,self.innercells[i].depth)
      self.depth += 1
#------------------------------------------------------------------------------  
  def set_level(self,level):
    self.level = level
    if not(self.is_leaf):
      for i in range(len(self.innercells)):
        self.innercells[i].set_level(level+1)
    return self
#------------------------------------------------------------------------------
  def reset_depth(self):
    if self.is_leaf:
      self.depth = 0
    else:
      self.depth = 0
      for i in range(len(self.innercells)):
        self.innercells[i].reset_depth()
        self.depth = max(self.depth,self.innercells[i].depth)
      self.depth += 1
#------------------------------------------------------------------------------
  def stdout(self,*vector):
    #~~~~~~~~~~~~~~~~~~~
    stdout_1 = "."*self.level+"["+str(self.level)+"] -> Cell["+str(int(self.cell.residual))+"]{"
    stdout_2 = '\033[95m' + stdout_1
    #~~~~~~~~~~~~~~~~~~~
    stdout_3 = "\033[94mo"*len(self.cell.organelles)+"\033[0m}"
    stdout_4 =         "o"*len(self.cell.organelles)+"}"
    #~~~~~~~~~~~~~~~~~~~
    display_float = lambda x :  str(int(10**10*x)/float(10**10))
    #~~~~~~~~~~~~~~~~~~~
    if not(self.compose_state):
      sys.stdout.write(stdout_1)
    else:
      sys.stdout.write(stdout_2)
    #~~~~~~~~~~~~~~~~~~~  
    if len(vector) != 0 and self.is_leaf:
      sys.stdout.write(stdout_3)
      for i in range(len(self.cell.organelles)):
        agreement = display_float(self.cell.agreement(i,vector[0]))
        if not(self.compose_state):
          sys.stdout.write("  \033[94m["+str(i)+"]\033[0m:"+agreement+"")
        else:
          sys.stdout.write("  \033[94m["+str(i)+"]\033[95m:"+agreement+"")
    else:
      sys.stdout.write(stdout_4)
    #~~~~~~~~~~~~~~~~~~~
    sys.stdout.write("\n")
    sys.stdout.write('\033[0m')
    #~~~~~~~~~~~~~~~~~~~  
    if not(self.is_leaf):
      for i in range(len(self.innercells)):
        self.innercells[i].stdout(*vector)
    return self
#------------------------------------------------------------------------------      
  def compress(self):
    the_cell = self.cell
    if not(self.is_leaf):
      composition_index = 0
      for k in range(len(self.innercells)):      
        r = self.innercells[k].compress()
        if  r!= None:
          the_cell = the_cell.compose(composition_index,r)
          composition_index += len(r.organelles)
        else:
          composition_index +=1
    return the_cell 
#------------------------------------------------------------------------------
  def action(self,vector,identity,option = "expected"):
    self.pre_action = list()
    if self.is_leaf:
      for i in range(len(self.cell.organelles)):
        if option == "specialized":
          self.pre_action.append(identity(self.cell.organelles[i]).action([vector]))
        else:
          self.pre_action.append(vector)
    else:
      for i in range(len(self.innercells)):
        self.pre_action.append(self.innercells[i].action(vector,identity,option))
    return self.cell.action(self.pre_action)
#------------------------------------------------------------------------------
  def base(self,identity):
    children = list()
    if self.is_leaf:
      for i in range(len(self.cell.organelles)):
        children.append(SuperCell(identity(self.cell.organelles[i])))
    else:
      for i in range(len(self.innercells)):
        children.append(SuperCell(self.innercells[i].cell))
    return SuperCell(self.cell,children)
#------------------------------------------------------------------------------      
  def compute_variables(self,gradient_descent):
    mu_var = list()
    alpha_var = list()
    kappa_var = list()
    lambda_var = list()
    #~~~~~~~~~~~~~~~~~~~
    def lambda_(i,u,modif):
      return self.cell.organelles[i][u] - modif
    #~~~~~~~~~~~~~~~~~~~
    def kappa_(i,u,modif):
      return  lambda_(i,u,modif) - \
      self.innercells[i].cell.K[u] + \
      self.innercells[i].cell.cytosol[u]
    #~~~~~~~~~~~~~~~~~~~
    for i in range(len(self.cell.organelles)):
      kappa_i = list()
      lambda_i = list()
      #~~~~~~~~~~~~~~~~~~~
      for u in range(self.cell.dimension):    
        if lambda_(i,u,gradient_descent[i][u]) < 0:
          gradient_descent[i][u] = 0
        kappa_i.append(kappa_(i,u,gradient_descent[i][u]))
        lambda_i.append(lambda_(i,u,gradient_descent[i][u]))
      #~~~~~~~~~~~~~~~~~~~   
      kappa_var.append(kappa_i)
      lambda_var.append(lambda_i)
      alpha_var.append(self.innercells[i].cell.residual)
      mu_var.append(len(self.innercells[i].cell.organelles))
         
    return [mu_var, alpha_var, kappa_var, lambda_var]
#------------------------------------------------------------------------------
  def homeostasis(self,gradient_descent,identity):
    cd = self.base(identity)
    c = cd.cell
    d = cd.innercells
    e = cd.compress()
    #~~~~~~~~~~~~~~~~~~~
    still_updating = True
    #~~~~~~~~~~~~~~~~~~~
    while still_updating:
      #~~~~~~~~~~~~~~~~~~~
      still_updating = False
      parameters_found = True
      #~~~~~~~~~~~~~~~~~~~
      mu_var, alpha_var, kappa_var, lambda_var = cd.compute_variables(gradient_descent)
      #~~~~~~~~~~~~~~~~~~~
      new_c = e.left(alpha_var,kappa_var,lambda_var)
      new_d = e.right(mu_var,alpha_var,kappa_var)
      #~~~~~~~~~~~~~~~~~~~
      innercells_well_defined = True
      for i in range(len(new_d)):
        innercells_well_defined = innercells_well_defined and new_d[i].well_defined
      #~~~~~~~~~~~~~~~~~~~
      if not(new_c.well_defined) or (not(innercells_well_defined) and not(self.is_leaf)): 
        parameters_found = False
        for i in range(len(self.cell.organelles)):
          for u in range(self.cell.dimension):
            #~~~~~~~~~~~~~~~~~~~
            if new_d[i].well_defined:
              if not(new_c.well_defined) and not(new_c.well_defined_cytosol[u]):
                #case where we try to remove too much from the cytosol of the parent new_c
                if gradient_descent[i][u] < 0:
                  gradient_descent[i][u] = 0
                  still_updating = True
            #~~~~~~~~~~~~~~~~~~~      
            elif not(self.is_leaf):
              if not(new_d[i].well_defined_cytosol[u]):
                #case where we try to remove too much from the cytosol of the child new_d[i]
                if gradient_descent[i][u] > 0:
                  gradient_descent[i][u] = 0
                  still_updating = True
            #~~~~~~~~~~~~~~~~~~~ 
      if parameters_found:
        self.cell = new_c
        if not(self.is_leaf):
          for i in range(len(c.organelles)):
            self.innercells[i].cell = new_d[i]
    
    return self
#-----------------------------------------------------------------
  def allostasis(self,vector,identity,gamma):
    self.pre_action = list()
    def display_average_change(s,p,n): 
      print "\033[92m["+s+"][\033[95m"+ \
            str(n)+"\033[92m]\033[0m\n" + \
            "\n".join(map(lambda x : str(sum(x)/len(x)),p))
    #~~~~~~~~~~~~~~~~~~~
    if self.is_leaf:
      for i in range(len(self.cell.organelles)):
        self.pre_action.append(vector)
    else:
      for i in range(len(self.innercells)):
        d_i = self.innercells[i].cell.copy()
        self.innercells[i].allostasis(vector,identity,gamma)
        self.pre_action.append(d_i.action(self.innercells[i].pre_action))
        #self.pre_action.append(self.innercells[i].cell.action(self.innercells[i].pre_action))
    #~~~~~~~~~~~~~~~~~~~ 
    algebra_operator = self.cell.algebra_operator(self.pre_action)
    #~~~~~~~~~~~~~~~~~~~
    parameter_gamma = gamma(self.cell,self.pre_action)
    display_average_change("Gamma parameters",parameter_gamma,0)
    #~~~~~~~~~~~~~~~~~~~      
    gradient_descent =list()    
    for i in range(len(self.cell.organelles)):
      gradient_descent_i = list()
      #~~~~~~~~~~~~~~~~~~~
      for u in range(self.cell.dimension):
        # Smart specialization
        if parameter_gamma[i][u] > 0:
          #The derivative of U^2: takes time to compute
          dU_de = self.cell.allostasis(self.pre_action,algebra_operator,i,u)
          gradient_descent_i.append(parameter_gamma[i][u] * dU_de)
        else:
          gradient_descent_i.append(0)
      #~~~~~~~~~~~~~~~~~~~  
      gradient_descent.append(gradient_descent_i)
    display_average_change("Allostasis matrix",gradient_descent,1)
    #~~~~~~~~~~~~~~~~~~~
    self.homeostasis(gradient_descent,identity)
    display_average_change("Allostasis matrix",gradient_descent,2)
    return self
#------------------------------------------------------------------------------  
  def spontaneous_reaction(self,vector,identity):
    self.pre_action = list()
    self.cell.spontaneous_reaction()
    if self.is_leaf:
      for i in range(len(self.cell.organelles)):
        self.pre_action.append(identity(self.cell.organelles[i]).action([vector]))
    else:
      for i in range(len(self.innercells)):
        self.innercells[i].spontaneous_reaction(vector,identity)
        self.pre_action.append(self.innercells[i].cell.action(self.innercells[i].pre_action))
        #Makes sure that the cell is in a homeostatic state
        self.cell.organelles[i] = self.innercells[i].cell.K
        self.cell.Sorg[i] = self.innercells[i].cell.SK
      self.cell.K = self.cell.content()
      self.cell.SK = sum(self.cell.K)
    return self
#------------------------------------------------------------------------------
  def tensor_pre_action(self):
    barycenters = [0] * self.cell.dimension
    #Computing prop increases the running time; instead we use self.Sorg[i]/float(self.SK)
    #prop = self.cell.content_proportions()
    if len(self.pre_action) == len(self.cell.organelles) and self.cell.SK != 0:
      for i in range(len(self.cell.organelles)):
        if len(self.pre_action[i]) == self.cell.dimension:
          for u in range(self.cell.dimension):
            #barycenters[u] += self.pre_action[i][u] * prop[i]
            barycenters[u] += self.pre_action[i][u] * (self.cell.Sorg[i]/float(self.cell.SK))
    return barycenters
#------------------------------------------------------------------------------
  def merge_base(self,list_of_organelles,tensor,order="non-sorted"):
    #If the following condition is removed, organelles of leaves can be merged (dynamic memory).
    if not(self.is_leaf):
      #~~~~~~~~~~~~~~~~~~~
      if order !="sorted":
        list_of_organelles.sort()
      #~~~~~~~~~~~~~~~~~~~
      self.cell.merge(list_of_organelles,"sorted")
      #~~~~~~~~~~~~~~~~~~~
      new_innercells = list()
      merged_innercells = list()
      #~~~~~~~~~~~~~~~~~~~
      new_pre_action = list()
      merged_pre_action = list()
      #~~~~~~~~~~~~~~~~~~~
      #Let j be the min of list_of_organelles
      j = list_of_organelles[0]
      for i in range(len(self.innercells)):
        if i in list_of_organelles:
          merged_innercells.append(self.innercells[i])
          merged_pre_action.append(self.pre_action[i])
          if i == j:
            new_innercells.append(None)
            new_pre_action.append(None)
        else:
          new_innercells.append(self.innercells[i])
          new_pre_action.append(self.pre_action[i])
      #~~~~~~~~~~~~~~~~~~~    
      new_innercells[j] = tensor(merged_innercells)
      new_innercells[j].pre_action = merged_pre_action
      new_innercells[j].set_level(self.level+1)
      new_pre_action[j] = new_innercells[j].tensor_pre_action()
      #new_pre_action[j] = new_innercells[j].cell.action(new_innercells[j].pre_action)
      #~~~~~~~~~~~~~~~~~~~
      self.innercells = new_innercells
      self.pre_action = new_pre_action
    return self
#------------------------------------------------------------------------------
  def divide_base(self,list_of_organelles,tensor):
    #~~~~~~~~~~~~~~~~~~~
    c_1,c_2 = self.cell.divide(list_of_organelles)
    #~~~~~~~~~~~~~~~~~~~
    new_innercells_1 = list()
    new_innercells_2 = list()
    #~~~~~~~~~~~~~~~~~~~
    new_pre_action_1 = list()
    new_pre_action_2 = list()
    #~~~~~~~~~~~~~~~~~~~
    j = 0
    for i in range(len(self.cell.organelles)):
      if j < len(list_of_organelles) and i == list_of_organelles[j]:
        if not(self.is_leaf):
          new_innercells_1.append(self.innercells[i])
        new_pre_action_1.append(self.pre_action[i])
        j += 1
      else:
        if not(self.is_leaf):
          new_innercells_2.append(self.innercells[i])
        new_pre_action_2.append(self.pre_action[i])
    #~~~~~~~~~~~~~~~~~~~
    sc_1 = SuperCell(c_1,new_innercells_1)
    sc_2 = SuperCell(c_2,new_innercells_2)
    #~~~~~~~~~~~~~~~~~~~
    sc_1.pre_action = new_pre_action_1 
    sc_2.pre_action = new_pre_action_2
    #~~~~~~~~~~~~~~~~~~~
    #here, we use pseudo_tensor
    t = tensor([sc_1,sc_2])
    t.set_level(self.level)
    self.depth          = t.depth
    self.level          = t.level
    self.cell           = t.cell
    self.is_leaf        = t.is_leaf
    self.compose_state  = t.compose_state
    self.innercells     = t.innercells
    self.pre_action     = [sc_1.tensor_pre_action(),sc_2.tensor_pre_action()]
    #self.pre_action = [sc_1.cell.action(sc_1.pre_action),sc_2.cell.action(sc_2.pre_action)]
    #~~~~~~~~~~~~~~~~~~~
    return self
#------------------------------------------------------------------------------
  def _fusion(self,tensor,filtering = 1.5):
    if not(self.is_leaf):
      #~~~~~~~~~~~~~~~~~~~
      for i in range(len(self.innercells)):
        self.innercells[i]._fusion(tensor,filtering)
      #~~~~~~~~~~~~~~~~~~~
      mer = self.cell.proposed_clustering(self.pre_action,"merging",filtering)
      if len(mer) in [0,len(self.cell.organelles)]:
        print "Warning: in SuperCell._fusion: total fusion canceled at level", self.level
      else:
        self.merge_base(mer,tensor)
        print "Merging occurring due to organelle(s):", mer, "at level =",self.level
#------------------------------------------------------------------------------
  def fusion(self,vector,operad,filtering  = 1.5):
    #spontaneous_reaction is necessary to update pre_action
    self.spontaneous_reaction(vector,operad.identity)
    self._fusion(operad.merging_tensor,filtering)
    self.reset_depth()
    return self
#------------------------------------------------------------------------------ 
  def _fission(self,tensor,filtering = 1.5):
    if not(self.is_leaf):
      #~~~~~~~~~~~~~~~~~~~
      for i in range(len(self.innercells)):
        self.innercells[i]._fission(tensor,filtering)
      #~~~~~~~~~~~~~~~~~~~
    div = self.cell.proposed_clustering(self.pre_action,"division",filtering)
    if len(div) in [0,len(self.cell.organelles)]:
      print "Warning: in SuperCell._fission: total fission canceled at level", self.level
    else:
      self.divide_base(div,tensor)
      print "Division occurring due to organelle(s):", div, "at level =",self.level
#------------------------------------------------------------------------------
  def fission(self,vector,operad,filtering  = 1.5):
    #spontaneous_reaction is necessary to update pre_action
    self.spontaneous_reaction(vector,operad.identity)
    self._fission(operad.dividing_tensor,filtering)
    self.reset_depth()
    return self
#------------------------------------------------------------------------------      
  def _compose(self,identity):
    if self.is_leaf:
      return self
    else:
      new_innercells = list()
      new_pre_action = list()
      valuation = list()
      index = 0
      is_new_leaf = True
      #~~~~~~~~~~~~~~~~~~~
      for i in range(len(self.innercells)):
        is_new_leaf = is_new_leaf and self.innercells[i].is_leaf and self.innercells[i].compose_state
      #~~~~~~~~~~~~~~~~~~~
      if is_new_leaf:
        
        for i in range(len(self.innercells)):
            new_innercells.extend([])
            new_pre_action.extend(self.innercells[i].pre_action)
            valuation.append([index,self.innercells[i].cell])
            index += len(self.innercells[i].cell.organelles)
      #~~~~~~~~~~~~~~~~~~~      
      else:
      
        for i in range(len(self.innercells)):  
          #~~~~~~~~~~~~~~~~~~~
          comp_i = self.innercells[i]._compose(identity)
          #~~~~~~~~~~~~~~~~~~~
          if comp_i.compose_state and comp_i.is_leaf:
            for j in range(len(comp_i.cell.organelles)):
              new_innercells.append(SuperCell(identity(comp_i.cell.organelles[j])))
            new_pre_action.extend(comp_i.pre_action)
            valuation.append([index,comp_i.cell])
            index += len(comp_i.cell.organelles)
          #~~~~~~~~~~~~~~~~~~~  
          elif comp_i.compose_state and not(comp_i.is_leaf):
            new_innercells.extend(comp_i.innercells)
            new_pre_action.extend(comp_i.pre_action)
            valuation.append([index,comp_i.cell])
            index += len(comp_i.innercells)
          #~~~~~~~~~~~~~~~~~~~
          else:
            new_innercells.append(comp_i)
            new_pre_action.append(self.pre_action[i])
            valuation.append([index,identity(comp_i.cell.K)])
            index += 1
      #~~~~~~~~~~~~~~~~~~~
      for k,c in valuation:
        self.cell = self.cell.compose(k,c)
      self.innercells = new_innercells
      self.is_leaf = is_new_leaf
      self.pre_action = new_pre_action
    return self
#------------------------------------------------------------------------------   
  def compose(self,identity):
    self._compose(identity)
    self.reset_depth()
    self.set_level(0)
    return self
#------------------------------------------------------------------------------
