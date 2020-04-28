import sys
sys.path.insert(0, '../intcyt/useful/')
from useful import *
#------------------------------------------------------------------------------
def intcyt(operad,supercell,index,events,vector,gamma,filtering = [1.5]*2):
  debug_time.set()
  #~~~~~~~~~~~~~~~~~~~
  debug_time.set("Allostasis "+str(index))
  supercell.allostasis(vector,operad.identity,gamma)
  debug_time.call()
  #~~~~~~~~~~~~~~~~~~~
  if index >= events[0] and index % events[1] in events[2]:
    debug_time.set("Fission "+str(index))
    supercell.fission(vector,operad,filtering[0])
    debug_time.call() 
  #~~~~~~~~~~~~~~~~~~~
  elif index >= events[0] and index % events[1] in events[3]:
    debug_time.set("Fusion "+str(index))
    supercell.fusion(vector,operad,filtering[1])
    debug_time.call()
  #~~~~~~~~~~~~~~~~~~~
  elif index >= events[0] and index % events[1] in events[4]:
    debug_time.set("Compose "+str(index))
    supercell.compose(operad.identity)
    debug_time.call()
#------------------------------------------------------------------------------
