import time

class _DebugTime:
  #------------------------------------------------------------------------------
  def __init__(self):
    self.event_names = list()
    self.iteration = 1
    self.global_time = time.time() 
  #------------------------------------------------------------------------------  
  def set(self,*event_names):
    self.event_names = ["Set_Debug_Time"]
    self.iteration = 1
    if len(event_names) > 0:
      self.event_names = list(event_names)
      print "\n[\033[94m"+str(self.event_names[0])+"\033[0m]"
    self.global_time = time.time()
    return self
  #------------------------------------------------------------------------------  
  def call(self):
    call_time = time.time()
    event_name = self.event_names[self.iteration % len(self.event_names)]
    time_display = " | time: " + str(call_time - self.global_time)+" sec\033[0m]"
    print"[\033[94mCHECK"+str(self.iteration)+" "+ event_name + time_display
    self.iteration += 1
    return self
  #------------------------------------------------------------------------------  
debug_time = _DebugTime()
