import sys
import math

class _Useful:
  #------------------------------------------------------------------------------
  def contrast(self,brightness,profiles):
  
    def contrast_profiles(brightness,profiles,vector,*challenge):
      maximum = 0
      for u in range(len(vector)):
        if len(challenge) == 1 and challenge[0][0] <= u % challenge[0][2] <= challenge[0][1]:
          continue
        elif maximum < vector[u]:
          maximum = vector[u]

      counts = [0] * len(brightness)
      total = 0
      for u in range(len(vector)):
        if len(challenge) == 1 and challenge[0][0] <= u % challenge[0][2] <= challenge[0][1]:
          continue
        total += 1.
        pixel_u = vector[u]/float(maximum)
        for j in range(len(counts)):
          if pixel_u > brightness[j]:
            counts[j] += 1

      percents = map(lambda x: x/total,counts)
      contrast_tests = list()
      
      for i in range(len(profiles)):
        profile_value = True
        for j in range(len(profiles[i])):
          profile_value = profile_value and profiles[i][j][0] <= percents[j] <= profiles[i][j][1]
        contrast_tests.append(profile_value)
      
      percents = map(lambda x: str(int(10000.*x)/10000.),percents)
      print "[brightness]  " + " ".join(percents)
      print " ".join(map(str,contrast_tests))+"\n"
      return contrast_tests

    return (lambda vector,*challenge : contrast_profiles(brightness,profiles,vector,*challenge))
  #------------------------------------------------------------------------------
  def gamma(self,E,F,brightness,profiles,scores,*challenge):
  
    def gamma_parameter(c,a):
      parameters = list()
      
      maximum = 0
      for i in range(len(c.organelles)):
      
        contrast = self.contrast(brightness,profiles)(c.organelles[i],*challenge)
        agreement = c.agreement(i,a[i],*challenge)
        #if True in contrast:
          #agreement = c.agreement(i,a[i])
        #else:
          #agreement = c.agreement(i,a[i],*challenge)
        good_scenario = True
        
        for j in range(len(profiles)):
          good_scenario = good_scenario and not(contrast[j] and agreement < scores[j])
        good_agreement = good_scenario * agreement
        maximum = max(maximum,good_agreement)
        
        row = list()
        for u in range(len(c.organelles[i])):
          if len(challenge) == 1 and challenge[0][0] <= u % challenge[0][2] <= challenge[0][1]:
            row.append(0)
          else:
            row.append(good_agreement)
          #row.append([good_agreement,i,agreement])
          
        parameters.append(row)
        
      maximum = 1 if maximum == 0 else maximum
      for i in range(len(parameters)):
        for u in range(len(parameters[i])):
          parameters[i][u] = (10 **E) * ((parameters[i][u]/maximum)**F)
          #parameters[i][u][0] = (10 **E) * ((parameters[i][u][0]/maximum)**F)
          #parameters[i][u] = [parameters[i][u][0],parameters[i][u][2],parameters[i][u][2]/maximum]
      
      return parameters
      
    return (lambda c, a : gamma_parameter(c,a))
  #------------------------------------------------------------------------------
  def print_data(self,inputs,image_size,output_file,option = 'display'):
    
    if option == 'save':
      output_file.write("\t".join(map(lambda x: str(int(x)),inputs)))
    
    else:
      display = list()
      
      for i in range(image_size[0]):
        image_row = inputs[image_size[1]*image_size[2]*i:image_size[1]*image_size[2]*(i+1)]
        
        row = list()
        for j in range(image_size[1]):
          #black and white tranformation
          intensity = 0
          for k in range(image_size[2]):
            intensity += image_row[image_size[2]*j+k]
          row.append(intensity)
          
        display.append(row)

      for i in range(len(display)):
        for j in range(len(display[i])):
          if display[i][j] == 0:
            string = "."
          else:
            string = str(int(display[i][j]))
          r = len(string)
          output_file.write(string+" "*(4-r))
        output_file.write("\n")
        
    return self  
  #------------------------------------------------------------------------------
  def print_root(self,super_cell,image_size,output_file,option = 'display'):
    if super_cell.level == 0:
      self.print_data(super_cell.cell.K,image_size,output_file,option)
      output_file.write("\n")
    return self 
  #------------------------------------------------------------------------------
  def print_organelles(self,super_cell,image_size,output_file,option = 'display'):
    if super_cell.is_leaf:
      for i in range(len(super_cell.cell.organelles)):
        self.print_data(super_cell.cell.organelles[i],image_size,output_file,option)
        output_file.write("\n")
    else:
      for i in range(len(super_cell.innercells)):
        self.print_organelles(super_cell.innercells[i],image_size,output_file,option)
    return self 
  #------------------------------------------------------------------------------
  def print_tree(self,super_cell,image_size,output_file,option = 'display',total = 0,local = 0):

    if super_cell.level == 0:
      total = super_cell.level + super_cell.depth
      local = total
      
    diff = local - (super_cell.level + super_cell.depth)
    for i in range(diff):
      output_file.write("@"+str(total - (super_cell.depth + diff - i)))
      output_file.write(":"+str(super_cell.depth + diff - i))
      output_file.write(":1")
      output_file.write(":"+str(super_cell.cell.residual)) #RESIDUAL
      output_file.write("\n")
      self.print_data(super_cell.cell.K,image_size,output_file,option)
      output_file.write("\n")
        
    output_file.write("@"+str(total - super_cell.depth))
    output_file.write(":"+str(super_cell.depth))
    output_file.write(":"+str(len(super_cell.cell.organelles)))
    output_file.write(":"+str(super_cell.cell.residual)) #RESIDUAL
    output_file.write("\n")
    self.print_data(super_cell.cell.K,image_size,output_file,option)
    output_file.write("\n")
       
    if super_cell.is_leaf:
    
      output_file.write("@"+str(super_cell.level))
      output_file.write(":-1")
      output_file.write(":"+str(len(super_cell.cell.organelles)))
      output_file.write(":1")
      output_file.write("\n")
      
      for i in range(len(super_cell.cell.organelles)):
        self.print_data(super_cell.cell.organelles[i],image_size,output_file,option)
        output_file.write("\n")
        
    else:
      local = super_cell.level + super_cell.depth
      
      for i in range(len(super_cell.innercells)):
        self.print_tree(super_cell.innercells[i],image_size,output_file,option,total,local)
        
    return self 
  #------------------------------------------------------------------------------
  def parse_image(self,string,separators,EOL_symbols):
    words = list()
    i = 0
    while i < len(string): 
      read = string[i]
      if not(read in EOL_symbols+['']):
        word = ''
        while not(read in separators) and not(read in EOL_symbols+['']) and i < len(string):
          word = word + read
          i += 1
          read = string[i]
        if word != '':
          words.append(float(word))
      i += 1
    return words 
  #------------------------------------------------------------------------------
  def get_memory(self,fload,batch,cycles):
  
    output = list()
    memory = list()
    cycle = 0
    
    for i, s in enumerate(fload):

      if cycle in cycles:
        memory.append(self.parse_image(s,['\t'],['\n']))
      
      if i - batch * cycle == batch -1:
        cycle += 1
        
        if memory != []:
          output.append(memory)
          memory = list()

      if cycle > max(cycles):
        break

    return output
  #------------------------------------------------------------------------------
  def get_trees(self,fload,cycles):
  
    output = list()
    tree = list()
    cycle = -1
    info = list()
    
    for i, s in enumerate(fload):
      
      if s[0] == '@':
        info = map(int,self.parse_image(s,['@',':'],['\n']))
        #construct a new tree
        if info != [] and info[0] == 0 and info[1] >= 0:
          cycle += 1
          if tree!=[] and tree[0] != []:
            output.append(tree)
          tree = list()
          for i in range(info[1]+2):
            tree.append([])

      elif cycle in cycles:
        
        if len(info) == 3:
          depth    = info[1]
          children = info[2]
          if depth == -1:
            tree[1+depth].append(self.parse_image(s,['\t'],['\n']))
          else:
            tree[1+depth].append([children,self.parse_image(s,['\t'],['\n'])])
        
        elif len(info) >= 4: #RESIDUAL
          depth    = info[1] #RESIDUAL
          children = info[2] #RESIDUAL
          residual = info[3] #RESIDUAL
          if depth == -1: #RESIDUAL
            tree[1+depth].append(self.parse_image(s,['\t'],['\n'])) #RESIDUAL
          else: #RESIDUAL
            content = self.parse_image(s,['\t'],['\n']) #RESIDUAL
            tree[1+depth].append([children,content,1+residual/sum(content)]) #RESIDUAL
            print "[DEBUG] learning stress["+str(1+depth)+"]:"+ str(1+residual/sum(content))
      
      
      if cycle > max(cycles):
        break
        
    return output
  #------------------------------------------------------------------------------ 
  def get_last_cycle(self,fload,ary):
    count = -1
    for i,s in enumerate(fload):
      if i % ary == 0:
        count += 1
    return count-1
  #-------------------------------------------------------------------------------------------  
  def join_images(self,image1,image2,spacing,explicit = 1,center = False):
    white_color = [500,500,500]
    grey_color = [200,200,200]
    new_image = list()
    if center:
      down1 = (len(image2)-len(image1))/2 if len(image2)-len(image1)> 0 else 0
      down2 = (len(image1)-len(image2))/2 if len(image1)-len(image2)> 0 else 0
    else:
      down1 = 0
      down2 = 0
    m = max(down1+len(image1),down2+len(image2))
    len1 = len(image1[0])
    len2 = len(image2[0])
    for i in range(m):
    
      if down1 <= i < down1+len(image1):
        im1 = image1[i-down1]
      else:
        im1 = [white_color] * len1
        
      if down2 <= i < down2+len(image2):
        im2 = image2[i-down2]
      else:
        im2 = [white_color] * len2
      
      if (down1 <= i < down1+len(image1) and down2 <= i < down2+len(image2) and explicit != 0) \
      or explicit == 2:
        delimitation = [white_color]*(spacing/2)+[grey_color]*2+[white_color]*(spacing/2)
      else:
        delimitation = [white_color]*(2*(spacing/2)+2)
        
      new_image.append(im1+delimitation+im2)
    
    return new_image
  #------------------------------------------------------------------------------ 
  def rgb_colormap(self,maximum,grayscale = False,residual_heat = 1.0): #RESIDUAL
    rescale = lambda t: int(200*t/float(maximum))
    #RGB codes [50,*,100] give blue/yellow contrast
    contrast = lambda x : [200-rescale(x[0])]*3 if grayscale else [int(residual_heat*50),rescale(x[0]),100] #RESIDUAL
    return lambda x : map(rescale,x) if len(x) != 1 else contrast(x)
  #------------------------------------------------------------------------------
  def make_rgb_panel(self,memory,image_size,grayscale = False,residual_heat = []): #RESIDUAL
    black_color = [1,1,1]
    rgb_panel = list()
    
    for j in range(image_size[0]):
      row = list()
      for i in range(len(memory)):
      
        line = memory[i][image_size[1]*image_size[2]*j:image_size[1]*image_size[2]*(j+1)]
        maximum = max(memory[i])
        
        for k in range(image_size[1]):
          pixel_k = line[image_size[2]*k:image_size[2]*(k+1)]
          if residual_heat != []: #RESIDUAL
            row.append(self.rgb_colormap(maximum,grayscale,residual_heat[i])(pixel_k)) #RESIDUAL
          else: #RESIDUAL
            row.append(self.rgb_colormap(maximum,grayscale)(pixel_k)) #RESIDUAL
        if i != len(memory)-1:
          row.append(black_color)
          
      rgb_panel.append(row)
      
    return rgb_panel
  #------------------------------------------------------------------------------
  def make_rgb_grid(self,memory,image_size,grayscale = False):
    black_color = [1,1,1]
    rgb_grid = list()
    
    horizontal_line = 0
    for i in range(len(memory)):
      panel = self.make_rgb_panel(memory[i],image_size,grayscale)
      
      if horizontal_line == 0 and panel != []:
        horizontal_line = len(panel[0])
        
      rgb_grid.extend(panel)
      rgb_grid.append([black_color]*horizontal_line)
        
    return rgb_grid
  #------------------------------------------------------------------------------
  def make_rgb_tree(self,tree,image_size,grayscale = False):
    white_color = [500,500,500]
    black_color = [1,1,1]
    rgb_tree = list()
    
    leaves = tree[0]
    leaf_number = len(leaves)
    rgb_tree.extend(self.make_rgb_panel(leaves,image_size,grayscale))
    
    
    positions_images = list()
    for i in range(leaf_number):
      half = image_size[1]/2
      positions_images.append([half,image_size[1]-half])
    
    #~~~~~~~~~~~~~~
    
    for index in range(len(tree)-1):
    
      branching_part = list()
   
      vertical_branche = lambda i : [white_color]*(positions_images[i][0]-1) + \
                                    [black_color] + \
                                    [white_color]*(positions_images[i][1])  
      
      #~~~~~~~~~~~~~~
      
      #building branches: | | | | ...
      vertical_branches = list()                   
      for i in range(len(positions_images)):                   
        vertical_branches.extend(vertical_branche(i))
        if i != len(positions_images)-1:
          vertical_branches.append(white_color)
      
      #building branches: | | | | ... x  height
      height =  5
      branching_part.extend([vertical_branches] * height)
      
      #getting data for building branches: ---  ---
      horizontal_cover = list() 
      parent_images = list()
      residual_heat = list() #RESIDUAL
      for i in range(len(tree[1+index])):
        horizontal_cover.append(tree[1+index][i][0])
        parent_images.append(tree[1+index][i][1])
        if len(tree[1+index][i]) >= 3: #RESIDUAL
          residual_heat.append(tree[1+index][i][2]) #RESIDUAL


      #~~~~~~~~~~~~~~
      
      #building branches: ---  ---
      horizontal_bars = list()
      sum_covers = 0
      for i in range(len(horizontal_cover)):
      
        if i != 0:
          horizontal_bars.append(white_color)
          
        for j in range(horizontal_cover[i]):
          if j == 0 and horizontal_cover[i] > 1:
            horizontal_bar = [white_color]*(positions_images[sum_covers+j][0]-1) + \
                             [black_color] + \
                             [black_color]*(positions_images[sum_covers+j][1]+1)
                 
          elif j == 0 and horizontal_cover[i] == 1:
            horizontal_bar = [white_color]*(positions_images[sum_covers+j][0]-1) + \
                             [black_color] + \
                             [white_color]*(positions_images[sum_covers+j][1])                   

          elif 0 < j < horizontal_cover[i]-1:
            horizontal_bar = [black_color]*(positions_images[sum_covers+j][0]-1) + \
                             [black_color] + \
                             [black_color]*(positions_images[sum_covers+j][1]+1)
      
          elif j == horizontal_cover[i]-1:
            horizontal_bar = [black_color]*(positions_images[sum_covers+j][0]-1) + \
                             [black_color] + \
                             [white_color]*(positions_images[sum_covers+j][1])
              
          horizontal_bars.extend(horizontal_bar)
        sum_covers += horizontal_cover[i]
          
      branching_part.append(horizontal_bars)

      rgb_tree.extend(branching_part)
      
      #~~~~~~~~~~~~~~
        
      new_positions_images = list()
      shifting_positions = list()
      count = 0
      cover = 0
      j = 0
      for i in range(len(positions_images)):
        cover += positions_images[i][0]+positions_images[i][1]
        count += 1
        if count == horizontal_cover[j]:
          #~~~~~~~~~~~~~~
          half1 = cover/2
          half2 = cover - half1
          
          add_borders1 = (horizontal_cover[j]-1)/2
          add_borders2 = horizontal_cover[j]/2
          
          new_positions_images.append([half1 + add_borders1, half2 + add_borders2])
          
          border_image1 = half1 + add_borders1 - image_size[1]/2
          border_image2 = half2 + add_borders2 - (image_size[1] - image_size[1]/2)
          shifting_positions.append([border_image1,border_image2])
          
          cover = 0
          count = 0
          j += 1
          #~~~~~~~~~~~~~~
      positions_images = new_positions_images
      
      #~~~~~~~~~~~~~~
      rgb_panel = self.make_rgb_panel(parent_images,image_size,grayscale,residual_heat) #RESIDUAL
      shifted_rgb_panel = list()
      for k in range(len(rgb_panel)):
        row = list()
        for i in range(len(horizontal_cover)):
          shifted_image = [white_color] * shifting_positions[i][0] + \
                          rgb_panel[k][i+image_size[1]*i:i+image_size[1]*(i+1)] + \
                          [white_color] * shifting_positions[i][1]
          row.extend(shifted_image)
          if i != len(horizontal_cover)-1:
            row.append(white_color)
        shifted_rgb_panel.append(row)           
      
      #~~~~~~~~~~~~~~
      
      rgb_tree.extend(shifted_rgb_panel)

    return rgb_tree
#-------------------------------------------------------------------------------------------  
  def draw_arrow(self,panel_width,length):
    white_color = [500,500,500]
    black_color = [1,1,1]
    #draw the brackground
    line = list()
    x = panel_width/2
    for j in range(50):
      line.append([white_color]*panel_width)
    a = 25
    b = a+5
    length = min(a,length)
    #draw the tail of the arrow
    for i in range(length):
      line[a-length+i][x] = black_color
      line[a-length+i][x-1] = black_color
      line[a-length+i][x+1] = black_color
    #draw the head of the arrow
    for j in range(a,b): 
      line[j][x] = black_color
      line[j][x-1] = black_color
      line[j][x+1] = black_color
      
      line[j][x-(b-1-j)] = black_color
      line[j][x+(b-1-j)] = black_color
      line[j][x-1-(b-1-j)] = black_color
      line[j][x+1+(b-1-j)] = black_color
      if not(j in [a]):
        line[j][x-2-(b-1-j)] = black_color
        line[j][x+2+(b-1-j)] = black_color
      if not(j in [a,a+1]):
        line[j][x-3-(b-1-j)] = black_color
        line[j][x+3+(b-1-j)] = black_color
    line[b][x] = black_color
    line[b][x-1] = black_color
    line[b][x+1] = black_color
    line[b][x-2] = black_color
    line[b][x+2] = black_color
    line[b+1][x] = black_color
    line[b+1][x-1] = black_color
    line[b+1][x+1] = black_color
    line[b+2][x] = black_color
    
    return line
  #------------------------------------------------------------------------------
  def zero_matrix(self,dim):
    matrix = list()
    for i in range(dim):
      matrix_i = list()
      for j in range(dim):
        matrix_i.append(0)
      matrix.append(matrix_i)
    return matrix
  #------------------------------------------------------------------------------
  #CODE TAKEN FROM THE PEDIGRAD LIBRARY
  def join_fibers(self,fibers1,fibers2,speed = "REG"):
    tmp1 = list()
    tmp2 = list()
    
    for i in range(len(fibers1)):
      tmp1.append(list(set(fibers1[i])))
      
    for i in range(len(fibers2)):
      tmp2.append(list(set(fibers2[i])))
      
    for i1 in range(len(tmp1)):
      for j1 in range(len(tmp1[i1])):
        for i2 in range(len(tmp2)):
          flag = False
          for j2 in range(len(tmp2[i2])):
            if tmp1[i1][j1] == tmp2[i2][j2]:
              tmp1[i1].extend(tmp2[i2])
              tmp2[i2]=[]
              tmp1[i1] = list(set(tmp1[i1]))
              flag = True
              break
          #tmp1[i1][j1] no longer needs to be searched in tmp2.
          if speed == "FAST" and flag == True:
            break
      tmp2.append(tmp1[i1])
      tmp1[i1] = []
      
    the_join = list()
    for i in range(len(tmp2)):
      if tmp2[i]!=[]:
        the_join.append(tmp2[i])
        
    return the_join
  #------------------------------------------------------------------------------
  def cliques(self,a_matrix):
    m = 0
    fibers = list()
    print "Searching cliques in the following graph:"
    for i in range(len(a_matrix)):
      print a_matrix[i]
      for j in range(i,len(a_matrix[i])):
        if m < a_matrix[i][j]:
          m = a_matrix[i][j]
          fibers = [[i,j]]
        elif m == a_matrix[i][j]:
          fibers.append([i,j])
    print "Maximal weight = ", m
    print "Edges with maximal weight = ", fibers
    join = self.join_fibers(fibers,fibers,"REG")
    join.sort(key = lambda x :-len(x))
    return join
    
usf = _Useful()    
    
