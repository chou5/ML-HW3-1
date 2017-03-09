import sys
import numpy as np

def read_fasta(input_file):
  fh = open(input_file, 'r')
  line = fh.readline()
  obs = []
  hidden = []
  index = 0
  
  while line:
    line = line.strip('\n')
    #print list(line)
    if '<>' in line:
      new_list = []
      new_list1 = []
      obs.append(new_list)
      hidden.append(new_list1)
      #print obs
      #print hidden
    elif 'end' in line:
      index = index + 1
    else:
      line.replace(" ", "")
      tmp = list(line)
      #print tmp
      obs[index].append(tmp[0])
      hidden[index].append(tmp[2])
    line = fh.readline()
  obs = filter(None, obs)
  hidden = filter(None, hidden)
  return obs, hidden 
  
def initial_prob(hidden):
  initial = {}

  helix = 0
  sheet = 0
  coil = 0
  for element in hidden:
    if element[0] == 'h':
      helix = helix + 1
    elif element[0]  == 'e':
      sheet = sheet + 1
    else:
      coil = coil + 1

  if helix == 0:
    initial['h'] = 0.00000000000000000001
  else:
    initial['h'] = float(helix)/len(hidden)
  if sheet == 0:
    initial['e'] = 0.00000000000000000001
  else:
    initial['e'] = float(sheet)/len(hidden)
  if coil == 0:
    initial['_'] = 0.00000000000000000001
  else:
    initial['_'] = float(coil)/len(hidden)
  return initial
  
def trans_prob(hidden):
  trans = {}
  h_e = 0
  h_c = 0
  e_h = 0
  e_c = 0
  c_h = 0 
  c_e = 0 
  sum_h = 0
  sum_e = 0
  sum_c = 0

  for element in hidden:
    for i in range(-1, len(element)-1):
      if element[i+1] == 'h' and element[i] == 'e':
        h_e = h_e + 1
        sum_h = sum_h + 1
      elif element[i+1] == 'h' and element[i] == '_':
        h_c = h_c + 1
        sum_h = sum_h + 1
      elif element[i+1] == 'e' and element[i] == 'h':
        e_h = e_h + 1
        sum_e = sum_e + 1
      elif element[i+1] == 'e' and element[i] == '_':
        e_c = e_c + 1
        sum_e = sum_e + 1
      elif element[i+1] == '_' and element[i] == 'h':
        c_h = c_h + 1
        sum_c = sum_c + 1
      elif element[i+1] == '_' and element[i] == 'e':
        c_e = c_e + 1
        sum_c = sum_c + 1
  
  if h_e == 0:
    trans['he'] = 0.00000000000000000001
  else:
    trans['he'] = float(h_e)/sum_h
  if h_c == 0:
    trans['h_'] = 0.00000000000000000001
  else:
    trans['h_'] = float(h_c)/sum_h
  if e_h == 0:
    trans['eh'] = 0.00000000000000000001
  else:
    trans['eh'] = float(e_h)/sum_e
  if e_c == 0:
    trans['e_'] = 0.00000000000000000001
  else:
    trans['e_'] = float(e_c)/sum_e
  if c_h == 0:
    trans['_h'] = 0.00000000000000000001
  else:
    trans['_h'] = float(c_h)/sum_c
  if c_e == 0:
    trans['_e'] = 0.00000000000000000001
  else:
    trans['_e'] = float(c_e)/sum_c
	
  trans['hh'] = 0.00000000000000000001
  trans['ee'] = 0.00000000000000000001
  trans['__'] = 0.00000000000000000001
  return trans  
  
def emit_prob(obs, hidden):
  emit = {}
  helix = 0
  sheet = 0
  coil = 0
  
  for i in range(len(obs)):
    for j in range(len(obs[i])):
      state = hidden[i][j] + obs[i][j]
      if emit.has_key(state):
        emit[state] += 1
      else:
        emit[state] = 1
      if hidden[i][j] == 'h':
        helix += 1
      elif hidden[i][j] == 'e':
        sheet += 1
      else:
        coil += 1

  for k in emit:
    if k[0] == 'h':
      emit[k] = float(emit[k])/helix
    elif k[0] == 'e':
      emit[k] = float(emit[k])/sheet
    else:
      emit[k] = float(emit[k])/coil
  return emit
  
def duration(hidden):
  length = {}
  helix = []
  sheet = []
  coil = []
  
  for element in hidden:
    count_h = 0
    count_e = 0
    count_c = 0
    for i in range(len(element)-1):
      if element[i] == 'h' and element[i+1] == 'h':
        count_h += 1
      elif element[i] == 'h' and element[i+1] != 'h':
        count_h += 1
        helix.append(count_h)
        count_h = 0
      elif element[i] == 'e' and element[i+1] == 'e':
        count_e += 1
      elif element[i] == 'e' and element[i+1] != 'e':
        count_e += 1
        sheet.append(count_e)
        count_e = 0
      elif element[i] == '_' and element[i+1] == '_':
        count_c += 1
      elif element[i] == '_' and element[i+1] != '_':
        count_c += 1
        coil.append(count_c)
        count_c = 0

  #helix.sort()
  #sheet.sort()
  #coil.sort()
  #print helix
  #print sheet
  #print coil
  helix_max = max(helix)
  sheet_max = max(sheet)
  coil_max = max(coil)
  helix_prob = []
  sheet_prob = []
  coil_prob = []
  #print helix_max
  #print sheet_max
  #print coil_max
  
  for i in range(helix_max):
    num = helix.count(i+1)
    if num == 0:
      helix_prob.append(0.00000000000000000001)
    else:
      helix_prob.append(float(num)/len(helix))
  for i in range(sheet_max):
    num = sheet.count(i+1)
    if num == 0:
      sheet_prob.append(0.00000000000000000001)
    else:
      sheet_prob.append(float(num)/len(sheet))
  for i in range(coil_max):
    num = coil.count(i+1)
    if num == 0:
      coil_prob.append(0.00000000000000000001)
    else:
      coil_prob.append(float(num)/len(coil))

  helix_prob.insert(0,0.00000000000000000001)
  sheet_prob.insert(0,0.00000000000000000001)
  coil_prob.insert(0,0.00000000000000000001)
  

  norm_helix = normListSumTo(helix_prob, 1)
  norm_sheet = normListSumTo(sheet_prob, 1)
  norm_coil = normListSumTo(coil_prob, 1)
 
  length['h'] = norm_helix
  length['e'] = norm_sheet
  length['_'] = norm_coil

  return length
		
def normListSumTo(L, sumTo=1):
  sum = reduce(lambda x,y:x+y, L)
  return [x/(sum*1.0)*sumTo for x in L]
		
def ghmm(obs, state, init_p, trans_p, emit_p, length):
  ghmm = [[1 for x in range(len(obs))]for x in range(len(states))]
  path = [[1 for x in range(len(obs))]for x in range(len(states))]
  
  # Initialize base cases aa length = 0
  # row0 = h, row1 = e, row2 = _
  ghmm[0][0] = np.log10(init_p['h'])
  ghmm[1][0] = np.log10(init_p['e'])
  ghmm[2][0] = np.log10(init_p['_'])
  
  # Run ghmm for aa length > 0
  for ele in range(1, len(obs)):
    for i in range(len(states)):
      #print len(length[states[i]])
      if (ele+1) >= len(length[states[i]]):
        prob_all_one_state = ghmm[i][0] + np.log10(0.00000000000000000001)
      else:
        #print states[i]
        #print ele+1
        #print length[states[i]][ele+1]
        prob_all_one_state = ghmm[i][0] + np.log10(length[states[i]][ele+1])
      for a in obs[:ele+1]:
        emit = states[i]+a.upper()
        if emit not in emit_p.keys():
          prob_all_one_state += np.log(0.00000000000000000001)
        else:
          prob_all_one_state += np.log(emit_p[emit])
      ghmm[i][ele] = prob_all_one_state
      path[i][ele] = (i,0)
      for x in range(1,ele):
        for y in range(len(states)):
          check = 0
          if x != y:
            seq = obs[x+1:ele+1]
            if (ele-x) >= len(length[states[i]]):
              check += ghmm[y][x] + np.log10(trans_p[states[y]+states[i]]) + np.log10(0.000000001)
            else:
              check += ghmm[y][x] + np.log10(trans_p[states[y]+states[i]]) + np.log10(length[states[i]][ele-x])
            for z in seq:
              emit = states[y]+z.upper()
              if emit not in emit_p.keys():
                check += np.log10(0.00000000000000000001)
              else:
                check += np.log10(emit_p[emit])
            if ghmm[i][ele] < check:
              ghmm[i][ele] = check
              path[i][ele] = (y,x)
  return ghmm,path
		
def backtrack(test, emission_prob, pro_matrix, pos_matrix):
  length = len(test) - 1
  annotation = ''
  score = max(pro_matrix[0][length],pro_matrix[1][length],pro_matrix[2][length])
  if score == pro_matrix[0][length]:
    tmp = pos_matrix[0][length]
    state = 'h'
  elif score == pro_matrix[1][length]:
    tmp = pos_matrix[1][length]
    state = 'e'
  else:
    tmp = pos_matrix[2][length]
    state = '_'
  repeat = length - tmp[1]
  annotation += state*repeat
  pre_state = tmp[0]
  pre_pos = tmp[1]
  
  while True:
    tmp = pos_matrix[pre_state][pre_pos]
    if pre_state == 0:
      state = 'h'
    elif pre_state == 1:
      state = 'e'
    else:
      state = '_'
    repeat = pre_pos - tmp[1]
    annotation += state*repeat
    pre_state = tmp[0]
    pre_pos = tmp[1]
    if pre_pos == 0:
      annotation += state
      break
	  
  annotation = annotation[::-1]
  return annotation
		
def accuracy(real, predict):
  error_count = 0
  for i in range(len(real)):
    if predict[i] != real[i]:
      error_count += 1
  accuracy = 1-(float(error_count)/len(real))
  return accuracy
		
def print_outcome(obs_list, annotation_list, output_file):
  f = open(output_file, 'w')
  for j in range(len(obs_list)):
    list_ann = list(annotation_list[j])
    f.write("<>\n")
    for i in range(len(obs_list[j])):
      f.write(str(obs_list[j][i]))
      f.write(" ")
      f.write(str(list_ann[i]))
      f.write("\n")
    f.write("end\n")
	
# main function
# Set initial_prob, trans_prob, emit_prob and length_duration from trainning data
obs,  hidden = read_fasta("train.fasta")
init_p = initial_prob(hidden)
#print init_p
trans_p = trans_prob(hidden)
#print trans_p
emit_p = emit_prob(obs, hidden)
#print emit_p
length = duration(hidden)
#print length

if len(sys.argv) != 5:
		print "USAGE: python PreProStr_ghmm.py -i [INPUT_FILE] -o [OUTPUT_FILE]"
		exit(-1)

if sys.argv[1] == "-i":
		if sys.argv[2].endswith('.fasta'):
			input_file = sys.argv[2]
		else:
			print "The input file should be in fasta format"
			exit(-1)
			
if sys.argv[3] == "-o":
    output_file = sys.argv[4]  
else:
		print "USAGE: python PreProStr_ghmm.py -i [INPUT_FILE] -o [OUTPUT_FILE]"
		exit(-1)

# Get states and observation_seq from test data
states = ['h', 'e', '_']
test_obs, test_hidden = read_fasta(input_file)
accuracy_list = []
annotation_list = []

for i in range(len(test_obs)):
  ghmm_matrix, path_matrix = ghmm(test_obs[i], states, init_p, trans_p, emit_p, length)
  #print ghmm_matrix
  #print path_matrix
  annotation = backtrack(test_obs[i], emit_p, ghmm_matrix, path_matrix)
  #print annotation
  annotation_list.append(annotation)
  acc = accuracy(test_hidden[i], list(annotation))
  accuracy_list.append(acc)
#print annotation_list
#print test_obs
print_outcome(test_obs, annotation_list, output_file)
print accuracy_list
print "The average accuracy for the prediction is %0.4f" %np.average(accuracy_list)
print "The result of the prediction is in output file."

