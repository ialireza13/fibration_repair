#the model is as follows:
#Given a directed network and a coloring of the network, 
#find the minimum number of edges to add or remove
#to ensure the coloring is inbalanced 
import networkx as nx
import kamei_cock as kamei
import numpy as np
import gurobipy as gp
from gurobipy import GRB
from gurobipy import abs_
from gurobipy import quicksum
import random
import pandas as pd
import itertools as itools
from collections import defaultdict
import time
import os
from contextlib import redirect_stdout

##precision parameter
epsilon = .001

#add and remove flag constants -- not to be changed
ADDONLY = 1
RMONLY = 2
BOTHADDRM = 3

charsep='\t'

def suppress_gurobi_output():
	"""Temporarily suppress Gurobi output."""
	return open(os.devnull, 'w')

#this code calculates the minimal number of balanced partitions 
##gg should be an nx.DiGraph
def FindMP(gg):
	gm = np.transpose(nx.adjacency_matrix(gg))
	network = kamei.CoupledCellNetwork(gm.todense())
	grain = network.top_lattice_node() #most partitions
	return(len(set(grain[0]))), network


def read_data(fname,colorfile,xlinks=None):
	
	#load starting network
	GraphData = pd.read_csv(fname,sep=charsep,index_col=[0,1],header=None, \
							comment="#")
	
	#check for duplicate edges
	assert not max(GraphData.index.duplicated()), 'Duplicate edges found.'
	
	#directed graph!
	EdgeDict = GraphData.to_dict()[2]
	
	#get the edge set, Edges, and edge weights
	edges, edge_weights = gp.multidict(EdgeDict)
			
	#cdict keys are read in color set
	ctable=pd.read_csv(colorfile,index_col=0,sep=charsep,header=None,\
					   comment="#")
		
	cdict = ctable.to_dict()[1]
	
	# create a list of nodes
	nodes = list({item for item in cdict.keys()} | {item for tup in EdgeDict.keys() for item in tup})
	
	# print("Read graph: n="+str(len(nodes))+" m="+str(len(edges)))

	#set up a list of colorsets
	color_dict = defaultdict(list)

	for key, value in cdict.items():
		color_dict[value].append(key)

	color_sets = list(color_dict.values())
	
	color_pairs =[]
	for C in color_sets:
		for p,q in itools.combinations(C,2):
			color_pairs.append((p,q))            

	# print("Read colors")
	
	rev_color_dict = {string: key for key, values in color_dict.items() for string in values}
			
	nc_tuples = []
	outer_imbalance_dict = defaultdict(dict)
	inner_imbalance_dict = defaultdict(dict)
	support_num=0
	
	for C,D in itools.combinations(color_sets,2):
		for p in C:
			p_color = rev_color_dict[p]
			for q in D:
				support_num = support_num + 1
				q_color = rev_color_dict[q]
				inner_imbalance_dict[p][q]=[p_color,q_color]
				base_colors = list(color_dict.keys())
				base_colors.remove(p_color)
				base_colors.remove(q_color)
				outer_imbalance_dict[p][q]=base_colors
				
				for c in color_dict.keys():
					nc_tuples.append((p,q,c))
					nc_tuples.append((q,p,c))

	all_pairs = [(i, j) for i in nodes for j in nodes if i != j]

	if xlinks != None:
		# print(xlinks)
		prohibited = pd.read_csv(xlinks,sep=charsep,index_col=[0,1],header=None)
		
		#directed graph!
		non_existing_EdgeDict = prohibited.to_dict()[2]
		non_existing_EdgeDict.update(EdgeDict)
		
		edges_to_avoid = non_existing_EdgeDict.copy()
		avoid_edges, ae_edge_weights = gp.multidict(edges_to_avoid)
		
		not_e = {(p,q):1 for (p,q) in all_pairs if (p,q) not in avoid_edges}
		not_edges,ne_weights = gp.multidict(not_e)      
	else:
		not_e = {(p,q):1 for (p,q) in all_pairs if (p,q) not in edges}
		not_edges,ne_weights = gp.multidict(not_e)
	
	inputs = {}     
	inputs['nodes'] = nodes
	inputs['edges'] = edges
	inputs['color_pairs'] = color_pairs
	inputs['color_sets'] = color_sets
	inputs['not_edges'] = not_edges
	inputs['color_dict'] = color_dict
	inputs['nc_tuples'] = nc_tuples
	inputs['outer_imbalance_dict'] = outer_imbalance_dict
	inputs['inner_imbalance_dict'] = inner_imbalance_dict
	inputs['support_num'] = support_num
	
	return inputs


##create a MIP
def CreateRMIP(inputs,env,Imbalance,HardFlag,FixedEdges,FixedNonEdges,AddRemoveFlag,InDegOneFlag):
	
	rmip = gp.Model(name='RepairKnown-Directed',env=env)

	nodes=inputs['nodes']
	edges=inputs['edges']
	color_pairs=inputs['color_pairs']
	color_sets=inputs['color_sets']
	not_edges=inputs['not_edges']
	color_dict=inputs['color_dict']
	nc_tuples=inputs['nc_tuples']
	outer_imbalance=inputs['outer_imbalance_dict']
	inner_imbalance_dict=inputs['inner_imbalance_dict'] 
	support_num=inputs['support_num']
	
	#initialize edge variables
	remove_edge=rmip.addVars(edges,vtype=GRB.BINARY,name='remove_edge')
	node_balance_pos = rmip.addVars(color_pairs,lb=0.0,vtype=GRB.CONTINUOUS,name='node_balance_pos')
	node_balance_neg = rmip.addVars(color_pairs,lb=0.0,vtype=GRB.CONTINUOUS,name='node_balance_neg')
	max_nodebalance = rmip.addVar(lb=0.0,vtype=GRB.CONTINUOUS,name='max_nodebalance')
	add_edge=rmip.addVars(not_edges,vtype=GRB.BINARY,name='add_edge')
	strict_balance = rmip.addVars(nc_tuples,vtype=GRB.BINARY,name='strict_balance')
	
	if Imbalance=='v2':
		auxiliary_var_1 = rmip.addVars(support_num,lb=-2,ub=2,vtype=GRB.SEMIINT,name='out_imbalance_one')
		auxiliary_var_2 = rmip.addVars(support_num,lb=-2,ub=2,vtype=GRB.SEMIINT,name='out_imbalance_two')

	rvars = {'re':remove_edge,'nb_p':node_balance_pos,\
			 'nb_n':node_balance_neg,'m_nb':max_nodebalance,\
				 'ae':add_edge,'sb':strict_balance}

	#constraint: colors in-balanced
	color_balance = []
	color_imbalance = []
	one_imbalance = []
	atleast_one = []
	indeg_one = []
	n = len(nodes)

	if InDegOneFlag:
			indeg_one.append(rmip.addConstrs((sum((1-remove_edge[i,j]) for (i,j) in edges if j == p) \
												+ sum(add_edge[i,j] for (i,j) in not_edges if j == p) >= 1 for p in nodes), name='indeg_one'))
	
	if HardFlag:
#########   
		for D in color_sets:
			for (p,q) in color_pairs:
				A = list((1-remove_edge[i,j]) for (i,j) in edges if j == p and i in D)
				B = list(add_edge[i,j] for (i,j) in not_edges if j == p and i in D)
				a = list((1-remove_edge[i,j]) for (i,j) in edges if j == q and i in D)
				b = list(add_edge[i,j] for (i,j) in not_edges if j == q and i in D)
				color_balance.append(rmip.addConstr((quicksum(A) + quicksum(B) == quicksum(a) + quicksum(b)), name='color_balance'+str(p)+'_'+str(q)))
				
		counter=0
		
		for C,D in itools.combinations(color_sets,2):
			for p in C:
				for q in D:
					for c in color_dict.keys():
						A = list((1-remove_edge[i,j]) for (i,j) in edges if j == p and i in color_dict[c])
						B = list(add_edge[i,j] for (i,j) in not_edges if j == p and i in color_dict[c])
						a = list((1-remove_edge[i,j]) for (i,j) in edges if j == q and i in color_dict[c])
						b = list(add_edge[i,j] for (i,j) in not_edges if j == q and i in color_dict[c])
						
						color_imbalance.append(rmip.addConstr((quicksum(A) + quicksum(B) >= quicksum(a) + quicksum(b) + strict_balance[p,q,c] - n*strict_balance[q,p,c]), name='imbalance_'+str(p)+'_'+str(q)+'_'+str(c)))
						color_imbalance.append(rmip.addConstr((quicksum(a) + quicksum(b) >= quicksum(A) + quicksum(B) + strict_balance[q,p,c] - n*strict_balance[p,q,c]), name='imbalance_'+str(q)+'_'+str(p)+'_'+str(c)))
						
						one_imbalance.append(rmip.addConstr((1 >= strict_balance[p,q,c] + \
															 strict_balance[q,p,c]) ,name='one_imbalance_'+str(p)+'_'+str(q)+'_'+str(c)))
					
					if Imbalance=='v1':                         
						atleast_one.append(rmip.addConstr((sum(strict_balance[p,q,i] for i in color_dict.keys()) +\
													   sum(strict_balance[q,p,i] for i in color_dict.keys()) >= 1),name='atleast_one_'+str(p)+'_'+str(q)))
					elif Imbalance=='v2':
						
						A = list(strict_balance[p,q,i] - strict_balance[q,p,i] for i in inner_imbalance_dict[p][q])
						rmip.addConstr(auxiliary_var_1[counter] == quicksum(A))
						rmip.addConstr(auxiliary_var_2[counter] == abs_(auxiliary_var_1[counter]))
						
						B = list(strict_balance[p,q,i] + strict_balance[q,p,i] for i in outer_imbalance[p][q])
						
						atleast_one.append(rmip.addConstr((quicksum(B) +\
											auxiliary_var_2[counter] >= 1),name='atleast_one_'+str(p)+'_'+str(q)))
					
					counter += 1
########
	else:
		for D in color_sets:
			color_balance.append(rmip.addConstrs((sum((1-remove_edge[i,j]) for (i,j) in edges \
										 if j == p and i in D) \
									+ sum(add_edge[i,j] for (i,j) in not_edges \
										  if j == p and i in D) - \
									sum((1-remove_edge[i,j]) for (i,j) in edges \
										if j == q and i in D) - \
									sum(add_edge[i,j] for (i,j) in not_edges \
										if j == q and i in D)
										 == \
									node_balance_pos[p,q] - node_balance_neg[p,q]\
									) for (p,q) in color_pairs))

	FElist = []
	for (i,j) in FixedEdges:
		FElist.append(rmip.addConstr(remove_edge[i,j]==0))

	FNElist = []
	for (i,j) in FixedNonEdges:
		FNElist.append(rmip.addConstr(add_edge[i,j]==0))        

	if AddRemoveFlag == RMONLY:
		for (i,j) in not_edges:
			FElist.append(rmip.addConstr(add_edge[i,j]==0))

	if AddRemoveFlag == ADDONLY:
		for (i,j) in edges:
			FElist.append(rmip.addConstr(remove_edge[i,j]==0))

	
	#keep track of edges/potential edges that are perturbed
	nodebalance_bounds_p = rmip.addConstrs((node_balance_pos[p,q] <= max_nodebalance \
										  for (p,q) in color_pairs))
	nodebalance_bounds_n = rmip.addConstrs((node_balance_neg[p,q] <= max_nodebalance \
										  for (p,q) in color_pairs))
	
	rcons={'cb':color_balance,'nb_b_p':nodebalance_bounds_p,\
		   'nb_b_n':nodebalance_bounds_n,'FEl':FElist,'FNEl':FNElist,\
			   'indeg_one':indeg_one}
		
	return rmip,rcons,rvars,remove_edge,add_edge,node_balance_pos,node_balance_neg

def set_rmip(graphpath,colorpath,Imbalance,HardFlag,\
				 FixedEdges,FixedNonEdges,InDegOneFlag,AddRemoveFlag,prohibit,verbose):
	
	#create the inputs
	# print("Reading data from " + graphpath + " and " + colorpath)
	inputs = read_data(graphpath,colorpath,prohibit)
	
	nodes = inputs['nodes']
	edges = inputs['edges']
	color_pairs = inputs['color_pairs']
	not_edges = inputs['not_edges']
	color_dict = inputs['color_dict']
	color_sets = inputs['color_sets']
	
	#set dictionary
	setdict = {'N':nodes,'E':edges,'CP':color_pairs,'NE':not_edges,'cd':color_dict}
	
	#initialize an environment
	env = gp.Env()
	env.setParam('OutputFlag', 1 if verbose else 0)
	
	#create the model
	# print("Creating model")
	rmip,rcons,rvars,remove_edge,add_edge,node_balance_pos,node_balance_neg = \
		CreateRMIP(inputs,env,Imbalance,HardFlag,FixedEdges,FixedNonEdges,AddRemoveFlag,InDegOneFlag)

	return rmip,rcons,rvars,setdict,color_sets,remove_edge,add_edge,node_balance_pos,node_balance_neg

def rmip_optimize(rmip,rcons,rvars,remove_edge,add_edge,node_balance_pos,node_balance_neg,rm_weight,add_weight,HardFlag,Solu_type,bal_weight=1):
	
	#need objective
	if HardFlag:
		if Solu_type=="Abs":
			#abs version, experimental
			w = rmip.addVar(vtype=gp.GRB.INTEGER, name="w")  
			rmip.addConstr(w >= rm_weight * gp.quicksum(remove_edge.select('*','*')) - add_weight * gp.quicksum(add_edge.select('*','*')))
			rmip.addConstr(w >= - rm_weight * gp.quicksum(remove_edge.select('*','*')) + add_weight * gp.quicksum(add_edge.select('*','*')))
		elif Solu_type=="Linear": #normal version used in our papers
			obj = rm_weight*gp.quicksum(remove_edge.select('*','*')) + \
				add_weight*gp.quicksum(add_edge.select('*','*'))
	
	else:
		obj = (epsilon + rm_weight)*(gp.quicksum(remove_edge.select('*','*'))) + \
			(epsilon + add_weight)*(gp.quicksum(add_edge.select('*','*'))) + \
			bal_weight*(gp.quicksum(node_balance_pos.select('*','*')) + \
			   gp.quicksum(node_balance_neg.select('*','*'))) 
		
		
	if Solu_type=="Abs":
		#abs version
		rmip.setObjective(w,GRB.MINIMIZE)
	elif Solu_type=="Linear":
		rmip.setObjective(obj,GRB.MINIMIZE)
	
	rmip.setParam('Threads', 0)
	# rmip.setParam('Presolve', 2)
	# rmip.setParam('Heuristics', 0.3)
	#optimize
	startTime_Prime = time.time()
	rmip.optimize()
	executionTime = round(time.time() - startTime_Prime,5)

	return rmip,rcons,rvars,executionTime


def solve_and_write(graphpath,colorpath,rm_weight,add_weight,fname,rmip,rcons,\
					   rvars,setdict,colorsets,remove_edge,add_edge,node_balance_pos,node_balance_neg,Solu_type,\
					   HardFlag=True,FixedEdges=[],FixedNonEdges=[],InDegOneFlag=True,\
					   prohibit=None,Save_info=True,NetX=False):
	
	rmip,rcons,rvars,executionTime = rmip_optimize(rmip,rcons,rvars,remove_edge,\
												   add_edge,node_balance_pos,node_balance_neg,rm_weight,add_weight,HardFlag,Solu_type,bal_weight=1)
	
	
	#find the edge removes
	E = setdict['E']
	NE = setdict['NE']

	cd = setdict['cd']

	re = rvars['re']
	ae = rvars['ae']
	sb = rvars['sb']
	sumremovals = 0
	sumadds = 0
	idealnum=len(colorsets)
	feasible = (rmip.Status == GRB.OPTIMAL)
	
	G_result = nx.DiGraph()

	if NetX==True:
		if feasible:
			
			G_result.add_edges_from([edge for edge in E if abs(re[edge[0], edge[1]].x - 1) > epsilon])
			sumremovals += len([edge for edge in E if abs(re[edge[0], edge[1]].x - 1) < epsilon])
			
			G_result.add_edges_from([edge for edge in NE if abs(ae[edge[0], edge[1]].x - 1) < epsilon])
			sumadds += len([edge for edge in NE if abs(ae[edge[0], edge[1]].x - 1) < epsilon])
			
	if Save_info==True:
		outfname = fname+"directed.output.txt"
		f = open(outfname,"w")
		gname = fname+"directed.out.graph.txt"
		gf = open(gname,"w")
	
		print(f'Total edges removed\n{sumremovals}',file=f)
		print('Edges removed',file=f)
		EdgesRemoved = []
		if feasible:
			for (i,j) in E:
				if abs(re[i,j].x - 1) < epsilon:
					print(f'{i} {j}',file=f)
					EdgesRemoved.append((i,j))
				else:
					print(f'{i} {j}',file=gf)
	
		print(f'Total edges added\n{sumadds}',file=f)
		print('Edges added',file=f)
		EdgesAdded = []
		if feasible:
			for (i,j) in NE:
				if abs(ae[i,j].x - 1) < epsilon:
					print(f'{i} {j}',file=f)
					print(f'{i} {j}',file=gf)
					EdgesAdded.append((i,j))
	
	
		CP = setdict['CP']
		m_nb = rvars['m_nb']
		nb_p = rvars['nb_p']
		nb_n = rvars['nb_n']
		if feasible:
			print(f'Maximum imbalance\n{m_nb.x}',file=f)
			
		print('Nonzero imbalances',file = f)    
		
		if feasible:
			for (i,j) in CP:
				imbalance = nb_p[i,j].x - nb_n[i,j].x
				if abs(imbalance) > epsilon:
					print(f'{i} {j} {imbalance}',file=f)
				
		print('\nImbalances for each node and color',file=f)
		
		if feasible:        
			for C, D in itools.combinations(colorsets, 2):
				for p in C:
					for q in D:
						# Log the imbalance between `p` and `q`
						print(f'Imbalances between {p} and {q}', file=f)
						
						# Filter and print relevant colors directly
						imbalanced_colors = [i for i in cd if sb[p, q, i].x == 1 or sb[q, p, i].x == 1]
						for color in imbalanced_colors:
							print(f'Color {color}', file=f)
	
		print("\n\n",end="",file=f)
		print("Input graph",file=f)
		GraphData = pd.read_csv(graphpath,sep=charsep,index_col=[0,1],header=None, \
								comment="#")
		GraphData.to_csv(f,sep=' ')
		
	
		print("Input colors",file = f)
		ctable=pd.read_csv(colorpath,index_col=0,sep=charsep,header=None, \
								comment="#")    
		ctable.to_csv(f,sep=' ')
		
		if prohibit!=None:
			print("Prohibited edges",file = f)
			prohibited = pd.read_csv(prohibit,sep=charsep,index_col=[0,1],header=None)
			prohibited.to_csv(f,sep=' ')
		
		f.close()
		gf.close()
	
		
	else:
		gname=[]; EdgesRemoved=[]; EdgesAdded=[]; outfname=[];

	minp, net = FindMP(G_result)
	
	# print(f"Found {minp} colors, minimal is {idealnum}")

	return gname,idealnum,EdgesRemoved,EdgesAdded,sumremovals,sumadds,outfname,rmip,rcons,rvars,G_result,executionTime
 
 
def repair_network(color_file_path, instance_file_path, output_file_path, alpha, beta, prohibit_file_path=None, verbose=False):
	
	with suppress_gurobi_output() as f, redirect_stdout(f):
		from settings import param_data
		
		rm_add_flag = param_data["rm_add_flag"]        
		if rm_add_flag == 'add_only':
			RM_AD = ADDONLY
		elif rm_add_flag == "rm_only":
			RM_AD = RMONLY
		elif rm_add_flag == "both":
			RM_AD = BOTHADDRM
		rmip,B,C,D,E,F,G,H,I = set_rmip(instance_file_path,color_file_path,param_data["model_type"],\
													param_data["hard_flag"],[],[],param_data["InDegOneFlag"],\
													RM_AD,prohibit_file_path,verbose)
		
		rmip.setParam("MIPGap",param_data["mip_gap"])
		rmip.setParam("Seed", random.randint(1, 1000000))
		gname,idealnum,EdgesRemoved,EdgesAdded,sumremovals,sumadds,outfname,rmip,rcons,rvars,G_result,executionTime = solve_and_write(instance_file_path,\
										color_file_path,alpha,beta,output_file_path,rmip,B,C,D,E,F,G,H,I,\
										"Linear",param_data["hard_flag"],[],[],param_data["InDegOneFlag"],prohibit_file_path,\
										Save_info=param_data["save_output"],NetX=True)
		return EdgesRemoved,EdgesAdded, G_result