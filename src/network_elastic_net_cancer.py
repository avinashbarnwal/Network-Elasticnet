from snap import *
from cvxpy import *
import numpy as np
from numpy import linalg as LA
import math
from multiprocessing import Pool
import csv
import os
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
import matplotlib
matplotlib.use('Agg')
matplotlib.rc('text',usetex=True)
import matplotlib.pyplot as plt
from z_u_solvers_Elastic_net import solveZ, solveU
import sys
import time
start_time=time.time()
import random
import pandas as pd
import pp


def solveX(data):
	nrow = data.shape[0]
	ncol = data.shape[1]
	class_no = data.shape[0]
	inputs = int(data[0][ncol-1])
	alpha  = data[0][ncol-2]
	lamb = data[0][ncol-3]
	rho = data[0][ncol-4]
	sizeData = int(data[0][ncol-5])
	mu = data[0][ncol-6]
	x = data[:,0:inputs]
	a = data[:,inputs:(inputs + sizeData)]
	neighs = data[:,(inputs + sizeData):(ncol-6)]
	xnew = cvxpy.Variable(class_no,inputs)
	#Fill in objective function here! Params: Xnew (unknown), a (side data at node)

	response_var_pos = a.shape[1]-1
	response = a[0,response_var_pos]
	y_response = numpy.zeros((class_no,1))

	class_type = 9

	for i in range(1,class_type+1):
		if response == i:
			y_response[i]=1


	row_score = numpy.zeros((class_no,1))
	exp_row_score = numpy.zeros((class_no,1))
	pro_score = numpy.zeros((class_no,1))

	for i in range(class_no):
		for j in range(inputs-1):
			row_score[i] = row_score[i] + xnew[i,j]*a[i,j]
		row_score[i] = row_score[i] + xnew[i,inputs-1]
		
	for i in range(class_no):
		exp_row_score[i] = math.exp(row_score[i])

	total_sum = numpy.sum(exp_row_score)

	g = 0
	for i in range(class_no):
		pro_score[i] = exp_row_score[i]/total_sum
		g = g + y_response[i]*math.log(pro_score[i])
	
	g = -g
	h = 0

	for i in range((neighs.shape[1])/(2*inputs+1)):
		weight = neighs[0,i*(2*inputs+1)]
		if(weight != 0):
			u = neighs[:i*(2*inputs+1)+1:i*(2*inputs+1)+(inputs+1)]
			z = neighs[:i*(2*inputs+1)+(inputs+1):(i+1)*(2*inputs+1)]
			h = h + rho/2*cvxpy.square(cvxpy.norm(xnew - z + u,"fro"))
	objective = cvxpy.Minimize(50*g+50*h)
	constraints = []
	p = cvxpy.Problem(objective, constraints)
	result = p.solve(solver=cvxpy.SCS)
	if(result == None):
		#CVXOPT scaling issue. Rarely happens (but occasionally does when running thousands of tests)
		objective = cvxpy.Minimize(51*g+52*h)
		p = cvxpy.Problem(objective, constraints)
		result = p.solve(verbose=False)
		if(result == None):
			print "SCALING BUG"
			objective = cvxpy.Minimize(52*g+50*h)
			p = cvxpy.Problem(objective, constraints)
			result = p.solve(verbose=False)
	return xnew.value, g.value

def solveZ(data):
	nrow = data.shape[0]
	ncol = data.shape[1]
	class_no = data.shape[0]
	inputs = int(data[0][ncol-1])
	alpha  = data[0][ncol-2]
	lamb = data[0][ncol-3]
	rho = data[0][ncol-4]
	useConvex = data[0][ncol-5]
	epsilon = data[0][ncol-6]
	weight = data[0][ncol-7]

	x1 = data[:,0:inputs]
	x2 = data[:,inputs:2*inputs]
	u1 = data[:,2*inputs:3*inputs]
	u2 = data[:,3*inputs:4*inputs]
	a =  x1 + u1
	b =  x2 + u2

	z1 = numpy.zeros((class_no,inputs))
	z2 = numpy.zeros((class_no,inputs))
	
	#(z1,z2)=(a,b)
	c1=lamb*(1-alpha)*weight
	c2=lamb*alpha*weight
	mu1=numpy.linalg.norm((rho*(a-b)-2*c2)/(rho))-2*c1/(rho)
	mu2=numpy.linalg.norm((rho*(a-b)+2*c2)/(rho))-2*c1/(rho)
	gamma1=mu1*c2/(2*c1+mu1*rho)
	gamma2=mu2*c2/(2*c1+mu2*rho)
	theta1=0.5+mu1*rho/(4*c1+2*mu1*rho)
	theta2=0.5+mu2*rho/(4*c1+2*mu2*rho)
	
	
	for iter_row in range(class_no):
		for iter_col in range(ncol)
			temp1=theta1*a[iter_row][iter_col]+(1-theta1)*b[iter_row][iter_col]-gamma1
			temp2=(1-theta1)*a[iter_row][iter_col]+theta1*b[iter_row][iter_col]+gamma1
			temp3=theta2*a[iter_row][iter_col]+(1-theta2)*b[iter_row][iter_col]+gamma2
			temp4=(1-theta2)*a[iter_row][iter_col]+theta2*b[iter_row][iter_col]-gamma2
			if (temp1>temp2):
				z1[iter_row][iter_col] = temp1
				z2[iter_row][iter_col] = temp2
			elif (temp3<temp4):
				z1[iter_row][iter_col] = temp3
				z2[iter_row][iter_col] = temp4
			else:
				z1[iter_row][iter_col] = z2[iter_row][iter_col] = 0.5*(a[iter_row][iter_col]+b[iter_row][iter_col])
	

	return z1,z2



def solveU(data):
	ncol = data.shape[1]
	u = data[:,0:leng/3]
	x = data[:,leng/3:2*leng/3]
	z = data[:,(2*leng/3):leng]
	rho = data[0][ncol-1]
	return u + (x - z)



def runADMM(G1, sizeOptVar, sizeData, lamb, rho, numiters, x, u, z, a, edgeWeights, useConvex, epsilon, mu, alpha):
	

	nodes = G1.GetNodes()
	edges = G1.GetEdges()
	#Find max degree of graph; hash the nodes
	(maxdeg, counter) = (0, 0)
	node2mat = snap.TIntIntH()
	for NI in G1.Nodes():
		maxdeg = max(maxdeg, NI.GetDeg())
		node2mat.AddDat(NI.GetId(), counter)
		counter = counter + 1
	#Stopping criteria
	eabs = math.pow(10,-2)
	erel = math.pow(10,-3)
	(r, s, epri, edual, counter) = (1,1,0,0,0)
	A = numpy.zeros((2*edges,class_no, nodes))
	for EI in G1.Edges():
		A[2*counter,class_no,node2mat.GetDat(EI.GetSrcNId())] = 1
		A[2*counter+1,class_no,node2mat.GetDat(EI.GetDstNId())] = 1 
		counter = counter+1
	(sqn, sqp) = (math.sqrt(nodes*sizeOptVar), math.sqrt(2*sizeOptVar*edges))

	#Run ADMM
	iters = 0
	#flag=1
	while(iters <= numiters and (r > epri or s > edual)):
		neighs = numpy.zeros((nodes,class_no,(2*sizeOptVar+1)*maxdeg))
		edgenum = 0
		numSoFar = snap.TIntIntH()
		for EI in G1.Edges():
			if (not numSoFar.IsKeye(EI.GetSrcNId())):
				numSoFar.AddDat(EI.GetSrcNId(), 0)
			counter  = node2mat.GetDat(EI.GetSrcNId())
			counter2 = numSoFar.GetDat(EI.GetSrcNId())
 			neighs[counter,:,counter2*(2*sizeOptVar+1)] = edgeWeights.GetDat(snap.TIntPr(EI.GetSrcNId(), EI.GetDstNId()))


 			neighs[counter,:,counter2*(2*sizeOptVar+1)+1:counter2*(2*sizeOptVar+1)+(sizeOptVar+1)] = u[2*edgenum,:,:]
 			neighs[counter,:,counter2*(2*sizeOptVar+1)+(sizeOptVar+1):(counter2+1)*(2*sizeOptVar+1)] = z[2*edgenum,:,:]
			numSoFar.AddDat(EI.GetSrcNId(), counter2+1)

			if (not numSoFar.IsKey(EI.GetDstNId())):
				numSoFar.AddDat(EI.GetDstNId(), 0)
			counter  = node2mat.GetDat(EI.GetDstNId())
			counter2 = numSoFar.GetDat(EI.GetDstNId())
 			neighs[counter,:,counter2*(2*sizeOptVar+1)] = edgeWeights.GetDat(snap.TIntPr(EI.GetSrcNId(), EI.GetDstNId()))
 			neighs[counter,:,counter2*(2*sizeOptVar+1)+1:counter2*(2*sizeOptVar+1)+(sizeOptVar+1)] = u[2*edgenum+1,:,:]
 			neighs[counter,:,counter2*(2*sizeOptVar+1)+(sizeOptVar+1):(counter2+1)*(2*sizeOptVar+1)] = z[2*edgenum+1,:,:]
			numSoFar.AddDat(EI.GetDstNId(), counter2+1)

			edgenum = edgenum+1

		const_x = [mu, sizeData,rho,lamb,alpha,sizeOptVar]
		const_xtemp = np.zeros((nodes,class_no,6))
		for i in range(nodes):
			for j in range(class_no)
				const_xtemp[i,j,0] = cons_x

		temp = numpy.concatenate((x,a,neighs,const_xtemp), axis=2)

		newx = numpy.zeros((nodes,class_no,sizeOptVar))

		for i in range(temp.shape[0]):
			flag1,flag2 = solveX(temp[i,:,:])
			newx[i,:,:] = flag1
		x = newx



		#z-update
		ztemp1 = numpy.zeros((edges,class_no,2*sizeOptVar))
		utemp1 = numpy.zeros((edges,class_no,2*sizeOptVar))

		for i in range(edges):
			ztemp1[i,:,0:sizeOptVar] = z[2*i,:,:]
			ztemp1[i,:,(sizeOptVar:2*sizeOptVar)] = z[2*i+1,:,:]

		for i in range(edges):
			utemp1[i,:,0:sizeOptVar] = u[2*i,:,:]
			utemp1[i,:,(sizeOptVar:2*sizeOptVar)] = u[2*i+1,:,:]	
		
		xtemp = numpy.zeros((2*edges,class_no,sizeOptVar))
		counter = 0
		weightsList = numpy.zeros((edges,class_no,1))

		for EI in G1.Edges():
			xtemp[2*counter,:,:] = numpy.array(x[node2mat.GetDat(EI.GetSrcNId()),:,:])
			xtemp[2*counter+1,:,:] = numpy.array(x[node2mat.GetDat(EI.GetDstNId()),:,:])
			weightsList[counter,:,0] = edgeWeights.GetDat(snap.TIntPr(EI.GetSrcNId(), EI.GetDstNId()))
			counter = counter+1

		xtemp1 = numpy.zeros((edges,class_no,2*sizeOptVar))
		for i in range(edges):
			xtemp1[i,:,0:sizeOptVar] = xtemp[2*i,:,:]
			xtemp1[i,:,(sizeOptVar:2*sizeOptVar)] = xtemp[2*i+1,:,:]

		xtemp = xtemp1	
		const_z = [epsilon, useConvex, rho,lamb,alpha,sizeOptVar]
		const_ztemp = np.zeros((edges,class_no,1))
		for i in range(nodes):
			for j in range(class_no):
				const_ztemp[i,j,:] = cons_z


		temp = numpy.concatenate((xtemp,utemp,ztemp,weightsList,const_ztemp), axis=2)
		flag1=0
		newz=np.zeros((2*edges,class_no,sizeOptVar))
		for i in range(temp.shape[0]):
			flag1,flag2=solveZ(temp[i,:,:])
			newz[2*i,:,:] = flag1
			newz[2*i+1,:,:] =flag2
		ztemp = newz
		s = numpy.linalg.norm(rho*numpy.dot(A.transpose(),(ztemp - z).transpose())) #For dual residual
		z = ztemp

		#u-update
		(xtemp, counter) = (numpy.zeros((2*edges,class_no,sizeOptVar)), 0)
		for EI in G1.Edges():
			xtemp[2*counter,:,:] = numpy.array(x[node2mat.GetDat(EI.GetSrcNId()),:,:])
			xtemp[2*counter+1,:,:] = numpy.array(x[:,node2mat.GetDat(EI.GetDstNId()),:,:])
			counter = counter + 1

		
		const_utemp = np.zeros((2*edge,class_no,1))
		const_utemp = rho
		temp = numpy.concatenate((u, xtemp, z, const_utemp), axis=2)

  
		newu=np.zeros((2*edges,class_no,sizeOptVar))
		for i in range(temp.shape[0]):
			newu[i,:,:]=solveU(temp[i,:,:])
		
		u = newu


		#Stopping criterion - p19 of ADMM paper
		epri = sqp*eabs + erel*max(numpy.linalg.norm(numpy.dot(A,x.transpose()), 'fro'), numpy.linalg.norm(z, 'fro'))
		edual = sqn*eabs + erel*numpy.linalg.norm(numpy.dot(A.transpose(),u.transpose()), 'fro')
		r = numpy.linalg.norm(numpy.dot(A,x.transpose()) - z.transpose(),'fro')
		s = s

		#print r, epri, s, edual
		iters = iters + 1

	return (x,u,z,0,0)


def main(alpha):
	#Set parameters
	useConvex = 1
	rho = 0.001
	numiters = 50
	thresh = 1000#10000
	lamb = 0.0
	startVal = 0.5#0.01
	useMult = 1 #1 for mult, 0 for add
	addUpdateVal = 0.1
	multUpdateVal = 1.5#1.1
	mu = 0.5 #For LS regularization
	#Test/Validation Set Information
	numNeighs = 5 #For data we keep
	testSetSize = 20
	validationSetSize = 0
	#For test/validation nodes
	numNewNeighs = 5 
	#Size of x
	sizeOptVar = 66
	#Size of side information at each node
	sizeData = 67
	#Number of Classes
	class_no = 9
	epsilon = 0.01
	no_var = sizeOptVar - 1
	index_counter = 1 
	dependent =  1
	
	var_start = index_counter 
	var_end = var_start + no_var
	#Generate graph, edge weights
	print("var_end")
	print(var_end)
	
	file = open("test_predictive_cancer.csv", "rU")
	file.readline() #ignore first line
	G1 = snap.TUNGraph.New()
	gene_expression = snap.TIntFltVH()
	dataset = snap.TIntFltVH()
	counter = 0
	num_patients = sum(1 for line in open("test_predictive_cancer.csv", "rU"))-1#Remove header

	flag = 1
	
	for line in file:
		G1.AddNode(counter)
		tempData = snap.TFltV()
		for i in range(var_start,var_end):
			tempData.Add(float(line.split(",")[i]))
		tempData.Add(float(line.split(",")[0]))
		tempData.Add(float(line.split(",")[var_end])) #19 for normalized; 5 for raw
		dataset.AddDat(counter, tempData)
		counter = counter + 1
		
	#Remove random subset of nodes for test and validation sets
	correlation_Coeff = pd.read_csv("correlation.csv",sep=",")
	numpy.random.seed(10)
	testid=numpy.random.choice(counter-1,testSetSize,replace=False)
	#print(testid)
	testList = snap.TIntV()
	for i in range(testSetSize):
		temp = testid[i]
		G1.DelNode(temp)
		testList.Add(temp)

	validationList = snap.TIntV()
	for i in range(validationSetSize):
		temp = G1.GetRndNId()
		G1.DelNode(temp)
		validationList.Add(temp)

	#For each node, find closest neightbors and add edge, weight = 5/distance
	edgeWeights = snap.TIntPrFltH()
        flag_flag=0
	
	for NI in G1.Nodes():
		this_node = NI.GetId()
		correlation_temp = correlation_Coeff[NI.GetId()]
		correlation_arg = numpy.argsort(correlation_Coeff[NI.GetId()])[::-1]
		correlation_arg = numpy.delete(correlation_arg,this_node,axis=0)
		correlation_arg = [ua for ua in correlation_arg if ua not in testid]
		#print(correlation_arg)
		for j in range(numNeighs):
			if (not G1.IsEdge(NI.GetId(),correlation_arg[j])):
				G1.AddEdge(NI.GetId(),correlation_arg[j])
				temp = snap.TIntPr(min(NI.GetId(),correlation_arg[j]), max(NI.GetId(),correlation_arg[j]))
				edgeWeights.AddDat(temp,correlation_temp[correlation_arg[j]])

	nodes = G1.GetNodes()
	edges = G1.GetEdges()
	print nodes, edges
	print snap.GetBfsFullDiam(G1, 1000, False);

	#Get side information
	a = numpy.zeros((sizeData, nodes))
	temp_a = numpy.zeros((nodes,class_no,sizeData))
	counter = 0
	for NI in G1.Nodes():
		for i in range(sizeData):
			a[i,counter] = dataset.GetDat(NI.GetId())[i]
		counter = counter + 1

	for i in range(nodes):
		for j in range(class_no):
			temp_a[i,j,:] = a[:,i]

	a = temp_a
	#Initialize variables to 0
	x = numpy.zeros((nodes,class_no,sizeOptVar))
	u = numpy.zeros((2*G1.GetEdges(),class_no,sizeOptVar))
	z = numpy.zeros((2*G1.GetEdges(),class_no,sizeOptVar))
	

	#Run regularization path
	[plot1, plot2, plot3] = [snap.TFltV(), snap.TFltV(), snap.TFltV()]
	while(lamb <= thresh):
		(x, u, z, pl1, pl2) = runADMM(G1, sizeOptVar, sizeData, lamb, rho + math.sqrt(lamb), numiters, x, u ,z, a, edgeWeights, useConvex, epsilon, mu,alpha)
		
		logloss = 0
		
		#Calculate accuracy on test set
		for i in testList:
			
			correlation_temp = correlation_Coeff[i]
			correlation_temp = numpy.delete(correlation_temp,testid,axis=0)
			correlation_arg = numpy.argsort(correlation_temp)[::-1]
			#correlation_arg = correlation_arg[1:]
			xpred = cvxpy.Variable(class_no,sizeOptVar)
		
	
			g = 0
			for j in range(numNewNeighs):
				weight = correlation_temp[correlation_arg[j]]
				g = g + weight*cvxpy.norm(xpred - x[correlation_arg[j],:,:],"fro")
			objective = cvxpy.Minimize(g)
			constraints = []
			p = cvxpy.Problem(objective, constraints)
			result = p.solve(verbose=False)
			xpred = xpred.value

			regressors = dataset.GetDat(i)
			response = float(dataset.GetDat(i)[var_end-1]
			y_response = np.zeros((class_no,1))

			class_type = 9

			for k in range(1,(class_type+1)):
				if response == k:
					y_response[k]=1

			row_score = np.zeros((class_no,1))
			exp_row_score = np.zeros((class_no,1))
			pro_score = np.zeros((class_no,1))

			for l in range(class_no):
				for m in range(sizeOptVar-1):
					row_score[l] = row_score[l] + xpred[l,m]*regressors[m]
				row_score[l] = row_score[l] + xpred[l,sizeOptVar-1]
		
			for q in range(class_no):
				exp_row_score[q] = exp(row_score[q])

			total_sum = numpy.sum(exp_row_score)

			g_g = 0
			for r in range(class_no):
				pro_score[r] = exp_row_score[r]/total_sum
				g_g = g_g + y_response[i]*log(pro_score[r])
	
			g_g = -g_g
			#Find log loss

			logloss = logloss + g_g

			
		lambda_array = numpy.zeros((nodes,class_no,1))
		lambda_array = lamb
		x_temp = x

		x_coeff = numpy.concatenate((x_temp,a),axis=2)
		x_coeff = numpy.concatenate((x_coeff,lambda_array),axis=2)

		cons = 0

		for i in range(edges):
			if(numpy.all(z[2*i,:,:] == z[2*i+1,:,:])):
				cons = cons + 1
		consensus = cons / float(edges)
		print "Lambda =", round(lamb,3),"Alpha =",alpha,"MSE=",round(mse,3), "Consensus=",round(consensus,3)
		plot1.Add(lamb)
		plot2.Add(logloss)
		plot3.Add(consensus)
		if(lamb == 0):
			lamb = startVal
		elif(useMult == 1):
			lamb = lamb*multUpdateVal
		else:
			lamb = lamb + addUpdateVal

	#Print/Save plot of results
	if(thresh > 0):
		pl1 = numpy.array(plot1)
		pl2 = numpy.array(plot2)
		pl3 = numpy.array(plot3)
		#plt.plot(pl1, pl3)

		if(useConvex == 1):
			#plt.savefig('consensus_gene_convex',bbox_inches='tight')
			file_1 = "CANCER_" + "_" + str(alpha) + ".out"
			file_2 = "CANCER_" + "_" + str(alpha) + ".csv"
			numpy.savetxt(file_1, (pl1, pl2, pl3), delimiter=',', fmt='%1.4f')
			numpy.savetxt(file_2,x_coeff,delimiter=',')
			


# tuple of all parallel python servers to connect with
ppservers = ()
#ppservers = ("10.0.0.1",)

if len(sys.argv) > 1:
    ncpus = int(sys.argv[1])
    # Creates jobserver with ncpus workers
    job_server = pp.Server(ncpus, ppservers=ppservers)
else:
    # Creates jobserver with automatically detected number of workers
    job_server = pp.Server(ppservers=ppservers)

print "Starting pp with", job_server.get_ncpus(), "workers"


alpha = [0,0.2,0.4,0.6,0.8,1]
ls_time = []


jobs = [(input, job_server.submit(main,(input,), (runADMM,solveX,solveZ,solveU,), ("math","snap","cvxpy","numpy","multiprocessing","csv","os","tempfile","matplotlib","sys","time","random","pandas",))) for input in alpha]
for input, job in jobs:
    print "Sum of primes below", input, "is", job()

print "Time elapsed: ", time.time() - start_time, "s"
job_server.print_stats()


