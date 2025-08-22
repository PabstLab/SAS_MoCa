#!/usr/bin/python

import time

###########################################################################################
###########################################################################################
def TimeCount(t):	
	'Time counting'
	minutes=int(t/60)
	seconds=int(t-minutes*60)
	return (minutes,seconds)
	

###########################################################################################
###########################################################################################
#	
def FreePar (PAR,FIX,L_LIM_org,H_LIM_org):
	'Redefine all the sub-array related to free parameters'
	FREE = []
	L_LIM = []
	H_LIM = []

	for i in range (len(PAR)):
		if FIX[i]!="f":
			FREE.append(PAR[i])
			L_LIM.append(L_LIM_org[i])
			H_LIM.append(H_LIM_org[i])
			
	return ( FREE , L_LIM , H_LIM )		

###########################################################################################
###########################################################################################
def MergePar (PAR,FREE,FIX):
	'Define the merged parameter list to evaluate the function'
	CALC=[]
	TEMP_PAR = PAR.copy()
	TEMP_FREE = FREE.copy()
	for el in (FIX):
		if el=="f":
			CALC.append(TEMP_PAR.pop(0)) 
		else:
			CALC.append(TEMP_FREE.pop(0)) 
			TEMP_PAR.pop(0)
	return CALC

###########################################################################################
###########################################################################################
def PrintResults(NAME,PAR,FIX,FREE,X2):	
	'Print Results'

	if X2>1e+7 or X2<1e-4:
		print("Reduced X^2 ...... %.1e\n" %X2 )
	else:
		print("Reduced X^2 ...... %.7s\n" %X2 )
	
	shift=0
	TEMP_NAME=list(NAME)
	TEMP_PAR=list(PAR)
	for i in range (len(TEMP_PAR)):
		if FIX[i]=="f":
			del TEMP_NAME[i-shift]
			del TEMP_PAR[i-shift]
			shift=shift+1
	for i in range(len(FREE)):
		print(" %d.\t%4s ...... %.4e" %(i,TEMP_NAME[i],FREE[i]) )
	print ()
	return

###########################################################################################
###########################################################################################
def PrintResultsFile(title,folder,datafile,function,NAME,PAR_RES,FIX,X2):	
	'Print Results on a File'
	localtime = time.asctime( time.localtime(time.time()) )
	
	resfile="./"+folder+"/Results.dat"	
	with open(resfile, 'a') as fl:
		fl.writelines(localtime+"\n")
		fl.writelines("Title: %s\n" %title)
		fl.writelines("Data File: %s\n" %datafile)
		fl.writelines("Model: %s\n\n" %function)
		for i in range(len(NAME)) :
			fl.writelines( "%d\t%4s\t%.7f\n" %(i,NAME[i],PAR_RES[i]) ) 

###########################################################################################
###########################################################################################

