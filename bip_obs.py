#calculates energy, mag^2, Ovbs, and VBS exponent
#  python bip_obs.py JQvbs_example
import os
import time
import datetime
import sys
import numpy as np

#if (len(sys.argv) < 2):
#	sys.exit("please input: (1) dir name  ex: J1J2_example")

if(len(sys.argv) > 2):
	drop=sys.argv[2]
else:
	drop='0'


# for old directory structure, not used here..
#dir_name = str(sys.argv[1])
#os.system('rm '+dir_name+'/obs 2>/dev/null')
#os.system('find `pwd`/'+dir_name+' -name *L*B*dmin* >tmp')
#f1array=open('tmp', 'r').read().split()
#N1=len(f1array)
N1 = 1 #assume serial
f1array = [os.getcwd()]
dir_name = os.getcwd()


#read jobs in sorted deltamin order
dvals=np.zeros(N1)
for i in range(N1):
	deltamin=float(open(f1array[i]+'/param.dat','r').read().split()[4])
	dvals[i]=deltamin
ind=np.argsort(dvals)


for j in range(N1):
	i=ind[j]
	start_time = time.time()
	# get the parameters from param.dat
	Lx=int(open(f1array[i]+'/param.dat','r').read().split()[0])
	Ly=int(open(f1array[i]+'/param.dat','r').read().split()[1])
	J=float(open(f1array[i]+'/param.dat','r').read().split()[2])
	Q=float(open(f1array[i]+'/param.dat','r').read().split()[3])
	deltamin=float(open(f1array[i]+'/param.dat','r').read().split()[4])
	deltamax=float(open(f1array[i]+'/param.dat','r').read().split()[5])
	Ndelta=int(open(f1array[i]+'/param.dat','r').read().split()[6])
	beta=float(open(f1array[i]+'/param.dat','r').read().split()[7])
	suN=int(open(f1array[i]+'/param.dat','r').read().split()[8])#should be 2
	eql=int(open(f1array[i]+'/param.dat','r').read().split()[9])
	mcs=int(open(f1array[i]+'/param.dat','r').read().split()[10])
	
	Nsite=1.0*Lx*Ly
	fac= 1.0/(beta*Nsite/2.0)

	#get values of delta
	tmpd=deltamin
	dvec=np.zeros(Ndelta)
	for n in range(Ndelta):
		dvec[n]=tmpd
		if Ndelta != 1:
			tmpd=tmpd*np.power(deltamax/deltamin, 1.0/(1.0*Ndelta - 1.0))



	#now copy over the data from $DATA (scratch folder)
	#os.system('cp /gpfs/data/fs71452/jondemidio/JQvbs2/'+f1array[i].split('data/')[1]+'/*.data '+f1array[i]+' 2>/dev/null;')
    
	#Nproc=48
	Nproc = 1 #assume 1 processor

	for c in range(Ndelta):
		delta=dvec[c]

		#gather all data
		for proc in range(Nproc//Ndelta):
			if proc == 0:
				#os.system('cat '+f1array[i]+'/'+str(c)+'.data > tmp_data 2>/dev/null;')
				os.system("awk 'FNR>"+drop+"' "+f1array[i]+"/"+str(c)+".data > tmp_data 2>/dev/null;")
			else:
				#os.system('cat '+f1array[i]+'/'+str(c + proc*Ndelta)+'.data >> tmp_data 2>/dev/null;')
				os.system("awk 'FNR>"+drop+"' "+f1array[i]+"/"+str(c + proc*Ndelta)+".data >> tmp_data 2>/dev/null;")



		if np.size(np.fromfile('tmp_data'))==0:
			print('skipping ')
			continue #skip this value of i if data set is empty


		data=np.loadtxt('tmp_data',ndmin=2)	
		Nbin =len(data[:,0])
		print(Nbin, " bins")

		Nfk = Nbin #number of fake data sets
		Nfk = 100	
	
		fkenergy=np.zeros(Nfk)
		fkmagsq=np.zeros(Nfk)
		fkOvbs=np.zeros(Nfk)
		fkdOvbs=np.zeros(Nfk)
		fkexp=np.zeros(Nfk)
		
		for f in range(Nfk):
		
			mybins = np.random.randint(Nbin,size=Nbin)
			fkenergy[f]=np.mean(data[mybins,0])
			fkmagsq[f]=np.mean(data[mybins,1])
			#fkOvbs[f]=np.mean(data[mybins,2]
			nxe = np.mean(data[mybins,2])
			nxo = np.mean(data[mybins,3])
			nxenxe = np.mean(data[mybins,4])
			nxonxo = np.mean(data[mybins,5])
			nxenxo = np.mean(data[mybins,6])

			
			Pe = nxe/(J+delta)
			Po = nxo/(J)

			####### the definition....
			fkOvbs[f] = (Pe-Po)/(beta*Nsite/2.0)
			fkdOvbs[f] = ((1.0/((J+delta)**2))*(nxenxe - nxe*nxe - nxe ) - (1.0/(J*(J+delta)))*(nxenxo - nxe*nxo ))/(beta*Nsite/2.0)
			fkexp[f] = delta*(fkdOvbs[f]/fkOvbs[f])
	


		#the statistical errors 
		errenergy = np.sqrt(np.var(fkenergy))
		errmagsq = np.sqrt(np.var(fkmagsq))
		errOvbs = np.sqrt(np.var(fkOvbs))
		errdOvbs = np.sqrt(np.var(fkdOvbs))
		errexp = np.sqrt(np.var(fkexp))
		
		#the mean values
		energy = np.mean(data[:,0])
		magsq = np.mean(data[:,1])
		#Ovbs = np.mean(data[:,2])
		nxe = np.mean(data[:,2])
		nxo = np.mean(data[:,3])
		nxenxe = np.mean(data[:,4])
		nxonxo = np.mean(data[:,5])
		nxenxo = np.mean(data[:,6])
		
		
		Pe = nxe/(J+delta)
		Po = nxo/(J)
		
		#the definition
		Ovbs = (Pe-Po)/(beta*Nsite/2.0)
		dOvbs = ((1.0/((J+delta)**2))*(nxenxe - nxe*nxe - nxe ) - (1.0/(J*(J+delta)))*(nxenxo - nxe*nxo ) )/(beta*Nsite/2.0)
		exp = delta*(dOvbs/Ovbs)
		
		
		outfile = open(dir_name+'/obs','a')
		outfile.write('%d %.6f %.12f %.6f ' % (Lx,J,delta,beta))
		outfile.write('%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f ' % (energy,errenergy,magsq,errmagsq,Ovbs,errOvbs,exp,errexp))
		outfile.write('\n')
		outfile.close()
	
	
	print("--- %s seconds ---" % (time.time() - start_time))
	print("Finished "+f1array[i])


#end
