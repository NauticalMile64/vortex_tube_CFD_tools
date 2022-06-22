import numpy as np
from ReadCFXExport import ReadCFXExport
from numpy.linalg import norm

#http://stackoverflow.com/a/31427697/1542146
from numpy.lib.recfunctions import append_fields

# Input Params #
l = 0.314		#Tube length
r = 0.0062611	#Tube radius

fast = True	#Fast switch

Prt = 0.9	# Turbulent Prandtl Number used by CFX

numSLs = 6	# Number of Streamlines

axisName = 'Z'
basefName = r'ModelA1_SL_SASSST{} of {}.csv'

TrnStr = ''

axisDict = {'X':0, 'Z':2}
aInd = axisDict[axisName]
cInd = abs(2-aInd)

PS = np.array(['X','Y','Z'])
#PSC = np.array(['Radius','Theta','Axial Distance'])
Ustr = f'Velocity{TrnStr}'
RadialVel = f'Velocity{TrnStr}_Radial'
RadialVelGrad = f'Velocity{TrnStr}_RadialGradient'
CircVel = f'Velocity{TrnStr}_Circumferential'
AxialVelGrad = f'Velocity{TrnStr}_AxialGradient'
AxialVel = f'Velocity{TrnStr}_Axial'

Enthalpy = f'Static_Enthalpy{TrnStr}'
Egrad = f'Static_Enthalpy{TrnStr}Gradient'
Eddyv = f'Eddy_Viscosity{TrnStr}'
Dens = f'Density{TrnStr}'
TKE = f'Turbulence_Kinetic_Energy{TrnStr}'
TKEgrad = f'Turbulence_Kinetic_EnergyGradient{TrnStr}'
TEgrad = f'Total_Enthalpy{TrnStr}Gradient'

AVGrad = 'AngVelGradient'

RVel = RadialVel
AVel = AxialVel
RVelGrad = RadialVelGrad
AVelGrad = AxialVelGrad

Tstr = 'Temperature'
Tgrad =  'TemperatureGradient'
Tcond = 'Thermal_Conductivity'
aLen = 'Arc Length'
Pos = 'Position'
Norm = 'Normal'
Nnum = 'Node_Number'
VelDiv = 'VelocityDivergence'
Dvisc ='Dynamic_Viscosity'

totShearT = np.zeros(numSLs)
totShearA = np.zeros(numSLs)
totShearX = np.zeros(numSLs)
totTKEtfr = np.zeros(numSLs)
totTEtfr = np.zeros(numSLs)
totCond = np.zeros(numSLs)
totE = np.zeros(numSLs)

from matplotlib import pyplot as plt
import matplotlib.cm as cm

labelsize = 16
plt.rc('legend',fontsize=labelsize-2)

if not fast:
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

fig = plt.figure(num=1)
plt.xlabel('Streamline Co-ordinate [m]', fontsize=labelsize)
plt.ylabel('Energy Transfer\nTo Hot Stream [W/m]', fontsize=labelsize)
plt.tick_params(axis='both', labelsize=labelsize-2)

tubePts = [(l, 0.),(l, r),(0., r),(0., 0.),]
tubePts.append(tubePts[0])
tubePts = np.array(tubePts).T

fig = plt.figure(num=0, figsize=(12,3))
plt.plot(tubePts[0],tubePts[1],c='k',lw=2)

XC = np.linspace(0,l,4000)
YC = np.linspace(0,r,200)

ax = plt.gca()
ax.set_aspect('equal', 'datalim')

from scipy.interpolate import interp1d

axialVec = np.array([0,0,1])

label = None
for i in range(0,numSLs):
	
	fName = basefName.format(i,numSLs)
	#fName = basefName
	Fields = ReadCFXExport(fName)
	
	#Skip streamlines which do not extend to the inlet
	if np.min(Fields[PS[2]]) > 0.0015:
		#skips += 1
		continue
	
	#Radial co-ordinates
	RCrd = np.sqrt(Fields[PS[0]]**2 + Fields[PS[1]]**2)
	#print(RCrd)
	
	#Position Vector
	Pos = np.vstack((Fields[PS[0]],Fields[PS[1]],Fields[PS[2]])).T
	PosCyl = np.vstack((RCrd,np.arctan2(Fields[PS[1]],Fields[PS[0]]),Fields[PS[2]])).T
	
	#Arc length of curve at each point
	AL = np.zeros_like(Fields[Tstr])
	VecNorm = np.zeros_like(Fields[Norm])
	VecNormCyl = np.zeros_like(Fields[Norm])
	totLen = 0
	thetaVecCyl = np.zeros(3)
	
	for j in range(Fields[Tstr].shape[0]-1):
		disp = PosCyl[j+1]-PosCyl[j]
		disp[1] = 0
		AL[j] = totLen + norm(disp)
		totLen = AL[j]
		
		rUnitVec = Pos[j].copy()
		rUnitVec[2] = 0
		thetaVec = np.cross(axialVec,rUnitVec)
		
		ct = np.cos(PosCyl[j,1])
		st = np.sin(PosCyl[j,1])
		
		thetaVecCyl += np.array([ct*thetaVec[0]+st*thetaVec[1],ct*thetaVec[1]-st*thetaVec[0],thetaVec[2]])
		
		cartDisp = np.array([ct*disp[0],st*disp[0],disp[2]])
		normVec = np.cross(thetaVec,cartDisp)
		
		normVecCyl = np.array([ct*normVec[0]+st*normVec[1],0,normVec[2]])
		normVec = np.array([ct*normVecCyl[0]-st*normVecCyl[1],st*normVecCyl[0]+ct*normVecCyl[1],normVecCyl[2]])
		
		VecNorm[j] = -normVec/norm(normVec)
		VecNormCyl[j] = -normVecCyl/norm(normVecCyl)
	
	VecNorm = np.nan_to_num(VecNorm,copy=True)
	VecNorm[-1] = VecNorm[-2]
	VecNormCyl = np.nan_to_num(VecNormCyl,copy=True)
	AL[-1] = AL[-2] + (AL[-2] - AL[-3])
	AL = AL[-1] - AL
	
	if np.any(np.isnan(AL)):
		print('Nan in AL array')
	
	thetaVecCyl /= Fields[Tstr].size
	
	Fields[Norm] = VecNorm
	
	nPos = Pos+VecNorm*(0.01/200)
	nPosR = np.sqrt(nPos[:,0]**2 + nPos[:,1]**2)
	
	#Need to compute 4 Types of energy transfer mean and turbulent conduction, and mean and turbulent shear
	
	#Mean Conduction Energy Transfer
	res = np.einsum('ij,ij->i',Fields[Tgrad],VecNorm)
	CondHM = -np.einsum('i,i,i->i',RCrd,Fields[Tcond],res)
	
	#Turbulent Conduction Energy Transfer
	res = np.einsum('ij,ij->i',Fields[Egrad],VecNorm)
	CondHT = -np.einsum('i,i,i->i',RCrd,Fields[Eddyv]/Prt,res)
	
	ViscEff = Fields[Eddyv] + Fields[Dvisc]
	
	shear_T = np.einsum('i,i,ij,ij->i',RCrd,Fields[CircVel],Fields[Norm],Fields[AVGrad])
	shear_A = np.einsum('i,ij,ij->i',Fields[AxialVel],Fields[Norm],Fields[AVelGrad])
	
	ShearT = -np.einsum('i,i,i->i',RCrd,ViscEff,shear_T)
	ShearA = -np.einsum('i,i,i->i',RCrd,ViscEff,shear_A)
	
	#Net Energy transfer
	netE = CondHM + CondHT + ShearT + ShearA
	
	#Compute the total energy transfer
	print(i)
	# print(f'average theta vector = {thetaVecCyl}')
	# print(f'average cylindrical normal = {np.average(VecNormCyl,axis=0)}')
	totShearT[i] = np.trapz(ShearT,AL)*2*np.pi
	totShearA[i] = np.trapz(ShearA,AL)*2*np.pi
	totCond[i] = np.trapz(CondHM+CondHT,AL)*2*np.pi
	totE[i] = np.trapz(netE,AL)*2*np.pi

	#Plot
	plt.figure(num=1)
	clr = f'C{i}'

	plt.plot(AL,CondHT,linestyle='dotted',color=clr)
	plt.plot(AL,ShearT,linestyle='dashed',color=clr)
	plt.plot(AL,ShearA,linestyle='dashdot',color=clr)
	plt.plot(AL,netE,label = label,linestyle='solid',color=clr)

	plt.figure(num=0)
	plt.plot(Fields[PS[aInd]],RCrd,marker=None,color=clr)
	plt.plot(nPos[:,2],nPosR,marker=None,color=clr,ls='dashed')
		
	plt.figure(num=5)
	plt.plot(AL,Fields[Tstr])
	
	energytfr = np.einsum('i,i,i,ij,ij->i',RCrd,Fields[Dens],Fields[Enthalpy],Fields[Norm],Fields[Ustr])
	plt.figure(num=1)

removalInds = np.where(totE == 0)
totShearT = np.delete(totShearT,removalInds)
totShearA = np.delete(totShearA,removalInds)
totCond = np.delete(totCond,removalInds)
totE = np.delete(totE,removalInds)

print('Circumfrential Shear Work Transfer = {} [W] Average = {} [W]'.format(totShearT,np.average(totShearT)))
print('Axial Shear Work Transfer= {} [W] Average = {} [W]'.format(totShearA,np.average(totShearA)))
print('Conduction Heat Transfer = {} [W] Average = {} [W]'.format(totCond,np.average(totCond)))
print('Net Energy transfer = {} [W] Average = {} [W]'.format(totE,np.average(totE)))

plt.figure(num=1)
ax = plt.gca()
ax.grid(True)
ax.set_xlim(0,AL[0])
plt.legend(fancybox=False, framealpha=1.0, loc = 'upper right')
plt.show()