import numpy as np
from ReadCFXExport import ReadCFXExport
from numpy.linalg import norm

#http://stackoverflow.com/a/31427697/1542146
from numpy.lib.recfunctions import append_fields

# Input Params #
l = 0.314617	#Tube length
r = 0.0062611	#Tube radius

Prt = 0.9	#Turbulent Prandtl Number

fast = True

axisName = 'Z'

fNames = [r'ModelA1_SL_ke0 of 6.csv',
		r'ModelA1_SL_kw0 of 6.csv',
		r'ModelA1_SL_kwSST0 of 6.csv',
		r'ModelA1_SL_SASSST0 of 6.csv']

labels = ['$k$-$\epsilon$', '$k$-$\omega$', '$k$-$\omega$ SST', 'SAS SST']

#TrnStr = 'Trnavg'
TrnStr = ''

axisDict = {'X':0, 'Z':2}
aInd = axisDict[axisName]
cInd = abs(2-aInd)

PS = np.array(['X','Y','Z'])
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

totShearT = np.zeros(4)
totShearA = np.zeros(4)
totShearX = np.zeros(4)
totTKEtfr = np.zeros(4)
totTEtfr = np.zeros(4)
totCond = np.zeros(4)
totE = np.zeros(4)
totMasstfr = np.zeros(4)
totMassEtfr = np.zeros(4)

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

cone_angle = np.radians(30)
cone_slope = np.tan(cone_angle)
nose_drad = 2.5/1000
rear_drad = 7.407477/1000
cone_dist = 1/1000
cone_offset_point = l + cone_dist
cone_rear_z = cone_offset_point + (rear_drad - r)/cone_slope
cone_tip_z = cone_offset_point + (nose_drad - r)/cone_slope
coldTubeRad = 2.5/1000
coldTubeLen = 10.0/1000

def get_cone_z(r):
	return cone_tip_z + (r - nose_drad)/cone_slope

def get_cone_r(z):
	return nose_drad + (z - cone_tip_z)*cone_slope

tubePts = [(cone_tip_z, 0.),(cone_tip_z, nose_drad),(cone_offset_point, get_cone_r(cone_offset_point)),(l, r),(0., r),(0., coldTubeRad),(-coldTubeLen, coldTubeRad),(-coldTubeLen, 0),]
tubePts.append(tubePts[0])
tubePts = np.array(tubePts).T

fig = plt.figure(num=0, figsize=(12,3))
plt.plot(tubePts[0],tubePts[1],c='k',lw=2)

XC = np.linspace(0,l,4000)
YC = np.linspace(0,r,200)

ax = plt.gca()
#ax.set_aspect('equal', 'datalim')
ax.set_aspect(4, 'datalim')

linestyles = ['dotted','dashed','dashdot','solid']

i=0
label = None
for fName, label, i in zip(fNames, labels, np.arange(len(fNames))):
	
	Fields = ReadCFXExport(fName)

	#Radial co-ordinates
	RCrd = np.sqrt(Fields[PS[0]]**2 + Fields[PS[1]]**2)

	#Position Vector
	Pos = np.vstack((Fields[PS[0]],Fields[PS[1]],Fields[PS[2]])).T

	PosCyl = np.vstack((RCrd,np.arctan2(Fields[PS[0]],Fields[PS[1]]),Fields[PS[2]])).T
	axialVec = np.array([0,0,1])

	#Arc length of curve at each point
	AL = np.zeros_like(Fields[Tstr])
	VecNorm = np.zeros_like(Fields[Norm])
	for j in range(1,Fields[Tstr].shape[0]):
		disp = PosCyl[j]-PosCyl[j-1]
		disp[1] = 0
		AL[j] = AL[j-1] + norm(disp)

		rUnitVec = Pos[j].copy()
		rUnitVec[2] = 0
		perpVec = np.cross(axialVec,rUnitVec)
		normVec = np.cross(perpVec,Pos[j]-Pos[j-1])
		VecNorm[j] = -normVec/norm(normVec)

	VecNorm = np.nan_to_num(VecNorm,copy=True)
	AL = AL[-1] - AL
	if np.any(np.isnan(AL)):
		print('Nan in AL array')

	Fields[Norm] = VecNorm
	RFields = Fields.copy()

	nPos = Pos+VecNorm*(0.01/200)

	#Need to compute 4 Types of energy transfer mean and turbulent conduction, and mean and turbulent shear

	#Mean Conduction Energy Transfer
	res = np.einsum('ij,ij->i',Fields[Tgrad],VecNorm)
	CondHM = -np.einsum('i,i,i->i',RCrd,Fields[Tcond],res)
	#print(Fields[Norm],Fields[Tgrad])

	#Turbulent Conduction Energy Transfer
	res = np.einsum('ij,ij->i',Fields[Egrad],VecNorm)
	CondHT = -np.einsum('i,i,i->i',RCrd,Fields[Eddyv]/Prt,res)

	ViscEff = Fields[Eddyv] + Fields[Dvisc]

	shear_T = np.einsum('i,i,ij,ij->i',RCrd,Fields[CircVel],RFields[Norm],RFields[AVGrad])
	shear_A = np.einsum('i,ij,ij->i',Fields[AxialVel],RFields[Norm],RFields[AVelGrad])

	ShearT = -np.einsum('i,i,i->i',RCrd,ViscEff,shear_T)
	ShearA = -np.einsum('i,i,i->i',RCrd,ViscEff,shear_A)

	#Net Energy transfer
	netE = CondHM + CondHT + ShearT + ShearA

	#Compute the total energy transfer
	print(i)
	totShearT[i] = np.trapz(ShearT,AL)*2*np.pi
	totShearA[i] = np.trapz(ShearA,AL)*2*np.pi
	totCond[i] = np.trapz(CondHM+CondHT,AL)*2*np.pi
	totE[i] = np.trapz(netE,AL)*2*np.pi

	#Plot
	plt.figure(num=1)
	clr = f'k'

	plt.subplot(221)
	plt.plot(AL,CondHT, label = label, linestyle=linestyles[i],color=clr)
	plt.subplot(222)
	plt.plot(AL,ShearT,linestyle=linestyles[i],color=clr)
	plt.subplot(223)
	plt.plot(AL,ShearA,linestyle=linestyles[i],color=clr)
	plt.subplot(224)
	plt.plot(AL,netE,label = label, linestyle=linestyles[i],color=clr)
	
	plt.figure(num=0)
	plt.axis('off')
	plt.plot(Fields[PS[aInd]],RCrd,marker=None,linestyle=linestyles[i],color=clr)

	energytfr = np.einsum('i,i,i,ij,ij->i',RCrd,Fields[Dens],Fields[Enthalpy],Fields[Norm],Fields[Ustr])

print('Circumfrential Shear Work Transfer = {} [W]'.format(totShearT))
print('Axial Shear Work Transfer= {} [W]'.format(totShearA))
print('Conduction Heat Transfer = {} [W]'.format(totCond))
print('Net Energy transfer = {} [W]'.format(totE))

fig = plt.figure(num=1)
axislabels = ['Heat Conduction\n Transfer [W]', 'Circumferential Shear\n Work Transfer [W]','Axial Shear\n Work Transfer [W]', 'Net Energy Transfer [W]']
for i in range(4):
	plt.subplot(221+i)
	ax = plt.gca()
	ax.set_xlabel('Streamline Co-ordinate [m]', fontsize=labelsize)
	ax.set_ylabel(axislabels[i], fontsize=labelsize)
	ax.grid(True)
	ax.set_xlim(0,AL[0]*0.45)
	if i == 0:
		#ax.legend(fancybox=False, framealpha=1.0, loc = 'lower right')
		handles, labels = ax.get_legend_handles_labels()

fig.legend(handles, labels, loc='upper center', ncol=2, borderaxespad=0.25)
plt.show()