from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import numpy as np
from odbAccess import *
path="D:/abaqustemp/abaqus2plane/50/"
mdb.ModelFromInputFile(name='topo',inputFileName=path+'topo.inp')
a = mdb.models['topo'].rootAssembly
p1 = mdb.models['topo'].parts['PART-1']
mdb.models['topo'].parts['PART-1'].setValues(space=TWO_D_PLANAR, type=DEFORMABLE_BODY)
mdb.models['topo'].steps['Step-1'].setValues(timePeriod=50000.0, maxNumInc=100,\
 timeIncrementationMethod=FIXED,initialInc=500.0, noStop=OFF, nlgeom=ON)
#mdb.models['topo'].StaticStep[0].(previous='Initial',
#    timePeriod=50000.0, maxNumInc=100, timeIncrementationMethod=FIXED,
#    initialInc=500.0, noStop=OFF, nlgeom=ON)
#n1 = a.instances['PART-1-1'].nodes
#nodes1 = n1[100:101]
#region = a.Set(nodes=nodes1, name='Set-F')
#mdb.models['topo'].ConcentratedForce(name='Load-1',createStepName='Step-1',region=region,cf1=0,cf2=50000.0)
#nodes1 = n1[0:21]+n1[1708:1729]
#region = a.Set(nodes=nodes1, name='Set-SPC')
#mdb.models['topo'].DisplacementBC(name='BC-1', createStepName='Step-1',region=region, u1=0.0, u2=0.0)
########  read the iter.mac and retrieve the iteration number
import re
fp = open(path+'iter.txt')
s = fp.readline()
iter_str = re.findall('\d+',s)
iter_str
iter=int(iter_str[0])
fp.close()

NumDesVar=2880
i=1
pen=3
xi=np.ones(NumDesVar*1,dtype=float)
file_object = open(path+'xi_i'+str(iter)+'.mac')
i=0
import re
import string
for line in file_object:
    a=line
    xi[i]=string.atof(a)
    i=i+1
i=1
for x in xi:
    p = mdb.models['topo'].parts['PART-1']
    e= mdb.models['topo'].parts['PART-1'].elements
    node=mdb.models['topo'].parts['PART-1'].nodes
    mdb.models['topo'].Material(name='steel_'+str(i))
    mdb.models['topo'].materials['steel_'+str(i)].Elastic(table=((210000.0*x**pen, 0.3), ))
    mdb.models['topo'].HomogeneousSolidSection(name='sec_'+str(i), material='steel_'+str(i), thickness=100.0)
    elements = e[i-1:i]
    #Set up the node Set
    region = p.Set(elements=elements, name='Set-TEM'+str(i))
    p.SectionAssignment(region=region, sectionName='sec_'+str(i), offset=0.0, offsetType=MIDDLE_SURFACE,offsetField='',thicknessAssignment=FROM_SECTION)
    i=i+1
#the requset for output
mdb.models['topo'].fieldOutputRequests['F-Output-1'].setValues(variables=('S','MISES','E','U','NFORC','ENER','ELEN'))
###########################################
##run the analysis
###########################################
mdb.Job(name='Job-'+str(iter), model='topo', type=ANALYSIS,resultsFormat=ODB)
mdb.jobs['Job-'+str(iter)].submit(consistencyChecking=OFF)
###########################################
##Read the rersult and run the sensitivity analysis
###########################################
from odbAccess import *
myodb=openOdb(path+'Job-'+str(iter)+'.odb')
nset=myodb.rootAssembly.nodeSets[' ALL NODES']
U1=myodb.steps['Step-1'].frames[99].fieldOutputs['U'].getSubset(region=nset)
U2=myodb.steps['Step-1'].frames[100].fieldOutputs['U'].getSubset(region=nset)
DU=U2-U1
DDU=DU.values
U=U2.values
RF1=myodb.steps['Step-1'].frames[100].fieldOutputs['NFORC1'].getSubset(region=nset).values
RF2=myodb.steps['Step-1'].frames[100].fieldOutputs['NFORC2'].getSubset(region=nset).values
#RF1 RF2 FieldValueArray
#check the node set
#myodb.rootAssembly.nodeSets.keys()
#sensitivity analysis
c=np.ones(NumDesVar*1,dtype=float)
sen=np.ones(NumDesVar*1,dtype=float)
een=np.ones(NumDesVar*1,dtype=float)
flag=0
comp=0
for i in xi :
    f_f=(flag+1)*4
    nnum1=RF1[f_f-4].nodeLabel
    nnum2=RF1[f_f-3].nodeLabel
    nnum3=RF1[f_f-2].nodeLabel
    nnum4=RF1[f_f-1].nodeLabel
    c[flag]=\
    U[nnum1-1].data[0]*RF1[f_f-4].data+U[nnum1-1].data[1]*RF2[f_f-4].data+\
    U[nnum2-1].data[0]*RF1[f_f-3].data+U[nnum2-1].data[1]*RF2[f_f-3].data+\
    U[nnum3-1].data[0]*RF1[f_f-2].data+U[nnum3-1].data[1]*RF2[f_f-2].data+\
    U[nnum4-1].data[0]*RF1[f_f-1].data+U[nnum4-1].data[1]*RF2[f_f-1].data
    #sen[flag]=c[flag]*pen/i
    c[flag]=-c[flag]
    sen[flag]=DDU[nnum1-1].data[0]*RF1[f_f-4].data+DDU[nnum1-1].data[1]*RF2[f_f-4].data+\
    DDU[nnum2-1].data[0]*RF1[f_f-3].data+DDU[nnum2-1].data[1]*RF2[f_f-3].data+\
    DDU[nnum3-1].data[0]*RF1[f_f-2].data+DDU [nnum3-1].data[1]*RF2[f_f-2].data+\
    DDU[nnum4-1].data[0]*RF1[f_f-1].data+DDU[nnum4-1].data[1]*RF2[f_f-1].data
    sen[flag]=sen[flag]*pen*100/i
    flag=flag+1
#calculate the valuses of responses
volume=np.sum(xi)
comp=np.sum(c)
# Sensitivity filtering:
from numpy import linalg as LA
#filter radius
rmin=25
p = mdb.models['topo'].parts['PART-1']
e= mdb.models['topo'].parts['PART-1'].elements
n = mdb.models['topo'].parts['PART-1'].nodes
loccentre = [[0 for col in range(3)] for row in range(NumDesVar)]
flag=0
#store the coordinates of elements
for j in range(0,NumDesVar) :
    i=e[j]
    nodes=i.getNodes()
    nc1=np.array(nodes[0].coordinates)
    nodes=i.getNodes()
    nc2=np.array(nodes[1].coordinates)
    nodes=i.getNodes()
    nc3=np.array(nodes[2].coordinates)
    nodes=i.getNodes()
    nc4=np.array(nodes[3].coordinates)
    nc=0.25*(nc1+nc2+nc3+nc4)
    loccentre[j]=nc
#filtering
filter_e=1
filter_f=1
numerator=0
denominator=0
for i in loccentre :
    for j in loccentre:
        dis=LA.norm(j-i)
        if (dis<=rmin) :
            numerator=numerator+(rmin-dis)*xi[filter_f-1]*sen[filter_f-1]
            denominator=denominator+rmin-dis
        filter_f=filter_f+1
    een[filter_e-1]=numerator/denominator/xi[filter_e-1]
    filter_f=1
    filter_e=filter_e+1
    numerator=0
    denominator=0
sen=np.copy(een)
#print
session.viewports['Viewport: 1'].setValues(displayedObject=myodb)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.printToFile(fileName=path+'topo_i'+str(iter), format=PNG, canvasObjects=(
    session.viewports['Viewport: 1'], ))
###########################################
##Formatted output for cap file (BOSS version)
###########################################
cpFile=open(path+'topo_i'+str(iter)+'.cap','w')
cpFile.write(':INFO\n')
cpFile.write('   Topology Optimization Developped\n')
cpFile.write(':VAR_TYPE\n')
cpFile.write('   UNKNOWN\n')
cpFile.write(':VARIABLE\n')
flag=1
for i in c :
    cpFile.write('  ph'+'%d\n' % flag)
    flag=flag+1
cpFile.write(':FUNCTION\n')
cpFile.write('   compl\n')
cpFile.write('   compl\n')
cpFile.write('   1\n')
cpFile.write('%19.5E\n' % comp)
for i in sen :
    cpFile.write('%19.5E\n' % (float(i)))
cpFile.write(':FUNCTION\n')
cpFile.write('   Volume\n')
cpFile.write('   Volume\n')
cpFile.write('   1\n')
cpFile.write('%19.5F\n' % volume)
for i in range(0,NumDesVar) :
    cpFile.write('%19.5F\n' % 1.0)
cpFile.close()
###########################################
##Plotting
###########################################
