from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import numpy as np
from odbAccess import *

path="E:/Workspace/NL/ABAQUS/shell/"

mdb.ModelFromInputFile(name='topo',inputFileName=path+'topo.inp')
a = mdb.models['topo'].rootAssembly
p1 = mdb.models['topo'].parts['PART-1']
mdb.models['topo'].parts['PART-1'].setValues(space=THREE_D, type=DEFORMABLE_BODY)
mdb.models['topo'].StaticStep(name='Step-1', previous='Initial',
    timePeriod=96000.0, maxNumInc=100, timeIncrementationMethod=FIXED,
    initialInc=960.0, noStop=OFF, nlgeom=ON)
n1 = a.instances['PART-1-1'].nodes
nodes1 = n1[190:191]
region = a.Set(nodes=nodes1, name='Set-F')
mdb.models['topo'].ConcentratedForce(name='Load-1',createStepName='Step-1',region=region,cf2=-60000.0)
nodes1 = n1[80:101]
region = a.Set(nodes=nodes1, name='Set-SPC')
mdb.models['topo'].DisplacementBC(name='BC-1', createStepName='Step-1',region=region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0)
########  read the iter.mac and retrieve the iteration number
import re
fp = open(path+'iter.txt')
s = fp.readline()
iter_str = re.findall('\d+',s)
iter_str
iter=int(iter_str[0])
fp.close()

NumDesVar=1600
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
    mdb.models['topo'].materials['steel_'+str(i)].Elastic(table=((4000.0*x**pen, 0.4), ))
    mdb.models['topo'].HomogeneousShellSection(name='sec_'+str(i), material='steel_'+str(i), thickness=100.0)
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
UR1=myodb.steps['Step-1'].frames[99].fieldOutputs['UR'].getSubset(region=nset)
U2=myodb.steps['Step-1'].frames[100].fieldOutputs['U'].getSubset(region=nset)
UR2=myodb.steps['Step-1'].frames[100].fieldOutputs['UR'].getSubset(region=nset)
DU=U2-U1
DUR=UR2-UR1
DDU=DU.values
DDUR=DUR.values
UC1=U2.values
UCR1=UR2.values

UC1=U2.values
UCR1=UR2.values
RF1=myodb.steps['Step-1'].frames[100].fieldOutputs['NFORC1'].getSubset(region=nset).values
RF2=myodb.steps['Step-1'].frames[100].fieldOutputs['NFORC2'].getSubset(region=nset).values
RF3=myodb.steps['Step-1'].frames[100].fieldOutputs['NFORC3'].getSubset(region=nset).values
RF4=myodb.steps['Step-1'].frames[100].fieldOutputs['NFORC4'].getSubset(region=nset).values
RF5=myodb.steps['Step-1'].frames[100].fieldOutputs['NFORC5'].getSubset(region=nset).values
RF6=myodb.steps['Step-1'].frames[100].fieldOutputs['NFORC6'].getSubset(region=nset).values
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
    UC1[nnum1-1].data[0]*RF1[f_f-4].data+UC1[nnum1-1].data[1]*RF2[f_f-4].data+UC1[nnum1-1].data[2]*RF3[f_f-4].data+\
    UCR1[nnum1-1].data[0]*RF4[f_f-4].data+UCR1[nnum1-1].data[1]*RF5[f_f-4].data+UCR1[nnum1-1].data[2]*RF6[f_f-4].data+\
    UC1[nnum2-1].data[0]*RF1[f_f-3].data+UC1[nnum2-1].data[1]*RF2[f_f-3].data+UC1[nnum2-1].data[2]*RF3[f_f-3].data+\
    UCR1[nnum2-1].data[0]*RF4[f_f-3].data+UCR1[nnum2-1].data[1]*RF5[f_f-3].data+UCR1[nnum2-1].data[2]*RF6[f_f-3].data+\
    UC1[nnum3-1].data[0]*RF1[f_f-2].data+UC1[nnum3-1].data[1]*RF2[f_f-2].data+UC1[nnum3-1].data[2]*RF3[f_f-2].data+\
    UCR1[nnum3-1].data[0]*RF4[f_f-2].data+UCR1[nnum3-1].data[1]*RF5[f_f-2].data+UCR1[nnum3-1].data[2]*RF6[f_f-2].data+\
    UC1[nnum4-1].data[0]*RF1[f_f-1].data+UC1[nnum4-1].data[1]*RF2[f_f-1].data+UC1[nnum4-1].data[2]*RF3[f_f-1].data+\
    UCR1[nnum4-1].data[0]*RF4[f_f-1].data+UCR1[nnum4-1].data[1]*RF5[f_f-1].data+UCR1[nnum4-1].data[2]*RF6[f_f-1].data
    #sen[flag]=c[flag]*pen/i
    c[flag]=-c[flag]
    sen[flag]=\
    DDU[nnum1-1].data[0]*RF1[f_f-4].data+DDU[nnum1-1].data[1]*RF2[f_f-4].data+DDU[nnum1-1].data[2]*RF3[f_f-4].data+\
    DDUR[nnum1-1].data[0]*RF4[f_f-4].data+DDUR[nnum1-1].data[1]*RF5[f_f-4].data+DDUR[nnum1-1].data[2]*RF6[f_f-4].data+\
    DDU[nnum2-1].data[0]*RF1[f_f-3].data+DDU[nnum2-1].data[1]*RF2[f_f-3].data+DDU[nnum2-1].data[2]*RF3[f_f-3].data+\
    DDUR[nnum2-1].data[0]*RF4[f_f-3].data+DDUR[nnum2-1].data[1]*RF5[f_f-3].data+DDUR[nnum2-1].data[2]*RF6[f_f-3].data+\
    DDU[nnum3-1].data[0]*RF1[f_f-2].data+DDU[nnum3-1].data[1]*RF2[f_f-2].data+DDU[nnum3-1].data[2]*RF3[f_f-2].data+\
    DDUR[nnum3-1].data[0]*RF4[f_f-2].data+DDUR[nnum3-1].data[1]*RF5[f_f-2].data+DDUR[nnum3-1].data[2]*RF6[f_f-2].data+\
    DDU[nnum4-1].data[0]*RF1[f_f-1].data+DDU[nnum4-1].data[1]*RF2[f_f-1].data+DDU[nnum4-1].data[2]*RF3[f_f-1].data+\
    DDUR[nnum4-1].data[0]*RF4[f_f-1].data+DDUR[nnum4-1].data[1]*RF5[f_f-1].data+DDUR[nnum4-1].data[2]*RF6[f_f-1].data
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
for i in e :
    nodes=i.getNodes()
    nc1=np.array(nodes[0].coordinates)
    nodes=i.getNodes()
    nc2=np.array(nodes[1].coordinates)
    nodes=i.getNodes()
    nc3=np.array(nodes[2].coordinates)
    nodes=i.getNodes()
    nc4=np.array(nodes[3].coordinates)
    nc=0.25*(nc1+nc2+nc3+nc4)
    loccentre[flag]=nc
    flag=flag+1
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

#cpFile=open('force1.txt','w')
#for i in RF1 :
#    cpFile.write('%d %7.4f\n' % (i.nodeLabel,i.data))
#else:
#    cpFile.close()

#cpFile=open('force2.txt','w')
#for i in RF2 :
#    cpFile.write('%d %7.4f\n' % (i.nodeLabel,i.data))
#else:
#    cpFile.close()

#cpFile=open('disp1.txt','w')
#for i in RS1 :
#    cpFile.write('%d %7.4f %7.4f\n' % (i.nodeLabel,i.data[0],i.data[1]))
#else:
#    cpFile.close()

#cpFile=open('disp2.txt','w')
#for i in RS2 :
#    cpFile.write('%d %7.4f %7.4f\n' % (i.nodeLabel,i.data[0],i.data[1]))
#else:
#    cpFile.close()






##define the variables for output
#stepName = "Step-1"
#instanceNo = 0
#frameNo = 0
#fieldOutputName = 'RF'
#rfField = odb.steps[stepName].frames[frameNo].fieldOutputs[fieldOutputName]
#rfFieldSub = rfField.getSubset(position = ELEMENT_NODAL)
#rfValues = rfFieldSub.values
#FX = RForce.getSubset(region=regS1).values[0].data[0]
#change the .values[0] to .values[1] or higher to change the node
#RForce.getSubset(region=regS1).values[0].data[0]


#e= mdb.models['topo'].parts['PART-1'].elements[0]
#n=e.getNodes() # n is a tuple of meshnode object
#sqn=node[n[0].label-1:n[0].label-1]
#en=p.Set(nodes=sqn, name='test')
#p.SetFromNodeLabels(name='test',nodeLabels=n)
#nset=myodb.rootAssembly.nodeSets[' ALL NODES']
#RS1=myodb.steps['Step-1'].frames[1].fieldOutputs['U'].getSubset(region=nset).values
#RS2=myodb.steps['Step-2'].frames[1].fieldOutputs['U'].getSubset(region=nset).values
#RF1=myodb.steps['Step-1'].frames[1].fieldOutputs['NFORC1'].getSubset(region=nset).values
#RF2=myodb.steps['Step-1'].frames[1].fieldOutputs['NFORC2'].getSubset(region=nset).values
################session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)
#RS1.__methods__
#newField = odb.steps[stepName].frames[frameNo].FieldOutput(name='microStress',
#        description='stress in microscale', validInvariants=(MISES, MAX_PRINCIPAL, MID_PRINCIPAL,
#         MIN_PRINCIPAL, TRESCA, PRESS, INV3), type=TENSOR_3D_FULL)
#    newField.addData(field=microStressField)
#    odb.save()

#    n=e[i-1].getNodes() # n is a tuple of meshnode object
#    a=n[0].label
#    b=n[1].label
#    c=n[2].label
#    d=n[3].label
#    nodes=node[a:a+1]+node[b:b+1]+node[c:c+1]+node[d:d+1]
#    p.Set(nodes=nodes, name='S'+str(i))
#e= mdb.models['topo'].parts['PART-1'].elements
#for i in e :
#    n=i.getNodes()
#    U1=myodb.steps['Step-1'].frames[1].fieldOutputs['U'].getSubset(region=n).values
#    U2=myodb.steps['Step-2'].frames[1].fieldOutputs['U'].getSubset(region=n).values
#    du=U1-U2
#    RF1=myodb.steps['Step-1'].frames[1].fieldOutputs['NFORC1'].getSubset(region=i).values
#    RF2=myodb.steps['Step-1'].frames[1].fieldOutputs['NFORC2'].getSubset(region=i).values
