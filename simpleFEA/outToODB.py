import re
import sys
from odbAccess import *
from odbMaterial import *
from odbSection import *
from abaqusConstants import *


if len(sys.argv)<1: 
        print("Should be invoked with the fileName")
        sys.exit()

modelName =sys.argv[-1]
outFileName = modelName+"_Out.txt"
f=open(outFileName,'r')
allTxt = f.read()
#print allTxt


#Get Nodes
print "\n\nCreating Nodes\n"
s=r"\*NODE START\n(.*)\*NODE END"
reNode = re.compile(s,re.I|re.DOTALL|re.M)
reSrc = re.search(reNode,allTxt)

nodeDataStr = reSrc.group(1)
lnodeData = nodeDataStr.strip().split("\n")
nodeList =[]

cnt=0
for item in lnodeData:
    nodeData=[]
    cnt+=1
    nodeCoors = item.strip().split(",")
    nodeData.append(cnt)
    nodeData.append(float(nodeCoors[0].strip()))
    nodeData.append(float(nodeCoors[1].strip()))
    nodeList.append(nodeData)

#Get Elements
print "Creating Elements\n"
s=r"\*ELEMENT START\n(.*)\*ELEMENT END"
reElem = re.compile(s,re.I|re.DOTALL|re.M)
reSrc = re.search(reElem,allTxt)

elemDataStr = reSrc.group(1)
lelemData = elemDataStr.strip()
lelemData = elemDataStr.strip()[:-1]
lelemData = lelemData.split("\n")

elemList =[]



for item in lelemData:
    elemData=[]    
    elemConn = item.strip().split(",")    
    
    for nodeId in elemConn:
        
        try:
            n = int(nodeId.strip())
        except:
            continue
        elemData.append(n)

    elemList.append(elemData)


s=r"\*INTEGRATION POINTS =(.*?)\n"
reRes = re.compile(s,re.I|re.DOTALL|re.M)
reSrc = re.search(reRes,allTxt)
numIntegrationPoints = int(reSrc.group(1).strip())
#print numIntegrationPoints


#Result
print "Writing Results\n"
s=r"\*RESULT START\n(.*)\*RESULT END"
reRes = re.compile(s,re.I|re.DOTALL|re.M)
reSrc = re.search(reRes,allTxt)

resText = reSrc.group(1)

s=r"\*FRAME START\n(.*?)\*FRAME END"
reFrame= re.compile(s,re.I|re.DOTALL|re.M)
reSrc = re.findall(reFrame,resText)

timeData={}
UData={}
SData={}
frameId=0

for txt in reSrc:
    
    frameId+=1    
    print "\tWriting Results for - "+str(frameId)

    s=r"\*TIME =(.*?)\n"
    reRes = re.compile(s,re.I|re.DOTALL|re.M)
    reSrc = re.search(reRes,txt)
    timeData[frameId]=float(reSrc.group(1).strip())


    s=r"\*U START\n(.*)\*U END"
    reRes = re.compile(s,re.I|re.DOTALL|re.M)
    reSrc = re.search(reRes,txt)

    uString = reSrc.group(1)
    uString = uString.strip()[:-1]
    uDataVal = []

    for item in uString.split(","):
        uDataVal.append(float(item))

    UData[frameId] = uDataVal


    s=r"\*S START\n(.*)\*S END"
    reRes = re.compile(s,re.I|re.DOTALL|re.M)
    reSrc = re.search(reRes,txt)

    sString = reSrc.group(1)
    sString = sString.strip()[:-1]
    sDataVal = []

    for item in sString.split(","):
        sDataVal.append(float(item))
    SData[frameId] = sDataVal

numFrames = frameId

odb = Odb(name=modelName,path=modelName+'.odb')
part1 = odb.Part(name='part-1', embeddedSpace=TWO_D_PLANAR,
        type=DEFORMABLE_BODY)
part1.addNodes(nodeData=nodeList, nodeSetName='nset-1')

part1.addElements(elementData=elemList, type='CPE3',
        elementSetName='eset-1')

instance1 = odb.rootAssembly.Instance(name='part-1-1',
        object=part1)

nodeLabelData = []
for item in nodeList:
    nodeLabelData.append(item[0])

elementLabelData = []
for item in elemList:
    elementLabelData.append(item[0])
    


step1 = odb.Step(name='step-1',
        description='first analysis step',
        domain=TIME, timePeriod=1.0)

for i in range(numFrames):

    analysisTime= timeData[i+1]
    frame1 = step1.Frame(incrementNumber=1,
            frameValue=analysisTime,
            description=\
                'results frame for time '+str(analysisTime))

     #  Write nodal displacements.       
    uField = frame1.FieldOutput(name='U',
        description='Displacements', type=VECTOR)
    
    allDispData=UData[i+1]
    
    dispData= [allDispData[x:x+2] for x in xrange(0, len(allDispData), 2)]

    #print "\n"
    #print dispData,nodeLabelData

    uField.addData(position=NODAL, instance=instance1,
        labels=nodeLabelData,
        data=dispData)



    #Write S
    allSData = SData[i+1]
    stressData = [allSData[x:x+4] for x in xrange(0, len(allSData), 4)]
    sField = frame1.FieldOutput(name='S',
        description='Stress', type=TENSOR_3D_PLANAR,
        componentLabels=('S11', 'S22', 'S33','S12'),
        validInvariants=(MISES,))
    sField.addData(position=INTEGRATION_POINT,
       instance=instance1,labels=elementLabelData, data=stressData)

odb.save()
odb.close()

print "done"