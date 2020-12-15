import csv
import random
import os
os.chdir(".")

# Get an empirical distribution from Ch22 and read it.

def ReadEmpDist(Filename):
    EmpDist={30: 0,
        40: 0,
        50: 0,
        60: 0,
        70: 0,
        80: 0,
        90: 0,
        100: 0,
        110: 0,
        120: 0,
        130: 0,
        140: 0,
        150: 0,
        160: 0,
        170: 0}
    with open(Filename, "r") as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            line = row[0].split(" ") 
            line = list(filter(len,line))
            for x,y in EmpDist.items():
                if(int(line[1])-int(x)>=0 and int(line[1])-int(x)<=9 ):
                    EmpDist[x]=int(EmpDist[x])+int(line[0])
    return(EmpDist)

def TSVtoDict(fileTSV):
    TSVDist={20: [],
        30: [],
        40: [],
        50: [],
        60: [],
        70: [],
        80: [],
        90: [],
        100: [],
        110: [],
        120: [],
        130: [],
        140: [],
        150: [],
        160: [],
        170: []}
    with open(fileTSV, "r") as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            for x,y in TSVDist.items():
                if(int(row[1])-int(x)>=0 and int(row[1])-int(x)<=9 ):
                    TSVDist[x].append(row[0])
    return(TSVDist)

## Create a function that return a list of keys from dictionary which has the given value
def getKeysByValue(dictOfElements, valueToFind):
    listOfKeys = list()
    listOfItems = dictOfElements.items()
    for item  in listOfItems:
        if item[1] == valueToFind:
            listOfKeys.append(item[0])
    return  listOfKeys

## Corresponding to our empirical distribution, what is the actual numbers of reads we have to pick in each bin
def GetNumberToPick(TSVDict,EmpDict):
    PickDict={30: 0,
        40: 0,
        50: 0,
        60: 0,
        70: 0,
        80: 0,
        90: 0,
        100: 0,
        110: 0,
        120: 0,
        130: 0,
        140: 0,
        150: 0,
        160: 0,
        170: 0}
    MostBin=getKeysByValue(EmpDict,max(EmpDict.values()))
    MostBin=MostBin[0]
    MostBinValue=max(EmpDict.values())
    AllFromMostBin=(len(TSVDict[MostBin]))
    PickDict[MostBin]=AllFromMostBin
    for x,y in EmpDict.items():
        PickDict[x]=int(int(AllFromMostBin)*int(y)/int(MostBinValue))
    return(PickDict)

def PickReads(TSVDict,PickDict,Mini):
    Keep={30: [],
        40: [],
        50: [],
        60: [],
        70: [],
        80: [],
        90: [],
        100: [],
        110: [],
        120: [],
        130: [],
        140: [],
        150: [],
        160: [],
        170: []}
    
    for x,y in PickDict.items():
        if(int(x)<170): #cause we didnt simulate reads higher than 170 and this is the higher window of the empirical distribution
            if(Mini==True):
                if(y>len(TSVDict[x])):
                    Keep[x]=random.sample(TSVDict[x], k=len(TSVDict[x]))
                else:
                    Keep[x]=random.sample(TSVDict[x], k=y)
            else:
                Keep[x]=random.sample(TSVDict[x], k=y)
    return(Keep)

def CreateIDFile(KeepDict,Name):
    with open(Name, "w+") as f:
        for x in KeepDict.keys():
            for i in KeepDict[x]:
                #f.write(i.split("|")[0]+"\n")
                f.write(i+"\n")
    f.close()




EmpDist=ReadEmpDist("Motala1_read_length_dist_chr22.txt")
for i in os.listdir("."):
    if i.endswith(".len"):
        filename=i
        DictTSV=TSVtoDict(filename)
        ToPick=GetNumberToPick(DictTSV,EmpDist)
        if("mini" in filename):
            ToKeep=PickReads(DictTSV,ToPick,True)
        else:
            ToKeep=PickReads(DictTSV,ToPick,False)
        Name=i.split(".map")[0]
        CreateIDFile(ToKeep, "IDtokeep_"+Name+".txt")
        os.system("bash PicardRun.sh {} {} {}".format(Name, "SubSample"+Name,"IDtokeep_"+Name+".txt"))
