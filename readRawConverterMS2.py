# coding = utf-8

"""
读取rawconverter导出的ms2文件
"""



def addRankRelInt2Spec(specInfoDic):
    mzList = sorted(list(specInfoDic.keys()))
    intsList = [specInfoDic[mz][0] for mz in mzList]
    relIntsList = [round(ele/max(intsList), 5) for ele in intsList]
    sortedIntsList = sorted(intsList, reverse = True)
    rankList = [sortedIntsList.index(ele) + 1 for ele in intsList]
    for mz in mzList:
        idx = mzList.index(mz)
        specInfoDic[mz].extend([relIntsList[idx], rankList[idx]])


def getScanInfoDic(ms2Filepath):
    f = open(ms2Filepath, 'r').readlines()
    scanInfoDic = {}
    i = 0
    while i < len(f):
        if f[i][0] != "S":
            i += 1
        else:
            scan = int(f[i].split("\t")[1])
            mzPre = float(f[i].split("\t")[-1])
            preCharge = f[i + 9].split("\t")[1]
            nce = f[i + 6].split("\t")[2].split(" ")[7].split("@")[1][3:5]

            specInfo = {}
            p = i + 10
            while p < len(f):
                if f[p][0] == "S":
                    break
                else:
                    lineList = f[p].split(" ")
                    mz = float(lineList[0])
                    ints = float(lineList[1])
                    #charge = int(lineList[2])
                    specInfo[mz] = [ints]
                    p += 1
            addRankRelInt2Spec(specInfo)
            scanInfoDic[scan] = [mzPre, preCharge, nce, specInfo]
            i = p
            
    return scanInfoDic

#print(getScanInfoDic("./test.ms2"))


def generate_ion_mass_range(num, mstol):
    deta = num * mstol / 1000000
    return num - deta, num + deta


def isMatched(mz, ms2_dic, mstol):
    mz_list = sorted(list(ms2_dic.keys()))
    lowTgtMZ, upTgtMZ = generate_ion_mass_range(mz, mstol)
    if lowTgtMZ > mz_list[-1]:
        return False
    else:
        i = 0
        while i < len(mz_list):
            if mz_list[i] < lowTgtMZ:
                i += 1
            elif mz_list[i] <= upTgtMZ:
                return True, mz_list[i], ms2_dic[mz_list[i]][2], ms2_dic[mz_list[i]][3]
            else:
                return False


# 计算intensity 大于cutoff的谱峰的匹配度
def exlpainRatio(spec, setMatched, cutoff):
    setTarget = set()   # 将谱峰大于cutoff的mz放在一个set里
    for mz in spec:
        if spec[mz][2] >= cutoff:
            setTarget.add(mz)
    
    overlap = setTarget & setMatched
    numTarge = len(setTarget)
    expRatio = len(overlap) / numTarge
    return numTarge, expRatio


#计算topN mz的匹配程度
def exlpainRatioTopN(spec, setMatched, cutoff):
    setTarget = set()
    for mz in spec:
        if spec[mz][3] <= cutoff:
            setTarget.add(mz)
    overlap = setTarget & setMatched
    numTarge = len(setTarget)
    expRatio = len(overlap) / numTarge
    return expRatio