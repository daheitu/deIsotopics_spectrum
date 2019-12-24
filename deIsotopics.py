import os, copy

mstol = 20
ms2path = "./DSSO_CID_MS2_HCD_IT_MS3_SinPair_re23_T1.ms2"
#print(scan_info_dic)


def generate_mass_range(num, delta_ppm):
    delta = num * delta_ppm / 1000000
    return num - delta, num + delta


def compareTwoNum(realNum, TheroNum, tolPPM):
    if abs(TheroNum- realNum)/TheroNum * 1000000 < tolPPM:
        return True
    else:
        return False


def findOneMZ(mz, ms2MzList, begin_idx, tol_ppm):
    low_ms, up_ms = generate_mass_range(mz, tol_ppm)
    if mz >  ms2MzList[begin_idx]:
        k = begin_idx + 1
        
        if ms2MzList[-1] < low_ms:
            return False, -1
        else:
            while k < len(ms2MzList):
                if float(ms2MzList[k]) > up_ms:
                    return False, -1
                elif float(ms2MzList[k]) >= low_ms:
                    return True, ms2MzList[k]
                else:
                    k += 1
            return False, -1
    else:
        k = begin_idx - 1
        if ms2MzList[0] > up_ms:
            return False, -1
        else:
            while k > 0:
                if float(ms2MzList[k]) < low_ms:
                    return False, -1
                elif float(ms2MzList[k]) <= up_ms:
                    return True, ms2MzList[k]
                else:
                    k -= 1
            return False, -1


def addRankRelInt2Spec(specInfoDic):
    mzList = sorted(list(specInfoDic.keys()))
    intsList = [specInfoDic[mz][0] for mz in mzList]
    relIntsList = [round(ele/max(intsList), 5) for ele in intsList]
    sortedIntsList = sorted(intsList, reverse = True)
    rankList = [sortedIntsList.index(ele) + 1 for ele in intsList]
    for mz in mzList:
        idx = mzList.index(mz)
        specInfoDic[mz].extend([relIntsList[idx], rankList[idx]])


def forword_looking(top1_mz, interval, ms2MzList, startPos_to_find, iso_cluster):
    i = startPos_to_find
    m = 2
    while m < 8:
        matchMoreBool, matchMoreMZ = findOneMZ(top1_mz + interval*m, ms2MzList, i, 20)
        if matchMoreBool:
            iso_cluster.append(matchMoreMZ)
            m += 1
        else:
            break


def back_looking(top1_mz, interval, ms2MzList, startPos_to_find, iso_cluster, ms2_dic):
    for n in range(1,6):
        matchLessBool, matchLessMZ = findOneMZ(top1_mz - interval*n, ms2MzList, startPos_to_find, 20)
        if not matchLessBool:
            break
        else:
            if ms2_dic[matchLessMZ][0]/ms2_dic[iso_cluster[0]][0] > 0.3:
                iso_cluster.insert(0, matchLessMZ)
            else:
                break


def deter_min_isotpic(ms2_dic, undeter_iso_dic, chged_iso_dic, chged_2peaks_dic):
    chrgDmassDic = {1: 1, 2: 0.5, 3: 0.33333, 4: 0.25, 5: 0.2}
    chr_list = [1, 2, 3, 4, 5]
    dmassChargeList = list(chrgDmassDic.values())
    ms2MzList = sorted(list(ms2_dic.keys()))
    ms2_info_list = ms2_dic.items()
    min_rank = min([x[1][-1] for x in ms2_info_list])
    top1_mz = [x[0] for x in ms2_info_list if x[1][-1] == min_rank][0]
    top1_mz_ints = ms2_dic[top1_mz][0]
    startPos = ms2MzList.index(top1_mz)
    #print(top1_mz, startPos)
    if startPos >= len(ms2MzList) - 2:
        undeter_iso_dic[top1_mz] = [0, top1_mz]
        del ms2_dic[top1_mz]
        return ms2_dic
    else:
        i = startPos + 1
        while i < len(ms2_info_list):
            mz_cur = ms2MzList[i]
            ints_cur = ms2_dic[mz_cur][0]
            if mz_cur > generate_mass_range(top1_mz+1, 20)[1]:
                undeter_iso_dic[top1_mz] = [0, top1_mz]
                del ms2_dic[top1_mz]
                return ms2_dic
            else:
                delta_to_top1 = mz_cur - top1_mz
                delta_min_ther = [abs(x - delta_to_top1) for x in dmassChargeList]
                min_delta_idx = delta_min_ther.index(min(delta_min_ther))
                ther_plus1_mz = top1_mz + dmassChargeList[min_delta_idx]
                mzPlus1Bool = compareTwoNum(mz_cur, ther_plus1_mz, mstol)
                intensPlus1Bool = ints_cur/top1_mz_ints > 0.3
                if not mzPlus1Bool or not intensPlus1Bool:
                    i += 1
                else:
                    charge = chr_list[min_delta_idx]
                    interval = dmassChargeList[min_delta_idx]
                    iso_cluster = [top1_mz, mz_cur]
                    forword_looking(top1_mz, interval, ms2MzList, startPos, iso_cluster)
                    back_looking(top1_mz, interval, ms2MzList, startPos, iso_cluster, ms2_dic)
                    if len(iso_cluster) == 2:
                        chged_2peaks_dic[iso_cluster[0]] = [charge, iso_cluster[0], iso_cluster]
                    else:
                        chged_iso_dic[iso_cluster[0]] = [charge, iso_cluster[0], iso_cluster]
                    
                    for peak in iso_cluster:
                        del ms2_dic[peak]
                    return ms2_dic

                undeter_iso_dic[top1_mz] = [0, top1_mz]        
                del ms2_dic[top1_mz]
                return ms2_dic
    

def detectIsotopic(ms2_dic):
    undeter_iso_dic = {}
    chged_iso_dic = {}
    chged_2peaks_dic = {}
    while ms2_dic != {}:
        ms2_dic = deter_min_isotpic(ms2_dic, undeter_iso_dic, chged_iso_dic, chged_2peaks_dic)

    return undeter_iso_dic, chged_iso_dic, chged_2peaks_dic


def write_dic2fl(chged_iso_dic, b):
    for mono in chged_iso_dic:
        b.write(",".join([str(x) for x in chged_iso_dic[mono]])+"\n") 


def detect_32Da(chged_iso_dic):
    filter_iso_dic = {}
    info_list = []
    for mono in chged_iso_dic:
        info_list.append(chged_iso_dic[mono][:2])
    info_list = sorted(info_list, key = lambda x:x[1])
    i = 0
    while i < len(info_list) - 1:
        charge, mz = info_list[i]
        mz_list = [x[1] for x in info_list if x[0] == charge]
        idx = mz_list.index(mz)
        if idx == len(mz_list) - 1:
            findBool = False
        else:
            findBool, matchedmz = findOneMZ(mz + 31.996/charge, mz_list, idx+1, 30)
        if not findBool:
            i += 1
        else:
            filter_iso_dic[mz] = chged_iso_dic[mz]
            filter_iso_dic[matchedmz] = chged_iso_dic[matchedmz]
            info_list.remove([charge, mz])
            info_list.remove([charge, matchedmz])
    
    return filter_iso_dic


def merge2Dic(dic1, dic2):
    merge_dic = {}
    for k,v in dic1.items():
        merge_dic[k] = v
    for k,v in dic2.items():
        merge_dic[k] = v
    return merge_dic


def reorginize_spec(undeter_iso_dic, chged_iso_dic, chged_2peaks_dic, specInfo):
    wlist = []
    for mz in undeter_iso_dic:
        line_list = [mz, specInfo[mz][0], 0]
        wlist.append(" ".join([str(x) for x in line_list]))
    for mz in chged_2peaks_dic:
        intns = sum([specInfo[m][0] for m in chged_2peaks_dic[mz][-1]])
        charge = chged_2peaks_dic[mz][0]
        line_list = [mz, round(intns, 1), charge]
        wlist.append(" ".join([str(x) for x in line_list]))
    for mz in chged_iso_dic:
        intns = sum([specInfo[m][0] for m in chged_iso_dic[mz][-1]])
        charge = chged_iso_dic[mz][0]
        line_list = [mz, round(intns, 1), charge]
        wlist.append(" ".join([str(x) for x in line_list]))
    
    wlist = sorted(wlist, key = lambda x: float(x.split(" ")[0]))
    return wlist  


def deIsotopic(ms2Filepath):
    flname = os.path.basename(ms2Filepath)
    parePath = os.path.dirname(ms2Filepath)
    deSisoName = flname[:-4] + "_deisotopic" + flname[-4:]
    f = open(ms2Filepath, "r").readlines()
    b = open(os.path.join(parePath, "report_pair.csv"), 'w')
    b2 = open(os.path.join(parePath, deSisoName), 'w')
    scanInfoDic = {}
    i = 0
    while i < len(f):
        if f[i][0] != "S":
            b2.write(f[i])
            i += 1
        else:
            scan = int(f[i].split("\t")[1])
            print("the scan is %d\n" % scan)
            mzPre = float(f[i].split("\t")[-1])
            preCharge = f[i + 9].split("\t")[1]
            nce = f[i + 6].split("\t")[2].split(" ")[7].split("@")[1][3:5]
            for k in range(i, i + 10):
                b2.write(f[k])

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
            ms2_dic = copy.deepcopy(specInfo)
            undeter_iso_dic, chged_iso_dic, chged_2peaks_dic = detectIsotopic(ms2_dic)
            wlist = reorginize_spec(undeter_iso_dic, chged_iso_dic, chged_2peaks_dic, specInfo)
            for line in wlist:
                b2.write(line + "\n")
            filter_iso_dic = detect_32Da(merge2Dic(chged_iso_dic, chged_2peaks_dic))
            if filter_iso_dic != {}:
                b.write("the scan is %d\n" % scan)
                #print(filter_iso_dic)
                write_dic2fl(filter_iso_dic, b)
            i = p
            
    b.close()
    b2.close()


if __name__ == "__main__":
    deIsotopic(ms2path)