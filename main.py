from readRawConverterMS2 import getScanInfoDic, addRankRelInt2Spec
from detect_isotopic import findOneMZ, generate_mass_range, compareTwoNum
mstol = 20
scan_info_dic = getScanInfoDic("./DSSO_CID_MS2_HCD_IT_MS3_SinPair_re23_T1.ms2")
#print(scan_info_dic)



def deter_min_isotpic(ms2_dic, undeter_iso_dic, chged_iso_dic):
    #print(ms2_dic)
    chrgDmassDic = {1: 1, 2: 0.5, 3: 0.33333, 4: 0.25, 5: 0.2}
    chr_list = [1, 2, 3, 4, 5]
    dmassChargeList = list(chrgDmassDic.values())
    ms2MzList = sorted(list(ms2_dic.keys()))
    ms2_info_list = ms2_dic.items()
    top1_mz = [x[0] for x in ms2_info_list if x[1][-1] == 1][0]
    startPos = ms2MzList.index(top1_mz)
    #print(top1_mz, startPos)
    if startPos >= len(ms2MzList) - 2:
        undeter_iso_dic[top1_mz] = [0, top1_mz]
        del ms2_dic[top1_mz]
        return ms2_dic
    else:
        i = startPos + 1
        #print(i)
        while i < len(ms2_info_list):
            mz_cur = ms2MzList[i]
            if mz_cur > generate_mass_range(top1_mz+1, 20)[1]:
                undeter_iso_dic[top1_mz] = [0, top1_mz]
                del ms2_dic[top1_mz]
                return ms2_dic
            else:
                delta_to_top1 = mz_cur - top1_mz
                delta_min_ther = [abs(x - delta_to_top1) for x in dmassChargeList]
                min_delta_idx = delta_min_ther.index(min(delta_min_ther))
                ther_plus1_mz = top1_mz + dmassChargeList[min_delta_idx]
                mzPlus1Bool = compareTwoNum(ms2MzList[i], ther_plus1_mz, mstol)
                intensPlus1Bool = ms2_dic[ms2MzList[i]][0]/ms2_dic[top1_mz][0] > 0.3
                if not mzPlus1Bool or not intensPlus1Bool:
                    i += 1
                else:
                    charge = chr_list[min_delta_idx]
                    interval = dmassChargeList[min_delta_idx]
                    match_bool, matchedMZ = findOneMZ(top1_mz+interval*2, ms2MzList, i, 20)
                    if not match_bool:
                        undeter_iso_dic[top1_mz] = [0, top1_mz]
                        del ms2_dic[top1_mz]
                        return ms2_dic
                    else:
                        iso_cluster = [top1_mz, ms2MzList[i], matchedMZ]
                        m = 3
                        while m < 8:
                            [matchMoreBool, matchMoreMZ] = findOneMZ(float(top1_mz)+interval*m, ms2MzList, i, 20)
                            if matchMoreBool:
                                iso_cluster.append(matchMoreMZ)
                                m += 1
                            else:
                                break
                        if iso_cluster[0] * charge > 1600:                                
                            for n in range(1,4):
                                matchLessBool, matchLessMZ = findOneMZ(top1_mz - interval*n, ms2MzList, startPos, 20)
                                if not matchLessBool:
                                    break
                                else:
                                    if ms2_dic[matchLessMZ][0]/ms2_dic[iso_cluster[0]][0] > 0.3:
                                        iso_cluster.insert(0, matchLessMZ)
                                    else:
                                        break
                        chged_iso_dic[iso_cluster[0]] = [charge, iso_cluster[0], iso_cluster]
                        for peak in iso_cluster:
                            del ms2_dic[peak]
                        return ms2_dic
                del ms2_dic[top1_mz]
                return ms2_dic
    

def reform_ms2dic(ms2_dic): 
    #print(ms2_dic) 
    if ms2_dic == {}:
        return ms2_dic
    else:
        new_ms2_dic = {}
        for mz in ms2_dic:
            new_ms2_dic[mz] = [ms2_dic[mz][0]]
        addRankRelInt2Spec(new_ms2_dic)
        return new_ms2_dic
    

def detectIsotopic(ms2_dic):
    #orig_ms2_dic = ms2_dic.copy.deepcopy()
    undeter_iso_dic = {}
    chged_iso_dic = {}
    
    while ms2_dic != {}:
        new_ms2_dic = deter_min_isotpic(ms2_dic, undeter_iso_dic, chged_iso_dic)
        ms2_dic = reform_ms2dic(new_ms2_dic)
        #print(ms2_dic)
    return undeter_iso_dic, chged_iso_dic


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



b = open("report.csv", 'w')
for key in scan_info_dic:
    print("the scan is %d\n" % key)
    
    ms2_dic = scan_info_dic[key][3]
    undeter_iso_dic, chged_iso_dic = detectIsotopic(ms2_dic)
    filter_iso_dic = detect_32Da(chged_iso_dic)
    if filter_iso_dic != {}:
        b.write("the scan is %d\n" % key)
        print(filter_iso_dic)
        write_dic2fl(filter_iso_dic, b)

b.close()
    
                        
                    
