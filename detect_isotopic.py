import os
from copy import deepcopy
from read_spec_from_mgf import read_spec
from delta_ppm import generate_mass_range

#os.chdir(r"D:\Onedriver\OneDrive\github\pFind_data_post_analysis.github.io\DETECT_PETIDE_isotpic_cluster")


def compareTwoNum(realNum, TheroNum, tolPPM):
    if abs(TheroNum- realNum)/TheroNum * 1000000 < tolPPM:
        return True
    else:
        return False


def findOneMZ(mz, ms2MzList, begin_idx, tol_ppm):
    [low_ms, up_ms] = generate_mass_range(mz, tol_ppm)
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
            if k == len(ms2MzList):
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
            if k == -1:
                return False, -1


def detectIsotopic(ms2_dic):
    chrgDmassDic = {1: 1, 2: 0.5, 3: 0.33333, 4: 0.25, 5: 0.2}
    chr_list = [1, 2, 3, 4, 5]
    dmassChargeList = list(chrgDmassDic.values())
    ms2_info_list = sorted(ms2_dic.items(), key=lambda d:d[1], reverse = True)
    ms2MzList = list(ms2_dic.keys())
    ms2MzList.sort()
    chged_iso_dic = {}
    i = 0
    while i < len(ms2_info_list):
        mz = ms2_info_list[i][0]
        mz_idx = ms2MzList.index(mz)
        idx = mz_idx + 1
        # fileToWrite.write("\t".join([str(mz), str(mz_idx)])+"\n")
        if mz_idx > len(ms2MzList)-2:
            i += 1
        else:
            find_bool = False
            while idx <= len(ms2MzList)-2:
                if float(ms2MzList[idx]) > generate_mass_range(float(mz)+1, 20)[1]:
                    break
                else:
                    delta =  float(ms2MzList[idx]) - float(mz)
                    delta_min_ther = [abs(x - delta) for x in dmassChargeList]
                    min_delta_idx = delta_min_ther.index(min(delta_min_ther))
                    ther_plus1_mz = float(mz) + dmassChargeList[min_delta_idx]
                    mzPlus1Bool = compareTwoNum(ms2MzList[idx], ther_plus1_mz, 20)
                    intensPlus1Bool = ms2_dic[ms2MzList[idx]]/ms2_dic[mz] > 0.3
                    if mzPlus1Bool and intensPlus1Bool:
                        chrg = chr_list[min_delta_idx]
                        interval = dmassChargeList[min_delta_idx]
                        #print(interval)
                        match_bool, matchedMZ = findOneMZ(float(mz)+interval*2, ms2MzList, idx, 20)
                        if match_bool:
                            find_bool = True
                            iso_cluster = [mz, ms2MzList[idx], matchedMZ]
                            # print(iso_cluster)
                            m = 3
                            while m < 6:
                                [matchMoreBool, matchMoreMZ] = findOneMZ(float(mz)+interval*m, ms2MzList, idx, 20)
                                if matchMoreBool:
                                    iso_cluster.append(matchMoreMZ)
                                    m += 1
                                else:
                                    break
                            if iso_cluster[0] * chrg > 1600:                                
                                for n in range(1,4):
                                    [matchLessBool, matchLessMZ] = findOneMZ(iso_cluster[0]-interval, ms2MzList, idx, 20)
                                    #print(matchLessBool, matchLessMZ)
                                    if matchLessBool:
                                        if ms2_dic[matchLessMZ]/ms2_dic[iso_cluster[0]] > 0.3:
                                            iso_cluster.insert(0, matchLessMZ)
                                        else:
                                            break
                                    else:
                                        break
                            else:
                                pass
                            # print(iso_cluster)
                            chged_iso_dic[iso_cluster[0]] = [chrg, iso_cluster]
                            for peak in iso_cluster:
                                if peak in ms2MzList:
                                    ms2MzList.remove(peak)
                                    ms2_info_list.remove((peak, ms2_dic[peak]))
                                else:
                                    pass
                            # fileToWrite.writelines([str(ele) for ele in ms2_info_list])
                        else:
                            idx += 1
                    else:
                        idx += 1
            if find_bool == True:
                i = i
            else:
                i += 1 
    # fileToWrite.writelines([str(ele) for ele in ms2_info_list])                    
    return chged_iso_dic


def pairedwithDeltaMass(chged_iso_dic, tol_ppm, delta_mass):
    ft_pair_dic = {}
    chr_ms_left_dic = {}
    chr_ms_dic = {}
    for key in chged_iso_dic:
        charge = chged_iso_dic[key][0]
        if charge not in chr_ms_dic:
            chr_ms_dic[charge] = [key]
        else:
            chr_ms_dic[charge].append(key)
    for chrg in chr_ms_dic:
        chr_ms_dic[chrg].sort()
        delta_mz = delta_mass/chrg
        bms_list = deepcopy(chr_ms_dic[chrg])
        
        i = 0        
        while i < len(bms_list)-1:
            mz = bms_list[i]
            ther_pair_mz = mz + delta_mz

            [fd_bool, realMZ] = findOneMZ(ther_pair_mz, bms_list, i, tol_ppm)
            if fd_bool:
                if chrg == 1:
                    [fd50_bool, real50MZ] = findOneMZ(ther_pair_mz+18.010, bms_list, i+1, tol_ppm)
                    if fd50_bool:
                        if chrg not in ft_pair_dic:
                            ft_pair_dic[chrg] = [(bms_list[i], realMZ, real50MZ)]
                        else:
                            ft_pair_dic[chrg].append((bms_list[i],realMZ, real50MZ))
                        bms_list.remove(mz)
                        bms_list.remove(realMZ)
                        bms_list.remove(real50MZ)
                    else:
                        if chrg not in ft_pair_dic:
                            ft_pair_dic[chrg] = [(bms_list[i],realMZ)]
                        else:
                            ft_pair_dic[chrg].append((bms_list[i],realMZ))
                        bms_list.remove(mz)
                        bms_list.remove(realMZ)
                else:        
                    if chrg not in ft_pair_dic:
                        ft_pair_dic[chrg] = [(bms_list[i], realMZ)]
                    else:
                        ft_pair_dic[chrg].append((bms_list[i],realMZ))
                    bms_list.remove(mz)
                    bms_list.remove(realMZ)
                i = i
            else:
                i += 1
        if chrg not in chr_ms_left_dic:
            chr_ms_left_dic[chrg] = bms_list
        else:
            print("wrong")
    return ft_pair_dic, chr_ms_left_dic


def statistcPairNum(ft_pair_dic):
    print(ft_pair_dic)
    pairNum = 0
    for chrg in ft_pair_dic:
        pairNum += len(ft_pair_dic[chrg])
    
    return pairNum


def deal_deisotopiced_mgf(mgfFl):
    ft_pair_dic ={}
    if mgfFl[-14:] != "deisotopic.mgf":
        print('wrong')
        return ft_pair_dic
    else:
        fl = open(mgfFl, 'r').readlines()
        
        i = 0
        while i < len(fl):
            if fl[i].strip() != "BEGIN IONS":
                i += 1
            else:
                scan = fl[i+1].split(".")[1]
                pairDic = {}
                p = i + 5
                while p < len(fl):
                    if fl[p].strip() == "END IONS":
                        break
                    else:
                        lineList = fl[p].strip().split(' ')
                        if len(lineList) == 3:
                            chrg = int(lineList[-1])
                            if chrg not in pairDic:
                                pairDic[chrg] = [float(lineList[0])]
                            else:
                                pairDic[chrg].append(float(lineList[0]))
                        p += 1
                print(pairDic)
                pickedPairDic = pairedwithDeltaMass(pairDic, 30, 32)[0]
                pairedNum  = statistcPairNum(pickedPairDic)
                #print(pairedNum)
                ft_pair_dic[scan] = pairedNum
                i = p+1
    return ft_pair_dic


# print(deal_deisotopiced_mgf('rep_dsso_deisotopic.mgf'))


def main():
    mgf_info = read_spec("./test1.mgf")
    # print(mgf_info)
    b = open("report.txt", 'w')
    pairNumDic = {}
    for ttl in mgf_info:
        ms2_dic = mgf_info[ttl]
        chged_iso_dic = detectIsotopic(ms2_dic)
        print(chged_iso_dic)
        t_pair_dic, chr_ms_left_dic = pairedwithDeltaMass(chged_iso_dic, 30, 31.97)
        # print(t_pair_dic)
        pairNum = statistcPairNum(t_pair_dic)
        if pairNum not in pairNumDic:
            pairNumDic[pairNum] = 1
        else:
            pairNumDic[pairNum] += 1
    for num in sorted(list(pairNumDic.keys())):
        b.write(str(num) + '\t' + str(pairNumDic[num]) + '\n')
    b.close()


if __name__ == "__main__":
    main()
