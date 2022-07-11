import math
# from math import comb
from itertools import chain, combinations


# Function to find PowerSet of a set
def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


def probabilityPGeo(PA, N, k1, k2, k):
    sum1 = 0
    for i in range(0, 2 ** (N-1-k2)):
        sum2 = 0
        for j in range(0, 2 ** k1):
            sum2 = sum2 + (PA * (2 ** (k2+1)*i + 2 ** k1 * k + j ))
        sum1 = sum1 + sum2
        return sum1



# Probability of P (Prediction bit group) for uniform distribution


def probabilityP(n):
    return 1 / (2 ** n)


# Probability of G0 (Generation with 0 carry-in  bit group) for uniform distribution

def probabilityG0(n):
    return (2 ** n - 1) / (2 ** (n + 1))


# Probability of G1 (Generation with 1 carry-in  bit group) for uniform distribution

def probabilityG1(n):
    return (1 / (2 ** n)) + (2 ** n - 1) / (2 ** (n + 1))


def minimum(a, n):
    minpos = a.index(min(a))
    return minpos


def check(list, element):
    for i in list:
        if i == element:
            return True
    return False


##################################################################################################################
# Algorithm 2 of paper "Probabilistic Error Modeling for Approximate Adders
#                                      (IEEE TRANSACTIONS ON COMPUTERS, VOL. 66, # NO. 3,# MARCH 2017)"
# Joint probability (of intersection) of error events.
# Inputs are (sub-adder-wise) Prediction bit array, Generation bit array and index of sub-adder involved in intersection
#####################################################################################################################
def JointProbability(S, R, I):  # Statement 2 from Algorithm is direct input
    #print("---------------------------------------------------------------------------------------")
    r_i = []
    S_k = 0
    for i in range(0, len(R)):
        S_k = S_k + S[i]
        r_i.append([j for j in range(S_k - R[i], S_k)])
    #print("###################################################################")
    #print("I  :", I)
    #print("r_i ", r_i)
    #print("###################################################################")
    R_U = []  # P (Prediction) bits set (RU used in Statement 3)
    g_all = []
    g_i = []
    G_D = []

    for ele in I:
        R_U = list(set().union(R_U, r_i[ele - 2]))  # Statement 3 RU Combining all prediction bits of all intersection events
    PP = probabilityP(len(R_U))  # Statement 4
    #print("###################################################################")
    #print("I  :", I)
    #print("R_U ", R_U)
    #print("###################################################################")
    S_R = 0
    for i in range(0, len(R)):
        S_R = S_R + S[i]
        g_all.append([j for j in range(0, S_R - R[i])])
    #print("g_all jP", g_all)
    for ele in I:
        g_i.append(g_all[ele - 2])  # Statement 5 set of set g_i: Combining all generation bits of all intersection events
    #print("###################################################################")
    #print("I  :", I)
    #print("g_i ", g_i)
    #print("###################################################################")
    for i in range(0, len(g_i)):
        g_i[i] = list(set(g_i[i]) - set(R_U))  # Statement 6: updating set of set gi=gi-RU
    #print("###################################################################")
    #print("I  :", I)
    #print("g_i-R_U ", g_i)
    #print("###################################################################")
    G_I = list(g_i[0])
    for index in range(1, len(I)):
        G_I = list(set(G_I).intersection(g_i[index]))  # Statement  7 GI: intersection of all updated g_i
    #print("###################################################################")
    #print("I  :", I)
    #print("G0=G_I  :", G_I)
    #print("###################################################################")
    SG = [G_I]  # Statement 8 initialization of set of all generation events G0 as well as G1
    #print("----------SG", SG)
    PP = PP * probabilityG0(len(G_I))  # Statement 9
    # print("P*G0", PP)
    for i in range(0, len(g_i)):
        if len(list(set(g_i[i]) - set(G_I))) != 0:
            G_D.append(list(set(g_i[i]) - set(G_I)))  # Statement 10: Difference set G_D
    #print("###################################################################")
    #print("I  :", I)
    #print("G_D=g_i - G_I difference set  :", G_D)
    #print("###################################################################")
    #  Statement 11 : for loop begins
    for i in range(0, len(G_D)):
        UPGI = G_D[i]
        for j in range(0, len(G_D)):
            intse = list(set(UPGI).intersection(G_D[j]))
            if len(intse) != 0:
                UPGI = intse
            else:
                UPGI = G_D[i]
            G_D[j] = list(set(G_D[j]) - set(UPGI))
        if len(UPGI) != 0:
            PP = PP * probabilityG1(len(UPGI))
            SG.append(UPGI)
    #print("###################################################################")
    #print("I  :", I)
    #print("Updated SG  :", SG)
    #print("###################################################################")
    #for i in range(0, len(SG)):
    #    SG[i].sort()
    #  Statement 16 : for loop ends
    s_i = [S[0]]
    for i in range(1, len(S)):
        s_i.append(S[i] + s_i[i - 1])
    #print("###################################################################")
    #print("s_i :", s_i)
    #print("###################################################################")
    V = 0
    #if len(SG) == 1:
    #    V= V + (2 ** s_i[max(I)-2])
    #else:
    #    for k in range(0, len(SG)):
    #        V = V + (2 ** (max(SG[k])+1))
    for i in range(0, len(SG)):
        # print("maxSG", max(SG[i]))
        #d_i = []
        #for j in range(0, len(s_i)):
        #    d_i.insert(j, s_i[j] - max(SG[i]))
        #    #d_i.append(s_i[j] - max(SG[i]))
        #minpos = d_i.index(min([k for k in d_i if k > 0]))
        if len(SG[i]) != 0:
            s = max(SG[i]) + R[i] + 1
        else:
            s = R[i] + 1
        #s = s_i[minpos]
        #if len(I) == 1:
        #    V = V + (2 ** s_i[I[0] - 2])
        #else:
        V = V + (2 ** s)
        # print("VA", V)
        #print("###################################################################")
        #print("i   :", i)
        #print("SG :", SG[i])
        #print("d_i :", d_i)
    return PP, V  # Statement 17


def algo1(SSET, RSET, L):
    S = []
    SIZEWISEINDEXOFPS = []
    for i in range(2, L + 1):
        S.append(i)
    PS = list(powerset(S))
    #print("POWERSET", PS)
    P_r_err = 0
    STARTINDEX = 0
    for i in range(0, L):
        SIZEWISEINDEXOFPS.append([STARTINDEX, STARTINDEX + math.comb(L - 1, i) - 1])
        STARTINDEX = STARTINDEX + math.comb(L - 1, i)
    #print("SIZEWISEINDEXOFPS", SIZEWISEINDEXOFPS)
    for i in range(1, len(SIZEWISEINDEXOFPS)):
        for j in range(SIZEWISEINDEXOFPS[i][0], SIZEWISEINDEXOFPS[i][1] + 1):
            P_r_err = P_r_err + (((-1) ** (i+1)) * JointProbability(SSET, RSET, PS[j])[0])
            #print("P_r_err", P_r_err)
    return P_r_err


def algo4(L, SSET, RSET, I):
    S = []
    SIZEWISEINDEXOFPJ = []
    for i in range(2, L + 1):
        S.append(i)
    #print("S--------------------", S)
    #print("I", I)
    J = list(set(S) - set(I))
    #print("J------------------", J)
    P_1 = JointProbability(SSET, RSET, I)[0]
    #print("P_1", P_1)
    P_J = list(powerset(J))
    #print("P_J", P_J)
    P_2 = 0
    STARTINDEX = 0
    for i in range(0, len(J)+1):
        SIZEWISEINDEXOFPJ.append([STARTINDEX, STARTINDEX + math.comb(len(J), i)-1])
        STARTINDEX = STARTINDEX + math.comb(len(J), i)
    #print("SIZEWISEINDEXOFPJ", SIZEWISEINDEXOFPJ)
    for i in range(1, len(SIZEWISEINDEXOFPJ)):
        for j in range(SIZEWISEINDEXOFPJ[i][0], SIZEWISEINDEXOFPJ[i][1] + 1):
            PJUI = list(set().union(P_J[j], set(I)))
            P_2 = P_2 + (((-1) ** (i+1)) * JointProbability(SSET, RSET, PJUI)[0])
    P_r_IJ = P_1 - P_2
    return P_r_IJ


def algo3(SSET, RSET, L):
    S = []
    for i in range(2, L + 1):
        S.append(i)
    PS = list(powerset(S))
    V = []
    p_v = []
    SIZEWISEINDEXOFPS = []
    STARTINDEX = 0
    for i in range(0, L):
        SIZEWISEINDEXOFPS.append([STARTINDEX, STARTINDEX + math.comb(L - 1, i) - 1])
        STARTINDEX = STARTINDEX + math.comb(L - 1, i)
    for i in range(1, len(SIZEWISEINDEXOFPS)):
        for j in range(SIZEWISEINDEXOFPS[i][0], SIZEWISEINDEXOFPS[i][1] + 1):
            P_r_i = algo4(L, SSET, RSET, PS[j])
            #print("P_r_i algo3", P_r_i)
            v = JointProbability(SSET, RSET, PS[j])[1]
            #print("v algo3",v)
            #if len(V) == 0:
            #    V.append(v)
            #    p_v.append(P_r_i)
                #print("V, p_v before__________________________________", V, p_v)
            #else:
            #    m=0
            #    for k in range(0, len(V)):
                    #print("------*************-------", k, v, V)
            #        if V[k] == v:
                        #print(k)
            #            m=1
            #            p_v[k] = p_v[k] + P_r_i
            #            break
            #    if m==0:
            #        V.append(v)
            #        p_v.append(P_r_i)
            m=0
            for k in range(0, len(V)):
                    # print("------*************-------", k, v, V)
                if V[k] == v:
                        # print(k)
                    m = 1
                    p_v[k] = p_v[k] + P_r_i
                    #break
            if m == 0:
                V.append(v)
                p_v.append(P_r_i)
    sumP_v = 0
    for i in range(0, len(p_v)):
        sumP_v = sumP_v + p_v[i]
    p_v.insert(0, 1 - sumP_v)
    return V, p_v


SIZE = int(input("Enter the  Size of adder.\n"))
ES = int(input("Is sub-adder of equal size? Enter 1 for equal size, 2 for unequal size\n"))
ECG = int(input("Is carry generator of equal size? Enter 1 for equal size, 2 for unequal size\n"))
#ZEWISEINDEXOFPS = []  # Start index and end index of powerset elements of same size (all elements of size1,
# all elements of size 2....etc)
S = []
R = []
##########################################################################################
# Input when number of sum bits and prediction bits for all sub-adders is same
##########################################################################################
if ES == 1 and ECG == 1:
    SUMBITS = int(input("Enter the number of sum bits.\n"))
    PREDBITS = int(input("Enter the number of prediction bits.\n"))
    S.append(SUMBITS + PREDBITS)
    for i in range(1, int((SIZE - SUMBITS - PREDBITS) / SUMBITS + 1)):
        S.append(SUMBITS)
    for i in range(1, int((SIZE - SUMBITS - PREDBITS) / SUMBITS + 1)):
        R.append(PREDBITS)
#print("No. of sum bits", S)
#print("No. of prediction bits", R)
###########################################################################################
# Input when size of all sub-adders is same but number of sum bits and prediction bits
# are not same for each sub-adder
###########################################################################################
if ES == 1 and ECG == 2:
    SIZESUBADDR = int(input("Enter the size of sub-adder.\n"))
    S.append(SIZESUBADDR)
    i = 1
    total_sum_bits = SIZESUBADDR
    while total_sum_bits < SIZE:
        SUMBITS = int(input("Enter the number of sum bits  for adder %d." % (i + 1)))
        total_sum_bits = total_sum_bits + SUMBITS
        S.append(SUMBITS)
        R.append(SIZESUBADDR - SUMBITS)
        i = i + 1
# print("No. of sum bits", S)
# print("No. of prediction bits", R)
#########################################################################################################
# Input when number of prediction bits in each sub-adder are same but size of each sub-adder is not same
#########################################################################################################
if ES == 2 and ECG == 1:
    PREDBITS = int(input("Enter the number of prediction bits in all Sub-adders.\n"))
    SIZESUBADDR = int(input("Enter the size of sub-adder 1.\n"))
    S.append(SIZESUBADDR)
    i = 1
    total_sum_bits = SIZESUBADDR
    while total_sum_bits < SIZE:
        SUMBITS = int(input("Enter the number of sum bits  for adder %d." % (i + 1)))
        total_sum_bits = total_sum_bits + SUMBITS
        S.append(SUMBITS)
        R.append(PREDBITS)
        i = i + 1
# print("No. of sum bits", S)
# print("No. of prediction bits", R)
#########################################################################################################
# Input when number of prediction bits, sum bits in each sub-adder are different and sizes of all
# sub-adders are different
#########################################################################################################
if ES == 2 and ECG == 2:
    SIZESUBADDR = int(input("Enter the size of sub-adder 1.\n"))
    S.append(SIZESUBADDR)
    i = 1
    total_sum_bits = SIZESUBADDR
    while total_sum_bits < SIZE:
        SUMBITS = int(input("Enter the number of sum bits  for adder %d." % (i + 1)))
        PREDBITS = int(input("Enter the number of prediction bits for adder %d." % (i + 1)))
        S.append(SUMBITS)
        R.append(PREDBITS)
        i = i + 1
        total_sum_bits = total_sum_bits + SUMBITS
print("No. of sum bits", S)
print("No. of prediction bits", R)

L = len(S)
print("L", L)
print("Probability of error", algo1(S, R, L))
#print("algo4", algo4(L, S, R, 2, 3))
print("algo3", algo3(S, R, L))