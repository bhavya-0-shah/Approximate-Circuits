import os
import numpy as np
# from matplotlib import pyplot as plt
import pandas as pd
from pandas import Series
import random
import array
import math


def get_bin(x, n): return format(x, 'b').zfill(n)
#format(x, n)-->The number x is formatted as a binary number 
#zfill(n)-->The binary number x will be preceded by as many number of zeroes as required to make x an 'n' bit binary number.

def BIT_G(A, B): return A & B
#function to obtain the Bit Generate Signal
#&-->Biwise AND Opertor


def BIT_P(A, B): return bin(int(A, 2) ^ int(B, 2))
#function to obtain the Bit Propagate Signal 
#int(A,2)-->A is a string. 2 is the base of the number contained in the string, i.e A is a mumber of base 2 (binary number).
        #int(A,2)-->converts A which is a binary number into decimal.
#int(B,2)-->B is a string. 2 is the base of the number contained in the string, i.e B is a mumber of base 2 (binary number).
        #int(B,2)-->converts B which is a binary number into decimal.
#^-->Bitwise XOR Operator
#bin(n)-->returns the binary version of n.

def BC(g, p, valency):
    gg = g #gg[3,2,1,0] = G3G2G1G0(3,2,1,0)
    gp = p #gp[3,2,1,0] = P3P2P1P0(3,2,1,0)
    gg[0] = g[0] #gg[0] = G0
    gp[0] = p[0] #gp[0] = P0

    for i in range(0, valency - 1):
        gg[i + 1] = g[i + 1] | (gg[i] & p[i + 1])
        gp[i + 1] = p[i + 1] & gp[i]
    return gg[valency - 1], gp[valency - 1]


def GC(g, p, valency):
    gg = g #gg[3,2,1,0] = G3G2G1G0[3,2,1,0] ; g = gg = [0,G1] ; g = gg = [G1,G2]; g = gg = [G2 + P2G1, G3]; g = gg = [0,G[4:1]]
    gp = p #gp[3,2,1,0] = P3P2P1P0[3,2,1,0] ; p = gp [0,P1] ; p = gp = [0,P2]; p = gp = [0,P3] ; p = gp  = [0,P[4:1]]
    gg[0] = g[0] #gg[0] = G1 ; gg[0] = g[0] = 0 ; gg[0] = g[0] = G1; gg[0] = g[0] = G2 + P2G1
    gp[0] = p[0] #gp[0] = P1 ; gp[0] = p[0] = 0 ; gp[0] = p[0] = 0 ; gp[0] = p[0] = P3

    for i in range(0, valency - 1): #valency = 2, i = 0; 
        gg[i + 1] = g[i + 1] | (gg[i] & p[i + 1])
        #gg[1] = g[1] + gg[0]&p[1] --> gg[1] = G1 + 0&P1 --> gg[1] = G1;
        #gg[1] = g[1] + gg[0]&p[1] --> gg[1] = G2 + P2G1 --> gg[1] = G2 + P2G1
        #gg[1] = g[1] + gg[0]&p[1] --> gg[1] = G3 + P3gg[0] --> gg[1] = G3 + P3(G2 + P2G1)
        #gg[1] = g[1] + gg[0]&p[1] --> gg[1] = G[4:1] + P[4:1]&0 --> gg[1] = G[4:1]
    return gg[valency - 1] #returns gg[1] = G1; returns gg[1] = G2 + P2G1; returns gg[1] = G3 + P3(G2 + P2G1); returns gg[1] = G[4:1]


def MUX_2_1(A, B, sel):
    if sel == 1:
        return A
    else:
        return B
#function to model 2x1 Multiplexer 
#input to the Multiplexer are A,B and the select line (sel)
#if the sel = 1, A is the output of the multiplexer.
#if the sel = 0, B is the output of the multiplexer.


def CLA_p(A, B, Cin, SIZE, valency): #function to model the Carry Look Ahead Adder
#Carry Look Ahead Adder uses black cells of valency 4 and gray cells of valency 2.
#Let us consider a carry look ahead adder which is a 4 bit adder.
#Output of the Black Cell --> G[4:1] = G4 + P4G3 + P4P3G2 + P4P3P2G1 (Cout of the CLA Adder)
                         #--> P[4:1] = P4P3P2P1
#G0-->G0 is the carry-in to the adder
#G[4:0] = G[4:1] + P[4:1]G0 --> Cout 
#S[4:1]-->S4S3S2S1 is the 4 bit sum.
#S4 = P4^G[3:0] and G[3:0] = G3 + P3G[2:0]
#S3 = P3^G[2:0] and G[2:0] = G2 + P2G[1:0]
#S2 = P2^G[1:0] and G[1:0] = G1 + P1G0
#S1 = P1^G[0:0] and G[0:0] = G0
#A is an operand to the adder
#B is an operand to the adder
#Cin - input carry 
#Size - Size of the adder

    g = [] # g is being initialized to an empty list
    g.insert(0, 0) #insert function is used to insert the value 0 at position 0 of the list g
    p = [] # p is being initialized to an empty list
    p.insert(0, 0) #insert function is used to insert the value 0 at the position 0 of the list p
    GG1 = [] # GG1 is being initialized to an empty list
    GG1.insert(0, 0) #insert function is used to insert the value 0 at the position 0 of the list GG1
    GP1 = [] # GP1 is being intialized to an empty list
    GP1.insert(0, 0) #insert function is used to insert the value 0 at the Position 0 of the list GP1
    SUM = [] #Sum is being intialized to an empty list

    for i in range(0, SIZE):
        # a_i = int(A & (1 << i) > 0)
        # b_i = int(B & (1 << i) > 0)
        g = g + [int(A[i], 2) * int(B[i], 2)]
        p = p + [int(A[i], 2) ^ int(B[i], 2)]
#Let SIZE = 4
#when i = 0, A[0] is a binary number given as a string and it is converted into its decimal equivalent 
            #B[0] is a binary number given as a string and it is converted into its decimal equivalent
#when i = 1, A[1] is a binary number given as a string and it is converted into its decimal equivalent 
            #B[1] is a binary number given as a string and it is converted into its decimal equivalent
#when i = 2, A[2] is a binary number given as a string and it is converted into its decimal equivalent 
            #B[2] is a binary number given as a string and it is converted into its decimal equivalent
#when i = 3, A[3] is a binary number given as a string and it is converted into its decimal equivalent 
            #B[3] is a binary number given as a string and it is converted into its decimal equivalent
#At the end of the for loop, g = G3G2G1G0(3,2,1,0)
                            #p = P3P2P1P0(3,2,1,0)
    GG = [] #initializing GG to an empty list

    GG.insert(0, Cin) #the first index of the list is made 0

    for i in range(0, int(SIZE / valency)): #let size be = 4, valency = 4, then size/valency = 1 and i = 0 and not including 1
        X, Y = BC(g[valency * i + 1:(i + 1) * valency + 1],
                  p[valency * i + 1:(i + 1) * valency + 1], valency)
        
        # when i = 0, X = Black_Cell(g(1:5), p(1:5), 4) --> X = (G(4:1)) = G4 + P4[G(3:1)]
                    # Y = Black Cell(g(1:5), p(1:5), 4) --> Y = (P(4:1)) = P4[P(3:1)] = P4P3P2P1     
        GG1.insert(i + 1, X) 
        GP1.insert(i + 1, Y)
        #when i = 0, GG1.insert(1,X)-->(0,X = G4 + P4[G(3:1)])(0,1)--> (0,G[4:1])
                   # GP1.insert(1,y)-->(0,Y = P4P3P2P1)(0,1) --> (0,P[4:1])
        
        for j in range(0, valency - 1): #Let valency = 4 and j - 0,1,2 and i = 0
            L = GC(list(GG[i * valency + j:i * valency + j + 1]) + list(g[i * valency + j + 1:i * valency + j + 2]),
                   list([0]) + list(p[i * valency + j + 1:i * valency + j + 2]), 2)

            #L will have the output returned by the function GC --> Gray Cell
            #Gray Cell is a function which takes the inputs g, p, valency. Here we are having a Gray Cell of valency of 2
            #when i = 0
                    #when j = 0
                    #GC(list(GG[0:1]) + list(g[1:2]), list([0]) + list(p[1:2], 2)
                    #list(GG[0:1]) = 0 converted into a list 
                    #list(g[1:2]) = G1. Therefore the final list is [0,G1]
                    #list([0]) = 0 converted into a list
                    #list(p[1:2]) = P1. Therefore the final list is [0,P1]
                    #Arguments Passed are = [0,G1],[0,P1] and valency = 2

                    #when j = 1
                    #GC(list(GG[1:2]) + list(g[2:3]), list([0]) + list(p[2:3], 2)
                    #list(GG[1:2]) = G1 converted into a list 
                    #list(g[2:3]) = G2. Therefore the final list is [G1,G2]
                    #list([0]) = 0 converted into a list.
                    #list(p[2:3]) = P2. Therefore the final list is [0,P2]
                    #Arguments Passed are = [G1,G2],[0,P2] and valency = 2

                    #when j = 2
                    #GC(list(GG[2:3]) + list(g[3:4]), list([0]) + list(p[3:4], 2)
                    #list(GG[2:3]) = G2 + P2G1 converted into a list 
                    #list(g[3:4]) = G3. Therefore the final list is [G2 + P2G1,G3]
                    #list([0]) = 0 converted into a list.
                    #list(p[2:3]) = P3. Therefore the final list is [0,P3]
                    #Arguments Passed are = [G2 + P2G1,G3],[0,P3] and valency = 2
                    

            GG.insert(i * valency + j + 1, L) #GG = [0,G1]; [0,G1,(G2 + P2G1)]; [0,G1,(G2 + P2G1),((G3 + P3(G2 + P2G1)))]
            #GG = [0, G[1:0], G[2:0], G[3:0]]

        Z = GC(list(GG[i * valency:i * valency + 1]) +
               list(GG1[i + 1:i + 2]), list([0]) + list(GP1[i + 1:i + 2]), 2)

            #when i = 0
            #GC(list(GG[0:1]) + list(GG1[1:2])) , (list([0]) + list(GP1[1:2])), 2)
            #GG[0:1] = 0 --> Carry in converted to a list
            #GG1[1:2] = G[4:1] converted to a list. Therefore, final list is [0,G[4:1]]
            #list([0]) = 0 converted to a list
            #GP1[1:2] = P[4:1] converted to a list. Therefore the final list is [0.P[4:1]]
            #The arguments passed to the Gray Cell function are ([0,G[4:1], [0,P[4:1]], 2])

        GG.insert((i + 1) * valency, Z) #[0,G1,(G2 + P2G1),((G3 + P3(G2 + P2G1))), G[4:1]]
    for i in range(0, SIZE): #SIZE = 4
        S = p[i + 1] ^ GG[i]
        SUM.insert(i, S) 

        # SUMINTR = int("".join(str(x) for x in SUMR), 2)

    return SUM, GG[SIZE]


def ACAII_p(a, b, Cin, N, R, P, valency):
    AR = get_bin(a, N) #AR = a is converted to an N bit binary number (MSB to LSB)
    BR = get_bin(b, N) #BR = b is converted to an N bit binary number (MSB to LSB)
    # Cin = get_bin(C,1)
    # print(int(Cin[0],2))
    A = AR[::-1] #A = Reversed AR (LSB to MSB)
    B = BR[::-1] #B = Reversed BR (LSB to MSB)
    cout = [] #Empty list to store Cout 
    SUM = [] #Empty list to store the sum 

    x, y = CLA_p(A[0:R + P], B[0:R + P], Cin, R + P, valency) #/*to calculate the sum of the last sub-adder with 8 bits using CLA Adder.
                                                #The last sub-adder produces an (R+P) bit result but the consequent sub adders
                                                #produce an R bit result andf the remaining P bits in the sub-adder are the 
                                                #prediction bits which predict the carry in to the R bits of the sub-adder which 
                                                #determine the result.*/
    SUM.extend(x) #x has the sum calculated by the last sub-adder. This sum gets added to the list of SUM. 
    cout.append(y) #y has the Cout calculated by the last sub-adder. This Cout gets appended to the list cout.

    for j in range(0, (int(N/R) - 2)):
        z, l = CLA_p(A[N - ((j + 2)*R):N - (j*R)],
                     B[N - ((j + 2)*R):N - (j*R)], 0, R + P, valency)
        SUM.extend(z[P:R + P]) # z will have the sum calculated by the sub adders except for the last sub-adder. This sum is
                            #extended onto the list of SUM.
        cout.append(l) # z will be have the carry-out computed by the sub-adders except for the last sub-adder.
    SUM.insert(N, cout[1]) #to insert the carry-out into the final SUM list
    SUMR = SUM[::-1] #to reverse the SUM  to (MSB to LSB) from (LSB to MSB)
    SUMINTR = int("".join(str(x) for x in SUMR), 2)
    #The base of the number contained in the string is 2
    #x in SUMR is converted to string and successive bits are converted into strings and get appended to the empty string ""
    #The final N+1 bit sum in the string is converted into an integer.
    return SUMINTR #returns the integer equivalent sum.


# *************************************************************************************************


N = int(input("enter the total size of adder\n")) #to get the ACA II adder size from the user 
R = int(input("enter the  number of result bits from each subadder.\n")) #to get the number of result bits produced as sum (k) from the user.
P = int(input("enter the  number of prediction bits for each subadder.\n")) #to get the number of prediction bits in each sub-adder from the user
valency = int(input("Enter valency to be used in each subadder.\n")) #valency of the blac-cell used in each sub-adder which are modelled as CLA adders

print("Result is: {0}" #To print the hexadecimal equivalent of the sum of A5C3 and CA53
      .format(
          hex(
              ACAII_p(
                  int("A5C3", 16),
                  int("CA53", 16),
                  0,
                  N,
                  R,
                  P,
                  valency
              )
          )
      )
      )

iterations = 100000 #Number of random number generation iterations
errorcount = 0 #Number of times error occurs
Maximum_Error = 0 #Maximum deviation or error from the actual result 
Total_Error = 0 #Totall error by summing up all the errors occurring.
Average_Error = 0 # absolute error/exact sum for one particular addition
Average_Error_Max = 0
A_Max = 0 #Maximum value of input a for which absolute error/exact sum is greater than maximum average error
B_Max = 0 #Maximum value of input b for which absolute error/exact sum is greater than maximum average error
Acceptance_Probability = 0
BIT_Flip_Error_Max = 0
Total_BIT_Flip_Error = 0

for i in range(0, iterations):
    # generate a, b and c randomly
    a = random.randint(0, 2 ** N - 1) #Random Number between 0 and 2^N - 1
    b = random.randint(0, 2 ** N - 1) #Random Number between 0 and 2^N - 1
    c = random.randint(0, 1) #Random input carry for the addition
    # a = i + 1
    # b = math.floor(iterations / 2)
    # c = 0

    # calculate exact and approximate sum
    EXACT_SUM = a + b + c
    APPROX_SUM = ACAII_p(a, b, c, N, R, P, valency)

    # print(a,b,c,EXACT_SUM, APPROX_SUM)
    # if error :
    if (EXACT_SUM != APPROX_SUM):
        errorcount = errorcount + 1 #incrementing the errorcount if there is an occurence for error
        absolute_error = abs(EXACT_SUM - APPROX_SUM) #stores the absolute value of error = Exact Sum - Approximate Sum 
        if absolute_error > Maximum_Error: #initially Maximum_Error = 0
            Maximum_Error = absolute_error
        Total_Error = Total_Error + absolute_error #to add the current error to the total error
        Average_Error = Average_Error + (absolute_error / EXACT_SUM) # to update the average error after every occurrence of error
        if (absolute_error / EXACT_SUM) > Average_Error_Max: #if average error of this addition is > Maximum Average Error
            Average_Error_Max = (absolute_error / EXACT_SUM) #then maximum average error = (absolute error/exact sum)
            A_Max = a #will have input a for which the absolute error is maximum
            B_Max = b #will have input b for which the absolute error is maximum
    if (EXACT_SUM - APPROX_SUM) / EXACT_SUM < 0.01: #Acceptance ratio is 0.01 i.e (exact sum - approximate sum)/exact sum must be less than 0.01
        Acceptance_Probability = Acceptance_Probability + 1 #to calculate the number of instances where the approximate result will be in the acceptable range

    EXACT_SUM_BR = get_bin(EXACT_SUM, N + 1) #to get the binary equivalent of the exact sum (in integer initially)
    EXACT_SUM_B = EXACT_SUM_BR[::-1] #to reverse this N+1 bit binary exact sum (LSB to MSB) (in integer initially)

    APPROX_SUM_BR = get_bin(APPROX_SUM, N + 1) #to get the binary equivalent of the spproximate sum(in integer initially)
    APPROX_SUM_B = APPROX_SUM_BR[::-1] #to get the binary equivalent of the approximate sum (in integer initially)

    BIT_Flip_Error = 0 #to calculate the number of instances of Bit Flip Error
    for i in range(0, N + 1):
        if int(APPROX_SUM_B[i]) != int(EXACT_SUM_B[i]): #comparing exact sum and approximate sum bit by bit from LSB to MSB
            BIT_Flip_Error = BIT_Flip_Error + 1 #if two bits are not the same, then there is a bit flip error.
        if BIT_Flip_Error > BIT_Flip_Error_Max:
            BIT_Flip_Error_Max = BIT_Flip_Error #to store the maximum number of times Bit Flip Error has occurred across all addition results
    Total_BIT_Flip_Error = Total_BIT_Flip_Error + BIT_Flip_Error #to calculate the total number of instances of Bit Flip

# *************************************************************************
# Results to be written in file
# *************************************************************************

curr_dict = os.getcwd()
dir_name = "ACAII_ERROR_ANALYSIS"
path_name = os.path.join(curr_dict, dir_name)

if dir_name not in os.listdir():
    os.mkdir(path_name)

# open file for write access
filename = "N{0}_R{1}_P{2}_valency{3}.txt".format(N, R, P, valency)
file = open(os.path.join(path_name, filename), "w")

# write into file : heading and parameters
file.write("ACAII ERROR ANALYSIS\n\n")
file.write("Parameters :\n")
file.write("\tN = {0}\n".format(N))
file.write("\tR = {0}\n".format(R))
file.write("\tP = {0}\n".format(P))
file.write("\tvalency = {0}\n".format(valency))

# write stats into file for chosen error attributes

file.write("\nNumber of iterations\t: {0}\n".format(iterations))
file.write("Number of error cases\t: {0}\n".format(errorcount))
file.write("Error Rate = errorcount/no. of iterations\n\t")
file.write("= {0:.3f}%\n\n".format((100 * errorcount / iterations)))

file.write("Maximum Hamming Distance = {0}\n".format(BIT_Flip_Error_Max))
file.write("Average Hamming Distance = Total_BIT_Flip_Error / iterations\n\t")
file.write("= {0:.3f} (Over all iterations)\n\t".format(
    Total_BIT_Flip_Error / iterations))
if errorcount != 0:
    file.write("= {0:.3f} (Over erroneous iterations)\n\n".format(
        Total_BIT_Flip_Error / errorcount))

file.write("Maximum Error = {0} (Absolute value)\n".format(Maximum_Error))
file.write(
    "Average error(sum of errors/iterations) = {0} \n\t".format(Average_Error))
if errorcount != 0:
    file.write("= {0:.3f} (Over erroneous iterations)\n\n".format(
        Total_Error / errorcount))
file.write("Average Maximum Error = {0} \n\n".format(Average_Error_Max))

file.write("Acceptance Probability Over iterations\n\t")
file.write("= {0:.3f}  (in %) \n\n".format(
    100 * Acceptance_Probability / iterations))