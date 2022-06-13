"""
Applied Programming Lab(EE2703) - Assignment 2
Name: Nipun Suresh
Roll No: EE20B091

In this program, we take a circuit as input and solve it using Modified Nodal Analysis
Remarks:
1. Branch current direction has been taken opposite to that of conventional algorithm, to ensure
it is positive for single source case
"""

from sys import argv, exit
import numpy as np
import cmath as cm
 # to avoid hardcoding netlist syntax
CKT = '.circuit' ; END = '.end' ; AC = ".ac"

class Component:    #class that represent generic electrical 2-terminal component
    def __init__(self,lst,omega=0):
        self.type = lst[0][0]
        self.nodeA = lst[1]
        self.nodeB = lst[2]
        val = 0
        val = float(lst[4]) if self.type == "V" or self.type == "I" else float(lst[3])
        if omega==0:    #DC case
            self.val = val
        else:
            type = self.type
            compVal = complex(0,0)  #compVal stores the complex impedance or voltage
            if type == "R":
                compVal = val
            elif type == "L":
                compVal = complex(0,omega*val)  # Xl = jwL
            elif type == "C":
                compVal = complex(0,-1/(omega*val)) # Xc = -j/wC
            elif type == "V" or type == "I":
                phase = float(lst[5])
                compVal = cm.rect(val/2,phase)
            self.val = compVal

#this function is the notation for branch current through voltage source
branch_i = lambda a,b: "i("+a+","+b+")"
circuit_list = [] ; rows_mna = {} ; omega = 0

try:    #try-except block to ensure file name entered is valid
    with open(argv[1]) as f: 
        lines = f.readlines()   #lines variable contains all the lines
        start = -1; end = -2
        for line in lines:  
            test_word = (line.split())[0]            
            if test_word == CKT:      #to know where circuit definition begins
                start = lines.index(line)
            elif test_word == END:    #to know where circuit definition ends
                end = lines.index(line)
            elif test_word == AC:
                omega = 2 * cm.pi * float( (line.split())[2] )  #value of angular frequency in case of AC

        count = 0
        for line in lines[start+1:end:1]:
            validLine = line.split("#") #to avoid including comments
            lst = validLine[0].split()
            ckt_comp = Component(lst,omega)
            if ckt_comp.nodeA not in rows_mna:  #rows_mna stores all nodes in circuit
                rows_mna[ckt_comp.nodeA] = count
                count+=1
            if ckt_comp.nodeB not in rows_mna:
                rows_mna[ckt_comp.nodeB] = count
                count+=1
            if ckt_comp.type == "V":
                temp_str = branch_i(ckt_comp.nodeA,ckt_comp.nodeB)
                rows_mna[temp_str] = count  #rows_mna also stores branch currents of voltage sources
                count+=1            
            circuit_list.append(ckt_comp)
except IOError:     #if some input given incorrectly, eg. file name
    print('Invalid file')
    exit()

#  the MNA matrix is of size dim-1 X dim-1, which will be implemented later
dim = len(rows_mna)
mna = np.zeros((dim,dim),dtype = complex)
b = np.zeros(dim, dtype = complex)

for ckt_comp in circuit_list:
    type = ckt_comp.type ; val = ckt_comp.val
    rowA = rows_mna[ckt_comp.nodeA] ; rowB = rows_mna[ckt_comp.nodeB]
    #every component gets an MNA stamp, which is added to the MNA matrix
    #the stamp represents the contribution of that component to the nodal equations
    if type == "R" or type == "L" or type == "C":
        mna[rowA][rowA]+= 1/val
        mna[rowA][rowB]+= -1/val
        mna[rowB][rowA]+= -1/val
        mna[rowB][rowB]+= 1/val
    elif type == "I":
        b[rowA]+= val
        b[rowB]+= -val
    elif type == "V":
        row_i = rows_mna[branch_i(ckt_comp.nodeA,ckt_comp.nodeB)]
        mna[rowA][row_i]+= -1
        mna[rowB][row_i]+= 1
        mna[row_i][rowA]+= 1
        mna[row_i][rowB]+= -1
        b[row_i]+= val

#we delete GND to ensure the rows are linearly independent
mna = np.delete(mna,rows_mna["GND"],axis = 0)
mna = np.delete(mna,rows_mna["GND"],axis = 1)
b = np.delete(b,rows_mna["GND"])
found = False
for key in rows_mna:
    if key=="GND":
        found = True
    if found:   #to ensure there is no jump in node numbers
        rows_mna[key] -= 1

x = np.linalg.solve(mna,b)  #numpy function to solve system of linear equations
try:
    i = 1
    while True:
        temp_num = rows_mna[str(i)]
        z = x[temp_num]
        if cm.phase(z)==0 or cm.phase(-z)==0:   #rounding off real values to ensure floating point error doesnt occur
            print("V",i," : ",round(z.real,5))  # i.e., 2.4999999 instead of 2.5
        else:   
            print("V",i," : ",z)
        i+=1
except: pass    #to ensure list traversal is stopped when index goes out of bounds

for ckt_comp in circuit_list:
    if ckt_comp.type == "V":
        temp_string = branch_i(ckt_comp.nodeA,ckt_comp.nodeB)
        row_i = rows_mna[temp_string]
        z = x[row_i]
        if cm.phase(z)==0 or cm.phase(-z)==0:
            print(temp_string," : ",round(z.real,5))
        else:
            print(temp_string," : ",z)