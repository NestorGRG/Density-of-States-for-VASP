import numpy as np
import matplotlib.pyplot as plt

#Reading firsts lines of DOSCAR file
with open("DOSCAR", "r") as doscar:
    lines = doscar.readlines()[0:6]
    index = 0
    natoms = int(lines[index].strip().split()[0])
    index = 5
    nedos = int(lines[index].strip().split()[2])
    efermi = float(lines[index].strip().split()[3])

#Reading POSCAR
with open("POSCAR") as poscar:
    lines = poscar.readlines()[5:7]
    index = 0
    types_atoms= lines[0].strip().split()
    num_per_atoms=lines[1].strip().split()
    num_per_atoms = np.array(num_per_atoms).astype(int)

#Reading whether the system is non-spin polarized or spin-polarized
with open("vasprun.xml","r") as vasprun:
    for count, line in enumerate(vasprun):
        if "ISPIN" in line:
            ispin = line.split()

if ispin[-1] == "1</i>": print("non-spin polarized")
else: print("spin polarized")

#Array for the number of atoms in the system
natomsindex =[]
for j in range(1,natoms+1):
    natomsindex.append(j)

#Creating a dictionary of type of atoms and number of them in the system
at=[]
for i in range(len(num_per_atoms)):
    at.append([types_atoms[i],num_per_atoms[i]])
print(at)
print("E.Fermi= ",efermi)

#Finding the total DOS
a = []
with open("vasprun.xml","r") as vasprun:
    for count, line in enumerate(vasprun):
        if "integrated" in line:
            start0 = count
            end0 = start0+nedos
        pass
with open("vasprun.xml","r") as vasprun:
        for count, line in enumerate(vasprun):
            if count in range(start0+3, end0+3): 
                a.append(line.split())
            pass    
if ispin[-1] != "1</i>":
    a_down=[]
    with open("vasprun.xml","r") as vasprun:
        for count, line in enumerate(vasprun):
            if count in range(start0+nedos+5, end0+nedos+5): 
                a_down.append(line.split())
            pass        

#Finding the line number where "ion j" appears and creating an array to store the lm-decomposed DOS
start1=0
end=0
b= []
b_down = []

for j in natomsindex:
    ion = "ion "+str(j)
    
    with open("vasprun.xml","r") as vasprun:
        for count, line in enumerate(vasprun):
            #print(count)
            if ion in line:
                start1 = count
                end1 = start1+nedos
            pass

    with open("vasprun.xml","r") as vasprun:
        for count, line in enumerate(vasprun):
            if count in range(start1+2, end1+2): 
                b.append(line.split())
            pass
    if ispin[-1] != "1</i>":
        with open("vasprun.xml","r") as vasprun:
            for count, line in enumerate(vasprun):
                if count in range(start1+nedos+4, end1+nedos+4): 
                    b_down.append(line.split())
                pass

#Converting into an numpy array and deleting <r> and </r> columns and refering the energy to the E.Fermi level
#Total DOS lm-decomposed
dost = np.array(a)
dost1 = np.delete(dost,4,1)
dost2 = np.delete(dost1,0,1)
dost3 = dost2.astype(float)

if ispin[-1] != "1</i>":
    dost_down = np.array(a_down)
    dost1_down = np.delete(dost_down,4,1)
    dost2_down = np.delete(dost1_down,0,1)
    dost3_down = dost2_down.astype(float)

#ions DOS lm-decomposed
dos=np.array(b)
dos1 = np.delete(dos,11,1)
dos2 = np.delete(dos1,0,1)
dos3 = dos2.astype(float)

if ispin[-1] != "1</i>":
    dos_down=np.array(b_down)
    dos1_down = np.delete(dos_down,11,1)
    dos2_down = np.delete(dos1_down,0,1)
    dos3_down = dos2_down.astype(float)

#refering to the E.Fermi level
dost3[:,0] = dost3[:,0] - efermi
dos3[:,0] = dos3[:,0] - efermi

if ispin[-1] != "1</i>":
    dost3_down[:,0] = dost3_down[:,0] - efermi
    dos3_down[:,0] = dos3_down[:,0] - efermi
    
    #Making the down spin states negative
    dost3_down[:,1] = 0.0 - dost3_down[:,1]
    for i in range(1,10):
        dos3_down[:,i] = 0.0 - dos3_down[:,i]

#Creating files with the lm-decomposed DOS for each atom and total DOS
#total dos lm-decomposed
np.savetxt("dos0.txt", dost3, fmt="%8.4f", delimiter=" ")
if ispin[-1] != "1</i>":
    np.savetxt("dos0_down.txt", dost3_down, fmt="%+8.4f", delimiter=" ")

#Saving ions DOS lm-decomposed
j=0
for i in natomsindex:
    name= "dos"+str(i)+".txt"
    np.savetxt(name, dos3[0+nedos*j:nedos*i,:], fmt="%+8.4f", delimiter=" ")

    if ispin[-1] != "1</i>":
        name2= "dos"+str(i)+"_down.txt"
        np.savetxt(name2, dos3_down[0+nedos*j:nedos*i,:], fmt="%+8.4f", delimiter=" ")
    j = j + 1
    
#Loading the dos-CAR files (dos"i".txt)
#Loading total lm-decomposed dos
dos_total = np.loadtxt("dos0.txt")

#Loading ions DOS lm-decomposed
dos_aux= np.zeros([natoms,nedos,10])
for i in range(1,natoms+1):
    name="dos"+str(i)+".txt"
    dos_aux[i-1] = np.loadtxt(name)

if ispin[-1] != "1</i>":
    dos_total_down = np.loadtxt("dos0_down.txt")
    dos_aux_down= np.zeros([natoms,nedos,10])
    for i in range(1,natoms+1):
        name="dos"+str(i)+"_down.txt"
        dos_aux_down[i-1] = np.loadtxt(name)


#Adding DOS of equal species
dos_aux2= np.zeros([len(at),nedos,10])

k = 0
i = 0
j = 0
for i in range(len(at)):
    for j in range(at[i][1]):
        dos_aux2[i]= dos_aux[k]+dos_aux2[i]
        k = k + 1
dos_aux2[:,:,0]=dos_aux[0,:,0]

if ispin[-1] != "1</i>":
    dos_aux2_down= np.zeros([len(at),nedos,10])
    k = 0
    i = 0
    j = 0
    for i in range(len(at)):
        for j in range(at[i][1]):
            dos_aux2_down[i]= dos_aux_down[k]+dos_aux2_down[i]
            k = k + 1    
    dos_aux2_down[:,:,0]=dos_aux_down[0,:,0]
    
#Saving the total lm decomposed for each species
i = 0
for i in range(len(at)):
    name3="dos_"+str(at[i][0])+"_lm-decomposed.txt"
    np.savetxt(name3,dos_aux2[i], fmt="%+8.4f", delimiter=" ")

if ispin[-1] != "1</i>":
    i = 0
    for i in range(len(at)):
        name3="dos_"+str(at[i][0])+"_lm-decomposed_down.txt"
        np.savetxt(name3,dos_aux2_down[i], fmt="%+8.4f", delimiter=" ")
        

#Adding DOS of equal l number for each species. From lm-decomposed to l-decomposed.
#Creating an array of l-decomposed DOS per species
dosls = np.zeros([len(at),nedos,4])
dosls[:,:,0] = dos_aux2[0,:,0]

if ispin[-1] != "1</i>":
    dosls_down = np.zeros([len(at),nedos,4])
    dosls_down[:,:,0] = dos_aux2_down[0,:,0]
    
#S orbitals
for i in range(len(at)):
    dosls[i,:,1] = dos_aux2[i,:,1]

if ispin[-1] != "1</i>":
    for i in range(len(at)):
        dosls_down[i,:,1] = dos_aux2_down[i,:,1]
    
#P orbitals
for i in range(len(at)):
    dosls[i,:,2] = dos_aux2[i,:,2] + dos_aux2[i,:,3] + dos_aux2[i,:,4]

if ispin[-1] != "1</i>":
    for i in range(len(at)):
        dosls_down[i,:,2] = dos_aux2_down[i,:,2] + dos_aux2_down[i,:,3] + dos_aux2_down[i,:,4]

#D orbitals
for i in range(len(at)):
    dosls[i,:,3] = dos_aux2[i,:,5] + dos_aux2[i,:,6] + dos_aux2[i,:,7] + dos_aux2[i,:,8] + dos_aux2[i,:,9]

if ispin[-1] != "1</i>":
    for i in range(len(at)):
        dosls_down[i,:,3] = dos_aux2_down[i,:,5] + dos_aux2_down[i,:,6] + dos_aux2_down[i,:,7] + dos_aux2_down[i,:,8] + dos_aux2_down[i,:,9] 

#Saving all the total PDOS per species and l- number in new files
for i in range(len(at)):
    name4 = str(at[i][0])+ "_total_dos.txt"
    np.savetxt (name4, dosls[i], fmt="%+8.4f", delimiter=" ")

if ispin[-1] != "1</i>":
    for i in range(len(at)):
        name4 = str(at[i][0])+ "_total_dos_down.txt"
        np.savetxt (name4, dosls_down[i], fmt="%+8.4f", delimiter=" ")
    
#Fitting to Legendre polynomial
#legendre= np.polynomial.legendre.Legendre.fit(dos_total[:,0], dos_total[:,1],deg=600)
#plt.plot(dos_total[:,0],legendre(dos_total[:,0]), label="total DOS", color="green")

#Plotting
plt.plot(dos_total[:,0], dos_total[:,1], label="Total DOS", color="black",linewidth=1.0)
plt.plot(dos_total[:,0], dosls[0,:,3], label="Ti (d)", color="#0099cc",linewidth=1.0)
plt.plot(dos_total[:,0], dosls[0,:,2], label="Ti (p)", color="crimson",linewidth=1.0)
plt.plot(dos_total[:,0], dosls[1,:,2], label="C (p)", color="limegreen",linewidth=1.0)

plt.xlim(-2,2)
plt.ylim(0,max(dos_total[:,1]))

#To view the Ef as a reference
plt.axvline(x=0, color="black", linestyle="--")

if ispin[-1] != "1</i>":
    plt.plot(dos_total[:,0], dos_total_down[:,1], color="black",linewidth=1.0)
    plt.plot(dos_total[:,0], dosls_down[0,:,3], color="#0099cc",linewidth=1.0)
    plt.plot(dos_total[:,0], dosls_down[0,:,2], color="crimson",linewidth=1.0)
    plt.plot(dos_total[:,0], dosls_down[1,:,2], color="limegreen",linewidth=1.0)

    plt.ylim(-15,15)

font= {"family": "normal"}
plt.rc("font",**font)
plt.rcParams["font.family"]="sans-serif"
plt.rcParams["font.sans-serif"]="Times New Roman"

plt.xlabel("∆E(E-Ef) (eV)")
plt.ylabel("DOS (states/eV)")
plt.legend(loc="center left", bbox_to_anchor=[1.01, 0.5])

#Saving the plot
name = "DOS"
plt.savefig(fname=name,dpi=1000,bbox_inches="tight")
plt.show()