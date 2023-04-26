from lxml import etree
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
pd.set_option('display.max_rows', None,
              'display.max_columns', None,
              'display.max_colwidth', None,
              'display.precision', 5)

# Load the XML file
tree = etree.parse("vasprun.xml")

# Get the root element
root = tree.getroot()

#Listing the atoms and creating list of lists with the species and number of atoms per species
# Find the array element with a name attribute of "atoms"
atoms_array = root.find(".//array[@name='atoms']")

# Get the set element
atoms_set = atoms_array.find("set")

# Loop over the rc elements and extract the data
list_atoms_raw = []
for rc in atoms_set.findall("rc"):
    element = rc.find("c[1]").text.strip()
    atomtype = int(rc.find("c[2]").text.strip())
    list_atoms_raw.append([element, atomtype])

# Group the atoms by species and count the number of atoms per species
species_counts = {}
for element, atomtype in list_atoms_raw:
    if element not in species_counts:
        species_counts[element] = 0
    species_counts[element] += 1

# Create a list of lists with the species and number of atoms per species
species_list = [[element, count] for element, count in species_counts.items()]

# Print the species list
print(species_list)

#Obtain the NEDOS value 
# Find the NEDOS tag
nedos_tag = tree.xpath('//i[@name="NEDOS"]')

# extract the NEDOS value
nedos_value = int(nedos_tag[0].text)

print("NEDOS= ",nedos_value)

#Is the system spin-polarized?
# Find the ISPIN tag
ispin_tag = tree.xpath('//i[@name="ISPIN"]')

# extract the ISPIN value
ispin = int(ispin_tag[0].text)

if ispin == 1:
    print("Non-spin polarized")
if ispin == 2:
    print("Spin polarized")

#Obtaining the Fermi Energy
# Find the efermi tag
efermi_tag = tree.xpath('//i[@name="efermi"]')

# extract the Efermi value
efermi = float(efermi_tag[0].text)
print("Fermi Energy: ", efermi, " eV")

#Finding the total DOS of the system
if ispin == 1:
        #Finding the total UP and DOWN DOS
        dos_total_tag = tree.xpath('//dos/total/array/set/set[@comment="spin 1"]')

        # Find the data element and extract the data
        data = []
        for r in dos_total_tag[0].iter('r'):
            data.append([float(x) for x in r.text.split()])

        # Convert the data into a DataFrame
        columns = ['Energy / eV', 'UP', 'integrated']
        df_total = pd.DataFrame(data, columns=columns).set_index('Energy / eV')
        df_total = df_total.drop('integrated',axis=1)

if ispin == 2:
    #Finding the total UP and DOWN DOS
    dos_total_up_tag = tree.xpath('//dos/total/array/set/set[@comment="spin 1"]')
    dos_total_dn_tag = tree.xpath('//dos/total/array/set/set[@comment="spin 2"]')

    # Find the data element and extract the data
    data_up = []
    for r in dos_total_up_tag[0].iter('r'):
        data_up.append([float(x) for x in r.text.split()])
    data_dn = []
    for r in dos_total_dn_tag[0].iter('r'):
        data_dn.append([float(x) for x in r.text.split()])

    # Convert the data into a DataFrame
    columns_up = ['Energy / eV', 'UP', 'integrated']
    columns_dn = ['Energy / eV', 'DOWN', 'integrated']
    df_total_up = pd.DataFrame(data_up, columns=columns_up).set_index('Energy / eV')
    df_total_dn = pd.DataFrame(data_dn, columns=columns_dn).set_index('Energy / eV')
    df_total_up = df_total_up.drop('integrated',axis=1)
    df_total_dn = df_total_dn.drop('integrated',axis=1)
    df_total = pd.merge(df_total_up,df_total_dn, on='Energy / eV')

    #Making the down DOS negative for plotting porpouses
    zero = df_total*0.0
    df_total["DOWN"] = zero["DOWN"].sub(df_total["DOWN"])

#Refering to the Fermi Energy
df_total.index = df_total.index - efermi

#Saving in a new file the Total DOS
with open("TDOS.txt", "w") as f:
    string = df_total.to_string(header =True, index=True)
    f.write(f"{string}\n")

#Finding the projected DOS per ion in the system
#Creating an empty list of DataFrames to store the UP and DOWN (in case ispin=2) PDOS for each atom
up_pdos = []
dn_pdos = []
for i in range(1,len(list_atoms_raw)+1):
    #Find the projected DOS per ion ´at´ in the list_atoms_raw
    path_up  = '//dos/partial/array/set/set[@comment="ion '+str(i)+'"]/set[@comment="spin 1"]'
    dos_pdos_up_tag = tree.xpath(path_up)

    # Find the data element and extract the data
    data = []
    for r in dos_pdos_up_tag[0].iter('r'):
        data.append([float(x) for x in r.text.split()])
    columns = ['Energy / eV','s','py','pz','px','dxy','dyz','dz2','dxz','d(x2-y2)']
    df_pdos_at_up = pd.DataFrame(data, columns=columns).set_index('Energy / eV')
    df_pdos_at_up.index = df_pdos_at_up.index - efermi
    up_pdos.append(df_pdos_at_up)

    if ispin == 2:
        #DOWN: Find the projected DOS per ion ´at´ in the list_atoms_raw
        path_dn  = '//dos/partial/array/set/set[@comment="ion '+str(i)+'"]/set[@comment="spin 2"]'
        dos_pdos_dn_tag = tree.xpath(path_dn)

        # DOWN: Find the data element and extract the data
        data_dn = []
        for r in dos_pdos_dn_tag[0].iter('r'):
            data_dn.append([float(x) for x in r.text.split()])
        df_pdos_at_dn = pd.DataFrame(data_dn, columns=columns).set_index('Energy / eV')
        df_pdos_at_dn.index = df_pdos_at_dn.index - efermi
        #Making the down DOS negative for plotting porpouses
        zero = df_pdos_at_dn*0.0
        df_pdos_at_dn = zero.sub(df_pdos_at_dn)
        dn_pdos.append(df_pdos_at_dn)

#From lm-decomposed to l-decomposed
up_pdos_l = []
dn_pdos_l = []
for df_lm in up_pdos:
    df_l = pd.DataFrame({'Energy / eV':[], 's':[],'p':[],'d':[]})
    df_l['Energy / eV'] = df_pdos_at_up.index
    df_l = df_l.set_index('Energy / eV')
    df_l['s'] = df_lm['s']
    df_l['p'] = df_lm['py'].add(df_lm['pz']).add(df_lm['px'])
    df_l['d'] = df_lm['dxy'].add(df_lm['dyz']).add(df_lm['dz2']).add(df_lm['dxz']).add(df_lm['d(x2-y2)']) 
    up_pdos_l.append(df_l)

if ispin == 2:
    for df_lm in dn_pdos:
        df_l = pd.DataFrame({'Energy / eV':[], 's':[],'p':[],'d':[]})
        df_l['Energy / eV'] = df_pdos_at_dn.index
        df_l = df_l.set_index('Energy / eV')
        df_l['s'] = df_lm['s']
        df_l['p'] = df_lm['py'].add(df_lm['pz']).add(df_lm['px'])
        df_l['d'] = df_lm['dxy'].add(df_lm['dyz']).add(df_lm['dz2']).add(df_lm['dxz']).add(df_lm['d(x2-y2)']) 
        dn_pdos_l.append(df_l)

#Saving the l-decomposed PDOS per atom
for i in range(1,len(up_pdos_l)+1):
    name_file = 'dos'+str(i)+'_l-decomposed.txt'
    with open(name_file,"w") as f:
        for item in up_pdos_l:
            string = item.to_string(header =True, index=True)
            f.write(f"{string}\n")
if ispin == 2:
    for i in range(1,len(dn_pdos_l)+1):
        name_file = 'dos'+str(i)+'_l-decomposed_down.txt'
        with open(name_file, "w") as f:
            for item in dn_pdos_l:
                string = item.to_string(header =True, index=True)
                f.write(f"{string}\n")

#Grouping DOS by species
#Creating a list with DataFrame = 0
dos_specie_up = []

#Looping over the list of species
for i in range(len(species_list)):
    dos_specie_up.append(up_pdos_l[0]*0)

#Looping over the number of species
k = 0
for i in range(len(species_list)):
    #Looping over the number of atoms of each specie
    for j in range(species_list[i][1]):
        #Adding the DOS to the auxiliar DOS tuple of DataFrame
        dos_specie_up[i] = dos_specie_up[i].add(up_pdos_l[k])
        k = k + 1

if ispin == 2:
    dos_specie_dn = []

    #Looping over the list of species
    for i in range(len(species_list)):
        dos_specie_dn.append(dn_pdos_l[0]*0)

    #Looping over the number of species
    k = 0
    for i in range(len(species_list)):
        #Looping over the number of atoms of each specie
        for j in range(species_list[i][1]):
            #Adding the DOS to the auxiliar DOS tuple of DataFrame
            dos_specie_dn[i] = dos_specie_dn[i].add(dn_pdos_l[k])
            k = k + 1

#Saving PDOS per species in files of UP and DOWN (in case ispin=2)
for i in range(len(dos_specie_up)):
    name_file = 'dos_'+str(species_list[i][0])+'_l-decomposed.txt'
    with open(name_file, "w") as f:
        for item in dos_specie_up:
            string = item.to_string(header =True, index=True)
            f.write(f"{string}\n")

if ispin == 2:
    for i in range(len(dos_specie_dn)):
        name_file = 'dos_'+str(species_list[i][0])+'_l-decomposed_down.txt'
        with open(name_file, "w") as f:
            for item in dos_specie_dn:
                string = item.to_string(header =True, index=True)
                f.write(f"{string}\n")

#Plotting Style
if ispin == 1:
    custom_cycler= (cycler(color=["black","#0099cc","crimson","limegreen"]) +
                    #cycler(marker=["o","o","o","","","","","",""]) +
                    #cycler(markersize=[8,8,8,0,0,0,0,0,0]) +
                    cycler(linestyle=["-","-","-","-"]) +
                    cycler(linewidth=[1,1,1,1]))
if ispin == 2:
    custom_cycler= (cycler(color=["black","black","#0099cc","#0099cc","crimson","crimson","limegreen","limegreen"]) +
                    #cycler(marker=["o","o","o","","","","","",""]) +
                    #cycler(markersize=[8,8,8,0,0,0,0,0,0]) +
                    cycler(linestyle=["-","-","-","-","-","-","-","-"]) +
                    cycler(linewidth=[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]))
font= {"family": "normal",
       "size"  : 14}
plt.rc("font",**font)
plt.rcParams["font.family"]="sans-serif"
plt.rcParams["font.sans-serif"]="Times New Roman"
plt.rcParams["axes.prop_cycle"] = custom_cycler

#Plotting
# Total DOS
plot = df_total.plot(legend="Total DOS")
plot.axvline(x=0, color="black", linestyle="--")

plt.xlim(-2,2)

# PDOS
if ispin == 1:
    plot1 = dos_specie_up[0]["d"].plot(ax=plot,label="Ti (d)")
    plot1 = dos_specie_up[0]["p"].plot(ax=plot,label="Ti (p)")
    plot1 = dos_specie_up[1]["p"].plot(ax=plot,label="C (p)")
    plt.ylim(0,10)
if ispin == 2:
    plot1 = dos_specie_up[0]["d"].plot(ax=plot,label="Ti (d)")
    plot1 = dos_specie_dn[0]["d"].plot(ax=plot,legend=None)
    plot1 = dos_specie_up[0]["p"].plot(ax=plot,label="Ti (p)")
    plot1 = dos_specie_dn[0]["p"].plot(ax=plot,legend=None)
    plot1 = dos_specie_up[1]["p"].plot(ax=plot,label="C (p)")
    plot1 = dos_specie_dn[1]["p"].plot(ax=plot,legend=None)
    plt.ylim(-15,15)

#Legend
handles, labels = plot.get_legend_handles_labels()
if ispin == 1:
    new_handles = [handles[0],handles[1],handles[2],handles[3]]
    new_labels = ["Total DOS","Ti (d)", "Ti (p)", "C (p)"]
    legend = plot.legend(handles=new_handles,labels=new_labels,loc="center left",bbox_to_anchor=[1.01, 0.5],frameon=False)
if ispin == 2:
    new_handles = [handles[0],handles[2],handles[4],handles[6]]
    new_labels = ["Total DOS","Ti (d)", "Ti (p)", "C (p)"]
    legend = plot.legend(handles=new_handles,labels=new_labels,loc="center left",bbox_to_anchor=[1.01, 0.5],frameon=False)
#                      handler_map={str: LegendTitle({'fontsize': 28})})

plt.xlabel("$\mathdefault{∆E(E-E_{f}}$) (eV)")
plt.ylabel("DOS (states/eV)")

#Saving the plot
name = "DOS-Ti3C2-afm1.png"
plt.savefig(fname=name,dpi=1000,bbox_inches="tight")
