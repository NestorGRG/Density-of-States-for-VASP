# Density-of-States-for-VASP
Python program to plot the Density of States (DOS) from VASP outputs

### The projection of the DOS must be calculated, in consequence, the `LORBIT =  11` keyword must be present in the `INCAR` file.

## Needed Files
The necessary input files are:
  1.  `vasprun.xml`

## Requirements  
The following Python modules must be installed:
  1. `Pandas`
  2. `lxml`
  3. `Cycler`
  2. `Matplotlib`

## Instructions
This program must be runned in the same directory as the input files.
The use of this program is very simple:
1.  Make sure that the input files have the aforementioned names.
2.  The above mentioned python modules are installed.
3.  In the Plotting section, the user is able to select the style of the plot. Including: color, linewidth, marker, font, etc.
4.  In the same section, the user must indicate which atom species and which orbital angular momentum (*l* number: s,p,d) is selected for the plot. The `dos_specie_up` variable  contains the DOS for each atom specie and orbital angular momentum for non-spin polarized calculations and for the spin-up channel. The first index states the atom specie following the sorting of the screen-displayed list of atom species and number of atoms of the same specie. e.g. index = 0: Ti, index = 1: C. The last index of the `dos_specie_up` variable corresponds to the orbital angular momentum, being index = 1: s, index = 2: p, index = 3: d, etc. In the case of spin-polarized system, the same must be done for the `dos_specie_dn` variable.
5.  The user must change the labels of the leyend and add or remove them according to the properties of the desired system.
6.  The user can change the name of the generated png file by changing the variable `name`.
