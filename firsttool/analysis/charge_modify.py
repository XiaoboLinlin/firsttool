
### find dictionary for electrode
fin = open("emim_data", "r")
fin.close
linelist = fin.readlines()
# find the line number for Bonds, which is 2 lines before the last line of 'Atoms'
lookup = 'Bonds'
with open('emim_data', "r") as myFile:
    for num, line in enumerate(myFile, 1):
        if lookup in line:
            bonds_num = num
# elec_R or elec_L is the mol number of right and left electrode, [-3] is correct.
elec_R = linelist[bonds_num-3].split()[1]
elec_R = int(elec_R)
elec_L = elec_R - 1 
print('elec_L is ', elec_L)
print('elec_R is ', elec_R)
d_data = dict()
start = False
for line in linelist:
    if start:
        if len(line.split()) >= 2 and (int(line.split()[1]) == elec_L or int(line.split()[1]) == elec_R): 
            d_data[line.split()[0]] = line.split()[3]
    try:
        if line.split()[0] == 'Atoms':
            start = True
        if line.split()[0] == 'Bonds':
            break
    except:
        continue
###

### making new file for modified charge
print('making new file for modified charge')
# print(d_data)
fin = open("customize_conp.lammpstrj", "r")
fout = open("customize_conp_modified.lammpstrj", "w")
fin.close
linelist = fin.readlines()
for line in linelist:
    if len(line.split()) >= 7 and line.split()[0] != 'ITEM:':
        if int(line.split()[1]) == elec_L or int(line.split()[1]) == elec_R:
#             print(float(line.split()[-1]))
#             print(line.split())
#             print(d_data[line.split()[0]])
            new_charge = float(line.split()[-1]) - float(d_data[line.split()[0]])
            line = line.rsplit(' ', 1)[0] + ' ' + str(round(new_charge,5))+'\n'
    fout.write(line)    
fout.close()
print('yes, done')