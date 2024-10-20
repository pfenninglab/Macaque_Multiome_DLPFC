# import required modules
import os
import re
import pandas as pd
# assign directory
directory = '/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/code/raw_code/label_multiomeATAC_cells/logs'

celltypes = []
lifted_macaque = 0
total_macaque = []
lifted_human = 0
total_human = []
lifted_rat = 0
total_rat = []
# iterate over files in
# that directory
for filename in os.listdir(directory):
    #print(f"filename: {filename}")
    f = os.path.join(directory, filename)
    # checking if it is the correct file
    if filename.startswith("halliftover_312"):
        # open the sample file used
        file = open(f)
        #print(f"file: {file}")
        # read the content of the file opened
        content = file.readlines()
        #print(f"content: {content}")
        celltype = content[3]
        start = celltype.find('Alzheimers.') + 11
        end = celltype.find('.mm10', start)
        celltype = celltype[start:end]
        celltypes.append(celltype)
        #print(f"celltype: {celltype}")
        macaque = content[13]
        human = content[20]
        rat = content[27]
        macaque_peaks = re.findall('\d+', macaque)
        macaque_peaks = [eval(i) for i in macaque_peaks]
        #print(f"macaque_peaks: {macaque_peaks}")
        total_macaque.append(macaque_peaks[0]/macaque_peaks[1])
        print(f"Macaque peaks lifted: {macaque_peaks}")
        #print(f"Macaque peaks lifted: {macaque_peaks[0]/macaque_peaks[1]}")
        human_peaks = re.findall('\d+', human)
        human_peaks = [eval(i) for i in human_peaks]
        total_human.append(human_peaks[0]/human_peaks[1])
        #print(f"Human peaks lifted: {human_peaks}")
        #print(f"Mouse peaks lifted: {mouse_peaks[0]/mouse_peaks[1]}")
        rat_peaks = re.findall('\d+', rat)
        rat_peaks = [eval(i) for i in rat_peaks]
        total_rat.append(rat_peaks[0]/rat_peaks[1])
        #print(f"Rat peaks lifted: {rat_peaks}")
        #print(f"Rat peaks lifted: {rat_peaks[0]/rat_peaks[1]}")
print(f"Macaque peaks lifted: {total_macaque}")
print(f"Human peaks lifted: {total_human}")
print(f"Rat peaks lifted: {total_rat}")

df = pd.DataFrame(list(zip(celltypes, total_rat, total_macaque, total_human)),
               columns =['Celltype', 'Mapped to Rat', 'Mapped to Macaque', 'Mapped to Human'])
df = df.sort_values('Celltype')

df

df.to_csv("10X_Mapped.csv", index=False)
