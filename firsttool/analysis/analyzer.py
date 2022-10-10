import mbuild as mb
import numpy as np
import parmed as pmd
from foyer import Forcefield
import MDAnalysis as mda


class analysis():
    def file_reader(self, file):
        with open(file,"r") as fi:
            lines = []
            copy = False
            for ln in fi:
                if ln.startswith("Step"):
                    copy = True
                    continue
                if ln.startswith("Loop") or ln.startswith("Lost"):
                    copy = False
                if ln.startswith("if \"$r == 0\" then \"log log.2\""):
                    copy = False
                if ln.endswith("if \"$r == 0\" then \"log log.2\""):
                    copy = False
                # if ln.startswith(""):
                #     copy = False
                if ln.startswith("WARN"):
                    continue
                if copy:
                    lines.append(ln)
        data = np.array(lines)  
        df = np.loadtxt(data)
        return df
    def profile_reader(self, file):
        with open(file,"r") as fi:
            lines = []
            for ln in fi:
                if ln.startswith(" "):
                    print('yes')
                    lines.append(ln)
        data = np.array(lines)  
        df = np.loadtxt(data)
        return df

# data = file_reader("log.2")
# plt.plot(data[:,0]*0.001, data[:,2], c= 'b', label = 'y+z+')
# plt.plot(data[:,0]*0.001, data[:,3], label = 'y-z+')
# plt.plot(data[:,0]*0.001, data[:,4], label = 'y+z-')
# plt.plot(data[:,0]*0.001, data[:,5], c= 'r',label = 'y-z-')
# plt.xlabel('time (ps)')
# plt.ylabel('q (e)')
# plt.legend()
# plt.ylim(-13,13)