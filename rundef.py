import os
import subprocess
mesadir   = '/Users/danielegana/Physicsapps/mesa-r23.05.1'

def runstar(x,y):
     Zdir=mesadir+'/Z'+str(x)+'M'+str(y)
     os.chdir(Zdir)
     output = subprocess.check_output('sh '+Zdir+'/mk && '+'sh '+Zdir+'/rn' , shell=True, text=True)
     return output
