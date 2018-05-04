#!/usr/bin/env python

#save as TargetOrtho_App.command file, sudo chmod +x TargetOrtho_App.command
#this make this script executable by double clicking it.

#need to add output directory name after job start, do this by generating job id and output dir name within this script then passing it to tartetortho.py
#add output dir path after job started

from tkinter import *
from tkinter.ttk import *
from tkinter import filedialog
import subprocess
import os
import datetime as dt

def system_call(command):
    p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
    return p.stdout.read()

window = Tk()
window.title("TargetOrtho 2.0")
window.geometry('1000x500')


#use first option for non-ipynb 
sys_path=os.path.dirname(os.path.abspath(__file__))
#sys_path=os.path.abspath(".")
now = dt.datetime.now()

jobID="j%s%s%s%s%s%s" %(now.year,now.month,now.day,now.hour,now.minute,now.second) 
#jobID="j20185116512"
output_dir="%s/%s_TargetOrtho2.0_Results" %(sys_path,jobID)

def clicked2():
    global file
    file=filedialog.askopenfilename()
    lbl2.configure(text= file)
    
def clicked():
    p_value =  txt1.get()
    distance=txt3.get()
    ref_species=selected.get()
    if ref_species==1:
        ref_species="C.elegans"
    if ref_species==2:
        ref_species="P.pacificus"
    if distance=='None':
        command="python \'%s/targetortho.py' -f \'%s\' -p %s -r %s -j %s" %(sys_path,file,p_value,ref_species,jobID)
    else:
        command="python '%s/targetortho.py' -f \'%s\' -p %s -d %s -r %s -j %s" %(sys_path,file,p_value,distance,ref_species,jobID)
    print(command)
    system_call(command)
    lbl5 = Label(window, text="results at: %s" %(output_dir),font=("Arial Bold", 12))
    lbl5.grid(column=0, row=25,sticky=W)
    system_call("open \'%s\'" %output_dir)

    
    
lbl2 = Label(window, text="TargetOrtho 2.0",font=("Arial Bold", 18))
lbl2.grid(column=0, row=0,sticky=W)

#############row
btn2 = Button(window, text="select motif file", command=clicked2)
btn2.grid(column=0, row=1,sticky=W)

lbl2 = Label(window, text="")
lbl2.grid(column=0, row=2,sticky=W)





##############row 
#get P value
lbl = Label(window, text="p value threshold for FIMO")
lbl.grid(column=0, row=3,sticky=W)

#get p value input from text box
x = StringVar(window, value='0.0001')
txt1 = Entry(window,width=10,textvariable=x)
txt1.pack()
txt1.grid(column=0, row=4,sticky=W)



#############row 
#get distance

lbl = Label(window, text="max distance from gene start (default is full intergenic distance)")
lbl.grid(column=0, row=5,sticky=W)

v = StringVar(window, value='None')
txt3 = Entry(window,width=10,textvariable=v)
txt3.pack()
txt3.grid(column=0, row=6,sticky=W)

##############row 
lbl = Label(window, text="Reference genome",font=("Arial Bold", 14))
lbl.grid(column=0, row=7,sticky=W)

selected = IntVar()

rad1 = Radiobutton(window,text='Caenorhabditis elegans', value=1, variable=selected)
rad2 = Radiobutton(window,text='Pristionchus pacficus', value=2, variable=selected)
selected.set(1)
rad1.grid(column=0, row=8,sticky=W)
rad2.grid(column=0, row=9,sticky=W)



###############row 
lbl = Label(window, text="Genomes",font=("Arial Bold", 14))
lbl.grid(column=0, row=10,sticky=W)


chk_state = BooleanVar()
chk_state.set(True) #set check state

chk_state2 = BooleanVar()
chk_state2.set(True) #set check state

i=11

speciesList=["Caenorhabditis elegans","Caenorhabditis briggsae","Caenorhabditis brenneri","Caenorhabditis remanei",
            "Caenorhabditis japonica","Pristionchus pacificus","Pristionchus exspectatus","Ascaris lumbricoides"]
for species in speciesList:
    lbl = Label(window, text=species)
    lbl.grid(column=0, row=i,sticky=W)
    i+=1


    
##################row 
    
btn = Button(window, text="Start job", command=clicked)
btn.grid(column=0, row=i+9,sticky=W)

window.mainloop()



# In[ ]:




# In[ ]:



