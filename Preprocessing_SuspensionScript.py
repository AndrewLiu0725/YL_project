# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 02/05/2021
# ===============================================================================
import os
import subprocess

'''
This program is a script for preprocessing the data for suspension system
'''

path = "/raid6/ctliao/Data/HI_ordering/"

def SendEmail(progress):
    f = open('email.txt','w')
    f.write('Subject: Current progress for preprocessing: {}%\n'.format(progress))
    f.close()
    command_line = 'sendmail e5n41kb9n@gmail.com  < email.txt'
    subprocess.Popen(command_line, shell = True)
    subprocess.run(["sendmail", "e5n41kb9n@gmail.com", "<", "email.txt"])

# parser
def parser(string):
    # output [phi, Ca, ensemble id]
    string = string.split("-")
    ensemble_id = int(string[1])
    # format: h24_phi4.4989_Re0.1_Ca0.06_WCA1_zero0.8-8
    if "_" in string[0]:
        string = string[0].split("_")
        phi = float(string[1][3:])
        Ca = float(string[3][2:])
    # format: h24phi4.9488Re0.1Ca0.06WCA1zero0.8-4
    else:
        phi = float(string[0][string[0].find('i')+1:string[0].find('R')])
        Ca = float(string[0][string[0].find('a')+1:string[0].find('W')])
    return [phi, Ca, ensemble_id]


parameter_set = {}
# create parameter_set (two layer dict)
for fn in os.listdir(path):
    if (fn[0] == "h") and (os.path.isdir(path+fn)):
        result = parser(fn)
        [phi, Ca, ensemble_id] = result
        if phi in parameter_set.keys():
            if Ca in parameter_set[phi].keys():
                parameter_set[phi][Ca].append(ensemble_id)
            else:
                parameter_set[phi][Ca] = [ensemble_id]
        else:
            parameter_set[phi] = {}
            parameter_set[phi][Ca] = [ensemble_id]


phis = []

number_of_parameter_set = 0
for phi in parameter_set.keys():
    if (len(parameter_set[phi].keys()) > 5):
        for Ca in parameter_set[phi].keys():
            number_of_parameter_set += len(parameter_set[phi][Ca])
        phis.append(phi)

# Preprocessing
count = 0
for fn in os.listdir(path):
    if (fn[0] == "h") and (os.path.isdir(path+fn)):
        result = parser(fn)
        if result[0] in phis:
            subprocess.call(["python3", "Preprocessing.py", fn, "system=1", "verbose=0"])
            count += 1
            if (count % int(number_of_parameter_set/10) == 0):
                print("Current progress is {}%".format(int(count*100/number_of_parameter_set)))
                #SendEmail(int(count*100/number_of_parameter_set))

#SendEmail(int(count*100/number_of_parameter_set))