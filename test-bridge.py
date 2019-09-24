import automate
import glob
import os
import shutil

site='bridge'
angles=[30.0,50.0,70.0]
dirs=[]

for angle in angles:
    dirname="PbPc-"+site+"-"+str(angle)
    dirs.append(dirname)
    automate.createPOSCAR('CONTCAR', [7,7,4], 15.0, angle, site, 4.6)

    os.makedirs(dirname,exist_ok=True)

    templete=os.listdir('./templete-vasp-imr')
    for file_name in templete:
        full_file_name = os.path.join('./templete-vasp-imr', file_name)
        if (os.path.isfile(full_file_name)):
            target_file_name=os.path.join(dirname, file_name)
            shutil.copy(full_file_name, target_file_name)

    target_file_name=os.path.join(dirname, 'POSCAR')
    shutil.copy('POSCAR',target_file_name)
    vaspfiles=os.listdir(dirname)

automate.sendDirsToIMR(dirs,vaspfiles)


