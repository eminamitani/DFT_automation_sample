from ase.build import fcc111
from ase import Atoms
from ase.io import read, write
from ase.constraints import FixAtoms
import paramiko
import os
from getpass import getpass
import math


def sendToIMR(files,dir):
    print("type your passwd for IMR supercomputer node")
    passwd = getpass()
    config_file = os.path.join(os.getenv('HOME'), '.ssh/config')
    ssh_config = paramiko.SSHConfig()
    ssh_config.parse(open(config_file, 'r'))
    lkup = ssh_config.lookup('super.imr')

    # ProxyCommand setup
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.load_system_host_keys()
    print("what found in config file")
    print(lkup)

    ssh.connect(
        lkup['hostname'],
        username=lkup['user'],
        sock=paramiko.ProxyCommand(lkup['proxycommand']),
        password=passwd
    )

    mkdircommand='mkdir '+dir
    stdin, stdout, stderr = ssh.exec_command(mkdircommand)

    #TODO catch stderr

    sftp = ssh.open_sftp()

    for f in files:
        target=dir+"/"+f
        sftp.put(f, target)

    sftp.close()
    ssh.close()

#send several directlies at one time, avoid typing passphrase so many times...
def sendDirsToIMR(dirs, vaspfiles):

    print("type your passwd for IMR supercomputer node")
    passwd = getpass()
    config_file = os.path.join(os.getenv('HOME'), '.ssh/config')
    ssh_config = paramiko.SSHConfig()
    ssh_config.parse(open(config_file, 'r'))
    lkup = ssh_config.lookup('super.imr')

    # ProxyCommand setup
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.load_system_host_keys()
    print("what found in config file")
    print(lkup)

    ssh.connect(
        lkup['hostname'],
        username=lkup['user'],
        sock=paramiko.ProxyCommand(lkup['proxycommand']),
        password=passwd
    )

    sftp = ssh.open_sftp()
    workingdir='/home/emi0716/work/'
    for i in dirs:
        mkdircommand='mkdir '+ workingdir+i
        stdin, stdout, stderr = ssh.exec_command(mkdircommand)
        outlines = stdout.readlines()
        result = ''.join(outlines)
        print(result)

    for i in dirs:

        #TODO catch stderr
        for f in vaspfiles:
            origin="./"+i+"/"+f
            print("origin:"+origin)
            #print(os.path.exists(origin))
            target=workingdir+i+"/"+f
            print("target:"+target)
            sftp.put(origin, target)

    sftp.close()
    ssh.close()


#Tb complex on Pb111 surface
def createPOSCAR(molecule, slabsize, vacuum, angle, site, dist):
    '''

    :param molecule: VASP POSCAR/CONTCAR format for molecule adsorbed on surface
    :param slabsize: [x,y,z] information to create slab
    :param angle: angle of molecular rotation in x-y plane
    :param site: ontop or bridge or hollow
    :return: None, create POSCAR file
    '''
    latt=4.97

    mol=read(molecule,format='vasp')
    mol.center()
    slab = fcc111('Pb', a=latt, size=(slabsize[0], slabsize[1], slabsize[2]), vacuum=vacuum)
    topLayers = []
    top = max(slab.positions[:, 2])
    for _ in slab.positions:
        if (abs(_[2] - top)) < 0.1:
            topLayers.append(_)

    ontop = topLayers[0]
    bridge = [ontop[0]+math.sqrt(2)/4.0*latt, ontop[1], ontop[2]]
    hollow= [ontop[0]+math.sqrt(2)/4.0*latt, ontop[1]+math.sqrt(6)/12.0*latt, ontop[2]]

    target=[]
    if (site =='ontop'):
        target=ontop
    elif (site=='bridge'):
        target=bridge
    elif(site=='hollow'):
        target=hollow


    adsorbates=Atoms(mol.symbols, mol.positions)
    adsorbates.set_cell(slab.get_cell(), scale_atoms=False)
    Tbposition = []

    for atom in adsorbates:
        if (atom.symbol == 'Tb'):
            Tbposition.append(atom.position)

    positions = sorted(Tbposition, key=lambda x: x[2])
    adsorbates.rotate(angle,'z',center=positions)
    vector = target - positions[0]
    adsorbates.translate(vector)
    adsorbates.translate([0.0, 0.0, dist])

    slab.extend(adsorbates)
    slab.center(vacuum=vacuum, axis=2)

    bottom = min(slab.positions[:, 2])
    fix = FixAtoms(indices=[atom.index for atom in slab if abs(atom.position[2] - bottom) < 0.5])
    slab.set_constraint(fix)
    write("POSCAR", slab, format='vasp', vasp5=True, direct=True)




