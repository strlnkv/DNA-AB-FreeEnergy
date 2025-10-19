from openmm.app import *
from openmm import *
from openmm.unit import *
import gc
import numpy
from dcdsubsetreporter import DCDSubsetReporter

chr = CharmmPsfFile('../rst.psf')
par = CharmmParameterSet('../par_all36_na.prm','../top_all36_na.rtf','../par_all36_carb.prm','../toppar_all36_carb_model.str','../toppar_water_ions_jccufix.str')
chr.setBox(9.06*nanometer,9.06*nanometer,9.06*nanometer)
system = chr.createSystem(par, nonbondedMethod=PME, nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0001, constraints=HBonds)

nbforce = [system.getForce(i) for i in range(system.getNumForces()) if isinstance(system.getForce(i), openmm.NonbondedForce)][0]
nbfix = [system.getForce(i) for i in range(system.getNumForces()) if isinstance(system.getForce(i), openmm.CustomNonbondedForce)][0]
nbforce.setNonbondedMethod(openmm.NonbondedForce.PME)
nbforce.setEwaldErrorTolerance(0.0001)
nbforce.setCutoffDistance(1.2*nanometer)
nbforce.setUseSwitchingFunction(True)
nbforce.setSwitchingDistance(1.0*nanometer)
nbfix.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
nbfix.setCutoffDistance(1.2*nanometer)
nbfix.setUseSwitchingFunction(True)
nbfix.setSwitchingDistance(1.0*nanometer)

platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision','mixed')
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(chr.topology, system, integrator, platform)
with open('restart19.dat', 'r') as f:simulation.context.setState(openmm.XmlSerializer.deserialize(f.read()))
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))

HBrstr = CustomBondForce('4.184*(1.06395*(r*10+4.76673)^2-4*0.09591663*((0.18481695/r)^12-(0.18481695/r)^6)+33.206371*0.66*0.26/r)')
HBrstr.setUsesPeriodicBoundaryConditions(True)
rstrid = system.addForce(HBrstr)
HBrstr.addBond(615,646)
HBrstr.addBond(17,1244)

indices_list = []
for atom in simulation.topology.atoms():
    if atom.residue.name != 'HOH' and atom.residue.name != 'ETOH' and atom.residue.name != 'MEOH' and (not atom.name.startswith('D')) and (not atom.name.startswith('LP')):
        indices_list.append(atom.index)
print(f'Selected {len(indices_list)} / {simulation.topology.getNumAtoms()} atoms')

WHAMrstr = CustomCVForce('4.184*5*((w1+w2+w3+w4+w5+w6+w7+w8+w9+w10+w11+w12)/12*10-19)^2')
w1=CustomBondForce('r')
w1.addBond(657,373)
w1.setUsesPeriodicBoundaryConditions(True)
w2=CustomBondForce('r')
w2.addBond(687,343)
w2.setUsesPeriodicBoundaryConditions(True)
w3=CustomBondForce('r')
w3.addBond(720,313)
w3.setUsesPeriodicBoundaryConditions(True)
w4=CustomBondForce('r')
w4.addBond(753,280)
w4.setUsesPeriodicBoundaryConditions(True)
w5=CustomBondForce('r')
w5.addBond(786,247)
w5.setUsesPeriodicBoundaryConditions(True)
w6=CustomBondForce('r')
w6.addBond(816,217)
w6.setUsesPeriodicBoundaryConditions(True)
w7=CustomBondForce('r')
w7.addBond(846,187)
w7.setUsesPeriodicBoundaryConditions(True)
w8=CustomBondForce('r')
w8.addBond(876,157)
w8.setUsesPeriodicBoundaryConditions(True)
w9=CustomBondForce('r')
w9.addBond(909,124)
w9.setUsesPeriodicBoundaryConditions(True)
w10=CustomBondForce('r')
w10.addBond(942,91)
w10.setUsesPeriodicBoundaryConditions(True)
w11=CustomBondForce('r')
w11.addBond(972,58)
w11.setUsesPeriodicBoundaryConditions(True)
w12=CustomBondForce('r')
w12.addBond(1002,28)
w12.setUsesPeriodicBoundaryConditions(True)
WHAMrstr.addCollectiveVariable('w1', w1)
WHAMrstr.addCollectiveVariable('w2', w2)
WHAMrstr.addCollectiveVariable('w3', w3)
WHAMrstr.addCollectiveVariable('w4', w4)
WHAMrstr.addCollectiveVariable('w5', w5)
WHAMrstr.addCollectiveVariable('w6', w6)
WHAMrstr.addCollectiveVariable('w7', w7)
WHAMrstr.addCollectiveVariable('w8', w8)
WHAMrstr.addCollectiveVariable('w9', w9)
WHAMrstr.addCollectiveVariable('w10', w10)
WHAMrstr.addCollectiveVariable('w11', w11)
WHAMrstr.addCollectiveVariable('w12', w12)
whamrstrid = system.addForce(WHAMrstr)
simulation.context.reinitialize(preserveState=True)
simulation.reporters.append(DCDSubsetReporter('trj19.dcd', 5000, indices_list,append=True))#
simulation.reporters.append(StateDataReporter('log19.txt', 5000, step=True, potentialEnergy=True, temperature=True, volume=True, speed=True, separator='	',append=True))#
simulation.runForClockTime(1.4)
state = simulation.context.getState( getPositions=True, getVelocities=True )
with open('restart19.dat', 'w') as f:f.write(openmm.XmlSerializer.serialize(state))
