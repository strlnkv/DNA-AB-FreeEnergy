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
HBrstr.addBond(300,331)
HBrstr.addBond(17,614)

indices_list = []
for atom in simulation.topology.atoms():
    if atom.residue.name != 'HOH' and atom.residue.name != 'ETOH' and atom.residue.name != 'MEOH' and (not atom.name.startswith('D')) and (not atom.name.startswith('LP')):
        indices_list.append(atom.index)
print(f'Selected {len(indices_list)} / {simulation.topology.getNumAtoms()} atoms')

WHAMrstr = CustomCVForce('4.184*5*((w1+w2)/2*10-19)^2')
w1=CustomBondForce('r')
w1.addBond(28,372)
w1.setUsesPeriodicBoundaryConditions(True)
w2=CustomBondForce('r')
w2.addBond(58,342)
w2.setUsesPeriodicBoundaryConditions(True)
WHAMrstr.addCollectiveVariable('w1', w1)
WHAMrstr.addCollectiveVariable('w2', w2)
whamrstrid = system.addForce(WHAMrstr)
simulation.context.reinitialize(preserveState=True)
simulation.reporters.append(DCDSubsetReporter('trj19.dcd', 5000, indices_list,append=True))#
simulation.reporters.append(StateDataReporter('log19.txt', 5000, step=True, potentialEnergy=True, temperature=True, volume=True, speed=True, separator='	',append=True))#
simulation.runForClockTime(1.4)
state = simulation.context.getState( getPositions=True, getVelocities=True )
with open('restart19.dat', 'w') as f:f.write(openmm.XmlSerializer.serialize(state))
