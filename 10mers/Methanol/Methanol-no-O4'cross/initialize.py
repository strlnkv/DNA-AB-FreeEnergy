from openmm.app import *
from openmm import *
from openmm.unit import *
import gc
import numpy
from dcdsubsetreporter import DCDSubsetReporter

chr = CharmmPsfFile('rst.psf')
par = CharmmParameterSet('par_all36_na.prm','top_all36_na.rtf','par_all36_carb.prm','toppar_all36_carb_model.str','toppar_water_ions_jccufix.str')
chr.setBox(6.2*nanometer,6.2*nanometer,6.2*nanometer)
crd = PDBFile('rst.pdb')
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

#equilibrate solvent with DNA atoms fixed in place and ions having an elevated temerature of 500K
restrini = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
rstridini=system.addForce(restrini)
restrini.addGlobalParameter('k', 100.0*kilojoules_per_mole/nanometer)
restrini.addPerParticleParameter('x0')
restrini.addPerParticleParameter('y0')
restrini.addPerParticleParameter('z0')
fixat=["C1'","C2","C2'","C3'","C4","C4'","C5","C5M","C5'","C6","C8","N1","N2","N3","N4","N6","N7","N9","O1P","O2","O2P","O3'","O4","O4'","O5'","O6","P"]
for atom in chr.topology.atoms():
    if atom.name in fixat and (atom.residue.name=="A" or atom.residue.name=="T" or atom.residue.name=="G" or atom.residue.name=="C"):
        restrini.addParticle(atom.index, crd.positions[atom.index])
print('Nrestr=',restrini.getNumParticles())

ions=[]
rest=[]
for atom in chr.topology.atoms():
    if atom.name == 'SOD' or atom.name == 'CLA':
        ions=numpy.append(ions,atom.index)
    else:
        if system.getParticleMass(atom.index) != 0*amu: rest=numpy.append(rest,atom.index)
print('Nion=',len(ions),'Nrest=',len(rest))
#subsystems are only available for the NH integrator
integrator = NoseHooverIntegrator(0.002*picoseconds)
integrator.addSubsystemThermostat(rest, [], 300*kelvin, 1/picosecond, 300*kelvin, 1/picosecond)
integrator.addSubsystemThermostat(ions, [], 500*kelvin, 1/picosecond, 500*kelvin, 1/picosecond)
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision','mixed')
simulation = Simulation(chr.topology, system, integrator, platform)
simulation.context.setPositions(crd.positions)

simulation.minimizeEnergy(maxIterations=500)

system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
simulation.context.reinitialize(preserveState=True)

simulation.reporters.append(StateDataReporter('log.txt', 5000, step=True, potentialEnergy=True, temperature=True, volume=True, speed=True, separator='	'))
simulation.step(500000)

#slowly cool down ions
for i in range(1,21):
    integrator.getThermostat(1).setTemperature((500-10*i)*kelvin)
    integrator.step(25000)

#don't run equilibration because C form is unstable

state = simulation.context.getState( getPositions=True, getVelocities=True )
with open('restart.dat', 'w') as f:f.write(openmm.XmlSerializer.serialize(state))
