from .errorhandler import assert_msg_critical
import sys
from importlib.metadata import version

try:
    import openmm.app as mmapp
    import openmm as mm
except ImportError:
    pass


class EvbReporter():
    #todo do this with force groups instead of different systems
    def __init__(self, file, report_interval, reference_reactant, reference_product, run_reactant, run_product, topology, Lambda, append = False):

        assert_msg_critical('openmm' in sys.modules, 'openmm is required for EvbReporter.')

        # OpenMM HIP version is slighly older and uses a different format for reporters
        if version('openmm') < '8.2':
            self.use_tuple = True
        else:
            self.use_tuple = False

        self.out = open(file, 'a' if append else 'w')
        self.report_interval = report_interval
        
        self.reference_product = reference_product
        self.run_reactant = run_reactant
        self.run_product = run_product

        self.Lambda = Lambda

        self.simulations = []
        for system in [reference_reactant, reference_product, run_reactant, run_product]:
            integrator = mm.VerletIntegrator(1)
            self.simulations.append(mmapp.Simulation(topology, system,integrator))
        if not append:
            header = "Lambda, E_ref_reactant, E_ref_product, E_run_reactant, E_run_product, E_m\n"
            self.out.write(header)

    def __del__(self):
        self.out.close()

    def describeNextReport(self, simulation):
        steps = self.report_interval - simulation.currentStep%self.report_interval
        if self.use_tuple:
            return (steps, True, True, False, False, True)
        else:
            return {'steps': steps, 'periodic': True, 'include':['positions','energy']}
        

    def report(self, simulation, state):

        positions = state.getPositions(asNumpy=True)
        E = [state.getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)]
        for simulation in self.simulations:
            simulation.context.setPositions(positions)
            state = simulation.context.getState(getEnergy=True)
            E.append(state.getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole))
        line = f"{self.Lambda}, {E[1]}, {E[2]}, {E[3]}, {E[4]}, {E[0]}\n"
        self.out.write(line)
        
            
