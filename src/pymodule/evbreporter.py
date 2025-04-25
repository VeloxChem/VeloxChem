#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from importlib.metadata import version
import sys

from .errorhandler import assert_msg_critical

try:
    import openmm.app as mmapp
    import openmm as mm
except ImportError:
    pass


class EvbReporter():
    #todo do this with force groups instead of different systems
    def __init__(self, file, report_interval, reference_reactant, reference_product, run_reactant, run_product, topology, Lambda, outputstream, append = False):

        assert_msg_critical('openmm' in sys.modules and version('openmm') >= '8.2', 'openmm >=8.2 is required for EvbReporter.')

        # # OpenMM HIP version is slighly older and uses a different format for reporters
        # if version('openmm') < '8.2':
        #     outputstream.print_info('Older version of OpenMM detected. Using tuple format for returning reporter information.')
        #     outputstream.flush()
        #     self.use_tuple = True
        # else:
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
            return (steps, True, False, False, True, True) #steps, positions, velocities, forces, energy, pbc
        else:
            return {'steps': steps, 'periodic': True, 'include':['positions','energy']}
        
    def report(self, simulation, state):

        positions = state.getPositions(asNumpy=True)
        E = [state.getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)]
        for sim in self.simulations:
            sim.context.setPositions(positions)
            state = sim.context.getState(getEnergy=True)
            E.append(state.getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole))
        line = f"{self.Lambda}, {E[1]}, {E[2]}, {E[3]}, {E[4]}, {E[0]}\n"
        self.out.write(line)
        
            
