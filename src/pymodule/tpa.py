from .inputparser import parse_frequencies
from .cppsolver import ComplexResponse
from .tpadriver import tpa
from .veloxchemlib import ElectricDipoleIntegralsDriver
import numpy as np
from .veloxchemlib import mpi_master
import sys
from .linearsolver import LinearSolver


class twophotonabs(LinearSolver):

    def __init__(self, comm, ostream):
        np.set_printoptions(threshold=sys.maxsize)
        # tpa setup
        self.iso = True
        self.reduced = False

        # cpp setup
        self.frequencies = 0
        self.comp = "zzzz,0,0,0"
        self.damping = 0.004556335294880438
        self.eri_thresh = 1.0e-15
        self.conv_thresh = 1.0e-8
        self.lindep_thresh = 1.0e-10
        self.max_iter = 50

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream

    def update_settings(self, rsp_dict, method_dict={}):
        if 'frequencies' in rsp_dict:
            if 'unit' in rsp_dict:
                if rsp_dict['unit'] in 'eV':
                    ev_frequencies = rsp_dict['frequencies']
                    star = float(
                        ev_frequencies.split("-")[0]) / 27.218679420751958
                    en = float(
                        ev_frequencies.split("-")[1]) / 27.218679420751958
                    step = float(
                        ev_frequencies.split("-")[2]) / 27.218679420751958
                    self.frequencies = str(star) + "-" + str(en) + "-" + str(
                        step)
                    print(self.frequencies)
            else:
                self.frequencies = rsp_dict['frequencies']
        if 'damping' in rsp_dict:
            self.damping = float(rsp_dict['damping'])
        if 'reduced' in rsp_dict:
            self.reduced = rsp_dict['reduced']
            if self.reduced in 'True':
                self.reduced = True
            else:
                self.reduced = False
        if 'component' in rsp_dict:
            self.comp = rsp_dict['component']
            if self.comp in 'isotropic':
                self.iso = True
            else:
                self.iso = False
        if 'conv_thresh' in rsp_dict:
            self.conv_thresh = rsp_dict['conv_thresh']
        if 'lindep_thresh' in rsp_dict:
            self.lindep_thresh = rsp_dict['lindep_thresh']
        if 'max_iter' in rsp_dict:
            self.max_iter = rsp_dict['max_iter']

    def compute(self, molecule, ao_basis, scf_tensors):
        """
        This routine computes all the single index response vectors and then call all the relevent functions to compute γ
        """

        if self.rank == mpi_master():
            print("Number of freqs")
            num = parse_frequencies(self.frequencies)
            print(len(num))
            S = scf_tensors['S']  # The overlap Matrix in AO basis
            da = scf_tensors['D'][0]  # Alpha Density
            mo = scf_tensors['C']  #  MO coeff matrix
            nocc = molecule.number_of_alpha_electrons(
            )  # The number of occupied orbitals
            norb = mo.shape[1]  # The number of total orbitals
            w = parse_frequencies(self.frequencies)
            d_a_mo = tpa().ao2mo_inv(mo, da)
        else:
            S = None
            da = None
            mo = None
            d_a_mo = None
            nocc = None
            norb = None
            w = None

        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
        dipole_mats = dipole_drv.compute(molecule, ao_basis)

        rhs = {}
        operator = 'dipole'
        component = 'xyz'

        b_rhs = self.get_complex_rhs(operator, component, molecule,
                                     ao_basis, scf_tensors)

        if self.rank == mpi_master():
            v1 = {(op, wi): v for op, v in zip(component, b_rhs) for wi in w}
            X = {
                'x': 2 * tpa().ao2mo(mo, dipole_mats.x_to_numpy()),
                'y': 2 * tpa().ao2mo(mo, dipole_mats.y_to_numpy()),
                'z': 2 * tpa().ao2mo(mo, dipole_mats.z_to_numpy())
            }
            self.comp = self.get_comp(w)
        else:
            v1 = None
            X = None
            self.comp = None

        Nx = {}
        kX = {}
        Focks = {}

        # Compute the first set of response vectors
        Nb_drv = ComplexResponse(self.comm, self.ostream)
        Nb_drv.update_settings({
            'frequencies': self.frequencies,
            'damping': self.damping,
            'eri_thresh': self.eri_thresh,
            'conv_thresh': self.conv_thresh,
            'lindep_thresh': self.lindep_thresh,
            'max_iter': self.max_iter,
        })

        Nb_Drv = Nb_drv.compute(molecule, ao_basis, scf_tensors, v1)

        if self.rank == mpi_master():
            plot = []

            Nx.update({'Nb': Nb_Drv['solutions']})
            kX.update({'Nb': Nb_Drv['kappas']})
            Focks.update({'Fb': Nb_Drv['Focks']})

            for k in Nx['Nb'].keys():
                plot.append(np.max(np.abs(Nx['Nb'][k].imag)))

            Nx['Nc'] = {}
            kX['Nc'] = {}
            Focks['Fc'] = {}
            Focks['Fd'] = {}
            """ Obtain all the Nc('x',-w) and kC('x',-w) from the Kb and Nb"""
            for m in list(Nx['Nb'].keys()):
                Nx['Nc'].update({(m[0], -m[1]): tpa().FlipYZ(Nx['Nb'][m])})
                kX['Nc'].update({
                    (m[0], -m[1]):
                        LinearSolver.lrvec2mat(Nx['Nc'][(m[0], -m[1])].real,
                                               nocc, norb) +
                        LinearSolver.lrvec2mat(Nx['Nc'][
                            (m[0], -m[1])].imag, nocc, norb) * 1j
                })
                Focks['Fc'].update({(m[0], -m[1]): Focks['Fb'][m]})
                Focks['Fd'].update({
                    (m[0], m[1]): np.conjugate(Focks['Fb'][m]).T
                })

            Focks['Fb'].update(Focks['Fd'])

            Nx['Na'] = Nx['Nb']
            kX['Na'] = kX['Nb']

            Nx['Nd'] = Nx['Nb']
            kX['Nd'] = kX['Nb']

        NaX3NyNz, NaA3NxNy, NaX2Nyz, NxA2Nyz, E3_dict, E4_dict, NaX2Nyz_red, NxA2Nyz_red, E3_dict_red = tpa(
        ).main(self.eri_thresh, self.conv_thresh, self.lindep_thresh,
               self.max_iter, Focks, self.iso, Nx, w, X, self.damping, d_a_mo,
               kX, self.comp, S, da, mo, nocc, norb, scf_tensors, molecule,
               ao_basis, self.comm, self.ostream)

        if self.rank == mpi_master():
            T4 = tpa().T4(self.damping, w, E4_dict, Nx, kX, self.comp, d_a_mo,
                          nocc, norb)
            T3 = tpa().T3(w, E3_dict, Nx, self.comp)
            T3_red = tpa().T3(w, E3_dict_red, Nx, self.comp)
        else:
            T4 = None
            T3 = None

        gamma = {}
        gamma_red = {}

        if self.rank == mpi_master():
            for i in range(len(w)):
                if self.iso is True:
                    gamma.update({
                        (w[i], -w[i], w[i]):
                            1 / 15 *
                            (T4[(w[i], -w[i], w[i])] + T3[(w[i], -w[i], w[i])] +
                             NaX3NyNz[(w[i], -w[i], w[i])] +
                             NaA3NxNy[(w[i], -w[i], w[i])] +
                             NaX2Nyz[(w[i], -w[i], w[i])] +
                             NxA2Nyz[(w[i], -w[i], w[i])])
                    })
                    gamma_red.update({
                        (w[i], -w[i], w[i]): 1 / 15 *
                                             (T3_red[(w[i], -w[i], w[i])] +
                                              NaX2Nyz_red[(w[i], -w[i], w[i])] +
                                              NxA2Nyz_red[(w[i], -w[i], w[i])])
                    })

            self.print_header()
            self.print_properties(self.iso, w, gamma, self.comp, T4, T3,
                                  NaX3NyNz, NaA3NxNy, NaX2Nyz, NxA2Nyz)

            self.make_plots(w, gamma, gamma_red, plot)
        else:
            gamma = {}
            gamma_red = {}

        if self.rank == mpi_master():
            return T4, T3, NaX3NyNz, NaA3NxNy, NaX2Nyz, NxA2Nyz, gamma, w
        else:
            return None, None, None, None, None, None, None, None

    def print_header(self):
        """
        Prints TPA setup header to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header("Two-Photon Absorbtion Driver Setup")
        self.ostream.print_header(31 * "=")
        self.ostream.print_blank()

        width = 50

        cur_str = "Frequencies : " + str(
            self.frequencies) + " (start-stop-step)"
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = "Damping     : " + \
            "{:.1e}".format(self.damping)
        self.ostream.print_header(cur_str.ljust(width))
        self.ostream.print_blank()
        self.ostream.print_header(('-' * len(cur_str)).ljust(width))
        cur_str = "Eri threshold for response solver : " + str(self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = "Convergance threshold for response solver : " + str(
            self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = "Linear-dep threshold for response solver : " + str(
            self.lindep_thresh)
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = "Max iterations for response solver : " + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(width))

        self.ostream.flush()

    def print_properties(self, iso, w1, gamma, comp, T4, T3, NaX3NyNz, NaA3NxNy,
                         NaX2Nyz, NxA2Nyz):
        """
        Prints properties.

        :param props:
            The dictionary of properties.
        """
        width = 50

        self.ostream.print_blank()
        self.ostream.print_blank()
        w_str = "Gamma tensor components computed per frequency"
        self.ostream.print_blank()
        self.ostream.print_header(w_str.ljust(width))
        self.ostream.print_blank()
        count = 1
        for a in range(int(len(comp) * (1 / len(w1)))):
            w_str = str(count) + '. ' + str(comp[a].split(",")[0])
            self.ostream.print_header(w_str.ljust(width))
            count += 1

        self.ostream.print_blank()

        for w in w1:
            if iso is False:
                w_str = "ΣT3term =  {:.8f}".format(T3[w, -w, w])
            else:
                w_str = "ΣT3term =  {:.8f}".format(T3[w, -w, w] / 15)
            self.ostream.print_header(w_str.ljust(width))

            if iso is False:
                w_str = "ΣT4term =  {:.8f}".format(T4[w, -w, w])
            else:
                w_str = "ΣT4term =  {:.8f}".format(T4[w, -w, w] / 15)
            self.ostream.print_header(w_str.ljust(width))

            if iso is False:
                w_str = "ΣNaX2Nyz =  {:.8f}".format(NaX2Nyz[w, -w, w])
            else:
                w_str = "ΣNaX2Nyz =  {:.8f}".format(NaX2Nyz[w, -w, w] / 15)
            self.ostream.print_header(w_str.ljust(width))

            if iso is False:
                w_str = "ΣNaX3NyNz =  {:.8f}".format(NaX3NyNz[w, -w, w])
            else:
                w_str = "ΣNaX3NyNz =  {:.8f}".format(NaX3NyNz[w, -w, w] / 15)

            self.ostream.print_header(w_str.ljust(width))

            if iso is False:
                w_str = "ΣNxA2Nyz =  {:.8f}".format(NxA2Nyz[w, -w, w])
            else:
                w_str = "ΣNxA2Nyz =  {:.8f}".format(NxA2Nyz[w, -w, w] / 15)
            self.ostream.print_header(w_str.ljust(width))

            if iso is False:
                w_str = "ΣNaA3NxNy =  {:.8f}".format(NaA3NxNy[w, -w, w])
            else:
                w_str = "ΣNaA3NxNy =  {:.8f}".format(NaA3NxNy[w, -w, w] / 15)
            self.ostream.print_header(w_str.ljust(width))

            w_str = 'Σ<<μ;μ,μ,μ>>= {:.8f}, ω=({:.4f},{:.4f},{:.4f}), hω =({:.4f} eV)'.format(
                gamma[w, -w, w], w, -w, w, w * 54.437358841503915 / 2)
            self.ostream.print_header(w_str.ljust(width))
            self.ostream.print_header(('-' * len(w_str)).ljust(width))
            self.ostream.print_blank()

    def print_properties_red(self, iso, w1, gamma, comp, T3, NaX2Nyz, NxA2Nyz):
        """
        Prints properties.

        :param props:
            The dictionary of properties.
        """
        width = 50

        self.ostream.print_blank()
        self.ostream.print_blank()
        w_str = "Gamma tensor components computed per frequency"
        self.ostream.print_blank()
        self.ostream.print_header(w_str.ljust(width))
        self.ostream.print_blank()
        count = 1
        for a in range(int(len(comp) * (1 / len(w1)))):
            w_str = str(count) + '. ' + str(comp[a].split(",")[0])
            self.ostream.print_header(w_str.ljust(width))
            count += 1

        self.ostream.print_blank()

        for w in w1:

            if iso is False:
                w_str = "T3term =  {:.8f}".format(T3[w, -w, w])
            else:
                w_str = "T3term =  {:.8f}".format(T3[w, -w, w] / 15)
            self.ostream.print_header(w_str.ljust(width))

            if iso is False:
                w_str = "NaX2Nyz =  {:.8f}".format(NaX2Nyz[w, -w, w])
            else:
                w_str = "NaX2Nyz =  {:.8f}".format(NaX2Nyz[w, -w, w] / 15)
            self.ostream.print_header(w_str.ljust(width))

            if iso is False:
                w_str = "NxA2Nyz =  {:.8f}".format(NxA2Nyz[w, -w, w])
            else:
                w_str = "NxA2Nyz =  {:.8f}".format(NxA2Nyz[w, -w, w] / 15)
            self.ostream.print_header(w_str.ljust(width))

            w_str = '<<A;B,C,D>>= {:.8f}, w=({:.4f},{:.4f},{:.4f}), hω =({:.4f} eV)'.format(
                gamma[w, -w, w], w, -w, w, w * 54.437358841503915 / 2)
            self.ostream.print_header(w_str.ljust(width))
            self.ostream.print_header(('-' * len(w_str)).ljust(width))
            self.ostream.print_blank()

    def get_comp(self, w):
        if self.iso is True:
            """ If the isotropic γ is to be computed this part of the code generates a list with all the tensor components to which are sent to the computation routine"""
            spat_A = ['x', 'y', 'z']
            comp_iso = []

            for a, b, e in ((a, b, e) for e in range(len(w))
                            for b in range(len(spat_A))
                            for a in range(len(spat_A))):
                comp_iso.append(
                    str(spat_A[a]) + str(spat_A[a]) + str(spat_A[b]) +
                    str(spat_A[b]) + ',' + str(w[e]) + ',' + str(-w[e]) + ',' +
                    str(w[e]))
                comp_iso.append(
                    str(spat_A[a]) + str(spat_A[b]) + str(spat_A[a]) +
                    str(spat_A[b]) + ',' + str(w[e]) + ',' + str(-w[e]) + ',' +
                    str(w[e]))
                comp_iso.append(
                    str(spat_A[a]) + str(spat_A[b]) + str(spat_A[b]) +
                    str(spat_A[a]) + ',' + str(w[e]) + ',' + str(-w[e]) + ',' +
                    str(w[e]))

            self.comp = sorted(comp_iso, key=comp_iso.index)

        else:

            Track = self.comp.split(",")
            count = 0
            comp = []
            for i in range(len(w)):
                for case in Track:
                    t = case[0] + case[1] + case[2] + case[3] + ',' + str(
                        w[i]) + ',' + str(-w[i]) + ',' + str(w[i])
                    comp.append(t)
                self.comp = comp
                count += 1
        return self.comp

    def make_plots(self, w, gamma, gamma_red, plot):
        g_comb = open('py_plot_comb.py', "w")
        add_stuff = 'import matplotlib \nimport matplotlib.pyplot as plt \n'
        omega_for_plot = 'ω = ['
        data_for_plot = 'σ = ['
        N = 'N = ['
        N1 = 4
        af = 1 / 137
        pi = 3.14159265359
        C = N1 * pi**(2) * af**(2)
        for i in range(len(w)):
            omega_for_plot += str(w[i] * 54.437358841503915 / 2) + ","
            data_for_plot += str(
                abs(C * w[i]**2 * gamma[(w[i], -w[i], w[i])].imag)) + ","
            N += str(plot[i]) + ","
        omega_for_plot += "] \n"
        data_for_plot += "]\n"
        N += "]\n"

        g_comb.write(add_stuff)
        g_comb.write(omega_for_plot)
        g_comb.write(data_for_plot)
        g_comb.write(N)

        data_for_plot = 'σ_red = ['
        for i in range(len(w)):
            data_for_plot += str(
                abs(C * w[i]**2 * gamma_red[(w[i], -w[i], w[i])].imag)) + ","
        omega_for_plot += "] \n"
        data_for_plot += "]\n"

        g_comb.write(data_for_plot)

        g_comb.write('\nplt.xlabel(\'ω eV\')')
        g_comb.write('\nplt.ylabel(\'σ TPA cross\')\nplt.legend()')
        g_comb.write('\nplt.plot(ω,σ,\'-Dk\',label=\'Full\')')
        g_comb.write('\nplt.plot(ω,σ_red,\'-or\',label=\'Reduced\')')
        g_comb.write('\nplt.show()')
        #g_comb.write('\nplt.savefig(\"filename.eps\", facecolor=\'w\', edgecolor=\'w\',orientation=\'portrait\', papertype=None, format=\'eps\',transparent=True,frameon=None, metadata=None)')
