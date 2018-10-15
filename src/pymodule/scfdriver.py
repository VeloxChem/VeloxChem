
class ScfDriver:

    def __init__(self):
        # scf type

        self.scf_type = "RHF"

        # initial guess

        self.gues_type = "SAD"

        # scf accelerator
        
        self.acc_type = "L2DIIS"
        self.max_vecs = 10

        # convergence

        self.max_iter = 50
        
        # screening scheme
        
        self.qq_type = "QQRDEN"
        self.qq_dyn = True

        # thresholds

        self.conv_thresh = 1.0e-6
        self.eri_thresh = 1.0e-12
        self.ovl_thresh = 1.0e-12
    
    def compute(self, comm, ostream):
        # set up MPI communicator data

        loc_rank = comm.Get_rank()
        loc_nodes = comm.Get_size()
        
        self.print_header(ostream)


    def print_header(self, ostream):
        ostream.new_line()
        ostream.put_header("Self Consistent Field Driver Setup")
        ostream.put_header(36 * "=")
        ostream.new_line()
        str_width = 80
        cur_str = "WaveFunction Model           : " + self.get_scf_type()
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "Initial Guess Model          : " + self.get_guess_type()
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "Convergence Accelerator      : " + self.get_acc_type()
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "Max. Number Of Iterations    : " + str(self.max_iter)
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "Max. Number Of Error Vectors : " + str(self.max_vecs)
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "Convergence Threshold        : " + \
            "{:.1e}".format(self.conv_thresh)
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "ERI screening scheme         : " + self.get_qq_type()
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "ERI screening mode           : " + self.get_qq_dyn()
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "ERI Screening Threshold      : " + \
            "{:.1e}".format(self.eri_thresh)
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "Linear Dependence Threshold  : " + \
            "{:.1e}".format(self.ovl_thresh)
        ostream.put_header(cur_str.ljust(str_width))
        ostream.new_line()

    def get_scf_type(self):
        if self.scf_type == "RHF":
            return "Spin Restricted Hatree-Fock"
        if self.scf_type == "UHF":
            return "Spin Unrestricted Hatree-Fock"
        return "Unknown Type"

    def get_guess_type(self):
        if self.gues_type == "SAD":
            return "Superposition of Atomic Densities"
        return "Unknown Type"

    def get_acc_type(self):
        if self.acc_type == "DIIS":
            return "Direct Inversion of Iterative Subspace"
        if self.acc_type == "L2DIIS":
            return "Two Level Direct Inversion of Iterative Subspace"
        return "Unknown Type"

    def get_qq_type(self):
        if self.qq_type == "QQ":
            return "Cauchy–Schwarz"
        if self.qq_type == "QQR":
            return "Distance Dependent Cauchy Schwarz"
        if self.qq_type == "QQDEN":
            return "Cauchy–Schwarz and Density"
        if self.qq_type == "QQRDEN":
            return "Distance Dependent Cauchy Schwarz and Density"
        return "Unknown Type"

    def get_qq_dyn(self):
        if self.qq_dyn:
            return "Dynamic"
        return "Static"
