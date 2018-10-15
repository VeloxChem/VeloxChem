class ScfDriver:
    def __init__(self):
        # scf type
        
        self.scf_type = "RHF"
        
        # initial guess
        
        self.gues_type = "SAD"
        
        # scf accelerator
        self.acc_type = "DIIS"
        self.max_vecs = 10
        
        # convergence
        
        self.max_iter = 50
        
        # thresholds
        
        self.conv_thresh = 1.0e-6
        self.eri_thresh = 1.0e-12
        self.ovl_thresh = 1.0e-12
    
    def compute(self, ostream):
        # print header to output stream
        
        self.print_header(ostream)

    def print_header(self, ostream):
        ostream.new_line()
        ostream.put_header("Self Consistent Field Driver Setup")
        ostream.put_header(36 * "=")
        ostream.new_line()
        cur_str = "WaveFunction Model           : " + self.get_scf_type()
        ostream.put_header(cur_str.ljust(70))
        cur_str = "Initial Guess Model          : " + self.get_guess_type()
        ostream.put_header(cur_str.ljust(70))
        cur_str = "Convergence Accelerator      : " + self.get_acc_type()
        ostream.put_header(cur_str.ljust(70))
        cur_str = "Max. Number Of Iterations    : " + str(self.max_iter)
        ostream.put_header(cur_str.ljust(70))
        cur_str = "Max. Number Of Error Vectors : " + str(self.max_vecs)
        ostream.put_header(cur_str.ljust(70))
        cur_str = "Convergence Threshold        : " + "{:.1e}".format(self.conv_thresh)
        ostream.put_header(cur_str.ljust(70))
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


