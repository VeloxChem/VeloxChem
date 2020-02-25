#!/usr/bin/env python3

import subprocess
import platform
import site
import sys
import os


def find_exe(executables):
    for exe in executables:
        for path in os.environ['PATH'].split(os.pathsep):
            fname = os.path.join(path, exe)
            if os.path.isfile(fname) and os.access(fname, os.X_OK):
                return exe, path
    return None, None


def get_command_output(command):
    try:
        output = subprocess.check_output(command)
    except subprocess.CalledProcessError:
        print()
        print('*** Error: Unable to execute \'{}\''.format(' '.join(command)))
        sys.exit(1)
    return output.decode('utf-8')


def find_avx_linux():
    cpuinfo = os.path.join(os.sep, 'proc', 'cpuinfo')
    if os.path.isfile(cpuinfo):
        for avx in ['avx512', 'avx2', 'avx']:
            with open(cpuinfo, 'r') as fh:
                for line in fh:
                    if line[:5] == 'flags' and avx in line:
                        return avx
    return None


def find_avx_macos():
    output = get_command_output(['sysctl', '-a'])
    lines = output.split(os.linesep)
    for line in lines:
        if 'machdep.cpu.leaf7_features' in line and 'AVX2' in line:
            return 'avx2'
    for line in lines:
        if 'machdep.cpu.features' in line and 'AVX1' in line:
            return 'avx'
    return None


def find_mkl_avx(is_linux, is_macos):

    print('*** Checking avx... ', end='')
    if is_linux:
        avx = find_avx_linux()
        if check_pdc_beskow():
            avx = 'avx2'
    elif is_macos:
        avx = find_avx_macos()
    print(avx)

    if avx is not None:
        mkl_avx = '-lmkl_{}'.format(avx)
    else:
        mkl_avx = '-lmkl_def'
    return mkl_avx


def check_ubuntu():
    for name in ['lsb-release', 'os-release']:
        fname = os.path.join(os.sep, 'etc', name)
        if os.path.isfile(fname):
            with open(fname, 'r') as fh:
                for line in fh:
                    if 'ubuntu' in line.lower():
                        return True
    return False


def check_pdc_beskow():
    is_pdc = ('SNIC_SITE' in os.environ and os.environ['SNIC_SITE'] == 'pdc')
    is_beskow = ('SNIC_RESOURCE' in os.environ and
                 os.environ['SNIC_RESOURCE'] == 'beskow')
    return (is_pdc and is_beskow)


def check_dir(dirname, label):
    if not os.path.isdir(dirname):
        print('*** Error: {} dir {} does not exist!'.format(label, dirname))
        sys.exit(1)


def check_file(filename, label):
    if not os.path.isfile(filename):
        print('*** Error: {} file {} does not exist!'.format(label, filename))
        sys.exit(1)


def generate_setup(template_file, setup_file):

    # OS information

    print('*** Checking operating system... ', end='')

    is_linux = ('Linux' == platform.system())
    is_macos = ('Darwin' == platform.system())
    if not (is_linux or is_macos):
        print()
        print('*** Error: Unsupported OS!')
        print('***        Only Linux and MacOS are supported.')
        sys.exit(1)

    is_ubuntu = (is_linux and check_ubuntu())

    if is_linux:
        print('Linux')
    elif is_macos:
        print('MacOS')

    # compiler information

    print('*** Checking c++ compiler... ', end='')

    if 'CXX' in os.environ:
        cxx, cxx_path = find_exe([os.environ['CXX']])
    else:
        cxx, cxx_path = find_exe(['mpiicpc', 'mpicxx', 'mpiCXX'])

    print(cxx)

    if cxx is None:
        print('*** Error: Unable to find c++ compiler!')
        if 'CXX' in os.environ:
            print('***        Please make sure that CXX is correctly set.')
        else:
            print('***        Please make sure that mpiicpc, mpicxx, or')
            print('***        mpiCXX is in your PATH.')
        sys.exit(1)

    if cxx in ['icpc', 'g++', 'clang++']:
        print('*** Error: {} is not a MPI compiler!'.format(cxx))
        sys.exit(1)

    if cxx in ['mpiicpc', 'mpicxx', 'mpiCXX']:
        cxxname = get_command_output([cxx, '-show'])
    else:
        cxxname = get_command_output([cxx, '--version'])
    if 'Cray clang' in cxxname:
        cxxname = cxxname.replace('Cray clang', 'Crayclang')
    cxxname = cxxname.split()[0]

    if cxxname in ['icc', 'gcc', 'clang']:
        print('*** Error: {} is not a c++ compiler!'.format(cxx))
        sys.exit(1)

    use_intel = (cxxname == 'icpc')
    use_gnu = cxxname in ['g++', 'x86_64-conda_cos6-linux-gnu-c++']
    use_clang = cxxname in ['clang++', 'Crayclang']

    if not (use_intel or use_gnu or use_clang):
        print('*** Error: Unrecognized c++ compiler!')
        print('***        Only Intel, GNU, and Clang compilers are supported.')
        sys.exit(1)

    # cxx and omp flags

    if use_intel:
        cxx_flags = '-xHost -qopenmp'
        if check_pdc_beskow():
            cxx_flags = '-qopenmp'
        omp_flag = '-liomp5'
    elif use_gnu:
        cxx_flags = '-fopenmp'
        omp_flag = '-lgomp'
    elif use_clang:
        cxx_flags = '-Xpreprocessor -fopenmp'
        omp_flag = '-lomp'

    # math library

    print('*** Checking math library... ', end='')

    use_mkl = 'MKLROOT' in os.environ
    use_openblas = 'OPENBLASROOT' in os.environ
    use_craylibsci = 'CRAY_LIBSCI_VERSION' in os.environ

    if not (use_mkl or use_openblas or use_craylibsci):
        print()
        print('*** Error: Unable to find math library!')
        print('***        Please make sure that you have set MKLROOT or')
        print('***        OPENBLASROOT. OpenBLAS can be downloaded from')
        print('***        https://github.com/xianyi/OpenBLAS')
        sys.exit(1)

    # mkl flags

    if use_mkl:
        print('MKL')

        mkl_inc = os.path.join(os.environ['MKLROOT'], 'include')
        check_dir(mkl_inc, 'mkl include')

        mkl_dir = os.path.join(os.environ['MKLROOT'], 'lib', 'intel64')
        if not os.path.isdir(mkl_dir):
            mkl_dir = os.path.join(os.environ['MKLROOT'], 'lib')
        check_dir(mkl_dir, 'mkl lib')

        if is_ubuntu and not use_intel:
            mkl_rt = '-lmkl_rt'
        else:
            mkl_avx = find_mkl_avx(is_linux, is_macos)
            mkl_rt = '-lmkl_intel_lp64 -lmkl_core {}'.format(mkl_avx)

        if use_intel or use_clang:
            mkl_thread = '-lmkl_intel_thread'
        elif use_gnu:
            mkl_thread = '-lmkl_gnu_thread'

        math_lib = 'MATH_INC := -I{}'.format(mkl_inc)
        math_lib += os.linesep + 'MATH_LIB := -L{}'.format(mkl_dir)
        math_lib += os.linesep + 'MATH_LIB += -Wl,-rpath,{}'.format(mkl_dir)
        math_lib += os.linesep + 'MATH_LIB += {} {}'.format(mkl_rt, mkl_thread)
        math_lib += os.linesep + 'MATH_LIB += {} -lpthread -lm -ldl'.format(
            omp_flag)

    # openblas flags

    elif use_openblas:
        print('OpenBLAS')

        openblas_inc = os.path.join(os.environ['OPENBLASROOT'], 'include')
        check_dir(openblas_inc, 'openblas include')

        openblas_dir = os.path.join(os.environ['OPENBLASROOT'], 'lib')
        check_dir(openblas_dir, 'openblas lib')

        math_lib = 'MATH_INC := -I{}'.format(openblas_inc)
        math_lib += os.linesep + 'MATH_LIB := -L{}'.format(openblas_dir)
        math_lib += os.linesep + 'MATH_LIB += -Wl,-rpath,{}'.format(
            openblas_dir)
        math_lib += os.linesep + 'MATH_LIB += -lopenblas {} {}'.format(
            omp_flag, '-lpthread -lm -ldl')

    # cray-libsci flags

    elif use_craylibsci:
        print('Cray LibSci')

        math_lib = 'MATH_INC := '
        math_lib += os.linesep + 'MATH_LIB := {} {}'.format(
            omp_flag, '-lpthread -lm -ldl')

    # extra flags for mac

    maclibs = ''
    if is_macos:
        maclibs = '-undefined dynamic_lookup'

    # lto flag

    lto_flag = ''
    if use_gnu:
        lto_flag = '-fno-lto'

    # pybind11

    try:
        import pybind11
    except ImportError:
        print()
        print('*** Error: Unable to find pybind11!')
        print('***        Please install via \"pip install pybind11 [--user]\"')
        sys.exit(1)

    python_user_base = site.getuserbase()

    # google test lib

    if 'GTESTROOT' in os.environ:
        gtest_root = os.environ['GTESTROOT']
        gtest_inc = os.path.join(gtest_root, 'include')
        check_dir(gtest_inc, 'GoogleTest include')
        if 'GTESTLIB' in os.environ:
            gtest_lib = os.environ['GTESTLIB']
        else:
            gtest_lib = os.path.join(gtest_root, 'lib', 'libgtest.a')
        check_file(gtest_lib, 'GoogleTest lib')
    else:
        gtest_root = None
        gtest_lib = None

    # print Makefile.setup

    with open(template_file, 'r', encoding='utf-8') as f_temp:
        lines = f_temp.readlines()

    with open(setup_file, 'w', encoding='utf-8') as f_mkfile:
        for line in lines:
            if '====placeholder====' in line:
                print('# Automatically generated settings', file=f_mkfile)
                print('', file=f_mkfile)

                print('USE_MPI := true', file=f_mkfile)
                print('USE_MKL := {}'.format('true' if use_mkl else 'false'),
                      file=f_mkfile)
                print('', file=f_mkfile)

                print(math_lib, file=f_mkfile)
                print('', file=f_mkfile)

                print('PYTHON :=',
                      'python{}.{}'.format(sys.version_info[0],
                                           sys.version_info[1]),
                      file=f_mkfile)
                python_version = 'python{}.{}{}'.format(sys.version_info[0],
                                                        sys.version_info[1],
                                                        sys.abiflags)
                print('PYTHON_USER_INC :=',
                      os.path.join(python_user_base, 'include', python_version),
                      file=f_mkfile)
                print('', file=f_mkfile)

                print('CXX :=', cxx, file=f_mkfile)
                print('', file=f_mkfile)

                print('CXX_REL_FLG :=', cxx_flags, file=f_mkfile)
                print('CXX_DEB_FLG :=', cxx_flags, file=f_mkfile)
                print('', file=f_mkfile)

                print('MACLIBS :=', maclibs, file=f_mkfile)
                print('', file=f_mkfile)

                print('LTOFLAG :=', lto_flag, file=f_mkfile)
                print('', file=f_mkfile)

                if gtest_root is not None and gtest_lib is not None:
                    print('GST_ROOT :=', gtest_root, file=f_mkfile)
                    print('GST_LIB :=', gtest_lib, file=f_mkfile)
                    print('', file=f_mkfile)

                print('# Generic settings', file=f_mkfile)
            else:
                print(line, end='', file=f_mkfile)

    print('*** Successfully generated {}'.format(setup_file))


if __name__ == '__main__':

    template_file = 'Setup.template'
    setup_file = os.path.join(os.pardir, 'src', 'Makefile.setup')

    if not os.path.isfile(template_file):
        template_file = os.path.join('config', 'Setup.template')
        setup_file = os.path.join('src', 'Makefile.setup')

    if not os.path.isfile(template_file):
        print('*** Error: Cannot find template file {}'.format(template_file))
        sys.exit(1)

    generate_setup(template_file, setup_file)
