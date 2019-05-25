#!/usr/bin/env python3

import subprocess
import platform
import sys
import os


def find_exe(executables):
    for exe in executables:
        for path in os.environ['PATH'].split(os.pathsep):
            fname = os.path.join(path, exe)
            if os.path.isfile(fname) and os.access(fname, os.X_OK):
                return exe
    return None


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

    print('*** Checking compiler... ', end='')

    if 'CXX' in os.environ:
        cxx = find_exe([os.environ['CXX']])
    else:
        cxx = find_exe(['mpiicpc', 'mpicxx', 'mpiCXX'])

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
    cxxname = cxxname.split()[0]

    if cxxname in ['icc', 'gcc', 'clang']:
        print('*** Error: {} is not a c++ compiler!'.format(cxx))
        sys.exit(1)

    use_intel = (cxxname == 'icpc')
    use_gnu = (cxxname == 'g++')
    use_clang = (cxxname == 'clang++')

    if not (use_intel or use_gnu or use_clang):
        print('*** Error: Unrecognized c++ compiler!')
        print('***        Only Intel, GNU, and Clang compilers are supported.')
        sys.exit(1)

    # cxx and omp flags

    if use_intel:
        cxx_flags = '-xHost -qopenmp'
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
    use_openblas = not use_mkl and 'OPENBLASROOT' in os.environ

    if not (use_mkl or use_openblas):
        print()
        print('*** Error: Unable to find math library!')
        print('***        Please make sure that you have set MKLROOT or')
        print('***        OPENBLASROOT. OpenBLAS can be downloaded from')
        print('***        https://github.com/xianyi/OpenBLAS')
        sys.exit(1)

    # mkl flags

    if use_mkl:
        print('MKL')

        mkl_dir = os.path.join(os.environ['MKLROOT'], 'lib', 'intel64')
        if not os.path.isdir(mkl_dir):
            mkl_dir = os.path.join(os.environ['MKLROOT'], 'lib')
        if not os.path.isdir(mkl_dir):
            print('*** Error: mkl lib dir {} does not exist!'.format(mkl_dir))
            sys.exit(1)

        if is_ubuntu and not use_intel:
            mkl_rt = '-lmkl_rt'
        else:
            mkl_avx = find_mkl_avx(is_linux, is_macos)
            mkl_rt = '-lmkl_intel_lp64 -lmkl_core {}'.format(mkl_avx)

        if use_intel or use_clang:
            mkl_thread = '-lmkl_intel_thread'
        elif use_gnu:
            mkl_thread = '-lmkl_gnu_thread'

        mkl_libs = 'MKLLIBS := -L{}'.format(mkl_dir)
        mkl_libs += os.linesep + 'MKLLIBS += -Wl,-rpath,{}'.format(mkl_dir)
        mkl_libs += os.linesep + 'MKLLIBS += {} {}'.format(mkl_rt, mkl_thread)
        mkl_libs += os.linesep + 'MKLLIBS += {} -lpthread -lm -ldl'.format(
            omp_flag)

    # openblas flags

    if use_openblas:
        print('OpenBLAS')

        openblas_dir = os.path.join(os.environ['OPENBLASROOT'], 'lib')
        if not os.path.isdir(openblas_dir):
            print('*** Error: openblas lib dir {} does not exist!'.format(
                openblas_dir))
            sys.exit(1)

        openblas_libs = 'OPENBLASLIBS := -L{}'.format(openblas_dir)
        openblas_libs += os.linesep + 'OPENBLASLIBS += -Wl,-rpath,{}'.format(
            openblas_dir)
        openblas_libs += os.linesep + 'OPENBLASLIBS += -lopenblas {} {}'.format(
            omp_flag, '-lpthread -lm -ldl')

    # extra flags for mac

    maclibs = ''
    if is_macos:
        maclibs = '-undefined dynamic_lookup'

    # pybind11 rootdir

    print('*** Checking pybind11... ', end='')

    if 'PYBIND11ROOT' not in os.environ:
        print()
        print('*** Error: Unable to find pybind11!')
        print('***        Please make sure that you have set PYBIND11ROOT.')
        print('***        pybind11 can be downloaded from')
        print('***        https://github.com/pybind/pybind11')
        sys.exit(1)

    pybind11_root = os.environ['PYBIND11ROOT']
    print(pybind11_root)

    if not os.path.isdir(pybind11_root):
        print(
            '*** Error: pybind11 dir {} does not exist!'.format(pybind11_root))
        sys.exit(1)

    # TODO: add GPU detection

    # print Makefile.setup

    with open(template_file, 'r') as f_temp:
        lines = f_temp.readlines()

    with open(setup_file, 'w') as f_mkfile:
        for line in lines:
            if '====placeholder====' in line:
                print('# Automatically generated settings', file=f_mkfile)
                print('', file=f_mkfile)

                print('USE_MPI := true', file=f_mkfile)
                print('USE_GPU := false', file=f_mkfile)
                print('', file=f_mkfile)

                if use_mkl:
                    print('USE_MKL := true', file=f_mkfile)
                    print(mkl_libs, file=f_mkfile)
                else:
                    print('USE_MKL := false', file=f_mkfile)
                    print(openblas_libs, file=f_mkfile)
                print('', file=f_mkfile)

                print('PYTHON :=', 'python3', file=f_mkfile)
                print('PYBIND11_ROOT :=', pybind11_root, file=f_mkfile)
                print('', file=f_mkfile)

                print('CXX :=', cxx, file=f_mkfile)
                print('', file=f_mkfile)

                print('CXX_REL_FLG :=', cxx_flags, file=f_mkfile)
                print('CXX_DEB_FLG :=', cxx_flags, file=f_mkfile)
                print('', file=f_mkfile)

                print('MACLIBS :=', maclibs, file=f_mkfile)
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
