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


if __name__ == '__main__':

    # OS information

    is_linux = ('Linux' == platform.system())
    is_macos = ('Darwin' == platform.system())
    if not (is_linux or is_macos):
        print('*** Error: Unsupported OS!')
        print('***        Only Linux and MacOS are supported.')
        sys.exit(1)

    is_ubuntu = (is_linux and 'Ubuntu' in platform.uname()[3])

    # compiler information

    if 'CXX' in os.environ:
        cxx = find_exe([os.environ['CXX']])
    else:
        cxx = find_exe(['mpiicpc', 'mpicxx', 'mpiCXX'])

    if cxx is None:
        print('*** Error: Unable to find c++ compiler!')
        if 'CXX' in os.environ:
            print('***        Please make sure that CXX is correctly set.')
        else:
            print('***        Please make sure that mpiicpc, mpicxx, or')
            print('***        mpiCXX is in your PATH.')
        sys.exit(1)

    if cxx in ['mpiicpc', 'mpicxx', 'mpiCXX']:
        try:
            cxxname = subprocess.check_output([cxx, '-show'])
        except subprocess.CalledProcessError:
            print('*** Error: Unable to execute \'{} -show\'!'.format(cxx))
            sys.exit(1)
        cxxname = cxxname.split()[0].decode('utf-8')
        use_intel = (cxxname == 'icpc')
        use_gnu = (cxxname == 'g++')
        use_clang = (cxxname == 'clang++')
    else:
        print('*** It seems that you are using a compiler wrapper.')
        print('*** Please specify which compiler is being used')
        answer = input('(intel/gnu/clang): ').lower()
        use_intel = (answer == 'intel')
        use_gnu = (answer == 'gnu')
        use_clang = (answer == 'clang')

    if not (use_intel or use_gnu or use_clang):
        print('*** Error: Unsupported c++ compiler!')
        print('***        Only Intel, GNU, and Clang compilers are supported.')
        sys.exit(1)

    # math library

    use_mkl = 'MKLROOT' in os.environ
    use_openblas = not use_mkl and 'OPENBLASROOT' in os.environ

    if not (use_mkl or use_openblas):
        print('*** Error: Unable to find math library!')
        print('***        Please make sure that you have set MKLROOT or')
        print('***        OPENBLASROOT')
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

    # mkl flags

    if use_mkl:
        if use_intel or use_clang:
            mkl_thread = '-lmkl_intel_thread'
        elif use_gnu:
            mkl_thread = '-lmkl_gnu_thread'

        if is_ubuntu and not use_intel:
            mkl_rt = '-lmkl_rt'
        else:
            mkl_rt = '-lmkl_intel_lp64 -lmkl_core'

        if is_linux:
            mkl_dir = os.path.join('$(MKLROOT)', 'lib', 'intel64')
        elif is_macos:
            mkl_dir = os.path.join('$(MKLROOT)', 'lib')

        print('*** It seems that you are using mkl.')
        print('*** Please specify avx for mkl')
        answer = input('(avx512/avx2/avx): ').lower()
        if answer in ['avx512', 'avx2', 'avx']:
            mkl_avx = '-lmkl_{}'.format(answer)
            print('*** Using -lmkl_{}'.format(answer))
        else:
            mkl_avx = '-lmkl_def'
            print('*** {} not recognized; using -lmkl_def'.format(answer))

        mkl_libs = 'MKLLIBS := -L{} -Wl,-rpath,{}'.format(mkl_dir, mkl_dir)
        mkl_libs += os.linesep + 'MKLLIBS += {} {}'.format(mkl_rt, mkl_thread)
        mkl_libs += os.linesep + 'MKLLIBS += {} {} -lpthread -lm -ldl'.format(
            mkl_avx, omp_flag)

    # openblas flags

    if use_openblas:
        openblas_dir = os.path.join('$(OPENBLASROOT)', 'lib')
        openblas_libs = 'OPENBLASLIBS := -L{} '.format(openblas_dir)
        openblas_libs += '-Wl,-rpath,{} -lopenblas '.format(openblas_dir)
        openblas_libs += os.linesep + 'OPENBLASLIBS += {} '.format(omp_flag)
        openblas_libs += '-lpthread -lm -ldl'

    # extra flags for mac

    maclibs = ''
    if is_macos:
        maclibs = '-undefined dynamic_lookup'

    # pybind11 rootdir

    if 'PYBIND11ROOT' not in os.environ:
        print('*** Error: Unable to find pybind11!')
        print('***        Please make sure that you have set PYBIND11ROOT')
        sys.exit(1)

    pybind11_root = os.environ['PYBIND11ROOT']

    # print Makefile.setup

    template_setup = 'Setup.template'
    if not os.path.isfile(template_setup):
        template_setup = os.path.join('config', template_setup)
    makefile_setup = template_setup.replace('Setup.template', 'Makefile.setup')

    with open(template_setup, 'r') as f_temp:
        lines = f_temp.readlines()

    with open(makefile_setup, 'w') as f_mkfile:
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

    print('*** Successfully generated {}'.format(makefile_setup))
    print('*** Please copy {} to src'.format(makefile_setup))
