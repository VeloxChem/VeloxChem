# This file is from autocmake, which is licensed under the 3-Clause BSD license.
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Copyright (c) 2015-2017 by Radovan Bast, Roberto Di Remigio, Jonas Juselius, and contributors.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of Autocmake nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#.rst:
#
# Take care of updating the cache for fresh configurations.
#
# Variables modified (provided the corresponding language is enabled)::
#
#   DEFAULT_Fortran_FLAGS_SET
#   DEFAULT_C_FLAGS_SET
#   DEFAULT_CXX_FLAGS_SET

macro(save_compiler_flags lang)
    if (NOT DEFINED DEFAULT_${lang}_FLAGS_SET)
        mark_as_advanced(DEFAULT_${lang}_FLAGS_SET)

        set (DEFAULT_${lang}_FLAGS_SET ON
            CACHE INTERNAL
            "Flag that the default ${lang} compiler flags have been set.")

        set(CMAKE_${lang}_FLAGS "${CMAKE_${lang}_FLAGS}"
            CACHE STRING
            "Flags used by the compiler during all builds." FORCE)

        set(CMAKE_${lang}_FLAGS_DEBUG "${CMAKE_${lang}_FLAGS_DEBUG}"
            CACHE STRING
            "Flags used by the compiler during debug builds." FORCE)

        set(CMAKE_${lang}_FLAGS_RELEASE "${CMAKE_${lang}_FLAGS_RELEASE}"
            CACHE STRING
            "Flags used by the compiler during release builds." FORCE)
    endif()
endmacro()

if(DEFINED CMAKE_Fortran_COMPILER_ID)
    save_compiler_flags(Fortran)
endif()

if(DEFINED CMAKE_C_COMPILER_ID)
    save_compiler_flags(C)
endif()

if(DEFINED CMAKE_CXX_COMPILER_ID)
    save_compiler_flags(CXX)
endif()
