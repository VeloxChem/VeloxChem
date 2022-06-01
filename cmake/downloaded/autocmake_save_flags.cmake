# (c) https://github.com/dev-cafe/autocmake/blob/master/AUTHORS.md
# licensed under BSD-3: https://github.com/dev-cafe/autocmake/blob/master/LICENSE

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
