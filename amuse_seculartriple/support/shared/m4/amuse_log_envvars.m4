# Macro for logging important environment variables

# AMUSE_LOG_ENVVARS()
#
AC_DEFUN([AMUSE_LOG_ENVVARS], [
    amuse_le_vars="SHELL CC CXX FC F77 MPICC MPICXX MPIFC CFLAGS CXXFLAGS CPPFLAGS LDFLAGS LIBS OMPI_CC OMPI_CXX OMPI_FC OMPI_CPPFLAGS OMPI_LDFLAGS OMPI_LIBS"

    AS_ECHO(["=== Shell environment ==="])

    for amuse_le_var in ${amuse_le_vars}
    do
        eval "amuse_le_val=\${$amuse_le_var}"
        AS_ECHO(["${amuse_le_var}=${amuse_le_val}"])
    done
    AS_ECHO(["=== End shell environment ==="])
])

