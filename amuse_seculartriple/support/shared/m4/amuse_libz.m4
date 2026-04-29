# Helper macro for detecting libz
#
# AMUSE_LIBZ()
#
# Searches for libz and sets FOUND_LIBZ to "yes" if found.
#
# There exists an AX_CHECK_ZLIB on the Autoconf Archive, but it checks specific
# directories instead of at least trying without specifying a directory first. As a
# result, it'll find a libz on the system even if we're in a Conda environment for
# example. This only searches in the default locations, relying on the package or
# environment manager to add appropriate flags to CPPFLAGS and LDFLAGS if needed.
#
AC_DEFUN([AMUSE_LIBZ], [
    amuse_libz_save_libs="$LIBS"

    AC_MSG_CHECKING([for libz])

    AC_LANG_PUSH([C])

    LIBS="-lz"

    AC_LINK_IFELSE([
        AC_LANG_PROGRAM([
            #include <zlib.h>
        ], [
            z_stream strm;
            deflateInit(&strm, 0);
        ])
    ], [
        FOUND_LIBZ="yes"
        AC_MSG_RESULT([yes])
    ], [
        AC_MSG_RESULT([no])
    ])

    AC_LANG_POP([C])

    AC_SUBST([LIBZ_CFLAGS])
    AC_SUBST([LIBZ_LIBS])
])

