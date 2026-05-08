# Helper macros for detecting download tools
#
# Some of the community codes aren't included with AMUSE, but are downloaded at build
# time. We need a tool for that, and here is where we find one.
#
# AMUSE_DOWNLOAD()
#
# Searches for a download tool.
#
# This macro tries to find a downloader and sets DOWNLOAD to a command that will take a
# URL and download its contents to standard output. This makes it easier to write that
# output to a file with a known name, which you usually want in a Makefile.
#
# This used to support both wget and curl, but we had problems with wget with BitBucket
# and a MESA server, so that's been removed and it's curl-only now.
#
# To download a file, use $(DOWNLOAD) https://example.com >example.html
#
AC_DEFUN([AMUSE_DOWNLOAD], [
    AC_CHECK_TOOL(CURL, curl)

    AC_MSG_CHECKING([for a tool to download files with])
    if test "x$CURL" != "x"
    then
        DOWNLOAD="$CURL -L"
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
    fi

    AC_SUBST([DOWNLOAD])
])

