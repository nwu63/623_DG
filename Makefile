all:
	f2py -c -m --link-lapack_opt fortran *.f90

debug:
	# f2py -c -DF2PY_REPORT_ON_ARRAY_COPY=1 --link-lapack_opt --debug --noopt -m fortran *.f90
	f2py -c --link-lapack_opt --debug --noopt -m fortran *.f90

clean:
	rm -fr *.o
	rm -fr *.so