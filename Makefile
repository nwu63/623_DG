all:
	f2py -c -m --link-lapack_opt fortran *.f90

debug:
	f2py -c -m --link-lapack_opt --debug --noopt fortran *.f90

clean:
	rm -fr *.so