all:
	f2py -c -m --link-lapack_opt dg_solver *.f90

debug:
	# f2py -c -DF2PY_REPORT_ON_ARRAY_COPY=1 --link-lapack_opt --debug --noopt -m fortran *.f90
	f2py -c --f90flags='-g -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan' --link-lapack_opt --debug --noopt -m dg_solver *.f90
clean:
	rm -fr *.o
	rm -fr *.so