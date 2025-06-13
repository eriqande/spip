
BIN=nbtest-$(uname)

#gcc -g -o $BIN  \
gcc -O3 -o $BIN  \
  -Wno-incompatible-pointer-types-discards-qualifiers \
	src/nbtest.c \
	eca-shared/ranlib/src/com.c \
	eca-shared/ranlib/linpack/linpack.c \
	eca-shared/ranlib/src/ranlib.c \
	eca-shared/ecalibs/ECA_Opt2.c \
	eca-shared/ecalibs/ECA_MemAlloc.c \
	eca-shared/ecalibs/MathStatRand.c \
	eca-shared/ecalibs/MCTypesEtc.c \
	-Ieca-shared/ranlib/src  \
	-Ieca-shared/ecalibs \
	-lm
