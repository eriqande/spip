
# I used this script to compile for Windows on my Mac.
# Before using, I installed mingw via Homebrew, like:
#
# brew install mingw-w64
#


BIN=spip-Windows.exe

i686-w64-mingw32-gcc -O3 -o $BIN  \
  -Wno-incompatible-pointer-types-discards-qualifiers \
	src/spip.c \
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
