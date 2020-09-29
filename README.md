SPIP – simulate pedigrees in populations
================

  - [Installation and compilation](#installation-and-compilation)
  - [Quick documentation](#quick-documentation)
  - [Tutorial files](#tutorial-files)

## Installation and compilation

You need a C compiler and git. With those, you can clone this project
and compile the code into the spip binary like this in your Unix
terminal:

``` sh
# first clone it from GitHub, and get a submodule for it
git clone https://github.com/eriqande/spip.git
cd spip
git submodule init
git submodule update


# then compile it using the provided shell script
./CompileSpip.sh 
```

There will be a couple of warnings when you compile, but nothing fatal.

## Quick documentation

After that, you can get the short help for spip by typing `./spip
--help`. It will spit out a short listing of the options.

``` sh
./spip --help
```

    ## spip_m  ---  a program for simulating pedigrees within populations exchanging migrants
    ## 
    ##      --help                                        short listing of options
    ##      --help-full                                   long listing of options
    ##      --help-nroff                                  long listing of options in nroff man-page format
    ##      --help-xml                                    output all options in XML format (in development)
    ##      --version                                     prints program version information
    ##      --version-history                             prints history of different versions
    ##      --command-file       F                        inserts contents of file F into the command line
    ## 
    ##    ****  Options for dealing with multiple population simulations  ****
    ## 
    ##      --new-pop                                     signifies completion of options for the current population
    ##      --num-pops           J                        number of populations in the simulation
    ## 
    ##    ****  Age-structure and reproduction parameters  ****
    ## 
    ## -A , --max-age            J                        the maximum age the organism can reach
    ## -s , --survival-probs     R1 R2 ... R_MA           age specific probabilities of survival
    ##      --fem-surv-probs     R1 R2 ... R_MA           age-specific survival probs for females
    ##      --male-surv-probs    R1 R2 ... R_MA           age-specific survival probs for males
    ## -f , --fem-asrf           R1 R2 ... R_MA           female age specific relative fecundities
    ##      --fem-prob-repro     R1 R2 ... R_MA           female age specific probabilities of reproducing
    ##      --repro-inhib        J                        age of offspring below which reproduction is inhibited
    ## -m , --male-asrp          R1 R2 ... R_MA           male age specific reproductive potential
    ##      --male-prob-repro    R1 R2 ... R_MA           male age specific probabilities of reproducing
    ##      --offsp-dsn          C                        distribution of offspring number
    ##      --fem-rep-disp-par   R                        female dispersion of reprod. success
    ##      --male-rep-disp-par  R                        male dispersion of reprod. success
    ##      --male-sticky-rep-var R                        male variance in reproductive success maintained over years
    ##      --mate-fidelity      R                        mate fidelity parameter
    ##      --fem-postrep-die    R1 R2 ... R_MA           additional prob of death for a female that mates
    ##      --male-postrep-die   R1 R2 ... R_MA           additional prob of death for a male that mates
    ##      --sex-ratio          R                        probability an offspring is male
    ##      --year-var-switch                             not currently implemented
    ## 
    ##    ****  Spawning group restriction parameters  ****
    ## 
    ##      --spawn-grp-size     J                        Set spawning group sizes of J females and J males each year
    ##      --spawn-grp-cnt      J                        Randomly partition each year spawners into J randomly sized spawner groups
    ## 
    ##    ****  Migration rate and related parameters  ****
    ## 
    ##      --male-prob-mig-out  G1 G2 R                  time-and-age-specific prob. that a male migrates out of a population
    ##      --male-prob-mig-in   G1 G2 R0 R1 ... R_NP     time-and-age-specific prob of a male migrating to a population
    ##      --fem-prob-mig-out   G1 G2 R                  time-and-age-specific prob. that a female migrates out of a population
    ##      --fem-prob-mig-in    G1 G2 R0 R1 ... R_NP-1   time-and-age-specific prob of a female migrating to a population
    ##      --sneaky-males       R0 R1 ... R_NP-1         rates of sneaky-male matings
    ## 
    ##    ****  Initial-condition and duration-of-simulation parameters  ****
    ## 
    ## -T , --number-of-years    J                        number of years to simulate
    ##      --initial-males      J0 J1 ... J_MA-1         initial number of males of different ages
    ##      --initial-females    J0 J1 ... J_MA-1         initial number of females of different ages
    ## 
    ##    ****  Population Size Parameters and Option  ****
    ## 
    ##      --random-cohort-size                          makes the cohort sizes random variables
    ##      --fixed-cohort-size                           makes the cohort sizes fixed quantities
    ##      --cohort-size        C[...]                   cohort size determined by model C with pars ...
    ## 
    ##    ****  Genetic Sampling Parameters  ****
    ## 
    ##      --discard-all        J                        make all pedigree members born b/t J and J+MaxAge-1 the founders
    ##      --gtyp-ppn-fem-pre   G R1 R2 ... R_MA         genotype a proportion of females randomly prior to death
    ##      --gtyp-ppn-male-pre  G R1 R2 ... R_MA         genotype a proportion of males randomly prior to death
    ##      --gtyp-ppn-fem-post  G R1 R2 ... R_MA-1       genotype a proportion of females randomly after death episode
    ##      --gtyp-ppn-male-post G R1 R2 ... R_MA-1       genotype a proportion of males randomly after death episode
    ##      --gtyp-ppn-fem-dur   G R                      genotype a proportion R of reproducing females
    ##      --gtyp-ppn-male-dur  G R                      genotype a proportion R of reproducing males
    ##      --gtyps-for-all      G                        draw genotypes for all individuals born at times in range G
    ##      --lethal-sampling    R                        makes the -pre and -post sampling lethal with probability R
    ##      --locus-file         F                        F=pathname to file with locus information
    ## 
    ##    ****  Program Implementation Options  ****
    ## 
    ##      --alloc-extra        J                        Integer factor by which to increase memory allocation.

If you want more complete information about that, use

``` sh
./spip --help-full
```

You may wish to put `spip` in your PATH.

## Tutorial files

Now, to learn how to use spip, you can use the tutorial. Do like this
first:

``` sh
cp -r data/examples_for_dis ./
mv examples_for_dis examples
```
