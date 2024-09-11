# To evaluate statistical extrapolation errors for a given parametrization.
# Functional parameters are given in 'for005.010.s*' where 's00' contains
# the basic parameter set and 's01'... the variations according to the
# principle axes of the error ellipsoid.
# The shell script assumes that the mean-field code 'haforpa.exe' resides
# two directory levels below the actual subdirectory. 
# The error estimate is taken as the r.m.s. deviation from the base result,
# see the code 'analyze.f90'.
#
# This shell script produces the errors for the whole density distributions
# of charge, proton, neutron and mass. The output is written to the file
# 'collect-variances-rho.dat'. The first block of resulots contains the
# densities together with their extrapolation errors. The second block
# contains the covariances.
#
# Input to the script is the wanted nucleus as:
#  $1  = proton number Z
#  $2  = neutron number N
#
echo "Z,N= $1 $2 " > collect-error-rho.dat
for i in s00 s01 s02 s03 s04 s05 s06 s07 s08 s09 s10 s11 s12 s13
do
cat for005.010.$i > for005
echo "$1,$2,,,3,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
echo "0" >> for005
echo "-1,,    IWUNIT" >> for005

../../haforpa.exe < for005 > /dev/null
grep -A 61 rho-coul for011  >> collect-error-rho.dat

done

gfortran analyze-rho.f90
./a.out > collect-variances-rho.dat
