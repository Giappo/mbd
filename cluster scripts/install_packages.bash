#!/bin/bash

echo "home_dir = substring(getwd(),1,13)" > Rsetup.R
echo "libs_dir = paste0(home_dir,'/mbd_like/libs')" >> Rsetup.R
echo "lib_files = list.files(pattern=paste0('[.]tar'),path=libs_dir, full.names=TRUE)" >> Rsetup.R
echo "mylibrary = paste0(home_dir,'/R/x86_64-pc-linux-gnu-library/3.3/')" >> Rsetup.R
#echo "install.packages(lib_files, repos = NULL, lib = mylibrary, dependencies = TRUE)" >> Rsetup.R
echo "devtools::install_github("Giappo/MBD")" >> Rsetup.R

echo "#!/bin/bash" > Rsetup2
echo "#SBATCH --time=00:59:00" >> Rsetup2
#echo "#SBATCH --output=testinst.out" >> Rsetup2
echo "module load R/3.3.1-foss-2016a" >> Rsetup2
echo "Rscript Rsetup.R" >> Rsetup2
echo "rm Rsetup.R" >> Rsetup2
echo "rm Rsetup2" >> Rsetup2

sbatch Rsetup2


