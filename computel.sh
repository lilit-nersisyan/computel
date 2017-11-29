#!/bin/bash

usage="\nProgram:\tcomputel
\nVersion:\t0.4.1 
\n\nusage:\t./computel.sh [options] {-1 <fq1> -2 <fq2> -3 <fq3> -o <o>}
\n\nInput:
\n\n\t<fq1>\tfastq file (the first pair or the only fastq file (for single end reads)
\n\n\t<fq2>\tfastq file (optional: the second pair of fastq files, if exists)
\n\n\t<fq3>\tfastq file (optional: the third pair of fastq files, if exists)
\n\n\t<o>\toutput directory (optional: the default is computel_out)
\n\nOptions (advanced):
\n\n\t<-proc>\tnumber of processors to be used (default: 4)
\n\n\t<-sam>\tsamtools path (optional: if not supplied, Computel will use the samtools installed on the system)
\n\n\t<-bowal>\tbowtie2-align path (optional: the bowtie2-align is located at computel's bin directory by default.)
\n\n\t<-bowb>\tbowtie2-build path (optional: the bowtie2-build is located at computel's bin directory by default.)
\n\n\t<-nchr>\tnumber of chromosomes in a haploid set (the default is 23)
\n\n\t<-lgenome>\twhole genome length (the default is 3244610000)
\n\n\t<-pattern>\ttelomere repeat pattern (the default is 'TTAGGG'; change this if you're using Computel for a non-human organism)
\n\n\t<-minseed>\tthe min seed length (read length minus the number of flanking N's in the telomeric index; should be in the range [12-read.length]; This is a tested and carefully set parameter (defualt = 12); Change this only if you REALLY KNOW what you're doing!)
\n\n\n********      TEST
\n\nTo test how this works navigate to computel directory and run:\n./computel.sh -1 src/examples/tel_reads1.fq.gz -2 src/examples/tel_reads2.fq.gz -o mytest
\n\nA successful test run, should return telomere length of 10683991 bp. 
\n\n********      TEST
\n"


#######################################################################
############### 	argument check		 ######################
#######################################################################

if [[ $# == 0 ]]; then
  echo -e $usage
  exit 1	
fi 



i=1
let n=$#+1
out="computel_out"
while [ $i -lt $n ]; do
	declare arg=${!i}
	case $arg in 
		"-1")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <fq1> argument not specified after \"-1\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare fq1=${!j}			
			fi		
		;;	

		"-2")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <fq2> argument not specified after \"-2\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare fq2=${!j}			
			fi
		;;	

		"-3")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <fq3> argument not specified after \"-3\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare fq3=${!j}			
			fi
		;;	
		"-o")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <o> argument not specified after \"o\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare out=${!j}							
			fi
		;;
		"-proc")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-proc> argument not specified after \"proc\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare proc=${!j}	
			fi
		;;	
		"-sam")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-sam> argument not specified after \"sam\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare samtools=${!j}	
			fi
		;;	
		"-bowal")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-bowal> argument not specified after \"bowal\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare bowtie_align=${!j}	
			fi
		;;	
		"-bowb")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-bowb> argument not specified after \"bowb\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare bowtie_build=${!j}	
			fi
		;;
		
		"-nchr")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-nchr> argument not specified after \"nchr\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare nchr=${!j}	
			fi
		;;
		"-lgenome")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-lgenome> argument not specified after \"lgenome\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare lgenome=${!j}	
			fi
		;;
		"-pattern")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-pattern> argument not specified after \"pattern\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare pattern=${!j}	
			fi
		;;
		"-minseed")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-minseed> argument not specified after \"minseed\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare minseed=${!j}	
			fi
		;;
	esac

	let i+=1
done


############### end of argument check ######################



#######################################################################
############### 	setup test 		#######################
#######################################################################

# find the directory of the script
computel_dir=`dirname $0`

bin_dir="$computel_dir/bin"
if [ ! -d $bin_dir ]; then
	echo "$(tput setaf 1;) \nerror: the bin directory is not found in the $computel_dir directory, nor it is specified. \"-1\"$(tput sgr0)"
	exit 1
fi

if [ -z $bowtie_build ]; then
	bowtie_build="$computel_dir/bin/bowtie2-build"
fi

bowtie_build=`readlink -m $bowtie_build`


if [ -z $bowtie_align ]; then
	bowtie_align="$computel_dir/bin/bowtie2-align"
fi
bowtie_align=`readlink -m $bowtie_align`


if [ -z $samtools ]; then		
	samtools=`which samtools`
	if [ -z $samtools ]; then 
		echo -e "$(tput setaf 1;) \nerror: samtools is not installed on your system.\nEither specify it with the -sam option or make sure it's in system path $(tput sgr0)"
	exit 1
	else		
		samtools=`readlink -m $samtools`
	fi
fi



samtofastq="$computel_dir/bin/SamToFastq.jar"

echo "Computel: testing setup:"

# testing bowtie2-build

if [ ! -f $bowtie_build ]; then
	echo -e "$(tput setaf 1;) \nError: bowtie2-build was not found at $bin_dir, nor was it specified. $(tput sgr0)"
	exit 1
fi
if [ ! -x $bowtie_build ]; then
	`chmod +x $bowtie_build`
	if [ ! -x $bowtie_build ]; then
		echo -e "$(tput setaf 1;) \nError: $bowtie_build is not executable. Could not change the permissions. Please run \"sudo chmod +x $bowtie_build\"  manually."
		exit 1
	else 
		echo -e "\t$bowtie_build set as executable"	
	fi
fi

`$bowtie_build 2>&1 | grep -A1 "Usage" | grep -q "bowtie2-build"`
if [ $? -eq 0 ]; then 
	echo -e "\t$bowtie_build is set properly"
else
	echo -e "$(tput setaf 1;) \nError: There are problems with $bowtie_build executable. Please, check if $bowtie_build works properly in your system. $(tput sgr0)"
	exit 1
fi

# testing bowtie2-align

if [ ! -f $bowtie_align ]; then
	echo -e "$(tput setaf 1;) \nError: bowtie2-align ($bowtie_align) does not exist. $(tput sgr0)"
	exit 1
fi

if [ ! -x $bowtie_align ]; then
	`chmod +x $bowtie_align`
	if [ ! -x $bowtie_align ]; then
		echo -e "$(tput setaf 1;) \nError: bowtie2-align ($bowtie_align) is not executable. Could not change the permissions. Please run \"sudo chmod +x $bowtie_align\"  manually."
		exit 1
	else 
		echo -e "\tbowtie2-align ($bowtie_align) set as executable"	
	fi
fi

`$bowtie_align 2>&1 | grep -A1 "Usage" | grep -q "bowtie2-align"`
if [ $? -eq 0 ]; then 
	echo -e "\t$bowtie_align is set properly"
else
	echo -e "$(tput setaf 1;) \nError: There are problems with bowtie2-align ($bowtie_align) executable. Please, check if it works properly in your system. $(tput sgr0)"
	exit 1
fi

# testing samtools
if [ ! -f $samtools ]; then
	echo -e "$(tput setaf 1;) \nError: the samtools executable $samtools does not exist. $(tput sgr0)"
	exit 1
fi
if [ ! -x $samtools ]; then
	`chmod +x $samtools`
	if [ ! -x $samtools ]; then
		echo -e "$(tput setaf 1;) \nError: $samtools is not executable. Could not change the permissions. Please run \"sudo chmod +x $samtools\"  manually. $(tput sgr0)"
		exit 1
	else 
		echo -e "\t$samtools set as executable"	
	fi
fi


`$samtools 2>&1 | grep -q "Usage.*samtools";`
if [ $? -eq 0 ]; then 
	echo -e "\t$samtools is set properly"
else
	echo -e "$(tput setaf 1;) \nerror: There are problems with samtools executable ($samtools). Please, check if it works properly in your system. $(tput sgr0)"
	exit 1	
fi

#check if -8000 is set properly for samtools

#first check if there is the "-d" option

`$samtools "depth" 2>&1  | grep -q "maximum coverage depth";`
if [ $? -eq 0 ]; then 
	echo -e "\t$samtools does have maximum coverage depth option"
else
	echo -e "$(tput setaf 1;) \nerror: $samtools does not have a -d option. This version of Computel works with samtools 1.3 or higher. Please use Computel.v0.3 for older samtools versions, or install a newer version of samtools. $(tput sgr0)"
	exit 1	
fi

#then check if it works
if [ ! -d "$computel_dir/setup_test" ]; then
	echo -e "$(tput setaf 3;)\nWarning: $(tput sgr0;) could not find directory $computel_dir/setup_test for testing samtools setup. \nPlease, locate the folder in the directory of computel.sh. \n"
fi

setup_dir="$computel_dir/setup_test"
test_bam="$setup_dir/test.bam"
if [ ! -f $test_bam ]; then
	echo -e "$(tput setaf 1;)\nError: $(tput sgr0;) could not find file $test_bam for testing samtools setup. \nPlease, locate the file from the originally downloaded computel package to the directory $setup_dir."
	exit 1
fi


depth="$setup_dir/test.depth.txt"

`$samtools depth -d 100000000 $test_bam > $depth`

max_cov=`sort -nrk3 $depth | head -1 | cut -f 3`

if [ $max_cov -gt 15000 ]; then 
	echo -e "\tsamtools checked, the 8000 limit is set properly." 
else
	echo -e "$(tput setaf 1;) \nError: samtools at $samtools is not setup to account for the limit of 8000.\nPlease, contact the application support team. $(tput sgr0;)"
	exit 1
fi




echo -e "Computel is setup properly.\n"
############### end of setup test ######################




#######################################################################
############### 	input files check		###############
#######################################################################

if [ -z $fq1 ]; then
	echo -e "$(tput setaf 1;) \nError:\tthe <fq1> argument is required but is not specified $(tput sgr0)"
	echo -e $usage
	exit 1	
fi

#check validity of fastq files
if [ ! -f $fq1 ]; then
	echo -e "$(tput setaf 1;) \nError:\tthe file <fq1> $fq1 does not exist $(tput sgr0)"
	echo -e $usage
	exit 1	
fi

if [ ! -z $fq2 ]; then
	if [ ! -f $fq2 ]; then
		echo -e "$(tput setaf 1;) \nError:\tthe file <fq2> $fq2 does not exist $(tput sgr0)"
		echo -e $usage
		exit 1	
	fi
fi

if [ ! -z $fq3 ]; then
	if [ ! -f $fq3 ]; then
		echo -e "$(tput setaf 1;) \nError:\tthe file <fq3> $fq3 does not exist $(tput sgr0)"
		echo -e $usage
		exit 1	
	fi
fi

# try not to overwrite the output file
if [ -f $out ]; then
	echo -e "$(tput setaf 1;)\nError:$(tput sgr0;)\tthe output directory $out is a regular file. Please provide a valid directory."
	exit 1
fi


if [ -d $out ]; then
	echo -e "$(tput setaf 3;)\nWarning:$(tput sgr0;)\tthe output directory $out already exists. Do you want to replace it? (y/n) "
	ans=""
	read ans
	while [ ! "$ans" == "y" ]; do
		if [ "$ans" == "n" ]; then
			echo "computel cancelled"
			exit 1
		fi
		read ans	
	done	
else
	mkdir $out
fi

#check the compression format of fastq files

function getCompression {
	mime=`file --mime-type $1`
#	echo $mime

	if file --mime-type $1 | grep -q gzip$ ; then
		compr="gz"
#		echo $1 is gzipped
	fi

	if file --mime-type $1 | grep -q bzip2$; then
		compr="bz2"
#		echo $1 is bzipped
	fi	

	if file --mime-type $1 | grep -q text\/plain$; then
		compr="none"
#		echo $1 is not compressed
	fi

	if [ -z $compr ]; then
		compr="unknown"
#		echo "could not determine mime type for file $1"
	fi
	echo $compr
}

compr1=`getCompression $fq1`

if [ ! -z $fq2 ]; then 
	compr2=`getCompression $fq2`
	if [ $compr1 != $compr2 ]; then
		echo -e "$(tput setaf 1;) \nerror:\tfiles $fq1 and $fq2 are of different compression types, should be the same. $(tput sgr0)"
		exit 1	
	fi
fi
if [ ! -z $fq3 ]; then 
	compr3=`getCompression $fq3`
	if [ $compr1 != $compr3 ]; then
		echo -e "$(tput setaf 1;) \nerror:\tfiles $fq1 and $fq3 are of different compression types, should be the same. $(tput sgr0)"
		exit 1	
	fi
fi

############### end of input files check ###################




#######################################################################
############	 setting configuration parameters	###############
#######################################################################

# determine the read length
if [ $compr1 == "none" ] || [ $compr1 == "unknown" ]; then
	first_line=`cat $fq1 | head -1`	
	if [[ $first_line != @* ]]; then 
		echo -e "$(tput setaf 1;) \nError:\tthe file $fq1 is not a valid fastq file (fastq files should start with '@') $(tput sgr0)"
		exit 1; 
	fi
	read_length=`cat $fq1 | head -4 | sed '2q;d' | wc -c`
	((read_length += -1))
fi 

# determine the read length
if [ $compr1 == "gz" ]; then
	first_line=`cat $fq1 | gunzip - | head -1`	
	if [[ $first_line != @* ]]; then 
		echo -e "$(tput setaf 1;) \nerror:\tthe file $fq1 is not a valid fastq file (fastq files should start with '@') $(tput sgr0)"
		exit 1; 
	fi
	read_length=`cat $fq1 | gunzip - | head -4 | sed '2q;d' | wc -c`
	((read_length += -1))
fi 

# if the files are bzipped, determine the read length
if [ $compr1 == "bz2" ]; then
	first_line=`cat $fq1 | bzip2 -dk - | head -1`	
	if [[ $first_line != @* ]]; then 
		echo -e "$(tput setaf 1;) \nerror:\tthe file $fq1 is not a valid fastq file (fastq files should start with '@') $(tput sgr0)"
		exit 1; 
	fi
	read_length=`cat $fq1 | bzip2 -dk - | head -4 | sed '2q;d' | wc -c`
	((read_length += -1))
fi 



if [ ! -z $fq2 ]; then 
	if [ ! -z $fq3 ]; then 
		fastq="$fq1,$fq2,$fq3" 
	else fastq="$fq1,$fq2" 
	fi 
else 
	fastq="$fq1" 
fi

if [ $compr1 == "none" ]; then
	file_compression="F"
else 
	file_compression=$compr1
fi

if [ -z $nchr ]; then
	nchr="23"
fi

if [ -z $lgenome ]; then
	lgenome="3244610000"
fi


if [ -z $pattern ]; then
	pattern="TTAGGG"
fi


if [ -z $minseed ]; then
	minseed="12"
fi

if [ -z $proc ]; then
	proc="4"
fi
############### end of setting configuration parameters ######################



#######################################################################
##########	generating config file and running	###############
#######################################################################


# generate config file
cat > "$out/config_file" << bzez
scripts.dir	$computel_dir/src/scripts
bowtie.build.path	$bowtie_build
bowtie.align.path	$bowtie_align
samtools.path	$samtools
picard.samtofastq.jar	$samtofastq


single	T
files.with.prefix	F
file.compression	$file_compression
fastq	$fastq
read.length	$read_length

pattern	$pattern
num.haploid.chr	$nchr
min.seed	$minseed
mode.local	F

compute.base.cov	F
##if files are compressed and the base coverage needs to be estimated during unzipping, set "estimate.base.cov" to T
estimate.base.cov	T
genome.length	$lgenome

output.dir	$out	

num.proc	$proc
ignore.err	F
quals	--phred33

bzez

Rscript $computel_dir/src/scripts/computel.cmd.R $out/config_file

########## 		the 	###############
