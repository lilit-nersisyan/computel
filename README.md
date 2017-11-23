# Thank you for visiting our page! 

Computel is designed for measuring mean telomere length from Illumina Whole Genome NGS Sequencing data.

### Operating Systems
Computel v0.4 works with Unix operating system (tested for Ubuntu Linux). Computel version v0.2 works with Windows and Unix operating systems, however this version is harder to apply.

### Releases
Version v0.4 works with shell script and has a straightforward usage. This version of Computel uses the samtools (version 1.3 or higher) installed on the user's system, instead of the precompiled version coming with the previous releases. 
Version v0.3 works with shell script and has a straightforward usage.  
Version v0.2 works also with compressed fastq files.
Version v0.1 works both Windows and Unix type operating systems, works with configuration files and Rscript. Works also with compressed fastq files.

### System requirements
You have to have R (version 3.0.3 or higher) and samtools version 1.3 or higher (for Computel v0.4) installed in your system. 

### Installation
For installation download and uncompress the Computel package in a local directory. The required binaries and files for setup configuration are set in the package.
Make computel.sh executable by running 'chmod +x computel.sh'. 

### Usage 
The basic usage is:

./computel.sh [options] -1 \<fq1\> -2 \<fq2\> -o \<outputpath\> 


#### Binaries
Note, Computel works with only some versions of Bowtie, and with Samtools version 1.3 or higher. You can specify your own binaries with respective options, however, you may first want to refer to the Computel User Manual. 

### Citation
Use the following citation to read about our software, and cite it in your research:

Nersisyan L, Arakelyan A (2015) Computel: Computation of Mean Telomere Length from Whole-Genome Next-Generation Sequencing Data. PLoS ONE 10(4): e0125201. doi:10.1371/journal.pone.0125201


### Feedback 
This is the betta version of the software: your feedback at this point is critical for its further development! 
Please, join the Computel discussion group at https://groups.google.com/forum/#!forum/computel-discussion-forum, and add your comments. Please, do not hesitate to share any inconvenience you face with our software, any bug you notice and any suggestion you have. 

If you're happy with Computel, let us know - make us happy too! 
