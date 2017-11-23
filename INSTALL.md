### Installation

For installation download and uncompress the Computel package in a local directory. The required binaries and files for setup configuration are set in the package.
Make computel.sh executable by running 'chmod +x computel.sh'. 



### Usage 
The basic usage is:

./computel.sh [options] -1 <fq1> -2 <fq2> -o <outputpath> 



### Test run
To test how this works navigate to computel directory and run:

./computel.sh -1 src/examples/tel_reads1.fq.gz -2 src/examples/tel_reads2.fq.gz -o mytest


In case of successful run, check if the content of mytest/tel.length.xls is the same as that in computel_out/tel.length.xls
You should get a warning that the coverage is low: this is ok for the test run, as the example files are of low coverage. 