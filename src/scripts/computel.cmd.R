
Introduction
============

What is Bowtie 2?
-----------------

[Bowtie 2] is an ultrafast and memory-efficient tool for aligning sequencing
reads to long reference sequences.  It is particularly good at aligning reads of
about 50 up to 100s or 1,000s of characters to relatively long (e.g. mammalian)
genomes.  Bowtie 2 indexes the genome with an [FM Index] (based on the
[Burrows-Wheeler Transform] or [BWT]) to keep its memory footprint small: for
the human genome, its memory footprint is typically around 3.2 gigabytes of RAM.
 Bowtie 2 supports gapped, local, and paired-end alignment modes.  Multiple
processors can be used simultaneously to achieve greater alignment speed. 
Bowtie 2 outputs alignments in [SAM] format, enabling interoperation with a
large number of other tools (e.g. [SAMtools], [GATK]) that use SAM.  Bowtie 2 is
distributed under the [GPLv3 license], and it runs on the command line under
Windows, Mac OS X and Linux.

[Bowtie 2] is often the first step in pipelines for comparative genomics,
including for variation calling, ChIP-seq, RNA-seq, BS-seq.  [Bowtie 2] and
[Bowtie] (also called "[Bowtie 1]" here) are also tightly integrated into some
tools, including [TopHat]: a fast splice junction mapper for RNA-seq reads,
[Cufflinks]: a tool for transcriptome assembly and isoform quantitiation from
RNA-seq reads, [Crossbow]: a cloud-enabled software tool for analyzing
reseuqncing data, and [Myrna]: a cloud-enabled software tool for aligning
RNA-seq reads and measuring differential gene expression.

If you use [Bowtie 2] for your published research, please cite the [Bowtie
paper].  Thank you!

[Bowtie 2]:        http://bowtie-bio.sf.net/bowtie2
[Bowtie]:          http://bowtie-bio.sf.net
[Bowtie 1]:        http://bowtie-bio.sf.net
[Burrows-Wheeler Transform]: http://en.wikipedia.org/wiki/Burrows-Wheeler_transform
[BWT]:             http://en.wikipedia.org/wiki/Burrows-Wheeler_transform
[FM Index]:        http://en.wikipedia.org/wiki/FM-index
[SAM]:             http://samtools.sourceforge.net/SAM1.pdf
[SAMtools]:        http://samtools.sourceforge.net
[GATK]:            http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit
[TopHat]:          http://tophat.cbcb.umd.edu/
[Cufflinks]:       http://cufflinks.cbcb.umd.edu/
[Crossbow]:        http://bowtie-bio.sf.net/crossbow
[Myrna]:           http://bowtie-bio.sf.net/myrna
[Bowtie paper]:    http://genomebiology.com/2009/10/3/R25
[GPLv3 license]:   http://www.gnu.org/licenses/gpl-3.0.html

How is Bowtie 2 different from Bowtie 1?
----------------------------------------

Bowtie 1 was released in 2009 and was geared toward aligning the relatively
short sequencing reads (up to 25-50 nucleotides) prevalent at the time. Since
then, technology has improved both sequencing throughput (more nucleotides
produced per sequencer per day) and read length (more nucleotides per read).

The chief differences between Bowtie 1 and Bowtie 2 are:

1. For reads longer than about 50 bp Bowtie 2 is generally faster, more
sensitive, and uses less memory than Bowtie 1.  For relatively short reads (e.g.
less than 50 bp) Bowtie 1 is sometimes faster and/or more sensitive.

2. Bowtie 2 supports gapped alignment with affine gap penalties. Number of gaps
and gap lengths are not restricted, except by way of the configurable scoring
scheme.  Bowtie 1 finds just ungapped alignments.

3. Bowtie 2 supports [local alignment], which doesn't require reads to align
end-to-end.  Local alignments might be "trimmed" ("soft clipped") at one or both
extremes in a way that optimizes alignment score. Bowtie 2 also supports
[end-to-end alignment] which, like Bowtie 1, requires that the read align
entirely.

4. There is no upper limit on read length in Bowtie 2.  Bowtie 1 had an upper
limit of around 1000 bp.

5. Bowtie 2 allows alignments to [overlap ambiguous characters] (e.g. `N`s) in
the reference.  Bowtie 1 does not.

6. Bowtie 2 does away with Bowtie 1's notion of alignment "stratum", and its
distinction between "Maq-like" and "end-to-end" modes.  In Bowtie 2 all
alignments lie along a continuous spectrum of alignment scores where the
[scoring scheme], similar to [Needleman-Wunsch] and [Smith-Waterman].

7. Bowtie 2's [paired-end alignment] is more flexible.  E.g. for pairs that do
not align in a paired fashion, Bowtie 2 attempts to find unpaired alignments for
each mate.

8. Bowtie 2 reports a spectrum of mapping qualities, in contrast fo Bowtie 1
which reports either 0 or high.

9. Bowtie 2 does not align colorspace reads.

Bowtie 2 is not a "drop-in" replacement for Bowtie 1.  Bowtie 2's command-line
arguments and genome index format are both different from Bowtie 1's.

[Needleman-Wunsch]: http://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm
[Smith-Waterman]:   http://en.wikipedia.org/wiki/Smith_waterman

What isn't Bowtie 2?
--------------------

Bowtie 1 and Bowtie 2 are not general-purpose alignment tools like [MUMmer],
[BLAST] or [Vmatch].  Bowtie 2 works best when aligning to large genomes, though
it supports arbitrarily small reference sequences (e.g. amplicons).  It handles
very long reads (i.e. upwards of 10s or 100s of kilobases), but it is optimized
for the read lengths and error modes yielded by recent sequencers, such as the
Illumina HiSeq 2000, Roche 454, and Ion Torrent instruments.

If your goal is to align two very large sequences (e.g. two genomes), consider
using [MUMmer].  If your goal is very sensitive alignment to a relatively short
reference sequence (e.g. a bacterial genome), this can be done with Bowtie 2 but
you may want to consider using tools like [NUCmer], [BLAT], or [BLAST].  These
tools can be extremely slow when the reference genome is long, but are often
adequate when the reference is short.

Bowtie 2 does not support alignment of colorspace reads.  This might be
supported in future versions.

[MUMmer]: http://mummer.sourceforge.net/
[NUCmer]: http://mummer.sourceforge.net/manual/#nucmer
[BLAST]:  http://blast.ncbi.nlm.nih.gov/Blast.cgi
[BLAT]:   http://genome.ucsc.edu/cgi-bin/hgBlat?command=start
[Vmatch]: http://www.vmatch.de/

What does it mean that some older Bowtie 2 versions are "beta"?
--------------------------------------------------------------

We said those Bowtie 2 versions were in "beta" to convey that it was not as
polished as a tool that had been around for a while, and was still in flux.
Since version 2.0.1, we declared Bowtie 2 was no longer "beta".

Obtaining Bowtie 2
==================

Download Bowtie 2 sources and binaries from the [Download] section of the
Sourceforge site.  Binaries are available for Intel architectures (`i386` and
`x86_64`) running Linux, and Mac OS X.  A 32-bit version is available for
Windows.  If you plan to compile Bowtie 2 yourself, make sure to get the source
package, i.e., the filename that ends in "-source.zip".

Building from source
--------------------

Building Bowtie 2 from source requires a GNU-like environment with GCC, GNU Make
and other basics.  It should be possible to build Bowtie 2 on most vanilla Linux
installations or on a Mac installation with [Xcode] installed.  Bowtie 2 can
also be built on Windows using [Cygwin] or [MinGW] (MinGW recommended). For a 
MinGW build the choice of what compiler is to be used is important since this
will determine if a 32 or 64 bit code can be successfully compiled using it. If 
there is a need to generate both 32 and 64 bit on the same machine then a multilib 
MinGW has to be properly installed. [MSYS], the [zlib] library, and depending on 
architecture [pthreads] library are also required. We are recommending a 64 bit
build since it has some clear advantages in real life research problems. In order 
to simplify the MinGW setup it might be worth investigating popular MinGW personal 
builds since these are coming already prepared with most of the toolchains needed.

First, download the source package from the [sourceforge site].  Make sure
you're getting the source package; the file downloaded should end in
`-source.zip`. Unzip the file, change to the unzipped directory, and build the
Bowtie 2 tools by running GNU `make` (usually with the command `make`, but
sometimes with `gmake`) with no arguments.  If building with MinGW, run `make`
from the MSYS environment.

Bowtie 2 is using the multithreading software model in order to speed up 
execution times on SMP architectures where this is possible. On POSIX 
platforms (like linux, Mac OS, etc) it needs the pthread library. Although
it is possible to use pthread library on non-POSIX platform like Windows, due
to performance reasons bowtie 2 will try to use Windows native multithreading
if possible.

[Cygwin]:   http://www.cygwin.com/
[MinGW]:    http://www.mingw.org/
[MSYS]:     http://www.mingw.org/wiki/msys
[zlib]:     http://cygwin.com/packages/mingw-zlib/
[pthreads]: http://sourceware.org/pthreads-win32/
[GnuWin32]: http://gnuwin32.sf.net/packages/coreutils.htm
[Download]: https://sourceforge.net/projects/bowtie-bio/files/bowtie2/
[sourceforge site]: https://sourceforge.net/projects/bowtie-bio/files/bowtie2/
[Xcode]:    http://developer.apple.com/xcode/

Adding to PATH
--------------

By adding your new Bowtie 2 directory to your [PATH environment variable], you
ensure that whenever you run `bowtie2`, `bowtie2-build` or `bowtie2-inspect`
from the command line, you will get the version you just installed without
having to specify the entire path.  This is recommended for most users.  To do
this, follow your operating system's instructions for adding the directory to
your [PATH].

If you would like to install Bowtie 2 by copying the Bowtie 2 executable files
to an existing directory in your [PATH], make sure that you copy all the
executables, including `bowtie2`, `bowtie2-align`, `bowtie2-build` and
`bowtie2-inspect`.

[PATH environment variable]: http://en.wikipedia.org/wiki/PATH_(variable)
[PATH]: http://en.wikipedia.org/wiki/PATH_(variable)

The `bowtie2` aligner
=====================

`bowtie2` takes a Bowtie 2 index and a set of sequencing read files and outputs
a set of alignments in SAM format.

"Alignment" is the process by which we discover how and where the read sequences
are similar to the reference sequence.  An "alignment" is a result from this
process, specifically: an alignment is a way of "lining up" some or all of the
characters in the read with some characters from the reference in a way that
reveals how they're similar.  For example:

      Read:      GACTGGGCGATCTCGACTTCG
                 |||||  |||||||||| |||
      Reference: GACTG--CGATCTCGACATCG

Where dash symbols represent gaps and vertical bars show where aligned
characters match.

We use alignment to make an educated guess as to where a read originated with
respect to the reference genome.  It's not always possible to determine this
with certainty.  For instance, if the reference genome contains several long
stretches of As (`AAAAAAAAA` etc) and the read sequence is a short stretch of As
(`AAAAAAA`), we cannot know for certain exactly where in the sea of `A`s the
read originated.

End-to-end alignment versus local alignment
-------------------------------------------

By default, Bowtie 2 performs end-to-end read alignment.  That is, it searches
for alignments involving all of the read characters.  This is also called an
"untrimmed" or "unclipped" alignment.

When the --local option is specified, Bowtie 2 performs local read alignment. In
this mode, Bowtie 2 might "trim" or "clip" some read characters from one or both
ends of the alignment if doing so maximizes the alignment score.

### End-to-end alignment example

The following is an "end-to-end" alignment because it involves all the
characters in the read.  Such an alignment can be produced by Bowtie 2 in either
end-to-end mode or in local mode.

    Read:      GACTGGGCGATCTCGACTTCG
    Reference: GACTGCGATCTCGACATCG
    
    Alignment:
      Read:      GACTGGGCGATCTCGACTTCG
                 |||||  |||||||||| |||
      Reference: GACTG--CGATCTCGACATCG

### Local alignment example

The following is a "local" alignment because some of the characters at the ends
of the read do not participate.  In this case, 4 characters are omitted (or
"soft trimmed" or "soft clipped") from the beginning and 3 characters are
omitted from the end.  This sort of alignment can be produced by Bowtie 2 only
in local mode.

    Read:      ACGGTTGCGTTAATCCGCCACG
    Reference: TAACTTGCGTTAAATCCGCCTGG
    
    Alignment:
      Read:      ACGGTTGCGTTAA-TCCGCCACG
                     ||||||||| ||||||
      Reference: TAACTTGCGTTAAATCCGCCTGG

Scores: higher = more similar
-----------------------------

An alignment score quantifies how similar the read sequence is to the reference
sequence aligned to.  The higher the score, the more similar they are.  A score
is calculated by subtracting penalties for each difference (mismatch, gap, etc)
and, in local alignment mode, adding bonuses for each match.

The scores can be configured with the `--ma` (match bonus), `--mp` (mismatch
penalty), `--np` (penalty for having an N in either the read or the
reference), `--rdg` (affine read gap penalty) and `--rfg` (affine reference
gap penalty) options.

### End-to-end alignment score example

A mismatched base at a high-quality position in the read receives a penalty of
-6 by default.  A length-2 read gap receives a penalty of -11 by default (-5 for
the gap open, -3 for the first extension, -3 for the second extension).  Thus,
in end-to-end alignment mode, if the read is 50 bp long and it matches the
reference exactly except for one mismatch at a high-quality position and one
length-2 read gap, then the overall score is -(6 + 11) = -17.

The best possible alignment score in end-to-end mode is 0, which happens when
there are no differences between the read and the reference.

### Local alignment score example

A mismatched base at a high-quality position in the read receives a penalty of
-6 by default.  A length-2 read gap receives a penalty of -11 by default (-5 for
the gap open, -3 for the first extension, -3 for the second extension).  A base
that matches receives a bonus of +2 be default.  Thus, in local alignment mode,
if the read is 50 bp long and it matches the reference exactly except for one
mismatch at a high-quality position and one length-2 read gap, then the overall
score equals the total bonus, 2 * 49, minus the total penalty, 6 + 11, = 81.

The best possible score in local mode equals the match bonus times the length of
the read.  This happens when there are no differences between the read and the
reference.

### Valid alignments meet or exceed the minimum score threshold

For an alignment to be considered "valid" (i.e. "good enough") by Bowtie 2, it
must have an alignment score no less than the minimum score threshold.  The
threshold is configurable and is expressed as a function of the read length. In
end-to-end alignment mode, the default minimum score threhsold is `-0.6 + -0.6 *
L`, where `L` is the read length.  In local alignment mdoe, the default minimum
score threshold is `20 + 8.0 * ln(L)`, where L is the read length.  This can be
configured with the `--score-min` option.  For details on how to set options
like `--score-min` that correpond to functions, see the section on [setting
function options].

Mapping quality: higher = more unique
-------------------------------------

The aligner cannot always assign a read to its point of origin with high
confidence.  For instance, a read that originated inside a repeat element might
align equally well to many occurrences of the element throughout the genome,
leaving the aligner with no basis for preferring one over the others.

Aligners characterize their degree of confidence in the point of origin by
reporting a mapping quality: a non-negative integer Q = -10 log10 p, where p is
an estimate of the probability that the alignment does not correspond to the
read's true point of origin.  Mapping quality is sometimes abbreviated MAPQ, and
is recorded in the [SAM] `MAPQ` field.

Mapping quality is related to "uniqueness."  We say an alignment is unique if it
has a much higher alignment score than all the other possible alignments. The
bigger the gap between the best alignment's score and the second-best
alignment's score, the more unique the best alignment, and the higher its mapping
quality should be.

Accurate mapping qualities are useful for downstream tools like variant callers.
For instance, a variant caller might choose to ignore evidence from alignments
with mapping quality less than, say, 10.  A mapping quality of 10 or less
indicates that there is at least a 1 in 10 chance that the read truly originated
elsewhere.

[SAM]: http://samtools.sourceforge.net/SAM1.pdf

Aligning pairs
--------------

A "paired-end" or "mate-pair" read consists of pair of mates, called mate 1 and
mate 2.  Pairs come with a prior expectation about (a) the relative orientation
of the mates, and (b) the distance separating them on the original DNA molecule.
Exactly what expectations hold for a given dataset depends on the lab procedures
used to generate the data.  For example, a common lab procedure for producing
pairs is Illumina's Paired-end Sequencing Assay, which yields pairs with a
relative orientation of FR ("forward, reverse") meaning that if mate 1 came from
the Watson strand, mate 2 very likely came from the Crick strand and vice versa.
 Also, this protocol yields pairs where the expected genomic distance from end
to end is about 200-500 base pairs.

For simplicity, this manual uses the term "paired-end" to refer to any pair of
reads with some expected relative orientation and distance.  Depending on the
protocol, these might actually be referred to as "paired-end" or "mate-paired."
Also, we always refer to the individual sequences making up the pair as "mates."

### Paired inputs

Pairs are often stored in a pair of files, one file containing the mate 1s and
the other containing the mates 2s.  The first mate in the file for mate 1 forms
a pair with the first mate in the file for mate 2, the second with the second,
and so on.  When aligning pairs with Bowtie 2, specify the file with the mate 1s
mates using the `-1` argument and the file with the mate 2s using the `-2`
argument.  This causes Bowtie 2 to take the paired nature of the reads into
account when aligning them.

### Paired SAM output

When Bowtie 2 prints a SAM alignment for a pair, it prints two records (i.e. two
lines of output), one for each mate.  The first record describes the alignment
for mate 1 and the second record describes the alignment for mate 2.  In both
records, some of the fields of the SAM record describe various properties of the
alignment; for instance, the 7th and 8th fields (`RNEXT` and `PNEXT`
respectively) indicate the reference name and position where the other mate
aligned, and the 9th field indicates the inferred length of the DNA fragment
from which the two mates were sequenced.  See the [SAM specification] for more
details regarding these fields.

### Concordant pairs match pair expectations, discordant pairs don't

A pair that aligns with the expected relative mate orientation and with the
expected range of distances between mates is said to align "concordantly".  If
both mates have unique alignments, but the alignments do not match paired-end
expectations (i.e. the mates aren't in the expcted relative orientation, or
aren't within the expected disatance range, or both), the pair is said to align
"discordantly".  Discordant alignments may be of particular interest, for
instance, when seeking [structural variants].

The expected relative orientation of the mates is set using the `--ff`,
`--fr`, or `--rf` options.  The expected range of inter-mates distances (as
measured from the furthest extremes of the mates; also called "outer distance")
is set with the `-I` and `-X` options.

To declare that a pair aligns discordantly, Bowtie 2 requires that both mates
align uniquely.  This is a conservative threshold, but this is often desirable
when seeking structural variants.

By default, Bowtie 2 searches for both concordant and discordant alignments,
though searching for discordant alignments can be disabled with the
`--no-discordant` option.

[structural variants]: http://www.ncbi.nlm.nih.gov/dbvar/content/overview/

### Mixed mode: paired where possible, unpaired otherwise

If Bowtie 2 cannot find a paired-end alignment for a pair, by default it will go
on to look for unpaired alignments for the constituent mates.  This is called
"mixed mode."  To disable mixed mode, set the `--no-mixed` option.

Bowtie 2 runs a little faster in `--no-mixed` mode, but will only consider
alignment status of pairs per se, not individual mates.

### Some SAM FLAGS describe paired-end properties

The SAM `FLAGS` field, the second field in a SAM record, has multiple bits that
describe the paired-end nature of the read and alignment.  The first (least
significant) bit (1 in decimal, 0x1 in hexidecimal) is set if the read is part
of a pair.  The second bit (2 in decimal, 0x2 in hexidecimal) is set if the read
is part of a pair that aligned in a paired-end fashion.  The fourth bit (8 in
decimal, 0x8 in hexidecimal) is set if the read is part of a pair and the other
mate in the pair had at least one valid alignment.  The sixth bit (32 in
decimal, 0x20 in hexidecimal) is set if the read is part of a pair and the other
mate in the pair aligned to the Crick strand (or, equivalently, if the reverse
complement of the other mate aligned to the Watson strand).  The seventh bit (64
in decimal, 0x40 in hexidecimal) is set if the read is mate 1 in a pair.  The
eighth bit (128 in decimal, 0x80 in hexidecimal) is set if the read is mate 2 in
a pair.  See the [SAM specification] for a more detailed description of the
`FLAGS` field.

### Some SAM optional fields describe more paired-end properties

The last severeal fields of each SAM record usually contain SAM optional fields,
which are simply tab-separated strings conveying additional information about
the reads and alignments.  A SAM optional field is formatted like this: "XP:i:1"
where "XP" is the `TAG`, "i" is the `TYPE` ("integer" in this case), and "1" is
the `VALUE`.  See the [SAM specification] for details regarding SAM optional
fields.

### Mates can overlap, contain, or dovetail each other

The fragment and read lengths might be such that alignments for the two mates
from a pair overlap each other.  Consider this example:

(For these examples, assume we expect mate 1 to align to the left of mate 2.)

    Mate 1:    GCAGATTATATGAGTCAGCTACGATATTGTT
    Mate 2:                               TGTTTGGGGTGACACATTACGCGTCTTTGAC
    Reference: GCAGATTATATGAGTCAGCTACGATATTGTTTGGGGTGACACATTACGCGTCTTTGAC

It's also possible, though unusual, for one mate alignment to contain the other,
as in these examples:

    Mate 1:    GCAGATTATATGAGTCAGCTACGATATTGTTTGGGGTGACACATTACGC
    Mate 2:                               TGTTTGGGGTGACACATTACGC
    Reference: GCAGATTATATGAGTCAGCTACGATATTGTTTGGGGTGACACATTACGCGTCTTTGAC

    Mate 1:                   CAGCTACGATATTGTTTGGGGTGACACATTACGC
    Mate 2:                      CTACGATATTGTTTGGGGTGAC
    Reference: GCAGATTATATGAGTCAGCTACGATATTGTTTGGGGTGACACATTACGCGTCTTTGAC

And it's also possible, though unusual, for the mates to "dovetail", with the
mates seemingly extending "past" each other as in this example:

    Mate 1:                 GTCAGCTACGATATTGTTTGGGGTGACACATTACGC
    Mate 2:            TATGAGTCAGCTACGATATTGTTTGGGGTGACACAT                   
    Reference: GCAGATTATATGAGTCAGCTACGATATTGTTTGGGGTGACACATTACGCGTCTTTGAC

In some situations, it's desirable for the aligner to consider all these cases
as "concordant" as long as other paired-end constraints are not violated. Bowtie
2's default behavior is to consider overlapping and containing as being
consistent with concordant alignment.  By default, dovetailing is considered
inconsistent with concordant alignment.

These defaults can be overridden.  Setting `--no-overlap` causes Bowtie 2 to
consider overlapping mates as non-concordant.  Setting `--no-contain` causes
Bowtie 2 to consider cases where one mate alignment contains the other as
non-concordant. Setting `--dovetail` causes Bowtie 2 to consider cases where
the mate alignments dovetail as concordant.

Reporting
---------

The reporting mode governs how many alignments Bowtie 2 looks for, and how to
report them.  Bowtie 2 has three distinct reporting modes.  The default
reporting mode is similar to the default reporting mode of many other read
alignment tools, including [BWA].  It is also similar to Bowtie 1's `-M`
alignment mode.

In general, when we say that a read has an alignment, we mean that it has a
[valid alignment].  When we say that a read has multiple alignments, we mean
that it has multiple alignments that are valid and distinct from one another. 

[BWA]: http://bio-bwa.sourceforge.net/

### Distinct alignments map a read to different places

Two alignments for the same individual read are "distinct" if they map the same
read to different places.  Specifically, we say that two alignments are distinct
if there are no alignment positions where a particular read offset is aligned
opposite a particular reference offset in both alignments with the same
orientation.  E.g. if the first alignment is in the forward orientation and
aligns the read character at read offset 10 to the reference character at
chromosome 3, offset 3,445,245, and the second alignment is also in the forward
orientation and also aligns the read character at read offset 10 to the
reference character at chromosome 3, offset 3,445,245, they are not distinct
alignments.

Two alignments for the same pair are distinct if either the mate 1s in the two
paired-end alignments are distinct or the mate 2s in the two alignments are
distinct or both.

### Default mode: search for multiple alignments, report the best one

By default, Bowtie 2 searches for distinct, valid alignments for each read. When
it finds a valid alignment, it generally will continue to look for alignments
that are nearly as good or better.  It will eventually stop looking, either
because it exceeded a limit placed on search effort (see `-D` and `-R`) or
because it already knows all it needs to know to report an alignment.
Information from the best alignments are used to estimate mapping quality (the
`MAPQ` [SAM] field) and to set SAM optional fields, such as `AS:i` and
`XS:i`.  Bowtie 2 does not gaurantee that the alignment reported is the best
possible in terms of alignment score.

See also: `-D`, which puts an upper limit on the number of dynamic programming
problems (i.e. seed extensions) that can "fail" in a row before Bowtie 2 stops
searching.  Increasing `-D` makes Bowtie 2 slower, but increases the
likelihood that it will report the correct alignment for a read that aligns many
places.

See also: `-R`, which sets the maximum number of times Bowtie 2 will "re-seed"
when attempting to align a read with repetitive seeds.  Increasing `-R` makes
Bowtie 2 slower, but increases the likelihood that it will report the correct
alignment for a read that aligns many places.

### -k mode: search for one or more alignments, report each

In `-k` mode, Bowtie 2 searches for up to N distinct, valid alignments for
each read, where N equals the integer specified with the `-k` parameter.  That
is, if `-k 2` is specified, Bowtie 2 will search for at most 2 distinct
alignments.  It reports all alignments found, in descending order by alignment
score.  The alignment score for a paired-end alignment equals the sum of the
alignment scores of the individual mates.  Each reported read or pair alignment
beyond the first has the SAM 'secondary' bit (which equals 256) set in its FLAGS
field.  See the [SAM specification] for details.

Bowtie 2 does not "find" alignments in any specific order, so for reads that
have more than N distinct, valid alignments, Bowtie 2 does not gaurantee that
the N alignments reported are the best possible in terms of alignment score.
Still, this mode can be effective and fast in situations where the user cares
more about whether a read aligns (or aligns a certain number of times) than
where exactly it originated.

[SAM specification]: http://samtools.sourceforge.net/SAM1.pdf

### -a mode: search for and report all alignments

`-a` mode is similar to `-k` mode except that there is no upper limit on the
number of alignments Bowtie 2 should report.  Alignments are reported in
descending order by alignment score.  The alignment score for a paired-end
alignment equals the sum of the alignment scores of the individual mates.  Each
reported read or pair alignment beyond the first has the SAM 'secondary' bit
(which equals 256) set in its FLAGS field.  See the [SAM specification] for
details.

Some tools are designed with this reporting mode in mind.  Bowtie 2 is not!  For
very large genomes, this mode is very slow.

[SAM specification]: http://samtools.sourceforge.net/SAM1.pdf

### Randomness in Bowtie 2

Bowtie 2's search for alignments for a given read is "randomized."  That is,
when Bowtie 2 encouters a set of equally-good choices, it uses a pseudo-random
number to choose.  For example, if Bowtie 2 discovers a set of 3 equally-good
alignments and wants to decide which to report, it picks a pseudo-random integer
0, 1 or 2 and reports the corresponding alignment.  Abitrary choices can crop up
at various points during alignment.

The pseudo-random number generator is re-initialized for every read, and the
seed used to initialize it is a function of the read name, nucleotide string,
quality string, and the value specified with `--seed`.  If you run the same
version of Bowtie 2 on two reads with identical names, nucleotide strings, and
quality strings, and if `--seed` is set the same for both runs, Bowtie 2 will
produce the same output; i.e., it will align the read to the same place, even if
there are multiple equally good alignments.  This is intuitive and desirable in
most cases.  Most users expect Bowtie to produce the same output when run twice
on the same input.

However, when the user specifies the `--non-deterministic` option, Bowtie 2
will use the current time to re-intiailize the pseud-random number generator.
When this is specified, Bowtie 2 might report different alignments for identical
reads.  This is counter-intuitive for some users, but might be more appropriate
in situations where the input consists of many identical reads.

Multiseed heuristic
-------------------

To rapidly narrow the number of possible alignments that must be considered,
Bowtie 2 begins by extracting substrings ("seeds") from the read and its reverse
complement and aligning them in an ungapped fashion with the help of the [FM
Index].  This is "multiseed alignment" and it is similar to what [Bowtie 1
does], except Bowtie 1 attempts to align the entire read this way.

This initial step makes Bowtie 2 much faster than it would be without such a
filter, but at the expense of missing some valid alignments.  For instance, it
is possible for a read to have a valid overall alignment but to have no valid
seed alignments because each potential seed alignment is interruped by too many
mismatches or gaps.

The tradeoff between speed and sensitivity/accuracy can be adjusted by setting
the seed length (`-L`), the interval between extracted seeds (`-i`), and the
number of mismatches permitted per seed (`-N`).  For more sensitive alignment,
set these parameters to (a) make the seeds closer together, (b) make the seeds
shorter, and/or (c) allow more mismatches.  You can adjust these options
one-by-one, though Bowtie 2 comes with some useful combinations of options
pre-packaged as "[preset options]."

`-D` and `-R` are also options that adjust the tradeoff between speed and
sensitivity/accuracy.

### FM Index memory footprint

Bowtie 2 uses the [FM Index] to find ungapped alignments for seeds.  This step
accounts for the bulk of Bowtie 2's memory footprint, as the [FM Index] itself
is typically the largest data structure used.  For instance, the memory
footprint of the [FM Index] for the human genome is about 3.2 gigabytes of RAM.

[Bowtie 1 does]: http://genomebiology.com/2009/10/3/R25
[Bowtie 1 paper]: http://genomebiology.com/2009/10/3/R25
[FM Index]: http://portal.acm.org/citation.cfm?id=796543
[bi-directional BWT approach]: http://www.computer.org/portal/web/csdl/doi/10.1109/BIBM.2009.42

Ambiguous characters
--------------------

Non-whitespace characters besides A, C, G or T are considered "ambiguous."  N is
a common ambiguous character that appears in reference sequences.  Bowtie 2
considers all ambiguous characters in the reference (including [IUPAC nucleotide
codes]) to be Ns.

Bowtie 2 allows alignments to overlap ambiguous characters in the reference. An
alignment position that contains an ambiguous character in the read, reference,
or both, is penalized according to `--np`.  `--n-ceil` sets an upper limit
on the number of positions that may contain ambiguous reference characters in a
valid alignment.  The optional field `XN:i` reports the number of ambiguous
reference characters overlapped by an alignment.

Note that the [multiseed heuristic] cannot find *seed* alignments that overlap
ambiguous reference characters.  For an alignment overlapping an ambiguous
reference character to be found, it must have one or more seed alignments that
do not overlap ambiguous reference characters.

[IUPAC nucleotide codes]: http://www.bioinformatics.org/sms/iupac.html

Presets: setting many settings at once
--------------------------------------

Bowtie 2 comes with some useful combinations of parameters packaged into shorter
"preset" parameters.  For example, running Bowtie 2 with the `--very-sensitive`
option is the same as running with options: `-D 20 -R 3 -N 0 -L 20 -i S,1,0.50`.
 The preset options that come with Bowtie 2 are designed to cover a wide area of
the speed/sensitivity/accuracy tradeoff space, with the presets ending in `fast`
generally being faster but less sensitive and less accurate, and the presets
ending in `sensitive` generally being slower but more sensitive and more
accurate.  See the [documentation for the preset options] for details.

Filtering
---------

Some reads are skipped or "filtered out" by Bowtie 2.  For example, reads may be
filtered out because they are extremely short or have a high proportion of
ambiguous nucleotides.  Bowtie 2 will still print a SAM record for such a read,
but no alignment will be reported and and the `YF:i` SAM optional field will be
set to indicate the reason the read was filtered.

* `YF:Z:LN`: the read was filtered becuase it had length less than or equal to
the number of seed mismatches set with the `-N` option.
* `YF:Z:NS`: the read was filtered because it contains a number of ambiguous
characters (usually `N` or `.`) greater than the ceiling specified with
`--n-ceil`.
* `YF:Z:SC`: the read was filtered because the read length and the match bonus
(set with `--ma`) are such that the read can't possibly earn an alignment
score gr