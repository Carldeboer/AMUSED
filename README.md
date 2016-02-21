# AMUSED
Auditing Motifs Using Statistical Enrichment &amp; Depletion

## Installation
Installation on *nix (including macs) should be straightforward.  Please contact me if you run into trouble.

### Requirements:
* Ruby (https://www.ruby-lang.org/en/documentation/installation/)
* gzip (if on *nix, and if you want to input/output to .gz compressed files)
* zdiff (for running unit tests)

### Unix (and mac)
Download the requisite files where you want them stored:
```
cd /path/to/where/you/want/the/files/to/be
git clone https://github.com/Carldeboer/AMUSED
git clone https://github.com/Carldeboer/Ruby-DNA-Tools
```
If you don't have git, you can download these repositories manually by clicking "Download ZIP" at the top right for both AMUSED and Ruby-DNA-Tools and saving these to an appropriate place.

Add  `Ruby-DNA-Tools` and`AMUSED` to your ruby path so that Ruby can find the required libraries:
```
export RUBYLIB=$RUBYLIB:/path/to/where/you/want/the/files/to/be/AMUSED:/path/to/where/you/want/the/files/to/be/Ruby-DNA-Tools
```
To make this persist between sessions, add the above line to your login script (e.g. ~/.bashrc, ~/.bash_profile, etc).

Now do one of the following two options:

(1) link the executable files into your bin:
```bash
ln -s /path/to/where/you/want/the/files/to/be/AMUSED/AMUSED ~/bin/AMUSED
ln -s /path/to/where/you/want/the/files/to/be/AMUSED/AMUSED-KS ~/bin/AMUSED-KS
ln -s /path/to/where/you/want/the/files/to/be/AMUSED/alignKMers ~/bin/alignKMers
```
or (2) add AMUSED to your PATH in your login script:
```bash
export PATH=$PATH:/path/to/where/you/want/the/files/to/be/AMUSED/
```
Note: if you choose the login-script option, you must re-login to re-run the script before your shell will be able to find AMUSED.

To test the installation, you can run the test script `runTests.bat` in the AMUSED directory:
```bash
cd /path/to/where/you/want/the/files/to/be/AMUSED/
./runTests.bat
```
which should output something along the lines of:
```
Scanning foreground sequences
Input All Data, All Trees Built.
Adding gaps...
Calculating statistics...
Sorting results...
Printing results...
All Done!
test1:  0  differences (should be 0)

...

Inputing All Data...
Making BG...
Calculating statistics...
Sorting results...
Calculating P-values...
Printing results...
Done!
test5:  0  differences (should be 0)
```
All tests should show 0 differences, e.g.
```
test2:  0  differences (should be 0)
```
`No such file or directory` exceptions indicate that AMUSED did not complete execution for some reason (e.g. missing libraries).  The relevant error message should be printed. Note also that these tests require `zdiff`.


### Windows
Windows is not officially supported and I have never tried it.  Theoretically, you should be able to either put all the ruby libraries (*.rb) and AMUSED executables in one directory and run it from there.  Let me know if it is especially problematic for whatever reason and I will look into it.  Just don't try to use the gzip reading/writing functionality.

## Usage
Run any of the included programs with no parameters to see usage information (also below).
### AMUSED
This program identifies enriched and depleted (gapped) k-mers in your set of query sequences, optionally comparing to a set of background sequences.
```
Usage: AMUSED -q <inFPQuery>  -o <outFP> [-b <inFPBG> | -r <randomizeNMer>] [-s <maxTreeSize>]
  -q <inFPQuery> = query sequences
  -b <inFPBG> = compare seqs to these background seqs
  -bp <bgPseudo> = pseudocount to add to background [default=0.5]
  -o <outFP> = output file
  -s <maxTreeSize> = max n-mer to consider [default=8]
  -z <subZCutoff> = minimum absolute Sub-Z-score [default = 0; print all]
  -t <numThreads> = number of CPU threads to use [default=1]
  -1p = sequences not in fasta format: each line is a full sequence
  -ng = no inserting gaps
  -nu = no changing to upper case before scan (non ATGC bases are discarded)
  -ds = double stranded (reverse complement sequences too)
  -ns = don't sort
  -do = descriptive output: lots of intermediate values also output (but many columns)
  -bc = add lines to output for base content
  -nsz = don't calculate super Zs
```
* The following parameters are required: `-q -o`.
* `-q` is used to specify the file containing the foreground (query) set of sequences. These are generally in fasta format (or gzipped fasta), although they can also be one sequence per line if `-1p` is also given.
* `-b` is used to specify the file containing the background set of sequences as above.
* `-bp` is the pseudocount added to the background k-mer counts to avoid cases where a k-mer is observed in the foreground but absent from the background. I do not recommend altering this parameter and for large amounts of sequence it should make very little difference.
* `-o` is used to specify the output file. If it ends in `.gz` the output will be compressed (gzip).  The output file is tab delimited and has a header (example below).
* `-s` specifies the maximum k-mer length k to consider.  The gap addition step is very slow for large k, so k<10 is recommended. If super-Z scores are being calculated, (k+1)-mers are also quantified and used to calculate Super-Z scores.
* `-z` you can specify with this parameter how significant the Sub-Z scores must be before they are reported.  By default, all are reported.
* `-t` specifies the number of threads to use. In general, AMUSED will only make use of two threads and this will only speed up the scanning step, and only if using a background.
* `-1p` use this if your input sequences are not in fasta format, but simply have one sequence per line.
* `-ng` use this if you want ungapped k-mers only
* `-nu` use this if your input sequences contain lower case letters and you don't want them to be scanned (e.g. repeat masked)
* `-ds` use this if you want both strands to be scanned so that strand-specific enrichments will not be looked for. **AMUSED is by default single-stranded and so will only scan one strand without this parameter set.**
* `-ns` don't sort the output file - it will instead be presented in k-mer order.
* `-do` more verbose output file including several of the intermediate statistics calculated but not usually provided.  This can be useful for debugging.
* `-bc` include base content (e.g. A, T, G, C frequencies) to output.  This can be useful for completeness, but no other statistics will be calculated for these as Sub-Z scores are meaningless for a k-mer of length 1.
* `-nsz` don't calculate super-Z scores

#### Example output
```
k       kMer    kMerRC  Num_FG  Size_FG Pct_FG  SubEnr_FG       SubDepl_FG      SubZ_FG Num_BG  Size_BG Pct_BG  logOR   SubEnr_BG       SubDepl_BG      SubZ_BG SubEnr_Rel      SubDepl_Rel     SubZ_Rel        SuperZ_FG       DeltaZ_FG   SuperZ_BG        SuperZ_Rel      DeltaZ_Rel
2       CT      AG      175     2446    7.1545380212592 1.53300948370107        1.53300948370107        5.69484017223451        710     10250   6.92682926829268        0.0497639818145644      0.987227374958513       0.987227374958513   -0.342531678449607       1.55282941211577        1.55282941211577        5.83130689640249        5.35514590667512        0.339694265559392       -3.73300144891046       3.71423500185696        2.11707189454553
4       TGGA    TCCA    35      2438    1.43560295324036        2.48727815039435        2.50118488903371        5.57909867359366        41      10242   0.400312438976762       1.84543512128644        0.946629886222897       1.44483902624608     0.0     2.62581905910962        1.74039457657653        5.57909867359366        3.01586623529439        2.56323243829927        0.0     1.70149810057736        3.87760057301631
2       TG      CA      208     2446    8.50367947669665        1.57550886038865        1.57550886038865        6.61261528784288        627     10250   6.11707317073171        0.477559525359313       1.09185846386172        1.09185846386172     2.20125206504686        1.44306632253758        1.44306632253759        5.55715858811527        5.07123185033956        1.54138343750332        3.62822045662428        4.31740146922982        1.23975711888545
4       TG-A    T-CA    77      2438    3.15832649712879        1.80187578059251        2.39905941796657        5.24191573322697        173     10242   1.68912321812146        0.908059985824659       0.952900803167581       0.98144315814809     -0.246373811467486      1.89068075228696        2.44428941331719        5.36322305241131        5.10303465536105        0.138881077865918       -1.48725044211292       4.83543013245826        0.527792919953058
4       T--A    T--A    123     2438    5.04511894995898        1.5455408031606 1.5455408031606 4.8667558196567 832     10242   8.12341339582113        -0.682211797245265      0.945877785584276       0.945877785584276       -1.60516509651097    1.6339220432293 1.6339220432293 5.34957873675918        5.41803467649389        -0.551278856837189      -2.88431542163064       4.55747805698061        0.792100679778571
...
```
Note that normally not all of these are present - to see them all, use the `-do` option
* `k` is length of the k-mer
* `kMer` is the k-mer
* `kMerRC` is the reverse complement of the k-mer
* `Num_FG` is the number of occurrences of the k-mer in the foreground
* `Size_FG` is the number of k-mers of this length in the foreground
* `Pct_FG` is the percent of k-ners of this length that are this k-mer (as a percent)
* `SubEnr_FG` is enrichment of this k-mer relative to k-1-mers
* `SubDepl_FG` is the depletion of this k-mer relative to k-1-mers
* `SubZ_FG` is the Sub-Z score of this k-mer
* `Num_BG` as above for the background
* `Size_BG` as above for the background 
* `Pct_BG` as above for the background 
* `logOR` is the log2 odds ratio of seeing the k-mer in the foreground vs. the background
* `SubEnr_BG` as above for the background 
* `SubDepl_BG` as above for the background 
* `SubZ_BG` as above for the background 
* `SubEnr_Rel` is the relative enrichment of the k-mer in the foreground vs. the background.
* `SubDepl_Rel` is the relative depletion of the k-mer in the foreground vs. the background.
* `SubZ_Rel` is the relative Sub-Z score of the k-mer in the foreground vs.  the background.
* `SuperZ_FG` is the Super-Z score of the k-mer in the foreground. 
* `DeltaZ_FG` is the difference in Z scores in the foreground (Sub-Z - Super-Z if Sub-Z>0, and Super-Z - Sub-Z if Sub-Z<0) such that positive values indicate that the k-mer is more enriched (or depleted) than its k+1-mers  
* `SuperZ_BG` is as above for the background 
* `SuperZ_Rel` is the relative super-Z score for the foreground vs. background.
* `DeltaZ_Rel` is the difference in relative Z scores (Sub-Z - Super-Z if Sub-Z>0, and Super-Z - Sub-Z if Sub-Z<0) such that positive values indicate that the k-mer is more enriched (or depleted) than its k+1-mers, relative to the background.
 
#####Tips:
Look for k-mers with extreme Sub-Z scores and positive DeltaZ scores.  Extreme positive Sub-Z scores indicate that the k-mer is enriched while extreme negative Sub-Z scores indicate that the k-mer is depleted (relative to the background for SubZ_Rel, or internally for SubZ_FG). Positive DeltaZ scores indicate that the k-mer is more enriched (or depleted) than the k+1-mers containing it.  That is to say, adding another base of information does not add significantly to the enrichment.  You should also consider those with negative DeltaZ scores, but whose k is the maximum considered since the k+1-mer that contains more information is not represented in your results. For example, if k=8 and ATGGATAA is enriched, but not as much as ATGGATAAC the DeltaZ score will be negative, but ATGGATAAC will not be in the presented results because it is a 9-mer.

### AMUSED-KS
This program uses the K-S test to identify k-mers that are positionally biased in a set of sequences.
```
Usage: AMUSED-KS -q <inFPQuery>  -o <outFP> [-b <backgroundSeqs>] [-s <maxTreeSize>] [-m <smoothBy>] [-t <numThreads] [-p] [-1p]
  -q <inFPQuery> = query sequences
  -b <inFPBG> = compare seqs to these background seqs
  -o <outFP> = output file
  -s <maxTreeSize> = max n-mer to consider [default=8]
  -m <smoothBy> = smooth data over how many bases to generate expected distribution for one-sample KS test (when no background set is included) [default=5]
  -t <numThreads> = number of CPU threads to use [default=1]
  -dsk = double stranded k-mer scan
  -rc = also scan reverse complement input sequences (overrides -dsk)
  -1p = sequences not in fasta format: each line is a full sequence
  -nu = no changing to upper case before scan (non ATGC bases are discarded)
```
* The following parameters are required: `-q -o`
* `-q` is used to specify the input file containing query (foreground) sequences.  These are by default in fasta format but can be just sequences, one per line if `-1p` is also specified. **All sequences must be the same length.**
* `-b` similarly is used to specify the background sequences.
* `-o` the file to which results will be output. If it ends in `.gz` the output will be compressed (gzip).  The output file is tab delimited and has a header (example below).
* `-s` specifies the maximum length of k for k-mers to consider
* `-m` specifies the number of bases to smooth across when creating the positional base content used to generate the expected k-mer occurrence. This is unused of `-b` is provided.
* `-t` specifies how many threads should be used. This is only useful at present for scanning background sequences, and only two threads are helpful.
* `-dsk` means that reverse complements of k-mers are considered as well.  Note that this is NOT the same thing as including the reverse complements of input sequences as done with `-rc`.
* `-rc` means that the sequences that are being scanned will also be reverse complemented.  This un-sets the `-dsk` parameter.
* `-1p` use this if your input sequences are not in fasta format, but simply have one sequence per line.
* `-nu` use this if your input sequences contain lower case letters and you don't want them to be scanned (e.g. repeat masked)


#### Example output
```
k       kMer    kMerRC  NumFG   KS-D    KS-K    MaxDiffPos
6       GACTCA  TGAGTC  408.0   0.262594482658521       5.30414854879273        191
6       GAGTCA  TGACTC  437.0   0.239704534747105       5.01091422382467        103
5       GACTC   GAGTC   501.0   0.214229619207197       4.79510784055751        191
4       ATGA    TCAT    1507.0  0.113765841627258       4.4164011082312 204
5       ACTCA   TGAGT   620.0   0.173545101359368       4.32123817529382        195
5       AGTCA   TGACT   679.0   0.164113076450645       4.27639756855729        95
3       TCA     TGA     4444.0  0.0637461152708169      4.24952852569111        191
...
```
* `k` is the length of the k-mer
* `kMer` is the sequence of the k-mer
* `kMerRC` is the reverse complement of the k-mer
* `NumFG` indicates how many instances of the k-mer were observed in the data
* `KS-D` is the maximum distance (from 0-1) between the two CDFs
* `KS-K` is the KS statistic K describing the positional bias (greater numbers are more biased)
* `MaxDiffPos` is the position within the sequence of maximum difference
* `NumBG`  same for the background set of sequences (only present when a background is provided)

##### Tips:
Note well the differences between default, double stranded k-mer mode (`-dsk`), and reverse complement mode (`-rc`). By default, the strands are not reverse complemented and k-mers are considered strand specific. Double stranded k-mers (`-dsk`) mean that ATGC and its reverse complement GCAT are considered the same k-mer and counts for the two are combined. For instance, when scanning promoters, if you don't care what strand the motif is present on (e.g. what orientation the TF is binding in) but you DO care that it is binding most around base 20 of 300 in your set of sequences, use `-dsk`. If you are scanning ChIP peaks which inherently have no orientation, use `-rc` and both strands will be scanned.  In short: if your sequences have an orientation (e.g. genes, promoters, intron-exon junctions, cleavage/poly-A sites), use `-dsk` (or neither).  If your sequences have no orientation (e.g. ChIP peaks, open chromatin sites, enhancers) use `-rc`.  

For now, AMUSED-KS does not support gapped k-mers.  This will likely be implemented at some future date.

### alignKMers
This program aligns a set of k-mers into motif models, taking the best matches first and stoping when no good groupings remain that would not break some other equally good grouping. 
```
Usage: alignKMers -i <inFile> -o <outFile> -l <logFile> -v <minOverlap> -m <mismatchPenalty>
 -i <inFile> = unaligned k-mers
 -v <minOverlap> = minimum overlap to align [default=3]
 -r = also reverse complement
 -m <mismatchPenalty> = penalty for mismatches [default=1]
 -o <outFile> = where to output results [default=stdout]
 -l <logFile> = where to output results [default=stderr]
```
* `-i` is used to specify a file of k-mers to align. The format is one k-mer per line. Gaps can be represented as "-". 
* `-v` specifies the minimum number of bases overlapping with an existing model required to attempt to merge the two k-mer groups.
* `-r` also consider reverse complements of the sequences
* `-m` how much penalty to deal when aligning k-mers.  Matches to consensus are worth +1, mismatches as specified here.
* `-o` the output file to yeild the k-mer models.
* `-l` where to document messages

#####Tips:
Neither `AMUSED` nore `AMUSED-KS` presently support ambiguous bases (other than gaps). Consequently, `alignKMers` can be useful if `AMUSED` or `AMUSED-KS` provide many k-mers enriched or depleted and you want to group these into a smaller set of motifs. In this case, you should take your top N k-mers, and put those that are enriched in one file, while those that are depleted in another file.  You can then align each set individually to see if they are comprised of larger motifs.

## Citation
If you find AMUSED to be useful, pleace cite:
TBD

