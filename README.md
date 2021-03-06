# VCFIXER
New updated version of dataset_fixer (discontinued), now it does not rely in any other program or library.
Checks vcf files and removes loci and samples with too many missing genotypes, it can also input genotypes in missing - various options for this.

It works smoothly with VCF files generated by `populations` [Stacks](http://catchenlab.life.illinois.edu/stacks/): populations.snps.vcf
 - I do not own or am related in any way to Stacks people. Is a good useful pipeline for SNP calling
 - When I say smoothly I mean in our lab server.
 - Included some options in order to adjust it to other VCF file formats, so it should work fine with most of them.

May be smoothly is an overstatement, **IT WILL TAKE LONG TIME** to run, specially in order to analyse huge files, run it on a server, still pending to implement multi-thread.
This was done long time ago, eventually this will be moved to python and the code should be transformed in a few functions that can be called instead of repeating so many lines.

Some population based analysis could fail if the sample names in the vcf file do not include the population names in their codename:
POP1_001; or popA002; or CoolPlace042; etc. 
You may solve that by using a popmap (file with one sample per row, first column with sample name, second column with population name). There is a `popmap_maker` [in other of my repositories](https://github.com/M-Osky/handle_edit_files).

`vcfixer` has many option and flags that can be used to change its behaviour or to adapt it to your file format.

### To see the full information, usage, arguments, and settings just use any of the usual help flags (-h --h help -help --help).

    vcfixer.pl --help

           vcfixer_1.5.pl   Help Information 30.VII.2020
        --------------------------------------------------

        This program will read a VCF file, first will delete empty or monomorphic SNP, then empty samples,
        then will delete loci with high missing rate, then samples with high missing rate,
        and finally  will input genotypes in the missing values left "./.".
        It is designed to work with the VCF file generated by the program "populations" (Stacks).
        Populations can be read either from sample names - saving the first x characters from each individual ID
        or from a popmap file - Stacks format: tab separated, one sample and one population per row.
        A log file, a list with deleted and kept SNPs, deleted and kept samples, and a new popmap will be saved in the output directory


        Command line arguments and defaults:

        --input / --vcf           name (or path) of the VCF file. Default: populations.snps.vcf

        --infocols                number of locus information columns before the first sample column. Default: 9

        --poplength               [int] how many characters of the sample name belong to the population code? Default: 2
        --popmap                  Alternatively provide a popmap to read which sample belongs to which population. No default

        --empty                   [float] missing rate from which a sample will be considered "empty" and deleted. Default: 0.8

        --miss_loci               [float] missing rate from which a loci should be deleted. Default: 0.3

        --miss_samples            [float] missing rate from which a sample should be deleted. Default: 0.3

        --minpop                  [int] minimum number of samples a population must have in order to keep it. Default: 8

        --gral_miss               How to replace the regular missing values? There are four options. Default: pop
                                  "pop" to replace it with the population mode* (most frequent genotype)
                                  "global" to replace the missing with the whole dataset mode*.
                                  "miss" to leave it as missing. "2/2" or any other value to input that value.

        --pop_miss                What to input if a SNP is missing in an entire population? Three options available. Default: global
                                  "global" to input the global mode*, "miss" to keep them as missing,
                                  "2/2", or "5/5", or any other value: to input a new genotype and remark its difference from the rest.

        --noquality               [flag] add this if you do not want vcfixer_1.5.pl to input quality/probability information in missing
                                  Some software cannot process genotypes without quality metadata, by default will be also inputed,
                                  quality information is not editable from command line and should be handled carefully, values are:
                                        0/0:1:1,0:12:-0.05,-0.82,-1.42
                                        0/1:1:1,1:12:-1.01,-0.05,-1.01
                                        1/1:1:0,1:12:-1.42,-0.82,-0.05

        --summary                 path/name for a summary table that will gather some details of populations and vcf_fixer outputs.
                                  if '--summary no' the file will not be created. Default: summary_table_vcf_fixer.txt

        --poplog                  [optional] path or name of the populations (Stacks) log file. Default: ref_map.log
                                  if '--poplog no' it will not look for a logfile. Only used for summary table.
                                  If no path (only name) provided will look for file in vcf file location.

        --out                     path or new name for the output file. By default will be "input name" + "tail" + ".vcf"

        --tail                    Something to the end of the input file name to generate the output name.
                                  If no tail provided, will add to the file name the sample number, SNP number, and missing handling:
                                   1p: regular missing replaced with population mode
                                   1g: regular missing replaced with global mode
                                   1m: regular missing not handled - left as missing.
                                   1x: regular missing replaced with custome input
                                   0m: if not regular missing found in the file
                                   2g: when a SNP is missing in the whole population, the global mode is input
                                   2m: when missing in the whole population is left as missing
                                   2x: custom input
                                   0p: there are not SNPs entirely missing in any population


        Command line call example:
        vcfixer.pl --input /usr/home/refmap/populations.snps.vcf --poplength 3 --gral_miss global --minpop 1 --summary no


                * When two or more alleles have the highest frequency, one of them will be picked randomly for each SNP








