# het_hap_phaser

This is program is written for phasing haplotypes associated with a 
heterozygous variant where at least one parent is available per sample of 
interest. It takes joint-called VCF files as input and produces a tab-delimited 
table of variants, phased genotypes and proposed haplotype associated with the 
allele of interest.

This is currently in alpha.
                       
## REQUIREMENTS

This is written in python and requires python 3 to run. In addition, the 
following python modules are required:

* VASE:         https://github.com/gantzgraf/vase.git
    
* parse_vcf:    installable with pip and should be installed automatically when
                following the VASE installation instructions.

## AUTHOR

Written by David A. Parry at the University of Edinburgh. 

## COPYRIGHT AND LICENSE

MIT License

Copyright (c) 2018 David A. Parry

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.



