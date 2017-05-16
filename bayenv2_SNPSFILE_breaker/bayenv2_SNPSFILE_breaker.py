'''
Copyright © 2017 Jaakko Tyrmi. All Rights Reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

3. The name of the author may not be used to endorse or promote products 
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS OR 
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.
'''


import optparse
import os
import shutil
import sys

VERSION = '17.01.30'
NAME = 'bayenv2_SNPSFILE_breaker'
descr = """
This program breaks a bayenv2 SNPSFILE into separate SNPFILEs containing
single SNPs for parallelization purposes. This operation could be achieved
with unix split command but using this script (optionally) renames the output
SNPSFILEs based on vcf file from which the variants originated. Vcf to
SNPSFILE conversion can be done with PGDSpider. There is one caveat in using
the vcf input: No sites should be omitted in the vcf to SNPSFILE conversion
and PGD spider omits monomorphic sites. VCF file should therefore be filtered
before conversion so sites will be omitted in the conversion. Output files
will have file extension ".SNPFILE". Bayenv2 requires that all input files
should be in the same directory with bayenv2 executable, so bayenv2 executable,
ENVIRONFILE file and MATRIXFILE will be copied to each SNPFILE output
directory to allow parallelized run. The files will be renamed to
MATRIXFILE.txt and ENVIRONFILE.txt for easy parallelization with
PipelineMaster1000. Refer to bayenv2 manual on how to generate ENVIRONFILE and
MATRIXFILE.
"""

print '\nRunning {0} v.{1}'.format(NAME, VERSION)

parser = optparse.OptionParser(description=descr)
parser.add_option('-i', '--input_file', help='path to a SNPSFILE')
parser.add_option('-o', '--output_dir', help='output directory')
parser.add_option('-v', '--vcf_file', help='original vcf file of which '
                                           'SNPSFILE was made of (optional)')
parser.add_option('-b', '--bayenv2_executable', help='path to bayenv2 executable file')
parser.add_option('-e', '--environfile', help='path to environfile')
parser.add_option('-m', '--matrixfile', help='path to matrixfile')

args = parser.parse_args()[0]


def bayenv2_SNPSFILE_breaker(in_file_path, out_dir_path, vcf_file_path, bayenv2_executable, environfile, matrixfile):
    print 'Using parameters:'
    print 'SNPSFILE', in_file_path
    print 'Output directory', out_dir_path
    print 'Bayenv2 executable', bayenv2_executable
    print 'Vcf file', vcf_file_path
    print 'ENVIRONFILE', environfile
    print 'MATRIXFILE', matrixfile

    if vcf_file_path is not None:
        print 'Reading vcf file...'
        in_handle = open(vcf_file_path)
        i = 0
        vcf_positions = []
        for line in in_handle:
            i += 1
            line = line.strip()
            if not line: continue
            if line.startswith('#'): continue
            line = line.split('\t')
            if len(line) < 9:
                print 'ERROR! Vcf file format should contain at least 9 columns!'
                print 'Line {0} had {1} columns!'.format(i, len(line))
                sys.exit(0)
            vcf_positions.append('{0}_{1}'.format(line[0], line[1]))
        in_handle.close()
        print 'Done. Found {0} sites.'.format(len(vcf_positions))

    print '\nBreaking SNPSFILE into parts...'
    in_handle = open(in_file_path)
    i = 0
    variants = []
    for line in in_handle:
        if not line.strip(): continue
        i += 1
        variants.append(line)
        if i%2 == 0:
            if vcf_file_path is not None:
                current_out_dir_path = os.path.join(out_dir_path, 'bayenv2_input_{0}.SNPFILEDIR'.format(vcf_positions[(i/2)-1]))
                try:
                    os.mkdir(current_out_dir_path)
                except OSError: # Directory exists already
                    pass
                output_path = os.path.join(current_out_dir_path, 'bayenv2_input_{0}.SNPFILE'.format(vcf_positions[(i/2)-1]))
                try:
                    out_handle = open(output_path, 'w')
                except IndexError:
                    print 'ERROR! The vcf file has fewer variants ({0}) than the input SNPFILE ({1})!'.format(len(vcf_positions), i/2)
                    sys.exit(0)
            else:
                current_out_dir_path = os.path.join(out_dir_path, 'bayenv2_input_{0}.SNPFILEDIR'.format(vcf_positions[(i / 2) - 1]))
                try:
                    os.mkdir(current_out_dir_path)
                except OSError: # Directory exists already
                    pass
                output_path = os.path.join(current_out_dir_path, 'bayenv2_input_{0}.SNPSFILE'.format(str([i/2]).rjust(10, '0')))
                out_handle = open(output_path)
            out_handle.write(''.join(variants))
            out_handle.close()
            shutil.copy(bayenv2_executable, current_out_dir_path)
            shutil.copy(environfile, os.path.join(current_out_dir_path, 'ENVIRONFILE.txt'))
            shutil.copy(matrixfile, os.path.join(current_out_dir_path, 'MATRIXFILE.txt'))
            variants = []

    if vcf_file_path is not None:
        if i/2 > len(vcf_positions):
            print 'ERROR! The vcf file has more variants ({0}) than the input SNPSFILE ({1})!'.format(len(vcf_positions), i/2)
            print 'The names of output files are wrong!'
            sys.exit(0)
    print 'Done. Created {0} output SNPFILEs!'.format(i/2)

    print '\nProgram run successful!'

bayenv2_SNPSFILE_breaker(args.input_file, args.output_dir, args.vcf_file,
                         args.bayenv2_executable, args.environfile, args.matrixfile)
