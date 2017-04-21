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


import operator
import optparse
import os
import sys

NAME = 'vcfcombineParallelizer'
VERSION = '17.04.11'
AUTHOR = 'Jaakko Tyrmi'
descr = """
Usage: (1) Run this script with appropriate parameters (2) edit and submit
SLURM_1.txt to SLURM (3) edit and submit SLURM_2.txt to SLURM (4) manually
concatenate and sort the output vcf files before downstream use.

VCF-files can be merged with vcftools vcf-merge utility or vcflib vcfcombine
tool. However vcflib vcfcombine consumes a large amount of RAM and is
computationally intensive making merging of large (or large number of) vcf
files impossible. vcf-merge does not consume much RAM but is even slower.

This program splits combine commands into n parts which are then run in
parallel. This is compatible with both vcftools vcftools and vcflib approach or
any other software with the following syntax:
path/to/executable input_1.vcf input_2.vcf ... > output.vcf
Path to appropriate executable is defined with -l parameter
(-l path/to/executable).

This program assumes that SLURM workload manager is available for use. The
program takes a directory containing multiple vcf files to be merged as an
input. User selects in how many parts the vcf files should be divided to.
Notice that reference genome scaffolds are not split, so each part consists
of one or more scaffolds. This means that when dealing with non-fragmented
"good quality" reference genomes only relatively low number of threads can be
generated. Notice that the input vcf files are expected to contain all
genotype calls (including monomorphic)!

Two SLURM array job files are generated: SLURM_1.txt and SLURM_2.txt. The
first one contains commands for splitting the input vcf files into n parts
(i.e. --process_count * [number of vcf files] processes will be created). The
second SLURM file then merges all samples for each genome part in parallel.
These files can be found in the output directory. The merged output files
should be manually concatenated and then sorted before downstream use!

When choosing appropriate process count (-p) notice that array jobs defined
in SLURM_1.txt creates process_count*number_of_samples parallel jobs.
SLURM_2.txt creates -p jobs. E.g. with 50 vcf files -p 5 would define 250
jobs in SLURM_1.txt and 5 jobs in SLURM_2.txt.

Also notice that some SLURM fields in SLURM_1.txt and SLURM_2.txt must be
manually edited to run the jobs (run time, memory usage etc.)
"""

print '\n\nRunning {0} v.{1}, by {2}'.format(NAME, VERSION, AUTHOR)

parser = optparse.OptionParser(description=descr)

parser.add_option('-j', '--job_name', help='SLURM job name')
parser.add_option('-v', '--input_vcf_file_dir', help='path to dir containing input vcf files')
parser.add_option('-r', '--reference_path', help='path to reference genome file (fasta)')
parser.add_option('-p', '--process_count', help='number of processes to use in parallelization', type=int)
parser.add_option('-l', '--vcf_merge_call', help='path to merge executable or a PATH variable (program parameters may be set also with triple underscore characters instead of single whitespaces e.g. vcf-merge___-s___-c___snps)')
parser.add_option('-t', '--temp_intermediate_dir', help='path to dir for intermediate vcf & bed files')
parser.add_option('-o', '--output_slurm_files_dir', help='path for slurm files to be created')
parser.add_option('-d', '--output_final_vcf_dir', help='output path to directory for merged vcf files')
parser.add_option('-m', '--load_environment_module', help='load an environment module (optional) e.g. -m module___load___mymodule (notice the use of triple underscore characters instead of single whitespace!)')

args = parser.parse_args()[0]


def vcfcombineParallelizer(output_final_vcf_dir, reference_path,
                           process_count, output_slurm_files_dir,
                           vcf_file_dir, temp_intermediate_dir, job_name,
                           vcf_merge_call, load_environment_module):

    if (output_final_vcf_dir is None
        or reference_path is None
        or process_count is None
        or output_slurm_files_dir is None
        or vcf_file_dir is None
        or temp_intermediate_dir is None
        or job_name is None
        or vcf_merge_call is None):
        print 'Error! At least one input parameter is missing!'
        sys.exit(0)

    if load_environment_module is None:
        load_environment_module = ''
    else:
        load_environment_module = load_environment_module.replace('___', ' ')

    vcf_merge_call = vcf_merge_call.replace('___', ' ')

    print '\nUsing parameters:'
    print '-r', reference_path
    print '-p', process_count
    print '-d', output_final_vcf_dir
    print '-o', output_slurm_files_dir
    print '-b', vcf_file_dir
    print '-t', temp_intermediate_dir
    print '-j', job_name
    print '-l', vcf_merge_call
    print '-m', load_environment_module

    # Open output file in good time to make sure it can be created
    try:
        out_handle_1 = open(os.path.join(output_slurm_files_dir, 'SLURM_1.txt'), 'w')
    except:
        print '\n\nERROR! Unable to open output file!\n'
        raise
    try:
        out_handle_2 = open(os.path.join(output_slurm_files_dir, 'SLURM_2.txt'), 'w')
    except:
        print '\n\nERROR! Unable to open output file!\n'
        raise

    # Identify vcf files
    vcf_file_names = []
    for file_name in os.listdir(vcf_file_dir):
        if not file_name.endswith('.vcf'): continue
        vcf_file_names.append(file_name)
    print '\n{0} vcf files found.'.format(len(vcf_file_names))

    # Read reference genome
    print '\nReading reference genome...'
    if not os.path.splitext(reference_path)[-1] in ['.fa', '.fas', '.fasta']:
        print 'Error! Reference file does not seem to be a fasta file!'
        sys.exit(0)

    try:
        in_handle = open(reference_path)
    except:
        print '\n\nError! Unable to open refence genome, reason\n\n'
        raise
    reference_sequences = {}
    reference_sequence_sizes = {}
    reference_sequence_order = {}
    first_line = True
    i = 0
    for line in in_handle:
        line = line.strip()
        if not line: continue
        if line.startswith('>'):
            if not first_line:
                i += 1
                seq = ''.join(seq)
                reference_sequences[current_scaffold] = seq
                reference_sequence_sizes[current_scaffold] = len(seq)
                reference_sequence_order[current_scaffold] = i
            current_scaffold = line.split(' ')[0]
            seq = []
            first_line = False
        else:
            try:
                seq.append(line)
            except NameError:
                print '\n\nError when reading reference file:\n{0}'.format(reference_path)
                print 'This does not seem like a fasta file, since first line does not contain > ' \
                      'character!'
                sys.exit(0)
    # Include the last scaffold
    i += 1
    seq = ''.join(seq)
    reference_sequences[current_scaffold] = seq
    reference_sequence_sizes[current_scaffold] = len(seq)
    reference_sequence_order[current_scaffold] = i
    in_handle.close()
    print 'Read {0} scaffolds'.format(len(reference_sequences))
    print 'Total sequence {0}b'.format(sum(reference_sequence_sizes.values()))
    print 'Largest scaffold contains {0}b'.format(max(reference_sequence_sizes.values()))

    # To balance load between processes (without splitting scaffolds)
    print '\nSplitting reference genome...'
    scaffolds_sorted_by_length = sorted(reference_sequence_sizes.items(),
                                        key=operator.itemgetter(1))
    scaffolds_sorted_by_length = list(reversed(scaffolds_sorted_by_length))
    output_scaffold_sets = []
    output_scaffold_set_sizes = []
    for scaffold in scaffolds_sorted_by_length:
        scaffold = scaffold[0]
        if len(output_scaffold_sets) < process_count:
            output_scaffold_sets.append([scaffold])
            output_scaffold_set_sizes.append(reference_sequence_sizes[scaffold])
        else:
            output_index = output_scaffold_set_sizes.index(min(output_scaffold_set_sizes))
            output_scaffold_sets[output_index].append(scaffold)
            output_scaffold_set_sizes[output_index] += reference_sequence_sizes[scaffold]
    assert len(output_scaffold_sets) == process_count

    # Generate bed file for each fragment
    print 'Generating bed files for reference genome parts...'
    bed_file_path_list = []
    i = 0
    for scaffold_set in output_scaffold_sets:
        i += 1
        bed_file_path = 'output_scaffold_set_{0}.bed'.format(i)
        bed_file_path = os.path.join(temp_intermediate_dir, bed_file_path)
        bed_file_path_list.append(bed_file_path)
        bed_out_handle = open(bed_file_path, 'w')
        bed_out_handle.write('#CHROM\tSTART\tEND\n')
        for scaffold in scaffold_set:
            bed_out_handle.write('\t'.join(map(str, [scaffold[1:], 0, reference_sequence_sizes[scaffold]])))
            bed_out_handle.write('\n')
        bed_out_handle.close()

    print '\nWriting shellscripts...'
    # Write array job file for splitting vcf files to parts based on bed files
    out_handle_1.write('#!/bin/sh\n')
    out_handle_1.write('#SBATCH --error={0}_%A_%a.err\n'.format(job_name))
    out_handle_1.write('#SBATCH --output={0}_%A_%a.out\n'.format(job_name))
    out_handle_1.write('#SBATCH --job-name={0}\n'.format(job_name))
    out_handle_1.write('#SBATCH --array={0}-{1}\n'.format(1, process_count*len(vcf_file_names)))
    out_handle_1.write('#SBATCH --time=\n')
    out_handle_1.write('#SBATCH --partition=\n')
    out_handle_1.write('#SBATCH --ntasks=\n')
    out_handle_1.write('#SBATCH --mem-per-cpu=\n')

    out_handle_1.write('\n\n')
    file_path = '{0}_split_subshell.sh'.format('"$SLURM_ARRAY_TASK_ID"')
    file_path = os.path.join(output_slurm_files_dir, file_path)
    out_handle_1.write('source {0}\n'.format(file_path))
    out_handle_1.close()

    # Write array job file for merging corresponding areas to single vcf file for all genome parts
    out_handle_2.write('#!/bin/sh\n')
    out_handle_2.write('#SBATCH --error={0}_%A_%a.err\n'.format(job_name))
    out_handle_2.write('#SBATCH --output={0}_%A_%a.out\n'.format(job_name))
    out_handle_2.write('#SBATCH --job-name={0}\n'.format(job_name))
    out_handle_2.write('#SBATCH --array={0}-{1}\n'.format(1, process_count))
    out_handle_2.write('#SBATCH --time=\n')
    out_handle_2.write('#SBATCH --partition=\n')
    out_handle_2.write('#SBATCH --ntasks=\n')
    out_handle_2.write('#SBATCH --mem-per-cpu=\n')

    out_handle_2.write('\n\n')
    file_path = '{0}_merge_subshell.sh'.format('"$SLURM_ARRAY_TASK_ID"')
    file_path = os.path.join(output_slurm_files_dir, file_path)
    out_handle_2.write('source {0}\n'.format(file_path))
    out_handle_2.close()

    current_split_file_number = 0
    current_merge_file_number = 0
    for bed_file_path in bed_file_path_list:
        current_merge_file_number += 1
        files_to_merge = []
        for vcf_file_name in vcf_file_names:
            vcf_file_path = os.path.join(vcf_file_dir, vcf_file_name)
            current_split_file_number += 1
            out_slurm_subscript_path = '{0}_split_subshell.sh'.format(current_split_file_number)
            out_slurm_subscript_path = os.path.join(output_slurm_files_dir, out_slurm_subscript_path)
            out_split_vcf_file_path = os.path.join(temp_intermediate_dir,
                                                   '{0}_{1}'.format(current_split_file_number, vcf_file_name))
            out_handle = open(out_slurm_subscript_path, 'w')
            out_handle.write('{0}\n'.format(load_environment_module))
            out_handle.write('vcftools --recode --vcf {0} --bed {1} --stdout > {2}\n'
                             .format(vcf_file_path,
                                     bed_file_path,
                                     out_split_vcf_file_path))
            out_handle.write('grep "#" {0} > {1}\n'.format(os.path.join(temp_intermediate_dir, out_split_vcf_file_path),
                                                         os.path.join(temp_intermediate_dir, out_split_vcf_file_path + '_sorted.vcf')))
            out_handle.write('grep -v "#" {0} | sort -k1,1V -k2,2n >> {1}\n'.format(os.path.join(temp_intermediate_dir, out_split_vcf_file_path),
                                                                      os.path.join(temp_intermediate_dir, out_split_vcf_file_path + '_sorted.vcf')))

            out_handle.write('bgzip -c {0} > {0}.gz\n'.format(os.path.join(temp_intermediate_dir, out_split_vcf_file_path + '_sorted.vcf')))
            out_handle.write('tabix {0}.gz\n'.format(os.path.join(temp_intermediate_dir, out_split_vcf_file_path + '_sorted.vcf')))

            out_handle.close()
            files_to_merge.append(os.path.join(temp_intermediate_dir, out_split_vcf_file_path + '_sorted.vcf.gz'))
        out_handle = open(os.path.join(output_slurm_files_dir,
                                       '{0}_merge_subshell.sh'.format(str(current_merge_file_number))),'w')
        out_handle.write('{0}\n'.format(load_environment_module))
        out_handle.write('{0} {1} > {2}\n'.format(vcf_merge_call,
                                                  ' '.join(files_to_merge),
                                                  os.path.join(output_final_vcf_dir, str(current_merge_file_number) + '.vcf')))
        out_handle.close()
    print '\nProgram run successful!'


vcfcombineParallelizer(args.output_final_vcf_dir,
                       args.reference_path,
                       args.process_count,
                       args.output_slurm_files_dir,
                       args.input_vcf_file_dir,
                       args.temp_intermediate_dir,
                       args.job_name,
                       args.vcf_merge_call,
                       args.load_environment_module)
