# Modified by Jaakko Tyrmi in accordance with the below licence agreement of
# the original script

'''
Copyright (c) 2015 Lex Flagel and Amanda Kenney
Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO
THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
'''

'''
This script calculates genome-wide average pi (or d) for all pairwise
comparisons of haploid samples in a data set.
The data format is:
Positions are rows,
Samples are columns,
Genotypes are represented by one letter from ATGC. Missing data can be
any simple character that is not one of those letters.
This script does not handle diploid genotypes.
The code is written for a data set in which the first two columns are
chromosome and position.
'''
from collections import defaultdict
import optparse
import sys

VERSION = '17.04.20'
NAME = 'Pimatic'
descr = """
This script calculates Pi for multiple. This is based on a script published in
Garner et al. 2016 New Phytologist. Haploid samples are assumed.

NOTICE that if the bed file (-b parameter) includes a non-standard fourth
column (-g parameter) indicating the gene in question (created with gff_to_bed_maker.py) the
contents of sequences with identical gene identifiers are concatenated!
"""

print '\nRunning {0} v.{1}'.format(NAME, VERSION)

parser = optparse.OptionParser(description=descr)
parser.add_option('-i', '--input', help='path to an input file (vcf)',
                  type='string')
parser.add_option('-o', '--output_file', help='output file path',
                  type='string')
parser.add_option('-m', '--min_genotypes', help='minimum number of genotypes',
                  type='int')
parser.add_option('-b', '--bed', help='path to bed file containing '
                                      'areas to include in the '
                                      'analysis')
parser.add_option('-g', '--bed_contains_gene_ids', help='bed file contains a '
                                                        'fourth column describing'
                                                        'gene id (generated '
                                                        'with '
                                                        'gff_to_bed_maker.py)'
                                                        '(def False)')
parser.add_option('-s', '--min_number_of_available_sites_per_gene',
                  help='defines a minimum number of available sites required '
                       'for calculating pi',
                  type='int')
args = parser.parse_args()[0]

def pairwise(li):
    """a convienience function that produces all pairwise comparisons from a list"""
    for i in range(len(li)):
        j = i+1
        while j < len(li):
            yield (li[i], li[j])
            j += 1

def Pimatic(in_path, out_path, min_genotypes, bed_path,
            min_sites_per_gene, bed_contains_gene_ids):
    # find files matching the string in quotes; * is wildcard
    # using glob is helpful if your data is separated into multiple files, e.g.,
    #  by chromosome
    # have all the files match this structure, with the * being their unique
    # identifier
    if bed_path is None:
        print 'Error! bed file path is required! Try the original pimatic if you wish to not use bed file.'
        sys.exit(0)
    if bed_contains_gene_ids.lower() in ('t', 'true'):
        bed_contains_gene_ids = True
    elif bed_contains_gene_ids.lower() in ('f', 'false', None):
        bed_contains_gene_ids = False
    gene_contents = split_vcf_to_genes(in_path, bed_path, bed_contains_gene_ids)
    number_of_gene_regions = len(gene_contents)
    print 'Read {0} gene regions.'.format(number_of_gene_regions)

    print '\nCalculating Pi...'
    allowed_genotypes = '01'

    gene_region_number = 0
    results = {}
    out_handle_header = []
    for gene_definition in sorted(gene_contents.keys()):
        num = defaultdict(float)
        den = defaultdict(float)
        gene_region_number += 1
        print 'Calculating Pi for gene region number {0}/{1}'.format(gene_region_number,
                                                                     number_of_gene_regions)
        i = 0
        genotype_list = gene_contents[gene_definition]
        genotypes_examined = 0
        too_few_genotypes = False
        for line in genotype_list:
            i += 1
            if not line.strip(): continue
            if line.startswith('##'): continue
            if line.startswith('#'):
                line = line.split('\t')
                sample_names = line[9:]
                continue
            line = line.split('\t')
            genotypes = dict(zip(sample_names, parse_genotypes(line[9:]))) #
            if genotypes.values().count('0') + genotypes.values().count('1') < min_genotypes: continue
            for sample_1, sample_2 in pairwise(sorted(genotypes.keys())):
                if sample_1 == sample_2: continue
                comparison_name = sample_1 + '/' + sample_2
                comparison_name = comparison_name.replace('\n', '')
                comparison_name = comparison_name.replace('\r', '')
                sample_1_genotype, sample_2_genotype = genotypes[sample_1], genotypes[sample_2]
                if sample_1_genotype in allowed_genotypes and sample_2_genotype in allowed_genotypes:
                    if sample_1_genotype == sample_2_genotype:
                        den[comparison_name]+=1 # if sequence is identical, only numerator is incremented
                        if comparison_name not in num:
                            num[comparison_name] = 0
                    else: # if there is a SNP, den and num are both incremented
                        den[comparison_name]+=1
                        num[comparison_name]+=1
                else: # make sure all comparisons are included in the output
                    if comparison_name not in den:
                        den[comparison_name] = 0
                    if comparison_name not in num:
                        num[comparison_name] = 0
            genotypes_examined += 1
        print 'analyzing gene', gene_definition

        assert len(num) == len(den) # for eachfile/chromosome, print file name
        # and number of pairwise samplecomparisons, should be the same for numerator and denominator
        print 'number of comparisons done', len(num)
        if len(num) != 0:
            out_handle_header.append(gene_definition)
        for comparison_name in num:
            #Store Pi into results
            if den[comparison_name] < min_sites_per_gene:
                try:
                    results[comparison_name].append('NA')
                except KeyError:
                    results[comparison_name] = ['NA']
            else:
                try:
                    results[comparison_name].append(num[comparison_name] / den[comparison_name])
                except KeyError:
                    results[comparison_name] = [(num[comparison_name] / den[comparison_name])]
            #print comparison_name, len(results[comparison_name])
            #raw_input('enter')
    out_handle = open(out_path, 'w')
    #Output file header line contains gene positional info:
    #out_handle.write('\t{0}\n'.format('\t'.join(gene_contents.keys())))
    out_handle.write('\t{0}\n'.format('\t'.join(out_handle_header)))
    for comparison_name in sorted(results.keys()):
        print comparison_name
        out_handle.write('{0}\t{1}\n'.format(comparison_name,
                                             '\t'.join(map(str, results[comparison_name]))))
    out_handle.close()

def parse_genotypes(gt_fields):
    return [gt.split(':')[0] for gt in gt_fields]

def split_vcf_to_genes(vcf_path, bed_path, bed_contains_gene_ids):
    #Read in the vcf file
    print '\nReading vcf file...'
    in_handle = open(vcf_path)
    header_lines = []
    genotype_lines = {}
    for line in in_handle:
        if line.startswith('#'):
            header_lines.append(line)
            continue
        line = line.strip()
        if not line: continue
        split_line = line.split()
        genotype_lines[(split_line[0], split_line[1])] = line
    in_handle.close()


    #Read bed file and select gene areas from vcf file
    print '\nReading bed file...'
    gene_contents = {}
    in_handle = open(bed_path)
    total_bed_area_size = 0
    sites_not_found = 0
    empty_gene_areas = 0
    first_line = True
    i = 0
    for line in in_handle:
        i += 1
        # Skip the possible header line
        if first_line:
            if line.startswith('#') or \
                    line.lower().startswith('browser') or \
                    line.lower().startswith('track'):
                continue
            first_line = False
            
        line = line.strip()
        if not line: continue
        line = line.split('\t')
        #If bed file line has four columns the last column is expected to contain gene name added by the gff_to_bed_maker.py
        if len(line) < 4 and bed_contains_gene_ids:
            print 'Error! Less than four columns found on line {i}, but four ' \
                  'expected as -g parameter was set True!'.format(i)
            sys.exit(0)
        if bed_contains_gene_ids:
            gene_definition = line[3]
        else:
            gene_definition = '_'.join(line) #ie. scaffold_startpos_endpos(_genename)
        gene_contents[gene_definition] = []
        total_bed_area_size += int(line[2]) - int(line[1])
        for position in xrange(int(line[1]), int(line[2])):
            try:
                gene_contents[gene_definition].append(genotype_lines[(line[0], str(position+1))])
            except KeyError:
                sites_not_found += 1
        if len(gene_contents[gene_definition]) == 0:
            empty_gene_areas += 1
            gene_contents.pop(gene_definition)
        else:
            gene_contents[gene_definition] = header_lines + gene_contents[gene_definition]

    in_handle.close()
    if sites_not_found == total_bed_area_size:
        print 'Error! No sites covered by bed file found in vcf-file!'
        sys.exit(0)
    print '{0} bed entries contain no sequence in vcf file.'.format(empty_gene_areas)
    if sites_not_found:
        print 'Warning! {0}/{1} sites ({2}%) covered in the bed file are ' \
              'not found in the vcf file!'.format(sites_not_found,
                                                  total_bed_area_size,
                                                  float(sites_not_found)/total_bed_area_size*100)
    return gene_contents

Pimatic(args.input, args.output_file, args.min_genotypes, args.bed,
        args.min_number_of_available_sites_per_gene, args.bed_contains_gene_ids)