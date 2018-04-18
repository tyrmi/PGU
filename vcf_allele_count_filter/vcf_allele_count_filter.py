import optparse
import sys

VERSION = '18.03.19'
NAME = 'vcf_allele_count_filter'
descr = """
This script can be used to filter vcf file sites by the number of detected
alleles. This is similar, but not identical to, vcftools --min-alleles,
--max-alleles, --mac and --max-mac options. The vcftools --min-alleles and
--max-alleles only look at the REF and ALT columns, but omit information from
sample columns. The --mac and --max-mac parameters only look at the numbers
of minor allele. It is therefore not possible with these parameters to filter
vcf file based on the number of alleles OBSERVED in the current site. This
script find the number of actual alleles on each site, and removes sites that
do not have proper number of alleles. This script works with data of any ploidy.
Notice that this script supports input either with -i parameter or STDIN and
output with either -o parameter of STDOUT.
"""

print '\nRunning {0} v.{1}'.format(NAME, VERSION)

parser = optparse.OptionParser(description=descr)
parser.add_option('-i', '--input_path', help='path to an input file (vcf)',
                  type='string')
parser.add_option('-o', '--output_path', help='output file path (vcf)',
                  type='string')
parser.add_option('-m', '--min_alleles', help='minimum number of different '
                                              'alleles observed on site',
                  type='int')
parser.add_option('-x', '--max_alleles', help='minimum number of different '
                                              'alleles observed on site',
                  type='int')
parser.add_option('-v', '--verbose', help='print progress during program run '
                                          '(true/false)')

args = parser.parse_args()[0]

def vcf_allele_count_filter(input_path, output_path, min_alleles,
                            max_alleles, verbose):
    # Parse verbose argument
    if verbose is None:
        verbose = False
    else:
        if verbose.lower() == 't' or verbose.lower == 'true':
            verbose = True
        else:
            verbose = False

    # Read input from stdin if no input file path is defined
    if input_path is None:
        in_handle = sys.stdin
    else:
        try:
            in_handle = open(input_path)
        except:
            sys.stderr.write('Unable to open the input file path. Reason:\n')
            raise

    # Write to stdout if no output file path is defined
    if output_path is None:
        out_handle = sys.stdout
    else:
        try:
            out_handle = open(output_path, 'w')
        except:
            sys.stderr.write('Unable to open the output file path. Reason:\n')
            raise

    sys.stderr.write('Starting to read vcf file...\n')
    i = 0
    sites_read = 0
    sites_kept = 0
    for line in in_handle:
        i += 1
        # Print progress if requested
        if verbose:
            if i % 1000000 == 0:
                sys.stderr.write('{0} M lines read...\n'.format(i / 1000000))

        # Parse vcf file
        if line.startswith('#'):
            out_handle.write(line)
            continue
        split_line = line.strip()
        if not split_line:
            out_handle.write(line)
            continue

        sites_read += 1

        split_line = split_line.split('\t')

        # If minimum allele number is 2, ALT field must contain a valid
        # genotype (i.e. not '.') for the line to pass filtering.
        if min_alleles > 1:
            if split_line[4] == '.': continue

        # Parse the FORMAT column to identify the index of genotype information
        # on sample columns
        GT_field_index = split_line[8].split(':').index('GT')

        # Identify the observed number of different alleles on current site
        allele_set = set([])
        for sample_field in split_line[9:]:
            try:
                genotypes = sample_field.split(':')[GT_field_index]
            except IndexError:
                # Sample field does not contain the GT field, assume as
                # missing data
                continue

            if len(genotypes) == 1:
                # Length of genotype field is 1, assume as haploid data
                genotypes = [genotypes]
            elif '/' in genotypes:
                # Unphased genotype field
                genotypes = genotypes.split('/')
            elif '|' in genotypes:
                # Phased genotype field
                genotypes = genotypes.split('|')

            for gt in genotypes:
                # Genotypes are always integers. If not int, assume
                # missing data.
                try:
                    gt = int(gt)
                except:
                    continue
                allele_set.add(gt)

        # Test if number of alleles falls below min_alleles limit
        if min_alleles is not None:
            if len(allele_set) < min_alleles:
                continue

        # Test if number of alleles exceeds max_alleles limit
        if max_alleles is not None:
            if len(allele_set) > max_alleles:
                continue

        out_handle.write(line)
        sites_kept += 1

    sys.stderr.write('Done. Kept {0} sites out of {1}.\n'.format(sites_kept, sites_read))


vcf_allele_count_filter(args.input_path, args.output_path, args.min_alleles,
                        args.max_alleles, args.verbose)
