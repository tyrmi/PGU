import optparse
import sys


VERSION = '17.04.27'
NAME = 'ParalogAreaBEDmatic'
AUTHOR = 'Jaakko Tyrmi'
DESCRIPTION = """
This program scans a vcf file for paralogous calls (heterozygous SNPs in
haploids) and creates an output BED file containing paralogous areas. By
default the site is considered as paralogous even if a single heterozygote call
is found. This behaviour can be modified with the -M parameter.

NOTICE! There may be some inaccuracy in the output BED file, about +-1 base!
Revise the code if this is important!
"""

parser = optparse.OptionParser(description=DESCRIPTION)
parser.add_option('-i', '--input', help='path to a .vcf file')
parser.add_option('-o', '--out_path', help='output file path')
parser.add_option('-M', '--max_number_of_heterozygotes', help='maximum number '
                                                              'of heterozygote '
                                                              'calls (def 0)',
                  type=int)
parser.add_option('-r', '--remove_range', help='radius of area where nearby '
                                               'snps are removed', type='int')
parser.add_option('-g', '--genome_size_file', help='genome size file')

args = parser.parse_args()[0]


def ParalogAreaBEDmatic(in_path, out_path, max_number_of_heterozygotes,
                        remove_range, genome_size_file):

    if max_number_of_heterozygotes is None:
        max_number_of_heterozygotes = 0

    print 'Running {0} v.{1}\n'.format(NAME, VERSION)
    print DESCRIPTION

    print '\nParameters:'
    print '--input', in_path
    print '--out_path', out_path
    print '--remove_range', remove_range
    print '--genome_size_file', genome_size_file

    print '\nReading genome size file...\n'
    chr_sizes = read_genome_size_to_dict(genome_size_file)

    print 'Identifying paralog areas...'
    try:
        in_handle = open(in_path)
    except IOError as err:
        print 'Unable to open input file, reason:\n'.format(err)
        sys.exit(1)

    BED_lines = []
    number_filtered = 0
    lines_read = 0
    for line in in_handle:
        lines_read += 1
        if lines_read % 1000000 == 0:
            print '{0} million lines read...'.format(lines_read/1000000)
        if line.startswith('#'): continue
        line = line.strip()
        if not line: continue

        line = line.split('\t')
        sample_fields = line[9:]
        heterozygous_genotypes = 0
        filter_this_site = False
        for field in sample_fields:
            genotype = field.split(':')[0]
            if '/' in genotype:
                genotype = genotype.split('/')
                if genotype[0] != genotype[1]:
                    heterozygous_genotypes += 1
                    if heterozygous_genotypes > max_number_of_heterozygotes:
                        filter_this_site = True
                        break

        if filter_this_site:
            number_filtered += 1
            scaffold = line[0]
            position = int(line[1])
            start = position - remove_range
            end = position + remove_range + 1
            if start < 0:
                start = 0
            if end > chr_sizes[scaffold]:
                end = chr_sizes[scaffold]
            BED_lines.append('{0}\t{1}\t{2}'.format(scaffold, start, end))

    in_handle.close()

    print '\nWriting output file...'
    try:
        out_handle = open(out_path, 'w')
    except IOError as err:
        print 'Unable to open output file, reason:\n'.format(err)
        sys.exit(1)
    out_handle.write('#CHR\tSTART\tEND\n')
    out_handle.write('\n'.join(BED_lines))

    print '\nProgram run successful!'
    print '\n{0} paralogous sites found.'.format(number_filtered)


def read_genome_size_to_dict(genome_path):
    try:
        handle = open(genome_path)
    except IOError as ex:
        print "Error! Unable to open genome file!\n{0}".format(ex)
        sys.exit(1)

    dictionary = {}
    for line in handle:
        line = line.strip()
        if not line: continue

        line = line.split('\t')
        if len(line) != 2: continue
        dictionary[line[0].strip('>')] = int(line[1].strip('>'))

    handle.close()
    return dictionary


ParalogAreaBEDmatic(args.input, args.out_path,
                    args.max_number_of_heterozygotes, args.remove_range,
                    args.genome_size_file)
