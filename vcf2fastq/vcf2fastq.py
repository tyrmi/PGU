import optparse
import sys

IUPAC_bases2dna = dict({"AA":"A","CC":"C","GG":"G","TT":"T",
               "AG":"R","CT":"Y","AC":"M","GT":"K",
               "AT":"W","CG":"S","CGT":"B","AGT":"D",
               "ACT":"H","ACG":"V","ACGT":"N"})
FASTQ_quality_string = '''!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'''

VERSION = '17.04.24'
NAME = 'vcf2fastq'
descr = """
This program converts single sample vcf-files into fastq files. Fastq quality
values are taken from the sample field's genotype quality. Input vcf files
should be sorted by position.

This program was created to produce suitable input files for psmc fq2psmcfa
script (Li & Durbin 2011 Inference of Human Population History From Whole
Genome Sequence of A Single Individual. Nature 475(7357):493-496)
"""

print '\nRunning {0} v.{1}'.format(NAME, VERSION)

parser = optparse.OptionParser(description=descr)
parser.add_option('-i', '--input', help='path to an input file (vcf)',
                  type='string')
parser.add_option('-o', '--output_file', help='output file path',
                  type='string')
args = parser.parse_args()[0]

def vcf2fastq(in_path, out_path):
    in_handle = open(in_path)
    out_handle = open(out_path, 'w')

    prev_scaffold = None
    lines_processed = 0
    for line in in_handle:
        lines_processed += 1
        if lines_processed % 1000000 == 0:
            print '{0} million vcf-file lines processed...'.format(lines_processed/1000000)
        # Skip empty and header lines
        line.strip()
        if line.startswith('##'): continue
        if not line: continue
        line = line.split('\t')
        if line[0].startswith('#'):
            if len(line) != 10:
                print 'Error! A single sample vcf file is expected as input ' \
                      'containing 10 columns, {0} columns found!'.format(len(line))
            continue
        #Parse and write output
        scaffold = line[0]
        position = int(line[1])
        if prev_scaffold is None: # i.e. the first row
            sequence = []
            sequence_quality = []
        elif scaffold != prev_scaffold: # Write output when scaffold changes
            out_handle.write('@{0}\n'.format(prev_scaffold))
            out_handle.write('{0}\n+\n'.format(''.join(sequence)))
            out_handle.write('{0}\n'.format(''.join(sequence_quality)))
            sequence = []
            sequence_quality = []
        else:
            if position != (prev_position + 1): # Write 'n' character for gaps in sequence
                sequence += ['n']*(position-prev_position-1)
                sequence_quality += ['!'*(position-prev_position-1)]
            sequence.append(parse_genotype(line))
            # If genotype is missing, use the lowest quality GQ value
            if sequence[-1] == 'n':
                sequence_quality.append('!')
            else:
                sequence_quality.append(parse_genotype_quality(line))
        prev_scaffold = scaffold
        prev_position = position
    # Write the last sequence
    out_handle.write('@{0}\n'.format(prev_scaffold))
    out_handle.write('{0}\n+\n'.format(''.join(sequence)))
    out_handle.write('{0}\n'.format(''.join(sequence_quality)))
    in_handle.close()
    out_handle.close()
    print 'File conversion successfull!'


def parse_genotype(line):
    """Extracts the genotype from vcf file line"""
    gt_field = line[9].split('\t')[0]
    gt_field = gt_field.split(':')[0]
    if gt_field in ('.', './.', '.|.'):
        return 'n'
    if '/' in gt_field:
        gt_field = gt_field.split('/')
    elif '|' in gt_field:
        gt_field.split('|')
    else:
        gt_field = [gt_field]

    genotypes = []
    for gt in gt_field:
        if gt == '0':
            genotypes.append(line[3])
        else:
            try:
                genotypes.append(line[4].split(',')[int(gt) - 1])
            except ValueError:
                print 'Error! Odd genotype value on line:'
                print '\t'.join(line)
                sys.exit(0)

    return IUPAC_bases2dna[''.join(sorted(genotypes))]


def parse_genotype_quality(line):
    """Extracts the genotype quality (GQ) from vcf file line"""
    sample_fields = line[9].split(':')
    try:
        GQ_index = line[8].split(':').index('GQ')
    except ValueError:
        print 'Error! GQ field was not found on line'
        sys.exit(0)
    try:
        try:
            # int(float()) to round down genotype quality
            return FASTQ_quality_string[int(float(sample_fields[GQ_index]))]
        except IndexError:
            return '~'
    except ValueError:
        # Genotype quality value is not integer on line
        print '\t'.join(line)
        sys.exit(0)


vcf2fastq(args.input, args.output_file)
