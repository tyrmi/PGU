import optparse
import sys

import numpy
from statsmodels.robust.scale import mad as statsmodels_mad

VERSION = '16.03.11'
NAME = 'MAD_MAX'
descr = """
This program can be used in detecting areas of "too" high and low coverage
in BAM files. The parameter -m sets the number of Median Absolute Deviations
that the coverage is allowed to differ from median (default 3). The input file
is a coverage file created with "bedtools genomecov" (1-based site numbering
is expected). The input file defines coverage data for each site of a bam
file. NOTICE that parameter -d should be included in the bedtools genomecov
command line to report depth at each genome position. The correct input file
format contains three columns: scaffold name, position and depth. By default
MAD_MAX includes coverage values of 0, but they can be omitted by using
defining "-z false". Output is a bed file containing filtered areas. NOTICE
that as reference genomes may contain areas that should be omitted from
analysis (repetitive sequence etc.), the -b parameter may be used to define
areas to omit from MAD analysis.
"""

print '\nRunning {0} v.{1}'.format(NAME, VERSION)

parser = optparse.OptionParser(description=descr)
parser.add_option('-i', '--input', help='path to a genomecov file',
                  type='string')
parser.add_option('-o', '--output', help='output BED file path', type='string')

parser.add_option('-m', '--MAD', help='max coverage in MAD (def 3.0)',
                  type='float')
parser.add_option('-z', '--include_zeros', help='include zero coverage '
                                                'values')
parser.add_option('-b', '--omit_bed', help='area (bed) to omit from analysis')
args = parser.parse_args()[0]


def BED_coverage_filter(in_path, out_path, MAD_max, include_zeros, omit_bed):

    if include_zeros is None or include_zeros.lower() == 'true':
        include_zeros = True
    elif include_zeros.lower() == 'false':
        include_zeros = False
    else:
        print 'ERROR! -z parameter value should be either "true" or "false"! ' \
              'Current value {0}'.format(include_zeros)
        sys.exit(1)

    if MAD_max is None:
        MAD_max = 3.0

    print '\nUsed parameters:'
    print '-i {0}'.format(in_path)
    print '-o {0}'.format(out_path)
    print '-m {0}'.format(MAD_max)
    print '-z {0}'.format(include_zeros)
    print '-b {0}'.format(omit_bed)

    omit_sites = []
    if omit_bed is not None:
        print '\nReading omit bed file...'
        in_handle = open(omit_bed)
        for line in in_handle:
            if line.startswith('#'): continue
            line = line.split('\t')
            if len(line) < 3: continue
            for i in range(int(line[1]), int(line[2])):
                line[0] = line[0][1:]
                line[0] = line[0].split(' ')[0]
                omit_sites.append(line[0] + '_' + str(i))
        in_handle.close()
        print 'Area to omit from analysis: {0}'.format(len(omit_sites))
        omit_sites = set(omit_sites)

    in_handle = open(in_path)
    print '\nReading depth information...'
    depth_values = []
    i = 0
    for line in in_handle:
        i += 1
        line = line.strip()
        if not line: continue
        line = line.split('\t')
        if len(line) < 3: continue
        if len(line) > 3:
            print 'Error! Input coverage file should have 3 columns, found {0}' \
                  ' on line {1}'.format(len(line), i)
            sys.exit(0)
        if line[2] == '0' and not include_zeros:
            continue
        if line[0] + '_' + str(int(line[1])-1) in omit_sites:
            continue
        depth_values.append(int(line[2]))
    if len(depth_values) == 0:
        print 'Error! No depth values found! Check the input file format!'
        sys.exit(0)
    print '{0} depth values found.'.format(len(depth_values))

    print '\nCalculating minimum & maximum coverage...'
    depth_values_array = numpy.array(depth_values)
    median_coverage = int(round(numpy.median(depth_values_array)))
    MAD = int(round(statsmodels_mad(depth_values_array, c=1)))
    min_coverage = median_coverage - MAD * MAD_max
    max_coverage = median_coverage + MAD * MAD_max

    print '\nAnalyzing input file coverage...'
    in_handle.seek(0)
    bed_out_lines = []
    ongoing_BED_interval = None
    current_scaffold = ''
    prev_site = None
    excluded_sites = 0
    for line in in_handle:
        line = line.strip()
        if not line: continue
        line = line.split('\t')
        if len(line) != 3: continue
        if line[2] == '0' and not include_zeros: continue
        #Change the site information to 0-based integer:
        if line[1] == '0':
            print 'Error! Position information should starts from 1, not 0! ' \
                  'Check file:\n{0}'.format(in_path)
            sys.exit(0)
        line[1] = int(line[1])-1

        #Is line in a new scaffold, in omit_bed or does it contain gaps
        if line[0] + '_' + str(line[1]) in omit_sites:
            if ongoing_BED_interval is not None:
                bed_out_lines.append('\t'.join(map(str, ongoing_BED_interval)))
                ongoing_BED_interval = None
        elif ongoing_BED_interval is not None and current_scaffold != line[0]:
            bed_out_lines.append('\t'.join(map(str, ongoing_BED_interval)))
            ongoing_BED_interval = None
        elif ongoing_BED_interval is not None and (prev_site + 1) != line[1]:
                print 'Error! There is a gap in the input coverage data just ' \
                      'before the following line:\n{0}'.format(line)
                print 'Check the input file format! Obsolete bedtools ' \
                      'versions contain bug that may omit areas when using ' \
                      'genomecov!'
                sys.exit(0)
                #bed_out_lines.append('\t'.join(map(str, ongoing_BED_interval)))
                #ongoing_BED_interval = None
        #Is the coverage within allowed limits or not:
        if int(line[2]) > max_coverage or int(line[2]) < min_coverage:
            excluded_sites += 1
            if ongoing_BED_interval is None:
                ongoing_BED_interval = [line[0], line[1], line[1]+1]
            else:
                ongoing_BED_interval[2] = line[1]+1
        elif ongoing_BED_interval is not None:
            bed_out_lines.append('\t'.join(map(str, ongoing_BED_interval)))
            ongoing_BED_interval = None
        current_scaffold = line[0]
        prev_site = line[1]
    #Add a possible last interval spanning the end of genomecov file
    if ongoing_BED_interval is not None:
        bed_out_lines.append('\t'.join(map(str, ongoing_BED_interval)))
    in_handle.close()

    print '\nWriting output file...'
    out_handle = open(out_path, 'w')
    out_handle.write('\n'.join(bed_out_lines))
    out_handle.write('\n')
    out_handle.close()

    print '\nProgram run successful!'
    print '\nMedian coverage\t{0}'.format(median_coverage)
    print 'MAD\t{0}'.format(MAD)
    print 'Min coverage\t{0}'.format(min_coverage)
    print 'Max coverage\t{0}'.format(max_coverage)
    print 'Sites marked to be excluded\t{0}'.format(excluded_sites)
    print 'Areas marked to be excluded\t{0}'.format(len(bed_out_lines))


BED_coverage_filter(args.input, args.output, args.MAD, args.include_zeros,
                    args.omit_bed)
