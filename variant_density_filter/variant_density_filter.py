import optparse
import sys


from VcfWalker import VCFTable
VERSION = '15.08.21'
NAME = 'variant_density_filter'
descr = """This tool removes areas containing "too much" variation.

NOTICE! This program reads vcf file in a given window sizes (increments of 1)
and removes areas where more than -m variants are detected. Each line of vcf
file is assumed to be a variant (i.e. even when there is no variation).
Filtering of the vcf file should therefore be performed before using this tool!

The --remove_range argument can be used to make sure that paralogous area is
completely removed by removing adjacent areas of a paralogous area.
"""

parser = optparse.OptionParser(description=descr)
parser.add_option('-i', '--input', help='path to a .vcf file')
parser.add_option('-o', '--output', help='output file path')
parser.add_option('-s', '--sliding_window_size', help='sliding window size',
                  type=int)
parser.add_option('-m', '--max_variants', help='maximum number of variants in window',
                  type = int)
parser.add_option('-r', '--remove_range', help='size of adjacent area which '
                                               'will be removed', type=int)

args = parser.parse_args()[0]


def variant_density_filter(in_path, out_path, sliding_window_size, max_variants,
                  remove_range):
    if remove_range is None:
        remove_range = 0
    else:
        try:
            remove_range = int(remove_range)
        except ValueError:
            print 'The --remove_range option required an integer as a value!'
            sys.exit(1)

    vcf_handle = open(in_path)
    vcf_file = VCFTable(vcf_handle)
    #remove_areas has sequence ID (scaffold/contig) as a key. Areas to remove
    #  in this ID is a set of positions (int).
    remove_areas = {}

    #First iteration of file - detect paralogous areas
    print '\nSliding the window...'
    sliding_windows = []
    size_of_removed_area = 0
    current_chr = ''
    for line in vcf_handle:
        if vcf_file.isHeader(line):
            continue

        line = line.strip()
        line = line.split('\t')
        print 'CURRENT LINE', vcf_file.seq(line), vcf_file.pos(line)
        if vcf_file.seq(line) == current_chr:
            print 'current chromosome'
            for window in sliding_windows:
                window.add_variant(vcf_file.pos(line))
            sliding_windows.append(sliding_window(max_variants,
                                                  sliding_window_size,
                                                  current_chr,
                                                  vcf_file.pos(line),
                                                  remove_range))
        else:
            print 'Initializing new chromosome'
            remove_areas[current_chr] = set()
            for window in sliding_windows:
                remove_areas[current_chr] = remove_areas[current_chr].union(window.remove_area)
            try:
                size_of_removed_area += len(remove_areas[current_chr])
            except KeyError:
                pass
            sliding_windows = [sliding_window(max_variants,
                                              sliding_window_size,
                                              current_chr,
                                              vcf_file.pos(line),
                                              remove_range)]
            current_chr = vcf_file.seq(line)

        #raw_input('enter')
    #Handle the last chromosome
    remove_areas[current_chr] = set()
    for window in sliding_windows:
        remove_areas[current_chr].union(window.remove_area)
        try:
            size_of_removed_area += len(remove_areas[current_chr])
        except KeyError:
            pass

    vcf_handle.close()


    #Second iteration of file - copy good areas to an output file
    print '\nWriting output file...'
    total_sites = 0
    sites_filtered = 0
    try:
        output_handle = open(out_path, 'w')
    except IOError as ex:
        print 'Error! Unable to create output file, reason: {0}'.format(str(ex))
        sys.exit(1)
    vcf_handle = open(in_path)
    for line in vcf_handle:
        if vcf_file.isHeader(line):
            output_handle.write(line)
            continue

        total_sites += 1
        split_line = line.strip()
        split_line = split_line.split('\t')

        if vcf_file.seq(split_line) in remove_areas:
            if vcf_file.pos(split_line) in \
                remove_areas[vcf_file.seq(split_line)]:
                sites_filtered += 1
                continue
        output_handle.write(line)
    output_handle.close()

    log_handle = open(out_path+'.log', 'w')
    log_handle.write('{0} {1}\n\n'.format(NAME, VERSION))
    log_handle.write('Size of the removed area: {0}\n'.format('size_of_removed_area'))
    log_handle.write('\nOriginal number of variants\tVariants removed\n')
    log_handle.write('{0}\t{1}'.format(str(total_sites), str(sites_filtered)))
    log_handle.close()

    print 'Done.'
    print 'Size of the removed area: {0}'.format(size_of_removed_area)
    print '\nOriginal number of variants\tVariants removed'
    print '{0}\t{1}'.format(str(total_sites), str(sites_filtered))


class sliding_window:
    def __init__(self, max_variants, window_size, chrom, start_pos,
                 remove_range):
        self.max_variants = max_variants
        self.window_size = window_size
        self.chrom = chrom
        self.start_pos = start_pos
        self.end_position = start_pos + window_size
        self.remove_range = remove_range

        self.positions = [start_pos]
        self.remove_area = set()
        print 'Created window:'
        print 'max_variants', max_variants
        print 'window_size', window_size
        print 'chrom', chrom
        print 'start_pos', start_pos
        print 'end_pos', self.end_position
        print 'remove_range', remove_range

    def add_variant(self, pos):
        print '\nadding pos {0} to window {1}, {2}'.format(pos,
                                                         self.chrom,
                                                         self.start_pos)
        if self.end_position < pos:
            print 'Not added, position over end point'
            return
        if self.remove_area:
            print 'not added, this area has already been removed'
            return

        self.positions.append(pos)
        print len(self.positions), self.max_variants
        print len(self.positions) > self.max_variants
        if len(self.positions) > self.max_variants:
            self.remove_area = range(self.start_pos-self.remove_range,
                                     self.end_position+self.remove_range+1)
            self.remove_area = set(self.remove_area)
            print '-----!!!!!Marking this window to be removed!!!!!-----'
            print self.remove_area
        print '\n\n'
        return


variant_density_filter(args.input, args.output, args.sliding_window_size,
                       args.max_variants, args.remove_range)
