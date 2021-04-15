import os
import optparse
import subprocess
import sys

VERSION = '13.04.2018'
NAME = 'find_EIBs'
descr = """This program prints positions of intron-exon boundaries for
probe design reference genome.
"""
parser = optparse.OptionParser(description=descr)
parser.add_option('-p', '--probe_target_area', help='probe_target_area_file')
parser.add_option('-r', '--probe_reference_transcriptome',
                  help='probe_target_file')
parser.add_option('-w', '--whole_genome_reference_blast_db',
                  help='probe_target_file')
parser.add_option('-o', '--output', help='output file name (.bed)')

args = parser.parse_args()[0]

#FUNCTIONS

def find_EIBs(probe_target_area, probe_reference_transcriptome,
              whole_genome_reference_blast_db, output):

    # Get the target coordinates
    print 'Reading probe text...'
    in_handle = open(probe_target_area)
    for line in in_handle:
        split_line = line.strip()
        split_line = split_line.split('\t')
        # Increment bait areas by 1 as bait file is 0-based and blast file 1-based
        break # as there is only one line
    in_handle.close()
    reference_name = split_line[0]
    bait_area_start = int(split_line[1])
    bait_area_end = int(split_line[2])
    print 'probe line is:'
    print split_line
    print 'bait area start', bait_area_start, 'end', bait_area_end

    # Get the target sequence
    print '\nExtracting target sequence...'
    in_handle = open(probe_reference_transcriptome)
    correct_sequence = False
    for line in in_handle:
        line = line.strip()
        if line.startswith('>'):
            line = line[1:]
            line = line.split(' ')
            if line[0] == reference_name:
                correct_sequence = True
        elif correct_sequence:
            bait_sequence = line[bait_area_start:bait_area_end]
            break
    print bait_sequence

    print '\nBlasting...'
    tmp_out_fasta_path = probe_target_area + '.temp.fasta'
    tmp_out_handle = open(tmp_out_fasta_path, 'w')
    tmp_out_handle.write('>query\n{0}\n'.format(bait_sequence))
    tmp_out_handle.close()
    blast_result = subprocess.check_output('module load biokit && '
                                           'blastn -db {0} '
                                           '-num_threads 1 '
                                           '-outfmt 7 '
                                           '-query {1}'
                                           ''.format(whole_genome_reference_blast_db,
                                                     tmp_out_fasta_path),
                                           shell=True)
    os.remove(tmp_out_fasta_path)
    blast_result = blast_result.strip()
    if '0 hits found' in blast_result:
        blast_result = ''
    EIBs = []
    if not blast_result:
        print 'No blast result found'
    else:
        print 'blast result found:'
        best_hit = ''
        for line in blast_result.split('\n'):
            line = line.strip()
            if not line: continue
            if line.startswith('#'):
                continue
            best_hit = line
            break
        print best_hit
        best_hit = best_hit.split('\t')

        if len(best_hit) != 12:
            print 'ERROR! Odd blast result'
            print blast_result
            print sys.exit(1)
        q_start = int(best_hit[6])
        q_end = int(best_hit[7])
        print 'quert start', q_start, 'end', q_end
        print q_start, q_end
        if q_start > 25:
            print q_start, 'is larger than 25'
            EIBs.append(bait_area_start + q_start - 1)
        if q_end < bait_area_end - bait_area_start - 25:
            print q_end, 'is smaller than', bait_area_end - bait_area_start - 25
            EIBs.append(bait_area_start + q_end - 1)
        print EIBs

    print '\nWriting output file'
    out_handle = open(output, 'w')
    if EIBs:
        for EIB in EIBs:
            out_handle.write('{0}\t{1}\t{2}\n'.format(reference_name, EIB, EIB))
    else:
        out_handle.write('NA')
    out_handle.close()

    print 'Script run successful!'


#MAIN
find_EIBs(args.probe_target_area,
          args.probe_reference_transcriptome,
          args.whole_genome_reference_blast_db,
          args.output)

