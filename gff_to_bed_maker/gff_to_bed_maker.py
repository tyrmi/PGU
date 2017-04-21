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
from collections import defaultdict
import sys

NAME = 'gff_to_bed_maker'
VERSION = '17.04.20'
AUTHOR = 'Jaakko Tyrmi'
descr = """This program converts areas defined in an input gff3 files to bed
file format. The -f parameter can be used to limit the output features to a
certain type (e.g exon, CDS...). The -I parameter takes a path to a file
containing list of ID names to include in the output. The ID is parsed from
9. column of the gff file. NOTICE that the ID field is defined as the string
between '=' and ':' or ';' characters. This makes it possible to list IDs
without detailing possible exon numbers etc. that may be included after ':'
character.
"""

print '\n\nRunning {0} v.{1}, by {2}'.format(NAME, VERSION, AUTHOR)

parser = optparse.OptionParser(description=descr)

parser.add_option('-i', '--input_path', help='path to a gff file')
parser.add_option('-o', '--output_path', help='path to output bed file')
parser.add_option('-f', '--feature_types', help='feature to extract (default '
                                                'all), comma delimited list '
                                                'may be used to include '
                                                'several types.')
parser.add_option('-I', '--id_list_path', help='list of id:s to include in '
                                               'the output')
parser.add_option('-g', '--include_gene_name_in_bed', help='creates an addional '
                                                           'column to bed file '
                                                           '(making it non-standard) '
                                                           'containing the gene name'
                                                           '(compatible with pimatic.py '
                                                           'script)')

args = parser.parse_args()[0]

def gff_to_bed_maker(input_path, output_path, feature_types, id_list_path,
                     include_gene_name_in_bed):
    if feature_types is None:
        feature_types = []
    elif feature_types.lower == 'all':
        feature_types = []
    else:
        feature_types = feature_types.lower()
        feature_types = feature_types.split(',')
    feature_types = set(feature_types)

    if include_gene_name_in_bed.lower() in ('true', 't'):
        include_gene_name_in_bed = True
    elif include_gene_name_in_bed.lower() in ('false', 'f', None):
        include_gene_name_in_bed = False
    if include_gene_name_in_bed and id_list_path is None:
        print 'Error! -g parameter can be set true only when id list is available (-I parameter)!'
        sys.exit(0)

    allowed_ids = set([])
    if id_list_path is not None:
        print 'Reading ID list...'
        in_handle = open(id_list_path)
        for line in in_handle:
            line = line.strip()
            if not line: continue
            allowed_ids.add(line)
        print '{0} IDs found.\n'.format(len(allowed_ids))
        in_handle.close()

    in_handle = open(input_path)
    out_handle = open(output_path, 'w')
    if include_gene_name_in_bed:
        out_handle.write('#CHR\tSTART\tEND\tGENE_NAME\n')
    else:
        out_handle.write('#CHR\tSTART\tEND\n')
    features_found = 0
    features_written = 0
    feature_types_included = defaultdict(int)
    feature_length_included = defaultdict(int)
    print 'Processing gff file...'
    for line in in_handle:
        if line.startswith('#'): continue
        line = line.strip()
        if not line: continue
        line = line.split('\t')
        if len(line) != 9: continue
        features_found += 1
        #Skip if incorrect ID
        if id_list_path is not None:
            attributes = line[8]
            attributes = attributes.replace(':', ';')
            attributes = attributes.split(';')
            ID_found = False
            for a in attributes:
                if a.startswith('ID='):
                    if a[3:] in allowed_ids:
                        ID_found = True
                        ID_string = a[3:]
                        break
            if not ID_found: continue
        #Skip if incorrect feature
        line[2] = line[2].lower()
        if line[2] in feature_types or not feature_types:
            features_written += 1
            feature_types_included[line[2]] += 1
            feature_length_included[line[2]] += int(line[4])-int(line[3])
            corrected_start_pos = int(line[3])
            corrected_start_pos = (corrected_start_pos-1)
            if include_gene_name_in_bed:
                out_handle.write('{0}\t{1}\t{2}\t{3}\n'.format(line[0],
                                                            corrected_start_pos,
                                                            line[4],
                                                            ID_string))
            else:
                out_handle.write('{0}\t{1}\t{2}\n'.format(line[0],
                                                          corrected_start_pos,
                                                          line[4]))
    in_handle.close()
    out_handle.close()
    print '\nDone.'
    print '{0}/{1} features kept.'.format(features_written, features_found)
    print '\nfeature\tnumber_included\tlength_of_sequence_included'
    for feature in feature_types_included:
        print '{0}\t{1}\t{2}'.format(feature,
                                     feature_types_included[feature],
                                     feature_length_included[feature])

gff_to_bed_maker(args.input_path, args.output_path, args.feature_types,
                 args.id_list_path, args.include_gene_name_in_bed)
