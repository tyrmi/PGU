import optparse
import sys
from copy import deepcopy

VERSION = '17.03.24'
NAME = 'pm_estimatic'
descr = """
This program estimates the average probability of misorientation per site,
Pm. This can be used to correct the shape of unfolded allele frequency spectrum. 
For theoretical background see (equation 3, in particular)
Baudry & Depaulis 2003 Effect of Misoriented Sites on Neutrality Tests With
Outgroup. Genetics 165: 1619-1622.

This program takes a vcf file as an input. By default it is assumed that the
reference genotype in the vcf file contains ancestral genotypes. If this is
not the case the ancestral genotypes may also be input as a separate file
using the -a parameter. The ancestral genotype file should consist of three
columns: chromosome, site number and ancestral genotype. First line is
assumed to be a header line and is therefore omitted. Only the sites defined in
ancestral genotype file are included in the analysis. Any ploidy level is
allowed in the vcf file.
"""

print '\nRunning {0} v.{1}'.format(NAME, VERSION)

parser = optparse.OptionParser(description=descr)
parser.add_option('-i', '--input', help='path to an input file (vcf)',
                  type='string')
parser.add_option('-a', '--ancestral_genotype_table', help='path to file '
                                                           'containing '
                                                           'ancestral '
                                                           'genotypes',
                  type='string')
parser.add_option('-o', '--output_file', help='output file path',
                  type='string')
args = parser.parse_args()[0]

ALLOWED_GENOTYPES = set(['A', 'C', 'G', 'T'])

def pm_estimatic(in_path, ancestral_genotypes_path, out_path):
    print '\nUsing options:'
    print '-i', in_path
    print '-a', ancestral_genotypes_path
    print '-o', out_path

    try:
        out_handle = open(out_path, 'w')
    except:
        print '\n\nERROR! Unable to open output file:\n'
        raise

    # Read ancestral genotypes
    ancestral_genotypes = {}
    if ancestral_genotypes_path is not None:
        print '\nReading ancestral genotypes...'
        try:
            in_handle = open(ancestral_genotypes_path)
        except:
            print 'ERROR! Unable to open ancestral genotype file, reason:\n'
            raise
        i = 0
        for line in in_handle:
            i += 1
            if i == 1:
                continue
            line = line.strip()
            if not line: continue
            line = line.split('\t')
            if line[2].upper() not in ALLOWED_GENOTYPES:
                print 'Warning! Ancestral genotype at chr {0} pos {1} is {2}, ' \
                      'and is not included in allowed genotype list: {3}'.format(line[0], line[1], line[2], ALLOWED_GENOTYPES)
                print 'Omitting...'
                continue
            ancestral_genotypes['{0}_{1}'.format(line[0], line[1])] = line[2]
        if not ancestral_genotypes:
            print '\nERROR! No ancestral genotypes found in file\n{0}'.format(ancestral_genotypes_path)
            sys.exit(0)
        in_handle.close()
        print 'Read {0} ancestral genotypes.'.format(len(
            ancestral_genotypes))

    # Read vcf file to count the proportions of Pd, alpha and beta
    print '\nReading vcf file...'
    try:
        in_handle = open(in_path)
    except:
        print '\nERROR! Unable to open vcf file, reason:\n'
        raise
    Pd_count = 0
    tri_tetrallelic_sites_omitted = 0
    monomorphic_sites = 0
    polymorphic_sites = 0
    polymorphic_sites_and_fixed_diff_no_doublemutations = 0
    polymorphic_or_fixed_transition_count = 0
    polymorphic_or_fixed_transversion_count = 0
    polymorphic_doublemutation_sites = 0
    doublemutation_transition_count = 0.0
    doublemutation_transversion_count = 0.0
    i = 1
    for line in in_handle:
        if line.startswith('#'): continue
        line = line.strip()
        if not line: continue
        line = line.split('\t')
        if ancestral_genotypes:
            try:
                ancestral_genotype = ancestral_genotypes['{0}_{1}'.format(line[0], line[1])]
            except KeyError: # Current site not found in ancestral genotype file
                continue
        else:
            ancestral_genotype = line[3]

        # Make a list of available genotypes for this site, first gt being
        # the REF genotype and the rest ALT genotypes
        if line[4] == '.':
            currently_available_genotypes = [line[3]]
        else:
            currently_available_genotypes = [line[3]] + line[4].split(',')

        # Sample field for current site
        sample_fields = line[9:]

        # Parse genotype for all samples
        genotypes = set([])
        for sample in sample_fields:
            genotype = sample.split(':')[0]
            if '.' in genotype: continue
            genotype = genotype.replace('/', '|')
            genotype = genotype.split('|')
            genotype = genotype_int2genotype_char(currently_available_genotypes, genotype)
            genotypes = genotypes.union(set(genotype))
        if not genotypes:
            print '\nERROR! No data found in chr {0} pos {1}!'.format(line[0], line[1])
            print 'The input vcf file should be filtered for missing data!'
            sys.exit(0)

        # Print progress
        i += 1
        if i % 1000000 == 0:
            print '{0} Mbp processed...'.format(i / 1000000)

        # Include only mono and biallelic sites
        if len(genotypes) > 2:
            tri_tetrallelic_sites_omitted += 1
            continue

        current_and_ancestral_genotypes = deepcopy(genotypes)
        current_and_ancestral_genotypes.add(ancestral_genotype)
        current_and_ancestral_genotypes = list(current_and_ancestral_genotypes)

        if len(genotypes) == 2:
            polymorphic_sites += 1
        # Count transitions and transversions
        if len(genotypes) == 1 and ancestral_genotype in genotypes:
            # Monomorphic site
            monomorphic_sites += 1
        elif len(genotypes) == 1 and ancestral_genotype not in genotypes or \
            len(genotypes) == 2 and ancestral_genotype in genotypes:
            # Fixed difference from ancestral or polymorpic difference from ancestral
            polymorphic_sites_and_fixed_diff_no_doublemutations += 1
            j = 0
            for gt1 in current_and_ancestral_genotypes:
                j += 1
                for gt2 in current_and_ancestral_genotypes[j:]:
                    if gt1 != gt2:
                        if is_transition([gt1, gt2]):
                            polymorphic_or_fixed_transition_count += 1
                        elif is_transversion([gt1, gt2]):
                            polymorphic_or_fixed_transversion_count += 1
                        else:
                            print 'Sanity check FAIL! Polymorphism detected, but it is not a transition or a trasversion!'
                            print 'VCF line\n', line
                            print 'Genotypes', genotypes
                            print 'Ancestral genotype', ancestral_genotype
                            sys.exit(0)

        elif len(genotypes) == 2 and ancestral_genotype not in genotypes:
            Pd_count += 1
            # Sites with two mutations in current line
            polymorphic_doublemutation_sites += 1
            j = 0
            for gt in genotypes:
                if is_transition([ancestral_genotype, gt]):
                    doublemutation_transition_count += 0.5
                elif is_transversion([ancestral_genotype, gt]):
                    doublemutation_transversion_count += 0.5
                else:
                    print 'Sanity check FAIL! Polymorphism detected, but it is not a transition or a trasversion!'
                    print 'VCF line\n', line
                    print 'Genotypes', genotypes
                    print 'Ancestral genotype', ancestral_genotype
                    sys.exit(0)
        else:
            print 'Sanity check FAIL! This should not be happening!'
            print 'VCF line\n', line
            print 'Genotypes', genotypes
            print 'Ancestral genotype', ancestral_genotype
            sys.exit(0)

    in_handle.close()
    print 'Done.'

    Pd = float(Pd_count) / (polymorphic_sites)
    alpha = float(polymorphic_or_fixed_transition_count)/polymorphic_sites_and_fixed_diff_no_doublemutations
    beta = float(polymorphic_or_fixed_transversion_count)/polymorphic_sites_and_fixed_diff_no_doublemutations
    dm_alpha = float(polymorphic_or_fixed_transition_count+doublemutation_transition_count) / \
               (polymorphic_sites_and_fixed_diff_no_doublemutations + polymorphic_doublemutation_sites)
    dm_beta = float(polymorphic_or_fixed_transversion_count+doublemutation_transversion_count) / \
              (polymorphic_sites_and_fixed_diff_no_doublemutations + polymorphic_doublemutation_sites)
    kappa = (2 * alpha) / beta
    dm_kappa = (2 * dm_alpha) / dm_beta
    Pm = ((alpha ** 2 + 2 * beta ** 2) / (2 * beta * (2 * alpha + beta))) * Pd
    dm_Pm = ((dm_alpha ** 2 + 2 * dm_beta ** 2) / (2 * dm_beta * (2 * dm_alpha + dm_beta))) * Pd
    Pm_kappa = (kappa ** 2 + 2)/(4 * kappa + 2) * Pd
    dm_Pm_kappa = (dm_kappa ** 2 + 2)/(4 * dm_kappa + 2) * Pd
    out_lines = ['{0} {1}\n'.format(NAME, VERSION)]
    out_lines.append('\nUsing options:')
    out_lines.append('-i ' + in_path)
    out_lines.append('-a' + str(ancestral_genotypes_path))
    out_lines.append('-o' + out_path)
    out_lines.append('\nSites included in analysis:\t{0}'.format(monomorphic_sites +
                                                                 polymorphic_sites_and_fixed_diff_no_doublemutations +
                                                                 polymorphic_doublemutation_sites))
    out_lines.append('Sites included in analysis (excluding double mutations):\t{0}'.format(monomorphic_sites +
                                                                                            polymorphic_sites_and_fixed_diff_no_doublemutations))
    out_lines.append('Tri- and tetra-allelic sites omitted:\t{0}'.format(tri_tetrallelic_sites_omitted))
    out_lines.append('Monomorphic sites:\t{0}'.format(monomorphic_sites))
    out_lines.append('Polymorphic sites:\t{0}'.format(polymorphic_sites))
    out_lines.append('Polymorphic sites & fixed differences (omitting doublemutations):\t{0}'.format(polymorphic_sites_and_fixed_diff_no_doublemutations))
    out_lines.append('Polymorphic sites & fixed difference (omitting doublemutations), polymorphic and fixed transitions:\t{0}'.format(polymorphic_or_fixed_transition_count))
    out_lines.append('Polymorphic sites & fixed difference (omitting doublemutations), polymorphic and fixed transversions:\t{0}'.format(polymorphic_or_fixed_transversion_count))
    out_lines.append('Polymorphic doublemutation sites:\t{0}'.format(polymorphic_doublemutation_sites))
    out_lines.append('Polymorphic doublemutation sites, transitions (x0.5):\t{0}'.format(doublemutation_transition_count))
    out_lines.append('Polymorphic doublemutation sites, transversions (x0.5):\t{0}'.format(doublemutation_transversion_count))

    out_lines.append('\nPd equation:\tPd_count/polymorphic_sites')
    out_lines.append('Pd count:\t{0}'.format(Pd_count))
    out_lines.append('Pd:\t{0}'.format(Pd))

    out_lines.append('\nAlpha equation:\tpolymorphic_or_fixed_transition_count)/polymorphic_sites_and_fixed_diff_no_doublemutations')
    out_lines.append('Alpha:\t{0}'.format(alpha))

    out_lines.append('\nBeta equation:\t(polymorphic_or_fixed_transversion_count)/polymorphic_sites_and_fixed_diff_no_doublemutations')
    out_lines.append('Beta:\t{0}'.format(beta))

    out_lines.append('\nKappa equation:\t(2 * {0}) / {1}'.format('alpha', 'beta'))
    out_lines.append('Kappa equation:\t(2 * {0}) / {1}'.format(alpha, beta))
    out_lines.append('Kappa:\t{0}'.format(kappa))

    out_lines.append('\nPm equation:\t(({0} ** 2 + 2 * {1} ** 2) / (2 * {1} * (2 * {0} + {1}))) * {2}'.format('alpha', 'beta', 'Pd'))
    out_lines.append('Pm equation:\t(({0} ** 2 + 2 * {1} ** 2) / (2 * {1} * (2 * {0} + {1}))) * {2}'.format(alpha, beta, Pd))
    out_lines.append('Pm:\t{0}'.format(Pm))

    out_lines.append('\nPm equation using kappa:\t({0} ** 2 + 2)/(4 * {0} + 2) * {1}'.format('kappa', 'Pd'))
    out_lines.append('Pm equation using kappa:\t({0} ** 2 + 2)/(4 * {0} + 2) * {1}'.format(kappa, Pd))
    out_lines.append('Pm:\t{0}'.format(Pm_kappa))

    out_lines.append('\n')
    out_lines.append('#'*80)
    out_lines.append('\n')

    out_lines.append('Pm calculation with double mutations included in alpha & beta estimation:')

    out_lines.append('\nAlpha equation:\t(polymorphic_or_fixed_transition_count+doublemutation_transition_count) / '
                     '(polymorphic_sites_and_fixed_diff_no_doublemutations + polymorphic_doublemutation_sites)')
    out_lines.append('Alpha:\t{0}'.format(dm_alpha))

    out_lines.append('\nBeta equation:\t(polymorphic_or_fixed_transversion_count+doublemutation_transversion_count) / \
              (polymorphic_sites_and_fixed_diff_no_doublemutations + polymorphic_doublemutation_sites)')
    out_lines.append('Beta:\t{0}'.format(dm_beta))

    out_lines.append('\nKappa equation:\t(2 * {0}) / {1}'.format('alpha', 'beta'))
    out_lines.append('Kappa equation:\t(2 * {0}) / {1}'.format(dm_alpha, dm_beta))
    out_lines.append('Kappa:\t{0}'.format(dm_kappa))

    out_lines.append('\nPm equation:\t(({0} ** 2 + 2 * {1} ** 2) / (2 * {1} * (2 * {0} + {1}))) * {2}'.format('alpha', 'beta', 'Pd'))
    out_lines.append('Pm equation:\t(({0} ** 2 + 2 * {1} ** 2) / (2 * {1} * (2 * {0} + {1}))) * {2}'.format(dm_alpha, dm_beta, Pd))
    out_lines.append('Pm:\t{0}'.format(dm_Pm))
    out_lines.append('')

    out_lines.append('\nPm equation using kappa:\t({0} ** 2 + 2)/(4 * {0} + 2) * {1}'.format('kappa', 'Pd'))
    out_lines.append('Pm equation using kappa:\t({0} ** 2 + 2)/(4 * {0} + 2) * {1}'.format(dm_kappa, Pd))
    out_lines.append('Pm:\t{0}'.format(dm_Pm_kappa))

    print '\n\n'
    print '\n'.join(out_lines)

    out_handle.write('\n'.join(out_lines))
    out_handle.close()

    print '\nProgram run successful.'


def genotype_int2genotype_char(current_possible_genotypes, genotype_list):
    """Convert vcf file's samples' numerical genotypes to genotype characters

    current_possible_genotypes: List of possible genotypes with first value
    being the REF genotype and the rest ALT genotypes.
    genotype list: List of genotype chars from sample's GT field.
    """
    gt_list = []
    for gt in genotype_list:
        gt_list.append(current_possible_genotypes[int(gt)])
    return gt_list


def is_transition(genotypes):
    if 'A' in genotypes and 'G' in genotypes or \
        'C' in genotypes and 'T' in genotypes:
        return True
    else:
        return False


def is_transversion(genotypes):
    if 'A' in genotypes and 'C' in genotypes or \
        'A' in genotypes and 'T' in genotypes or \
        'G' in genotypes and 'C' in genotypes or \
        'G' in genotypes and 'T' in genotypes:
        return True
    else:
        return False

