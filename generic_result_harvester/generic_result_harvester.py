import os
import sys
import optparse

VERSION = '17.03.01'
NAME = 'generic_result_harvester'
AUTHOR = 'Jaakko Tyrmi'
DESCR = """
This program can be used to process any number of input files and collect
specific lines into a single output file.

If you use several conditions for finding line, all of the defined conditions
must be met for the line to be a match!

If there are several input files, the output is written into singe output
file by default. This behaviour can be setting the -m parameter "true".

If --search_subdirectories is set true, subdirectories AND the current
input directory are both searched.

Usage example:
python generic_result_harvester.py -i in/path/dir -x .txt -o out/path/file.out -n 155 -s
string1 -e string2 -c string3 -d 1 -m true
"""

print '\nRunning {0} v.{1} by {2}\n'.format(NAME, VERSION, AUTHOR)

parser = optparse.OptionParser(description=DESCR)
parser.add_option('-i', '--input_path', help='path to input file or dir')
parser.add_option('-u', '--search_subdirectories', help='if input path is '
                                                        'directory, include '
                                                        'subdirectories into '
                                                        'search (true/false) '
                                                        '(def false)')
parser.add_option('-x', '--input_extension', help='open files with this '
                                                  'extension')
parser.add_option('-o', '--output_path', help='output dir or file path')
parser.add_option('-n', '--line_number', help='number of line to extract',
                  type='int')
parser.add_option('-s', '--line_start', help='beginning of line to extract')
parser.add_option('-e', '--line_end', help='end of line to extract')
parser.add_option('-c', '--line_contains', help='line to extract contains')
parser.add_option('-d', '--delay', help='extract n:th line after matching '
                                        'line', type=int)
parser.add_option('-m', '--multiple_output_files', help='If input is a dir '
                                                      'and this arg is set '
                                                      '"true", output is '
                                                      'written into multiple '
                                                      'files instead of one.')
parser.add_option('-f', '--include_file_name', help='includes the input file '
                                                    'name to the output '
                                                    'strings (true/false)')


args = parser.parse_args()[0]




def generic_result_harvester(input_path, input_extension, output_path, line_number,
                             line_start, line_end, line_contains, delay,
                             multiple_output_files, include_file_name, search_subdirectories):

    # Set -d parameter value
    if delay is None:
        delay = 0

    # Check -u parameter value
    if search_subdirectories is not None:
        if search_subdirectories.lower() in ['t', 'true']:
            search_subdirectories = True
        elif search_subdirectories.lower() in ['f', 'false']:
            search_subdirectories = False
        else:
            print 'Error! Odd value for -u parameter! Should be "true" or "false"!'
            sys.exit(0)
    else:
        search_subdirectories = False
    if search_subdirectories and not os.path.isdir(input_path):
        print 'Error! -u parameter is available only when input path is a ' \
              'directory!'
        sys.exit(0)

    # Check -f parameter value
    if include_file_name is not None:
        if include_file_name.lower() == 'true':
            include_file_name = True
        elif include_file_name.lower() == 'false':
            include_file_name = False
        else:
            print 'Error! The "--include_file_name" should have a value ' \
                  'true/false!'
            sys.exit(0)
    else:
        include_file_name = False

    # Find & list input files
    if os.path.isfile(input_path):
        input_file_list = [input_path]
        if os.path.isdir(output_path):
            print 'Error! Input path is a file, but output path is a directory!'
            sys.exit(1)
    elif os.path.isdir(input_path):
        input_file_list = extract_files_in_dir(input_path, search_subdirectories)
    else:
        print 'Error! Input path is not a file or a directory:\n{0}'\
            .format(input_path)
        sys.exit(0)

    # Check -m parameter value
    if multiple_output_files is None or multiple_output_files in ['false', 'f']:
        multiple_output_files = False
        if os.path.exists(output_path):
            print 'Error! The output file path should be a non existing file, ' \
                  'when the multiple_output_files argument is set as "false"!'
            sys.exit(0)
    elif multiple_output_files.lower() in ['true', 't']:
        multiple_output_files = True
        if not os.path.exists(output_path) or os.path.isfile(output_path):
            print 'Error! The output file should be an existing directory ' \
                  'when multiple_output_files argument is set as "true"!'
            sys.exit(0)
    else:
        print 'Error! The multiple_output_files argument value contains an ' \
              'unknown value "{}"! Allowed values are "true" and "false". If ' \
              'the argument is not defined, "false" will be used as default.'
        sys.exit(0)

    if multiple_output_files and search_subdirectories:
        print 'Error! -u and -m parameters can not both be true!'
        sys.exit(0)


    out_lines = []
    i = 0
    for file_path in input_file_list:
        if input_extension:
            if not file_path.endswith(input_extension):
                continue

        file_handle = open(file_path)
        if include_file_name:
                file_name_string = file_path + '\t'
        else:
            file_name_string = ''
        matching_lines_found = 0
        j = 0
        line_numbers_to_extract = set([])
        for line in file_handle:
            j += 1
            line = line.strip()

            # Assume every line is the correct line
            matching_line = True

            # Check if any condition is NOT met:
            if line_number:
                if j != int(line_number):
                    matching_line = False

            if line_start:
                if not line.startswith(line_start):
                    matching_line = False

            if line_end:
                if not line.endswith(line_end):
                    matching_line = False

            if line_contains:
                if line_contains in line:
                    matching_line = False

            # Set capture the line
            if matching_line:
                line_numbers_to_extract.add(j + delay)

            if j in line_numbers_to_extract:
                out_lines.append(file_name_string + line)
                matching_lines_found += 1

        if matching_lines_found != 1:
            print 'WARNING! Weird number of matching lines ({0}) found from ' \
                  'file:\n{1}'.format(matching_lines_found, file_path)

        file_handle.close()
        i += 1

        # Write output file if multiple_output_files is true
        if multiple_output_files:
            current_out_path = os.path.join(output_path, file_path)
            write_output(current_out_path, out_lines, 1)
            out_lines = []

    # Write output file if multiple_output_files is false
    if not multiple_output_files:
        write_output(output_path, out_lines, i)


def write_output(output_path, out_lines, scanned_files):
    """Writes the output file

    :param
    output_path: Path to output file.
    out_lines: List of lines to be written,
    scanned_files: Number of input files scanned for this output file.
    """
    try:
        output_handle = open(output_path, 'w')
    except OSError:
        print 'Error! Unable to open output file:\n{0}'.format(output_path)
        sys.exit(1)
    output_handle.write('\n'.join(out_lines))
    output_handle.close()

    print '{0} matching lines found from {1} files!'.format(len(out_lines),
                                                            scanned_files)


def extract_files_in_dir(dir_path, search_subdirectories):
    """ Returns all file paths of a directory and subdirectories

    :param
    dir_path: path to a directory
    :return: list of file paths
    """
    file_paths = []
    for file_name in os.listdir(dir_path):
        absolute_path = os.path.join(dir_path, file_name)
        if os.path.isfile(absolute_path):
            file_paths.append(absolute_path)
        elif os.path.isdir(absolute_path) and search_subdirectories:
            file_paths += extract_files_in_dir(absolute_path, search_subdirectories)
    return file_paths


generic_result_harvester(args.input_path,
                         args.input_extension,
                         args.output_path,
                         args.line_number,
                         args.line_start,
                         args.line_end,
                         args.line_contains,
                         args.delay,
                         args.multiple_output_files,
                         args.include_file_name,
                         args.search_subdirectories)
