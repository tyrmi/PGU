'''
Copyright 2017 Jaakko Tyrmi. All Rights Reserved.

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


import os
import sys
import optparse

VERSION = '18.01.03'
NAME = 'generic_result_harvester'
AUTHOR = 'Jaakko Tyrmi'
DESCR = """
This program can be used to process any number of input files and collect
specific lines into a single output file. Output is written to stdout if -o
parameter is not defined.

If you use several conditions for finding line, all of the defined conditions
must be met for the line to be a match!

By default each condition can take only single value. Multiple alternative 
values can be defined (with -n, -s, -e and -c options) by giving a path to
a text file containing an alternative value on each row. For example
-c /path/to/file.txt
where file.txt would contain rows:
im_looking_for_this_string
or_this
or_maybe_this

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
parser.add_option('-n', '--line_number', help='number of line to extract')
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
parser.add_option('-D', '--debug', help='print all sorts of info during runtime (true/false)')

args = parser.parse_args()[0]


def generic_result_harvester(input_path, input_extension, output_path, line_number,
                             line_start, line_end, line_contains, delay,
                             multiple_output_files, include_file_name, search_subdirectories,
                             debug):
    if debug is None:
        debug = False
    elif debug.lower() == 'false' or debug.lower() == 'f':
        debug = False
    elif debug.lower() == 'true' or debug.lower() == 't':
        debug == True
    else:
        sys.stderr.write('Error! Odd value for -D option. Allowed values are "true" and "false"!')
        sys.exit(0)

    if input_path is None:
        sys.stderr.write('Error! -i parameter is required!\n')
        sys.exit(0)

    # Make condition values as lists of strings (and split them if -M parameter has been defined)

    if line_number is not None:
        if os.path.isfile(line_number):
            line_number = set(map(int, read_multiple_input_parameters(line_number)))
        else:
            line_number = set([int(line_number)])
    if line_start is not None:
        if os.path.isfile(line_start):
            line_start = read_multiple_input_parameters(line_start)
        else:
            line_start = [line_start]
    if line_end is not None:
        if os.path.isfile(line_end):
            line_end = read_multiple_input_parameters(line_end)
        else:
            line_end = [line_end]
    if line_contains is not None:
        if os.path.isfile(line_contains):
            line_contains = read_multiple_input_parameters(line_contains)
        else:
            line_contains = [line_contains]

    if debug:
        sys.stderr.write('Using search parameters:')
        sys.stderr.write('\n')
        sys.stderr.write('-n')
        sys.stderr.write('\n')
        sys.stderr.write(str(line_number))
        sys.stderr.write('\n')
        sys.stderr.write('-s')
        sys.stderr.write('\n')
        sys.stderr.write(str(line_start))
        sys.stderr.write('\n')
        sys.stderr.write('-e')
        sys.stderr.write('\n')
        sys.stderr.write(str(line_end))
        sys.stderr.write('\n')
        sys.stderr.write('-c')
        sys.stderr.write('\n')
        sys.stderr.write(str(line_contains))
        sys.stderr.write('\n')

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
            sys.stderr.write('Error! Odd value for -u parameter! Should be "true" or "false"!')
            sys.exit(0)
    else:
        search_subdirectories = False
    if search_subdirectories and not os.path.isdir(input_path):
        sys.stderr.write('Error! -u parameter is available only when input path is a directory!')
        sys.exit(0)

    # Check -f parameter value
    if include_file_name is not None:
        if include_file_name.lower() == 'true':
            include_file_name = True
        elif include_file_name.lower() == 'false':
            include_file_name = False
        else:
            sys.stderr.write('Error! The "--include_file_name" should have a value true/false!')
            sys.exit(0)
    else:
        include_file_name = False

    # Find & list input files
    if os.path.isfile(input_path):
        input_file_list = [input_path]
        if output_path is not None and os.path.isdir(output_path):
            sys.stderr.write('Error! Input path is a file, but output path is a directory!')
            sys.exit(1)
    elif os.path.isdir(input_path):
        input_file_list = extract_files_in_dir(input_path, search_subdirectories)
    else:
        sys.stderr.write('Error! Input path is not a file or a directory:\n')
        sys.stderr.write(str(input_path))
        sys.exit(0)

    # Check -m parameter value
    if multiple_output_files is None or multiple_output_files in ['false', 'f']:
        multiple_output_files = False
        if output_path is not None and os.path.exists(output_path):
            sys.stderr.write('Error! The output file path should be a non existing file, '
                  'when the multiple_output_files argument is set as "false"!')
            sys.exit(0)
    elif multiple_output_files.lower() in ['true', 't']:
        if output_path is None:
            sys.stderr.write('Error! -m parameter can not be used if -o parameter is omitted!')
            sys.exit(0)
        multiple_output_files = True
        if not os.path.exists(output_path) or os.path.isfile(output_path):
            sys.stderr.write('Error! The output file should be an existing directory '
                  'when multiple_output_files argument is set as "true"!')
            sys.exit(0)
    else:
        sys.stderr.write('Error! The multiple_output_files argument value contains an unknown value:\n')
        sys.stderr.write(str(multiple_output_files))
        sys.stderr.write('Allowed values are "true" and "false". If the argument is not defined, '
                         '"false" will be used as default.')
        sys.exit(0)

    if multiple_output_files and search_subdirectories:
        sys.stderr.write('Error! -u and -m parameters can not both be true!')
        sys.exit(0)


    out_lines = []
    input_file_number = 0
    for file_path in input_file_list:
        if input_extension:
            if not file_path.endswith(input_extension):
                continue

        if debug:
            sys.stderr.write('Starting to read file ')
            sys.stderr.write(file_path)
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

            if debug:
                if j % 1000000 == 0:
                    sys.stderr.write('{0}M lines read...\n'.format(j/1000000))

            matching_line = False

            # Check if any condition is met:
            if line_number is not None:
                if j in line_number:
                    matching_line = True

            if not matching_line and line_start is not None:
                for ls in line_start:
                    if line.startswith(ls):
                        matching_line = True
                        break

            if not matching_line and line_end is not None:
                for le in line_end:
                    if line.endswith(le):
                        matching_line = True
                        break

            if not matching_line and line_contains is not None:
                for lc in line_contains:
                    if lc in line:
                        matching_line = True
                        break

            # Specify the line number to catch
            if matching_line:
                line_numbers_to_extract.add(j + delay)

            # Capture the line
            if j in line_numbers_to_extract:
                out_lines.append(file_name_string + line)
                matching_lines_found += 1

        file_handle.close()
        input_file_number += 1

        # Write output file if multiple_output_files is true
        if multiple_output_files:
            current_out_path = os.path.join(output_path, file_path)
            write_output(current_out_path, out_lines)
            out_lines = []

    # Write output file if multiple_output_files is false
    if not multiple_output_files:
        write_output(output_path, out_lines)

    if output_path is not None:
        sys.stderr.write('{0} matching lines found from {1} files!\n'.format(len(out_lines), input_file_number))


def read_multiple_input_parameters(in_path):
    """Reads a file containing a list of input parameters

    Exits program if unable to open input file or input file is empty.

    :param
    in_path: path to a file

    :return: list of parameters read from input file
    """
    try:
        in_handle = open(in_path)
    except:
        sys.stderr.write('\n\nError! Unable to open parameter file:\n{0}\n'.format(in_path))
        sys.stderr.write('Reason:\n\n')
        raise
    parameters = []
    for line in in_handle:
        line = line.strip()
        if not line: continue
        parameters.append(line)
    if len(parameters) == 0:
        sys.stderr.write('Error! Parameter file is empty:\n{0}\n'.format(in_path))
        sys.exit(0)

    return parameters


def write_output(output_path, out_lines):
    """Writes the output file

    :param
    output_path: Path to output file.
    out_lines: List of lines to be written,
    scanned_files: Number of input files scanned for this output file.
    """
    if output_path is not None:
        try:
            output_handle = open(output_path, 'w')
        except OSError:
            sys.stderr.write('Error! Unable to open output file:\n{0}\n'.format(output_path))
            sys.exit(1)
        output_handle.write('\n'.join(out_lines))
        output_handle.close()
    else:
        for line in out_lines:
            print line


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
                         args.search_subdirectories,
                         args.debug)
