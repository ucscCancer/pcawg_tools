#!/usr/bin/env python2.7

"""Tool to read a TCGA Variant Call Format (VCF) file and output an
equivalent file with a different header.

Returns exit code 1 for bad parameters and 2 for header errors detected.
"""


import argparse
import sys

import yaml


__version__ = '1.0.0'


def main():
    args = parse_args()
    with open(args.parameter_file_path) as yaml_file:
        args.parameter_map = yaml.load(yaml_file)
    # TODO: Configure logging
    errors = run(args)
    if errors:
        sys.exit(2)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_file_path', help='the VCF to read')
    parser.add_argument('output_file_path', help='the VCF to write')
    parser.add_argument('parameter_file_path', help='the YAML with details')
    args = parser.parse_args()
    return args


def run(args):
    """Main entry point for testing and higher-level automation"""
    CONFIG = args.parameter_map['config']
    fixed_headers = CONFIG['fixed_headers']
    with open(args.input_file_path) as fin:
        with open(args.output_file_path, 'w') as fout:
            write_fixed_headers(fout, fixed_headers)
            write_sample_lines(fout, CONFIG, args.parameter_map['samples'])
            errors = process_headers(fin, fout, fixed_headers)
            for raw_line in fin:
                fout.write(raw_line)
    return errors


def write_fixed_headers(fout, fixed_headers):
    for name, ignored, value in fixed_headers:
        write_meta_line(fout, name, value)


def write_sample_lines(fout, config, samples):
    SAMPLE_LINE_FORMAT = '##' + config['sample_line_format'].replace(' ', '')
    for id, params in samples.items():
        sample_line = SAMPLE_LINE_FORMAT.format(
            id=id, **dict(params, **config['fixed_sample_params'])
        )
        write_stripped_line(fout, sample_line)


def process_headers(fin, fout, fixed_headers):
    """Keep processing until we write the data header line."""
    filtered_headers = set(item[0] for item in fixed_headers)
    filtered_headers.add("SAMPLE")
    expected_values = {
        name: value for name, asserted, value in fixed_headers if asserted
    }
    errors = False
    for raw_line in fin:
        if raw_line.startswith('##'):
            # TODO: This will break if the metadata header is bad.
            name, value = raw_line[2:].rstrip().split('=', 1)
            if name in filtered_headers:
                if name in expected_values:
                    if value != expected_values[name]:
                        errors = True
                        # TODO: propper logging
                        sys.stderr.write(
                            'tcga-vcf-reheader: mismatch {}={}\n'.format(
                                name, value
                            )
                        )
            else:  # Just some other header...
                fout.write(raw_line)
        else:
            break
    fout.write(raw_line)  # raw_line should now be the data header line.
    return errors


def write_meta_line(fout, name, value):
    fout.write('##{}={}\n'.format(name, value))


def write_stripped_line(fout, line):
    """Just adds the newline."""
    fout.write(line)
    fout.write('\n')


if __name__ == '__main__':
    main()
