#!/usr/bin/env python3

import argparse
import re
import pdb


def match_id_to_name(control_file):
    with open(control_file, 'r') as fobj:
        header = fobj.readline()
        id_idx, name_idx = find_id_name_cols(header)
        county_dict = dict()

        line_num = 1
        for line in fobj:
            line = line.split('\t')
            county_id = line[id_idx]
            name = line[name_idx]
            if county_id not in county_dict:
                county_dict[county_id] = name
            elif county_dict[county_id] != name:
                raise RuntimeError('Problem in line {linenum} of {filename}: '
                                   'different name for county ID "{id}" than expected ({name})'.
                                   format(linenum=line_num, filename=control_file, id=county_id, name=county_dict[county_id]))
            line_num += 1

    return county_dict


def write_matched_table(csv_file, county_dict):
    new_csv_file = csv_file.replace('.csv', '_with_county_names.csv')
    with open(csv_file, 'r') as infile, open(new_csv_file, 'w') as outfile:
        header = infile.readline()
        header = header.split(',')
        moves_file_idx = header.index('moves_run_file')
        insert_idx = moves_file_idx + 1
        header.insert(insert_idx, 'county_description')
        outfile.write(','.join(header))

        for line in infile:
            line = line.split(',')
            run_file = line[moves_file_idx]
            county_id = re.search('(?<=moves_pan_)\d+(?=_)', run_file).group()
            county_name = county_dict[county_id]
            line.insert(insert_idx, county_name)
            outfile.write(','.join(line))


def find_id_name_cols(header):
    header = header.split('\t')
    pdb.set_trace()
    id_idx = header.index('CountyID')
    name_idx = header.index('County Description')
    return id_idx, name_idx


def parse_args():
    parser = argparse.ArgumentParser(description='Insert county names into a MOVES output .csv')
    parser.add_argument('output_csv', help='The .csv file containing MOVES output, preprocessed by '
                                           'my sum_emissions SQL script')
    parser.add_argument('control_file', help='The MOVES control file used to generate the county-specific '
                                             'run specs')

    return vars(parser.parse_args())


def driver(output_csv, control_file):
    county_dict = match_id_to_name(control_file)
    write_matched_table(output_csv, county_dict)


def main():
    args = parse_args()
    driver(**args)

if __name__ == '__main__':
    main()
