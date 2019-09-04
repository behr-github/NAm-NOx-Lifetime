#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pdb

def convert_value(value):
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def read_csv(csv_file, species):
    with open(csv_file, 'r') as infile:
        header = infile.readline()
        header = header.strip().split(',')
        species_idx = header.index('species_name')
        data = []
        for line in infile:
            line = line.strip().split(',')
            if line[species_idx] == species:
                data.append({k: convert_value(v) for k, v in zip(header, line)})

    return data


def county_to_index(data):
    counties = [row['county_description'] for row in data]
    counties = sorted(list(set(counties)))
    for row in data:
        row['county_index'] = counties.index(row['county_description'])

    # Format the county names a little nicer
    counties = [name.split('-')[1].strip() for name in counties]

    return data, counties


def make_data_array(data, counties):
    # ultimately want three columns: county index, emissions, and month. Will
    # average over years, so start with an n-by-12-by-3 array, with year in the third dim, then
    # average over the third dim
    n_counties = len(counties)
    data_array = np.full((3, 12, 3, n_counties), np.nan, dtype=np.float)
    data_array[:, :, 0, :] = np.resize(np.arange(n_counties), (3, 12, n_counties))

    months = np.arange(12)+1
    months_array = np.repeat(months, n_counties)
    data_array[:, :, 2, :] = np.resize(months_array, (3, 12, n_counties))

    years = np.array([2012, 2013, 2014])

    for row in data:
        month_idx = months == row['emis_month']
        year_idx = years == row['emis_year']
        county_idx = row['county_index']
        data_array[year_idx, month_idx, 1, county_idx] = row['total_emis_kg']

    # need the variable dimension in front before shaping to correctly get it in a
    # variable-by-(county+month) matrix
    data_array = np.mean(data_array, 0).transpose((1, 2, 0)).reshape((3, -1))
    return data_array


def plot_emiss(data_array, counties):
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.35, 0.9, 0.6])
    sc = ax.scatter(data_array[0, :], data_array[1, :], c=data_array[2, :])
    cb = fig.colorbar(sc)

    ax.set_ylabel('NO$_x$ emissions (kg)')
    ax.set_xticks(np.arange(len(counties)))
    ax.set_xticklabels(counties)
    ax.tick_params(axis='x', rotation=90)
    cb.set_ticks(np.arange(1, 13))
    cb.set_label('Month')
    return fig, ax, sc, cb


def driver(csv_file, species='NOx', figure_file='moves_emiss.png'):
    data = read_csv(csv_file, species=species)
    data, counties = county_to_index(data)
    data_array = make_data_array(data, counties)
    fig, _, _, _ = plot_emiss(data_array, counties)
    fig.savefig(figure_file)


def parse_args():
    parser = argparse.ArgumentParser(description='Plot MOVES emissions')
    parser.add_argument('csv_file', help='The .csv file to read from - must have county names inserted with '
                                         '"match_moves_county_ids.py"')
    parser.add_argument('-s', '--species', default='NOx', choices=['NOx','NO','NO2'], help='Which species emissions to plot')
    parser.add_argument('-f', '--figure-file', default='moves_emiss.png', help='The filename to give the figure')

    return vars(parser.parse_args())


def main():
    driver(**parse_args())


if __name__ == '__main__':
    main()
