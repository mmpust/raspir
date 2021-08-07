#!/bin/python

# Marie-Madlen Pust
# last updated: 05 August 2021

import os
import argparse
import pandas as pd


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass


desc = 'Clean your file up and add taxonomy information'
epi = """DESCRIPTION:
Input: merged .CSV file, taxonomy file
Output: final .CSV file
"""

parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)


parser.add_argument('csv_file', metavar='csv_file', type=str,
                    help='Input csv file')
parser.add_argument('out_prefix', metavar='out_prefix', type=str,
                    help='output file prefix')


def move_column_inplace(df, col, pos):
    col = df.pop(col)
    df.insert(pos, col.name, col)
    return df


def process_csv(file_name, out_prefix):
    with open(file_name, 'r') as csv_file, open('taxonomy_file.csv') as tax_list:
        csv_df = pd.read_csv(csv_file, sep=",")
        # remove in-text lines with header information
        csv_df = csv_df[~csv_df['organism'].isin(['organism'])]
        # convert all columns (except for organism names into numeric dtypes)
        cols = csv_df.columns.drop('organism')
        csv_df[cols] = csv_df[cols].apply(pd.to_numeric, errors='coerce')

        id_values_0 = csv_df["organism"]

        id_values_1 = []
        for items in id_values_0:
            if items.startswith('NZ'):
                id_values_1.append(items[14:])
            elif items.startswith('NC'):
                id_values_1.append(items[12:])
            elif items.startswith('CP'):
                id_values_1.append(items[11:])
            elif items.startswith('CR'):
                id_values_1.append(items[11:])
            elif items.startswith('AE'):
                id_values_1.append(items[11:])
            elif items.startswith('AP'):
                id_values_1.append(items[11:])
            elif items.startswith('BK'):
                id_values_1.append(items[20:])
            elif items.startswith('HG'):
                id_values_1.append(items[11:])
            elif items.startswith('FM'):
                id_values_1.append(items[11:])
            elif items.startswith('CF'):
                id_values_1.append(items[11:])
            elif items.startswith('Ca22chr1A_C_albicans_SC5314___organism_'):
                id_values_1.append(items[39:])
            elif items.startswith('CH476594___organism_'):
                id_values_1.append(items[20:])
            elif items.startswith('Chr1_A_fumigatus_Af293___organism_'):
                id_values_1.append(items[34:])
            elif items.startswith('ChrI_A_nidulans_FGSC_A4___organism_'):
                id_values_1.append(items[35:])
            else:
                id_values_1.append(items)

        id_values_2 = []
        for p in id_values_1:
            if p.startswith("_"):
                new_string = p.split("_", 4)[:4]
                id_values_2.append(new_string)
            elif p.startswith("__"):
                new_string = p.split("_", 5)[:4]
                id_values_2.append(new_string)
            elif p.startswith("___"):
                new_string = p.split("_", 6)[:4]
                id_values_2.append(new_string)
            elif p.startswith("____"):
                new_string = p.split("_", 7)[:4]
                id_values_2.append(new_string)
            else:
                new_string = p.split("_", 3)[:2]
                id_values_2.append(new_string)

        id_values_3 = []
        for i in id_values_2:
            i2 = ' '.join(i)
            id_values_3.append(i2)

        csv_df = csv_df.iloc[:, 1:]
        csv_df.insert(0, 'Species', id_values_3)
        csv_df.groupby(by=csv_df.columns, axis=1).sum()

        # remove all rows that sum to zero
        cols_to_sum = csv_df.columns[: csv_df.shape[1]-1]
        csv_df['sum_all'] = csv_df[cols_to_sum].sum(axis=1)
        csv_df = csv_df[csv_df.sum_all != 0]
        csv_df = csv_df.drop('sum_all', 1)

        # grep taxonomy file
        tax_df = pd.read_csv(tax_list, sep=";")
        tax_df = tax_df.astype(str)
        tax_df['Species_new'] = tax_df[['Genus', 'Species']].apply(lambda x: ' '.join(x), axis=1)
        tax_df = tax_df.drop('Species', 1)

        # subset taxonomy file based on samples in csv_df
        tax_df_sub_0 = tax_df[tax_df.set_index(['Species_new']).index.isin(csv_df.set_index(['Species']).index)]
        csv_tax_df_0 = csv_df.merge(tax_df_sub_0, how='left', left_on='Species', right_on='Species_new',)
        csv_tax_df_1 = csv_tax_df_0.drop('Species_new', 1)

        # move columns
        csv_tax_df_2 = move_column_inplace(csv_tax_df_1, 'Domain', 0)
        csv_tax_df_3 = move_column_inplace(csv_tax_df_2, 'Phylum', 1)
        csv_tax_df_4 = move_column_inplace(csv_tax_df_3, 'Class', 2)
        csv_tax_df_5 = move_column_inplace(csv_tax_df_4, 'Order', 3)
        csv_tax_df_6 = move_column_inplace(csv_tax_df_5, 'Family', 4)
        csv_tax_df_7 = move_column_inplace(csv_tax_df_6, 'Genus', 5)
        csv_tax_df_7.sort_values("Species", inplace=True)
        csv_tax_df_8 = csv_tax_df_7.drop_duplicates(keep=False)

        # export file
        outfile = '{}_clean.csv'.format(out_prefix)
        csv_tax_df_8.to_csv(outfile, index=False, header=True, sep=";")


def main():
    """
    Main interface
    """
    # output directory
    outdir = os.path.split(args.out_prefix)[0]
    if outdir != '' and not os.path.isdir(outdir):
        os.makedirs(args.outdir)
    # processing each file
    process_csv(args.csv_file, out_prefix=args.out_prefix)


if __name__ == "__main__":
    args = parser.parse_args()
    main()
