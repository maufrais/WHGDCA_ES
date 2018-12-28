'''
@author: corinne.maufrais@pasteur.fr
Bioinformatics and Biostatistics Hub C3BI - Institut Pasteur
Paris, France
'''
import sys
import argparse


def header(line):
    all_snps = {}
    fld = line.strip().split('\t')
    tag_strains = {}
    for i in range(len(fld)):
        if 'Tag#' in fld[i]:
            tag_strains[fld[i]]=i
            all_snps[fld[i]] = {}
    return tag_strains, all_snps


def analyse_matrix(fhin, tag_strains, all_snps):
    line = fhin.readline()
    while line:
        fld = line.strip().split('\t')
        chr_name = fld[0].split('-')[0]
        position = int(fld[1])
        for tag in all_snps.keys():
            col = tag_strains[tag]
            snp = fld[col]
            if snp not in ['0/0', "'2/2", 'H/H']:
                if chr_name in all_snps[tag]:
                    all_snps[tag][chr_name][position]={'snp':snp, 'FT':fld[5]}
                else:
                    all_snps[tag][chr_name] = {position:{'snp':snp, 'FT':fld[5]}}
        line = fhin.readline()
    return all_snps


def print_loh(fhout, all_snps, list_pos=False):
    if list_pos:
        print >>fhout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % ("Cple", "chr_name", "pos_up", "pos_down", "L", "pos_med_up", "pos_med_down", "l", "nb_snps", "pos_snps")
    else:
        print >>fhout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % ("Cple", "chr_name", "pos_up", "pos_down", "L", "pos_med_up", "pos_med_down", "l", "nb_snps")
    for tag in all_snps.keys():
        chromosomes = all_snps[tag].keys()
        chromosomes.sort()
        for chr_name in chromosomes:
            positions = all_snps[tag][chr_name].keys()
            positions.sort()
            pos_up, pos_down = None, None
            tag_up = None
            pos_med_up, pos_med_down = None, None
            tag_med_up = None
            cmpt = 0
            l_snps_h0 = []
            chr_begin = True
            for pos in positions:
                if all_snps[tag][chr_name][pos]['snp'] in ['h/h', '0/h']:
                    if not tag_up:
                        pos_up = pos
                        tag_up = True
                    else:
                        pos_down = pos
                        if cmpt >= 2:
                            if list_pos:
                                print >>fhout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (tag.replace('Tag#', ''), chr_name, pos_up, pos_down, pos_down-pos_up+1, pos_med_up, pos_med_down, pos_med_down-pos_med_up+1, cmpt, l_snps_h0)
                            else:
                                print >>fhout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (tag.replace('Tag#', ''), chr_name, pos_up, pos_down, pos_down-pos_up+1, pos_med_up, pos_med_down, pos_med_down-pos_med_up+1, cmpt)
                        pos_up = pos
                        tag_med_up, pos_med_up, pos_med_down = None, None, None
                        cmpt = 0
                        l_snps_h0 = []
                elif all_snps[tag][chr_name][pos]['snp'] in ['h/0']:
                    if chr_begin and not tag_up:
                        chr_begin = False
                        pos_up = pos
                        tag_up = True
                    if not tag_med_up:
                        pos_med_up = pos
                        tag_med_up = True
                    else:
                        pos_med_down = pos
                    cmpt += 1
                    l_snps_h0.append(pos)
                    strain1 = tag.replace('Tag#', '').split('/')[0]

def main(fhin, fhout, list_pos):
    first_line = fhin.readline()
    tag_strains, all_snps = header(first_line)
    all_snps = analyse_matrix(fhin, tag_strains, all_snps)
    print_loh(fhout, all_snps, list_pos)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='snp_extract_microLOH',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Micro LOH analysis between strains")

    parser.add_argument('-i', '--SNPs_matrix', metavar='file',
                        dest='fhin',
                        type=argparse.FileType('r'),
                        help="SNPs shared between strains in a tabulated file with heterozygous information. See example.",
                        required=True)
    parser.add_argument('-o', '--new_SNPs_matrix', metavar='file',
                        dest='fhout',
                        type=argparse.FileType('w'),
                        help="Output file. ",
                        required=True)
    parser.add_argument("-p", "--list_pos",
                        dest="list_pos",
                        help="Add microLOH positions in the output file. Default: False",
                        action='store_true',
                        default=False)
        
    args = parser.parse_args()
    main(args.fhin, args.fhout, args.list_pos)

