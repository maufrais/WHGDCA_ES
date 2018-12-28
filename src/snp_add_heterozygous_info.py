'''

@author: corinne.maufrais@pasteur.fr
Bioinformatics and Biostatistics Hub C3BI - Institut Pasteur
Paris, France
'''
import sys
import argparse


def read_matrix_header(line):
    fld = line.strip().split('\t')
    l = len(fld)
    strains = {}
    nbcol_info = 26
    if fld[nbcol_info] == 'ref':
        nbcol_info = 27
    for i in range(nbcol_info, l):
        strains[i]={'name':fld[i]}
    return strains, nbcol_info


def write_matrix_header(fhout, line, strains, nbcol_info ):
    """simple analysis without sam and gatk complements"""
    bidule = ''
    fld = line.strip().split('\t')
    l = len(fld)
    #l = len(fld) -1
    list_of_cple = []
    list_of_strain = []
    for i in range(nbcol_info, l):
        list_of_strain.append(strains[i]['name'])
        for j in range(nbcol_info, l):
            if i == j:
                continue
            bidule += strains[i]['name'] + '/'+ strains[j]['name'] + '\t' + 'Tag#' + strains[i]['name'] + '/'+ strains[j]['name'] + '\t'+ 'Position' +'\t'
            list_of_cple.append((strains[i]['name'], strains[j]['name']))
    ori = '\t'.join(fld[0:nbcol_info])
    print >>fhout, "%s\t%s" % (line.strip(), bidule)
    return list_of_strain, list_of_cple

def compare(a, b):
    if a == b:
        if a[0] == a[1]:
            return '0', '0'
        else:
            return 'h', 'h'
    elif a[0] == a[1] and b[0] != b[1]:
        return '0', 'h'
    elif a[0] != a[1] and b[0] == b[1]:
        return 'h', '0'
    else:
        if a[0] == a[1]: 
            return "'2", '2'
        else: 
            return 'H', 'H'  # trimorphic


def analyse_matrix(fhout, line, strains, list_of_strain, list_of_cple, nbcol_info):
    fld = line.strip().split('\t')
    ori = line.strip()
    l = len(fld)
    snp_pos = fld[1]
    chr_name = fld[0].split('-')[0]
    lne = ori

    x = 0
    for i in range(nbcol_info, l):
        for j in range(nbcol_info, l):
            if i == j:
                continue
            cpl = list_of_cple[x]
            x += 1
            res = compare(fld[i], fld[j])
            position = ''
            if res in [('0','h'), ('h', 'h')]:
                position = fld[1]
            lne += '\t' + fld[i] + '/' + fld[j]
            lne += '\t' + res[0] + '/' + res[1]
            lne += '\t' + position
    print >>fhout, lne

def main(fhin, fhout):
    line = fhin.readline()
    strains, nbcol_info = read_matrix_header(line)
    list_of_strain, list_of_cple = write_matrix_header(fhout, line, strains, nbcol_info)
    line = fhin.readline()
    while line:
        analyse_matrix(fhout, line, strains, list_of_strain, list_of_cple, nbcol_info)
        line = fhin.readline()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='snp_add_heterozygous_info',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Add the heteorzygous status of each SNPs.")

    parser.add_argument('-i', '--gtx_matrix', metavar='file',
                        dest='fhin',
                        type=argparse.FileType('r'),
                        help="SNPs shared between strains in a tabulated file. See example",
                        required=True)
    parser.add_argument('-o', '--out', metavar='file',
                        dest='fhout',
                        type=argparse.FileType('w'),
                        help="Output file. SNPs shared between strains in a tabulated file with heterozygous information.",
                        required=True)
        
    args = parser.parse_args()
    main(args.fhin, args.fhout)
    
