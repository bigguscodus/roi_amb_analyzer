from Bio import SeqIO
import itertools
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('locus', help="HLA prefix and gene. Example: HLA-A", type=str)
parser.add_argument('exons', help='Enter the numbers of exons, separated by space. Example: 2 3 4', nargs='+')
args = parser.parse_args()


def roi_bio_analyzer(locus_name, exons, ins_amb='insoluble_ambiguities', dec_uni=' decidable_uniqueness'):
    insoluble_ambiguities = []  # лист хранит тьюплы неразличимых аллелей
    decidable_uniqueness = []  # лист хранит тьюплы различимых аллелей
    allele_names = []  # лист имен аллелей
    total_exon_names_seqs = []  # лист для хранение отображений последовательностей на экзоны
    hla = SeqIO.parse('hla.dat', 'imgt')  # создает генератор
    for allele in hla:
        if allele.description.split(',')[0].startswith(locus_name) and len(
                allele.seq) > 1:  # сортировка по названию и существующей последовательности
            allele_names.append(
                allele.description.split(',')[0])  # захват классификационного названия аллели из описания
            exon_seqs = {f.qualifiers['number'][0]: allele.seq[f.location.start:f.location.end] for f in allele.features
                         if f.type == 'exon'}  # номер экзона: экзонная последовательность, если тип экзон
            total_exon_names_seqs.append(exon_seqs)
    if len(allele_names) == 0:
        return f'{locus_name} is not a part of the HLA system, please double check or contact authors'
    allele_exon = dict(zip(allele_names, total_exon_names_seqs))
    for name in itertools.combinations(allele_names, 2):  # все возможные комбинации 2 аллелей
        flags = []  # лист хранит экзонные флаги на отличие
        for exon in exons:
            if exon not in allele_exon[name[0]] and exon not in allele_exon[name[1]]:
                continue  # экзона нет нигде, на нет и суда нет
            if exon in allele_exon[name[0]] and exon in allele_exon[name[1]] and allele_exon[name[0]][exon] == \
                    allele_exon[name[1]][exon]:
                flags.append(False)  # экзоны есть и одинаковые, нельзя различить
            if exon in allele_exon[name[0]] and exon in allele_exon[name[1]] and allele_exon[name[0]][exon] != \
                    allele_exon[name[1]][exon]:
                flags.append(True)  # экзоны есть и они не одинаковые, можно различить
            if exon not in allele_exon[name[0]] or exon not in allele_exon[name[1]]:
                flags.append(True)  # экзон только в 1 аллели, можно различить
        if sum(flags) > 0:
            decidable_uniqueness.append(name)  # если есть хоть 1 флаг различия, то различить можно
        else:
            insoluble_ambiguities.append(name)  # если сумма равна нули, значит лист содержит только False
    with open(dec_uni, 'a') as out:  # запись файла различимых аллелей
        for elm in decidable_uniqueness:
            out.write(f'{elm}' + '\n')
    with open(ins_amb, 'a') as out:  # запись файла неразличимых аллелей
        for elm in insoluble_ambiguities:
            out.write(f'{elm}' + '\n')
    return "The program has finished it's work"


print(roi_bio_analyzer(locus_name=args.locus, exons=args.exons))
