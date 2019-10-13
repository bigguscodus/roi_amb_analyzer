from Bio import SeqIO
import pandas as pd
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--locus', default='HLA-A', type=str, help="HLA prefix and gene. Example: HLA-A", )
parser.add_argument('-i', '--imgt_file', default='hla.dat', type=str, help='Enter path for the database')
parser.add_argument('-e', '--exons', help='Enter the numbers of exons, separated by space. Example: 2 3 4', nargs='+')
parser.add_argument('-a', '--amb', type=str,
                    help='Enter the path for ambiguous alleles output. Example ~/Documents/amb')
parser.add_argument('-u', '--unq', type=str,
                    help='Enter the path for nonambiguous alleles output. Example ~/Documents/uni')
args = parser.parse_args()


def roi_bio_analyzer(locus_name, exons, ins_amb='insoluble_ambiguities', dec_uni='decidable_uniqueness',
                     base='hla.dat'):
    time_start = time.time()
    exons = [int(i) for i in exons]
    amb_alleles = []  # содержит множество неуникальных аллелей
    uniq_alleles = []  # содержит названия уникальных аллелей
    allele_names = []  # лист имен аллелей
    total_exon_names_seqs = []  # лист для хранение отображений последовательностей на экзоны
    hla = SeqIO.parse(base, 'imgt')  # создает генератор
    for allele in hla:
        if allele.description.split(',')[0].startswith(locus_name) and len(
                allele.seq) > 1:  # сортировка по названию и существующей последовательности
            allele_names.append(
                allele.description.split(',')[0])  # захват классификационного названnpия аллели из описания
            exon_seqs = {int(f.qualifiers['number'][0]): str(allele.seq[f.location.start:f.location.end]) for f in
                         allele.features if
                         f.type == 'exon'}  # номер экзона: экзонная последовательность, если тип экзон
            total_exon_names_seqs.append(exon_seqs)
    if len(allele_names) == 0:
        return f'{locus_name} is not a part of the HLA system, please double check or contact authors'
    allele_exon = dict(zip(allele_names, total_exon_names_seqs))
    df = pd.DataFrame(allele_exon)  # создаем dataframe,  где индексы - экзоны, колонки -аллели
    df = df.loc[exons,].reset_index().drop(labels='index', axis=1)  # ресет индексов для итерации
    solved_alleles = []  # лежат имена аллелей, по которым уже прошлись
    for allele_name in allele_names:
        if allele_name in solved_alleles:
            continue  # аллель уже решена, идем к другой
        dt = df.copy()  # выглядит страшно, но по времени занимает около 0.00004 секунды
        dt.loc[:, allele_name] = dt.loc[:, allele_name].fillna(value='press F to pay respect for Misha')  # пасхалочка
        for i in dt.index:
            dt.loc[i,] = dt.loc[i,].fillna(value=dt[allele_name][i])  # Nan = экзонная последовательность сравнения
            bool = dt.loc[i,] == dt.loc[i, allele_name]  # pandas.Series совпадает/не совпадает allele
            dt = dt.loc[:, bool]  # оставляем колонки, которые только совпали
        set_of_names = set(dt.columns)  # множество из имен совпавщих колонок (аллелей)
        solved_alleles.extend(set_of_names)  # все аллели из множества считаются решенными
        if len(set_of_names) > 1:
            amb_alleles.append(set_of_names)  # если множество больше 1, добавляем в неразрешимые неоднозначности
        else:
            uniq_alleles.append(set_of_names)  # аллель совпала сама с собой
    with open(ins_amb, 'a') as out:
        for elm in amb_alleles:
            out.write(f'{elm}' + '\n')  # записывет неразрешимые неоднозначности
    with open(dec_uni, 'a') as out:
        for elm in uniq_alleles:
            out.write(f'{elm}' + '\n')  # записывает  уникальные аллели
    return f"The program has finished it's work, it has been taking {time.time() - time_start} seconds for finishing"


if __name__ == '__main__':
    print(roi_bio_analyzer(locus_name=args.locus, exons=args.exons, ins_amb=args.amb, dec_uni=args.unq,
                           base=args.imgt_file))
