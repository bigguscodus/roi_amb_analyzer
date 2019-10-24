from Bio import SeqIO
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--locus', help="HLA prefix and gene. Example: HLA-A", )
parser.add_argument('-i', '--imgt_file', default='hla.dat', type=str, help='Enter path for the database')
parser.add_argument('-e', '--exons', help='Enter the numbers of exons, separated by space. Example: 2 3 4', nargs='+')
parser.add_argument('-a', '--amb',
                    help='Enter the path for ambiguous alleles output. Example ~/Documents/amb')
parser.add_argument('-u', '--unq',
                    help='Enter the path for nonambiguous alleles output. Example ~/Documents/uni')
args = parser.parse_args()


def roi_bio_analyzer(locus_name, exons, unsolvable_ambiguity, decidable_uniqueness, base):
    exons = [int(i) for i in exons]
    ambiguity_alleles = []  # содержит множество неуникальных аллелей
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
        dt = df.copy()
        for i in dt.index:  # для каждого ряда, а ряд это экзонная последовательность
            if pd.isna(dt.loc[i, allele_name]):
                continue  # аллель сравнения неотсеквенирована в этом экзоне, можно и не сравнивать, все равно там может быть что угодно
            dt.loc[i,] = dt.loc[i,].fillna(
                value=dt[allele_name][i])  # где нет сиквенса, вставляем сиквенс из аллели сравнения
            bools = dt.loc[i,] == dt.loc[i, allele_name]  # строим вектор True/False  по совпадающим аллелям
            dt = dt.loc[:, bools]  # оставляем колонки, которые только совпали
        set_of_names = set(dt.columns)  # множество из имен совпавщих колонок (аллелей)
        solved_alleles.extend(set_of_names)  # все аллели из множества считаются решенными
        if len(set_of_names) > 1:
            ambiguity_alleles.append(set_of_names)  # если множество больше 1, добавляем в неразрешимые неоднозначности
        else:
            uniq_alleles.append(set_of_names)  # аллель совпала только сама с собой, она уникальна
    with open(unsolvable_ambiguity, 'a') as out:
        for allele in ambiguity_alleles:
            out.write(f'{allele}' + '\n')  # записывет неразрешимые неоднозначности
    with open(decidable_uniqueness, 'a') as out:
        for allele in uniq_alleles:
            out.write(f'{allele}' + '\n')  # записывает  уникальные аллели
    return ambiguity_alleles  # возвращает неразрешимые аллели


if __name__ == '__main__':
    roi_bio_analyzer(locus_name=args.locus, base=args.imgt_file, exons=args.exons, unsolvable_ambiguity=args.amb,
                     decidable_uniqueness=args.unq)
