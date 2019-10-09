from Bio import SeqIO
import pandas as pd
from itertools import compress
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('locus', help="HLA prefix and gene. Example: HLA-A", type=str)
parser.add_argument('exons', help='Enter the numbers of exons, separated by space. Example: 2 3 4', nargs='+')
args = parser.parse_args()


def roi_bio_analyzer(locus_name, exons, ins_amb='insoluble_ambiguities', dec_uni='decidable_uniqueness'):
    exons = [int(i) for i in exons]
    uniq_alleles = set()  # содержит названия уникальных аллелей
    allele_names = []  # лист имен аллелей
    total_exon_names_seqs = []  # лист для хранение отображений последовательностей на экзоны
    hla = SeqIO.parse('hla.dat', 'imgt')  # создает генератор
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
    solved_alleles = []  # хранит имена неразрешимых неоднозначных аллелей, чтобы лишний раз по ним не проходить
    for allele_name in allele_names:
        if allele_name in solved_alleles:
            continue  # пропускаем, так как аллель уже обработан в прошлом запросе
        for i in df.index:
            df.loc[i,] = df.loc[i,].fillna(value=df[allele_name][i])  # Nan=последовательность аллели
        bool_array = df.apply(lambda x: sum(x == df[allele_name]) == len(df.index)).tolist()  # поэкзонное сравнение
        if sum(bool_array) == 1:  # аллель совпала только с собой
            for elm in compress(allele_name, bool_array):  # находим эту аллель по бинарному ключу
                uniq_alleles.add(elm)  # добавляем в множестно уникальных аллелей
        else:
            amb_alleles = set(compress(allele_names, bool_array))  # множество уникальных аллелей по бинарному ключу
            solved_alleles.extend(amb_alleles)
            with open(ins_amb, 'a') as out:
                out.write(f'{amb_alleles}' + '\n')  # записываем множество неуникальных аллелей
    with open(dec_uni, 'a') as out:
        for elm in uniq_alleles:
            out.write(f'{elm}' + '\n')  # записываем множество уникальных аллелей
    return "The program has finished it's work"


if __name__ == '__main__':
    print(roi_bio_analyzer(locus_name=args.locus, exons=args.exons))
