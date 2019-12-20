import pandas as pd
from Bio import SeqIO

amb_alleles = []  # содержит множество неуникальных аллелей
uniq_alleles = []  # содержит названия уникальных аллелей
allele_names = []  # лист имен аллелей
exons = [2, 3, 4]
total_exon_names_seqs = []  # лист для хранение отображений последовательностей на экзоны
hla = SeqIO.parse('hla.dat', 'imgt')  # создает генератор
for allele in hla:
    if allele.description.split(',')[0].startswith('HLA-DRB1') and len(
            allele.seq) > 1:  # сортировка по названию и существующей последовательности
        allele_names.append(allele.description.split(',')[0])  # захват классификационного названnpия аллели из описания
        exon_seqs = {int(f.qualifiers['number'][0]): str(allele.seq[f.location.start:f.location.end]) for f in
                     allele.features if f.type == 'exon'}  # номер экзона: экзонная последовательность, если тип экзон
        total_exon_names_seqs.append(exon_seqs)
allele_exon = dict(zip(allele_names, total_exon_names_seqs))
df = pd.DataFrame(allele_exon)  # создаем dataframe,  где индексы - экзоны, колонки -аллели
df = df.loc[exons,].reset_index().drop(labels='index', axis=1)  # ресет индексов для итерации
pd.set_option('display.expand_frame_repr', False)
wtf = df.loc[:, ['HLA-DRB1*09:01:02:02', 'HLA-DRB1*09:01:02:01', 'HLA-DRB1*09:31', 'HLA-DRB1*09:21']]
print(wtf.apply(lambda x: x == wtf['HLA-DRB1*09:01:02:02']))
