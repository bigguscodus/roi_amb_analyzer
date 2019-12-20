if __name__=='__main__':
    import pandas as pd
    from Bio import SeqIO

    amb_alleles = []  # содержит множество неуникальных аллелей
    uniq_alleles = []  # содержит названия уникальных аллелей
    allele_names = []  # лист имен аллелей
    exons = [2, 3, 4]
    total_exon_names_seqs = dict()  # словарь вида {region:{exon_number:sequence}
    hla = SeqIO.parse('hla.dat', 'imgt')  # парсит hla
    for allele in hla:
        if allele.description.split(',')[0].startswith('HLA-DRB1') and len(allele.seq) > 1:# сортировка по названию и существующей последовательности
            total_exon_names_seqs[allele.description.split(',')[0]] = {int(f.qualifiers['number'][0]): str(allele.seq[f.location.start:f.location.end]) for f in allele.features if f.type == 'exon'}  # наполнение словаря значениями
    df = pd.DataFrame(total_exon_names_seqs)  # создаем dataframe,  где индексы - экзоны, колонки -аллели
    df = df.loc[exons,].reset_index().drop(labels='index', axis=1)  # ресет индексов для итерации
    pd.set_option('display.expand_frame_repr', False)
    wtf = df.loc[:, ['HLA-DRB1*09:01:02:02', 'HLA-DRB1*09:01:02:01', 'HLA-DRB1*09:31', 'HLA-DRB1*09:21']]
    print(wtf)
    print(wtf.apply(lambda x: x == wtf['HLA-DRB1*09:21']))
