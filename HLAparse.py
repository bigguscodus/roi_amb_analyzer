from Bio import SeqIO
def hla_parser(region:str,data='hla.dat') -> dict:
    '''This function takes imgt file and returns a dictionary {allele_name:{exon_number:sequence}'''
    total_exon_names_seqs = dict()  # словарь вида {region:{exon_number:sequence}
    hla = SeqIO.parse(data, 'imgt') # парсит hla
    for allele in hla:
        if allele.description.split(',')[0].startswith(region) and len(
                allele.seq) > 1:  # сортировка по названию и существующей последовательности
            total_exon_names_seqs[allele.description.split(',')[0]] = {
                int(f.qualifiers['number'][0]): str(allele.seq[f.location.start:f.location.end]) for f in
                allele.features if f.type == 'exon'} # наполнение словаря значениями
    return total_exon_names_seqs

