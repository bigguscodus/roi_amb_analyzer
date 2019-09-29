from Bio import SeqIO
import re

allele_names = []  # лист имен аллелей
total_exon_names_seqs = []  # лист для хранение отображений последовательностей на экзоны
hla = SeqIO.parse('hla.dat', 'imgt')  # создает генератор
for allele in hla:
    allele_names.append(allele.description.split(',')[0])  # захват классификационного названия аллели из описания
    exon_seqs = []
    for feature in allele.features:
        if feature.type == 'exon':  # захват экзонов
            slices = re.findall('[0-9]*', str(feature.location))[1:4:2]  # захват концов экзонных отрезков
            exon_seqs.append(allele.seq[int(slices[0]):int(slices[1])])  # слайс аллельной последовательности по экзонам
    exon_names = [str(feature) + ' exon' for feature in
                  range(1, len(exon_seqs) + 1)]  # формирование one-based последовательности экзонов
    exon_names_seqs = dict(zip(exon_names, exon_seqs))  # отображение экзонных последовательностей на номера экзонов
    total_exon_names_seqs.append(exon_names_seqs)
allele_info = dict(zip(allele_names, total_exon_names_seqs))  # отображение экзонов на имена аллелей
#Вывести имена всех аллелей
#for allele_name in allele_info:
    #print(allele_name)
#Вывести имена и все экзонные последовательности
#for allele_name, value in allele_info.items():
    #print(allele_name, value)
#Адресный запрос
#print(allele_info['HLA-A*01:01:01:01']['1 exon'])
