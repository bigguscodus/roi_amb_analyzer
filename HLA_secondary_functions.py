from Bio import SeqIO
import pandas as pd
import re


def hla_parser(gene: str, data='hla.dat') -> dict:
    '''Input: gene in which we make a search
       data .dat file which contains allele information. Can be obtained from https://github.com/ANHIG/IMGTHLA
Action: parse data from .dat file and put in to the dictionary with following structure {allele_name : {exon_number : exon_sequence}}
Return: a dictionary
'''
    try:
        total_exon_names_seqs = dict()  # словарь вида {region:{exon_number:sequence}
        hla = SeqIO.parse(data, 'imgt')  # парсит hla
        for allele in hla:
            if allele.description.split(',')[0].startswith(gene) and len(
                    allele.seq) > 1:  # сортировка по названию и существующей последовательности
                total_exon_names_seqs[allele.description.split(',')[0]] = {
                    int(f.qualifiers['number'][0]): str(allele.seq[f.location.start:f.location.end]) for f in
                    allele.features if f.type == 'exon'}  # наполнение словаря значениями
        if len(total_exon_names_seqs) < 2:
            print('Check your inputs. There is an error in gene name')
            return None
        return total_exon_names_seqs
    except FileNotFoundError:
        print("Check your inputs. The program can't find a base")
        return None


def hla_filter_framer(allele_exon_dict: dict, exon_numbers: list) -> pd.DataFrame:
    '''Input: A dictionary with following structure {allele_name : {exon_number : exon_sequence}}
       List that contains exon numbers
Action: Convert dictionary to dataframe with allele_names as columns, exons as rows, sequence as value
Return: pd.Dataframe
'''
    if allele_exon_dict is None:
        return None
    exon_numbers = [int(exon_number) for exon_number in exon_numbers]
    df = pd.DataFrame(allele_exon_dict)  # создаем датафрэйм
    df = df.loc[exon_numbers,].reset_index().drop(labels='index', axis=1)  # сбрасываем индексы, для быстрой работы
    return df


def hla_allele_solver(data_frame: pd.DataFrame, resolution: str) -> tuple:
    '''Input: dataframe that contains allele names as columns, exons as rows and sequence as value.
       resolution is the level of detalization.
Action: finds  ambiguity_alleles and unique alleles in dataframe.
Return: tuple with two lists. Each list contains set of allele according to their ambiguity
    '''
    if data_frame is None:
        return None
    first_level = r'HLA-[A-Z1-9]*\*[0-9]*'
    second_level = r'HLA-[A-Z1-9]*\*[0-9]*:[0-9]*'
    third_level = r'HLA-[A-Z1-9]*\*[0-9]*:[0-9]*:[0-9]*'
    fourth_level = r'HLA-[A-Z1-9]*\*[0-9]*:[0-9]*:[0-9]*:[0-9]*'
    resolutiion_image = {'1': first_level, '2': second_level, '3': third_level, '4': fourth_level}
    ambiguity_alleles = []  # лист для хранения множеств неразрешимых аллелей
    uniq_alleles = []  # лист для хранения уникальных аллелей
    solved_alleles = []  # лист для хранения решенных аллелей, чтобы лишний раз их не сравнивать со всеми
    for allele_name in data_frame.columns:
        prefix_set = set()
        after_set = set()  # хранит множество неразрешимых аллелей после детализации по разрешению
        if allele_name in solved_alleles:
            continue  # если аллель мы решили, то дальше с ней работь
        data_temp = data_frame.copy()
        for i in data_temp.index:
            if pd.isna(data_temp.loc[i, allele_name]):
                continue
            data_temp.loc[i,] = data_temp.loc[i,].fillna(
                value=data_temp[allele_name][i])  # где нет сиквенса, вставляем сиквенс из аллели сравнения
            bools = data_temp.loc[i,] == data_temp.loc[
                i, allele_name]  # строим вектор True/False  по совпадающим аллелям
            data_temp = data_temp.loc[:, bools]  # оставляем колонки, которые только совпали
        set_of_names = list(data_temp.columns)
        solved_alleles.extend(set_of_names)
        if len(set_of_names) > 1:
            for allele in set_of_names:
                prefix = re.findall(resolutiion_image[resolution], allele)
                if prefix == []:  # если у аллели нет нужного разрешения, оставляем как есть
                    prefix = [allele]
                if prefix[0] not in prefix_set:
                    prefix_set.add(prefix[0])
                    after_set.add(prefix[0])
            if len(after_set)>1:
                ambiguity_alleles.append(after_set)  # если множество больше 1, добавляем в неразрешимые неоднозначности
        else:
            uniq_alleles.append(set_of_names[0])
    return (ambiguity_alleles, uniq_alleles)


def hla_allele_writer(allele_lists: tuple, amb_allele_file: str, uniq_allele_file: str, gene: str, exon_numbers: list,
                      resolution: str):
    '''Input: tuple, that contains two lists: one with sets of  ambiguity alleles another with sets of unique alleles.
Action: writes each list by element if user defined path.
Return:None
'''
    if allele_lists is None:
        return None
    if amb_allele_file is not None:  # если юзер что-то ввел пробуем туда записать
        with open(amb_allele_file, 'w+') as out:
            out.write(f'#Result for {gene} gene for {exon_numbers} exons with {resolution} resolution level\n')
            for allele_set in allele_lists[0]:  # из листа с неразрешимыми неозднозначностями тащим множества
                out.write(f'{allele_set}\n')
    if uniq_allele_file is not None:
        with open(uniq_allele_file, 'w+') as out:
            out.write(f'#Result for {gene} gene for {exon_numbers} exons with {resolution} resolution level\n')
            for allele_set in allele_lists[1]:  # из листа с уникальными аллелями тащим множества
                out.write(f'{allele_set}\n')


if __name__ == '__main__':
    help(hla_parser)
    help(hla_filter_framer)
    help(hla_allele_solver)
    help(hla_allele_writer)
