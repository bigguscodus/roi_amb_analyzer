import pandas as pd
import re
def hla_allele_solver(data_frame: pd.DataFrame, resolution:str) -> tuple:
    '''This function takes dataframe as an argument and returns tuple with two lists. Each list contains set of allele according to their ambiguity'''
    first_level = r'HLA-[A-Z]*\*[0-1]*'
    second_level = r'HLA-[A-Z]*\*[0-9]*:[0-9]*'
    third_level = r'HLA-[A-Z]*\*[0-9]*:[0-9]*:[0-9]*'
    fourth_level = r'HLA-[A-Z]*\*[0-9]*:[0-9]*:[0-9]*:[0-9]*'
    resolutiion_image = {'1':first_level,'2':second_level,'3':third_level,'4':fourth_level}
    ambiguity_alleles = [] # лист для хранения множеств неразрешимых аллелей
    uniq_alleles = [] # лист для хранения уникальных аллелей
    solved_alleles = [] # лист для хранения решенных аллелей, чтобы лишний раз их не сравнивать со всеми
    for allele_name in data_frame.columns:
        prefix_set = set()
        after_set = set() # хранит множество неразрешимых аллелей после детализации по разрешению
        if allele_name in solved_alleles:
            continue # если аллель мы решили, то дальше с ней работь
        data_temp = data_frame.copy()
        for i in data_temp.index:
            if pd.isna(data_temp.loc[i, allele_name]):
                continue
            data_temp.loc[i,] = data_temp.loc[i,].fillna(
                value=data_temp[allele_name][i])  # где нет сиквенса, вставляем сиквенс из аллели сравнения
            bools = data_temp.loc[i,] == data_temp.loc[
                i, allele_name]  # строим вектор True/False  по совпадающим аллелям
            data_temp = data_temp.loc[:, bools]  # оставляем колонки, которые только совпали
        set_of_names = set(data_temp.columns)
        solved_alleles.extend(set_of_names)
        if len(set_of_names) > 1:
            for allele in set_of_names:
                prefix = re.findall(resolutiion_image[resolution], allele)
                if prefix == []: # если у аллели нет нужного разрешения, оставляем как есть
                    prefix = [allele]
                if prefix[0] not in prefix_set:
                    prefix_set.add(prefix[0])
                    after_set.add(allele)
            ambiguity_alleles.append(after_set)  # если множество больше 1, добавляем в неразрешимые неоднозначности
        else:
            uniq_alleles.append(set_of_names)
    return (ambiguity_alleles,uniq_alleles)