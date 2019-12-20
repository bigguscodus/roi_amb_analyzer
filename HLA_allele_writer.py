def hla_allele_writer(allele_lists:tuple, amb_allele_file:str, uniq_allele_file:str):
    '''This function takes tuple, that contains two lists and  writes that list by element if user defines path '''
    if isinstance(amb_allele_file, str): # если юзер что-то ввел пробуем туда записать
        with open(amb_allele_file,'w+') as out:
            for allele_set in allele_lists[0]: # из листа с неразрешимыми неозднозначностями тащим множества
                out.write(f'{allele_set}' + '\n')
    if isinstance(uniq_allele_file, str):
        with open(uniq_allele_file, 'w+') as out:
            for allele_set in allele_lists[1]: # из листа с уникальными аллелями тащим множества
                out.write(f'{allele_set}' + '\n')