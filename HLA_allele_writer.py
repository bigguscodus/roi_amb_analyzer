def hla_allele_writer(allele_lists:tuple, amb_allele_file:str, uniq_allele_file:str):
    '''This function takes tuple, that contains two lists and  writes that list by element if user defines path '''
    if isinstance(amb_allele_file, str):
        with open(amb_allele_file,'w+') as out:
            for allele_set in allele_lists[0]:
                out.write(f'{allele_set}' + '\n')
    if isinstance(uniq_allele_file, str):
        with open(uniq_allele_file, 'w+') as out:
            for allele_set in allele_lists[0]:
                out.write(f'{allele_set}' + '\n')