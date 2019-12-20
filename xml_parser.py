from xml.etree import ElementTree
def xml_parser(xml_file = 'hla_ambigs.xml'):
    tree = ElementTree.parse(xml_file) # создаем дерево
    root = tree.getroot() # находим корень
    amb_dict = {el.attrib['name']:[] for el in root.iter('{http://www.example.org/ambig-aw}gene')} # {имя гена:х}
    for group in root.iter('{http://www.example.org/ambig-aw}gGroup'): # для каждой группы
        amb = {group.attrib['name']:set()} # {имя группы :[ неразлечимые аллели]}
        for allele in group:
            amb[group.attrib['name']].add(allele.attrib['name']) #  присоедияет имя неразлечимой аллели
        for amb_group in amb: # для каждой группы
            for prefix in amb_dict: # для каждого гена
                if amb_group.startswith(prefix) and len(amb[amb_group])>1: #  если имя группы начинается с имени гена и длиннее 1
                    amb_dict[prefix].append(amb[amb_group]) #   присоединяем лист с сгруппированными неразлечимымы аллелями к гену
    return amb_dict
