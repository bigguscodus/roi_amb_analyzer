# Name of the project: 
*unsolvable ambiguity generator*
## Description:
HLA typing is a necessary step in the selection of donors for organ and tissue transplantation. HLAs are proteins located on the surface of cells that provide the presentation of antigen fragments to immune cells. HLAs are located in the specialized [IMGT database **hla.dat**](https://github.com/ANHIG/IMGTHLA), currently numbering about 24,000 variants of HLA genes. HLA typing is done through NGS sequencing and subsequent read alignment on the database. Moreover, different test systems cover different regions of the HLA genes, as a result of which there are natural limitations of the method due to the fact that some of the alleles cannot be discriminated because the differences between them are in regions uncovered by the test system. Identification of such restrictions is a necessary step in the verification of test systems, which must be carried out every time after quarterly database updates.
### Files for work:
⋅⋅* roi_amb_analyzer.py
⋅⋅* HLA_secondary_functions.py
⋅⋅* hla.dat
⋅⋅* test_hla_secondary_functions.py
#### Examples:
`python3 roi_amb_analyzer.py -i hla.dat -g HLA-A -e 2 3 -r 2'
'python3 roi_amb_analyzer.py -i hla.dat -g HLA-A -e 2 3 -r 2 -a amb.txt -u uni.txt'
