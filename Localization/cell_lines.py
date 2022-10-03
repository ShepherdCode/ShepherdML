class Cell_Lines():
    def get_count():
        return 15
    def get_mapping():
        line_to_index={}
        line_to_index['A549']    =0
        line_to_index['H1.hESC'] =1
        line_to_index['HeLa.S3'] =2
        line_to_index['HepG2']   =3
        line_to_index['HT1080']  =4
        line_to_index['HUVEC']   =5
        line_to_index['MCF.7']   =6
        line_to_index['NCI.H460']=7
        line_to_index['NHEK']    =8
        line_to_index['SK.MEL.5']=9
        line_to_index['SK.N.DZ'] =10
        line_to_index['SK.N.SH'] =11
        line_to_index['GM12878'] =12
        line_to_index['K562']    =13
        line_to_index['IMR.90']  =14
        return line_to_index
    def get_ordered_list():
        ordered_list = \
        ['A549',\
         'H1.hESC',\
         'HeLa.S3',\
         'HepG2',\
         'HT1080',\
         'HUVEC',\
         'MCF.7',\
         'NCI.H460',\
         'NHEK',\
         'SK.MEL.5',\
         'SK.N.DZ',\
         'SK.N.SH',\
         'GM12878',\
         'K562',\
         'IMR.90']
        return ordered_list