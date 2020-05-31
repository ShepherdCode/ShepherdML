class Bucket():
    '''Utility class for Manber-Myers algorithm.'''

    def __init__(self,prefix,stringT):
        self.prefix = prefix # one or more letters
        self.stringT = stringT # needed for shortcut sort
        self.suffixes = []  # array of int

    def __str__(self):
        viz = ""
        viz = viz + str(self.prefix)
        viz = viz + " "
        viz = viz + str(self.suffixes)
        return (viz)

    def get_prefix(self):
        return self.prefix

    def add_suffix(self,index):
        self.suffixes.append(index)

    def get_suffixes(self):
        return self.suffixes

    def sort_suffixes_shortcut(self):
        self.suffixes.sort(key=self.get_suffix_string)

    def get_suffix_string(self,i):
        offset = i - 1
        suffix_string = self.stringT[offset:]
        return suffix_string
