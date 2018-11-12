
# coding: utf-8

# In[ ]:


# Dna class. Inherited from str class.
class Dna(str):
    def __init__(self, sequence):
        DNAcharlist = ('A', 'C', 'G', 'T', 'S', 'W', 'R', 'Y', 'K', 'M', 'B', 'D', 'H', 'V', 'N')
        badseq = False
        trueseq = ''
        for char in sequence:
            if char.upper() in DNAcharlist:
                trueseq += char
            else:
                badseq = True
                continue
        self.sequence = trueseq
        if badseq:
            print ('While initializing Dna object non-DNA characters were omitted')
        if len (sequence) == 0:
            print ('Empty DNA sequence was created. Using gc() function will return None')
            
    def __repr__(self):
        return (self.sequence)
    
    def __str__(self):
        return (self.sequence)
            
# transcribe function returns transcription in coordinates [transcrstart : transcrend] as an Rna class object
# whole sequence is transcribed by default
    def transcribe(self, transcr_start = 1, transcr_end = None): 
        if transcr_end == None:
            transcr_end = len(self)
        Rnatext = self.sequence[(transcr_start-1):transcr_end].replace('T', 'U').replace('t', 'u')
        return Rna(Rnatext)

# gc() function returns None if sequence length is 0 or all characters are ambiguous
    def gc(self):
        DNAgcdict = {'A': 0, 'C': 1, 'G': 1, 'T': 0, 'S': 1, 'W': 0}
        notforgc = ['R', 'Y', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
        gcs = 0
        definedbps = len(self.sequence)
        if len(self.sequence) == 0:
            return None
        for char in self.sequence:
            if char in notforgc:
                definedbps -= 1
            else:
                gcs += DNAgcdict[char.upper()]
        if definedbps == 0:
            return None
        else:
            return (round(gcs / definedbps, 2))

# reverse_complement() function returns Dna object of complementary strand in 5' - 3' orientation
    def reverse_complement(self):
        DNArcdict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'S': 'W', 'W': 'S', 'R': 'Y', 'Y': 'R', 
                     'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'}
        rcseq = ''
        for char in self.sequence:
            if char.isupper():
                rcseq = DNArcdict[char] + rcseq
            else:
                rcseq = DNArcdict[char.upper()].lower() + rcseq
        return Dna(rcseq)


# Rna class. Inherited from str class. Same as Dna, except transcribe() function is impossible
class Rna(str):
    def __init__(self, sequence):
        RNAcharlist = ('A', 'C', 'G', 'U', 'S', 'W', 'R', 'Y', 'K', 'M', 'B', 'D', 'H', 'V', 'N')
        badseq = False
        trueseq = ''
        for char in sequence:
            if char.upper() in RNAcharlist:
                trueseq += char
            else:
                badseq = True
                continue
        self.sequence = trueseq
        if badseq:
            print ('While initializing Rna object non-RNA characters were omitted')
        if len (sequence) == 0:
            print ('Empty RNA sequence was created. Using gc() function will return None')
            
    def __repr__(self):
        return (self.sequence)
    
    def __str__(self):
        return (self.sequence)
    
# gc() function returns None if sequence length is 0 or all characters are ambiguous
    def gc(self):
        RNAgcdict = {'A': 0, 'C': 1, 'G': 1, 'U': 0, 'S': 1, 'W': 0}
        notforgc = ['R', 'Y', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
        gcs = 0
        definedbps = len(self.sequence)
        if len(self.sequence) == 0:
            return None
        for char in self.sequence:
            if char in notforgc:
                definedbps -= 1
            else:
                gcs += RNAgcdict[char.upper()]
        if definedbps == 0:
            return None
        else:
            return (round(gcs / definedbps, 2))

# reverse_complement() function returns Rna object of complementary sequence in 5' - 3' orientation
    def reverse_complement(self):
        RNArcdict = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A', 'S': 'W', 'W': 'S', 'R': 'Y', 'Y': 'R', 
                     'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'}
        rcseq = ''
        for char in self.sequence:
            if char.isupper():
                rcseq = RNArcdict[char] + rcseq
            else:
                rcseq = RNArcdict[char.upper()].lower() + rcseq
        return Rna(rcseq)
    
