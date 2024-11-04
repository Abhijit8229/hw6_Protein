"""
15-110 Hw6 - Protein Sequencing Project
Name:
AndrewID:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    s = ""
    try:
        with open(filename, 'r') as f:
            for line in f:
                # print(line)
                s+=line.strip()
        return s
    except IOError:
        print("Error: could not read file " + filename)



'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    # print(dna,startIndex)
    rna = dna.replace("T","U")
    lt = []
 
    for i in range(startIndex,len(rna),3):
        # print(rna[i:i+3])
        codon = rna[i:i+3]
        if(len(codon)<3):
            break
        lt.append(codon)
        if codon in {'UAA', 'UAG', 'UGA'}:
            break
        
    
    # print(lt)
    return lt


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):

    import json
    with open(filename, 'r') as f:
        am_ac_dt = json.load(f)
    
    codon_dt = {}
    
   
    for i, j in am_ac_dt.items():
        for k in j:
            codon_u = k.replace('T', 'U')
            codon_dt[codon_u] = i
            
    return codon_dt
   


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    
    protein = []
    
    for i, codon in enumerate(codons):
       
        if i == 0 and codon == "AUG":
            protein.append("Start")
        elif codon in codonD:
            protein.append(f"{codonD[codon]}")
        
        if codon in {'UAA', 'UAG', 'UGA'}:
            break
    
    return protein



'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    with open(dnaFilename, 'r') as dnafile:
        dna_sequence = readFile(dnaFilename).strip()
    codonD = makeCodonDictionary(codonFilename)    
    i = 0
    tp = []
    while i < len(dna_sequence):
        if dna_sequence[i:i+3] == "ATG":
            codonS = dnaToRna(dna_sequence, i)
            p = generateProtein(codonS, codonD)
            tp.append(p)
            i += 3 * len(codonS)
        else:
            i += 1
            
  
    return tp

    


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    set1 = set(tuple(protein) for protein in proteinList1)
    set2 = set(tuple(protein) for protein in proteinList2)
    common = set1.intersection(set2)
    return [list(protein) for protein in common]


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    combined = []
    for protein in proteinList:
        combined.extend(protein)
    return combined



'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    from collections import defaultdict
    aa_dict = defaultdict(int)
    for aa in aaList:
        aa_dict[aa] += 1
    return dict(aa_dict)
   


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    combined1 = combineProteins(proteinList1)
    combined2 = combineProteins(proteinList2)
    aa_dict1 = aminoAcidDictionary(combined1)
    aa_dict2 = aminoAcidDictionary(combined2)
    
    total1 = len(combined1)
    total2 = len(combined2)
    
    differences = []
    
    all_amino_acids = set(aa_dict1.keys()).union(aa_dict2.keys())
    
    for aa in all_amino_acids:
        if aa == "Start" or aa == "Stop":
            continue
        freq1 = aa_dict1.get(aa, 0) / total1
        freq2 = aa_dict2.get(aa, 0) / total2
        if abs(freq1 - freq2) > cutoff:
            differences.append([aa, freq1, freq2])
    
    return differences



'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("Common Proteins:")
    for protein in commonalities:
        if protein != ["Start", "Stop"]:
            print(f"- {' '.join(protein)}")
    
    print("\nAmino Acids with Biggest Differences in Frequency:")
    for aa, freq1, freq2 in differences:
        print(f"Amino Acid: {aa}")
        print(f"  Gene 1 Frequency: {round(freq1 * 100, 2)}%")
        print(f"  Gene 2 Frequency: {round(freq2 * 100, 2)}%")
        print()


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    combined1 = combineProteins(proteinList1)
    combined2 = combineProteins(proteinList2)
    unique_amino_acids = set(combined1 + combined2)
    return sorted(unique_amino_acids)


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    combined = combineProteins(proteinList)
    aa_dict = aminoAcidDictionary(combined)
    total = len(combined)
    return [aa_dict.get(aa, 0) / total for aa in labels]


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    import numpy as np
    x = np.arange(len(xLabels))
    bar_width = 0.4

    plt.bar(x - bar_width / 2, freqList1, width=bar_width, label=label1, edgecolor=edgeList)
    plt.bar(x + bar_width / 2, freqList2, width=bar_width, label=label2, edgecolor=edgeList)

    plt.xticks(ticks=x, labels=xLabels, rotation=90)
    plt.xlabel("Amino Acids")
    plt.ylabel("Frequency")
    plt.legend()
    plt.show()


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    diff_amino_acids = {diff[0] for diff in biggestDiffs}
    return ["black" if aa in diff_amino_acids else "white" for aa in labels]


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():

    
    proteinList1 = [["A", "T", "G"], ["C", "G", "A"], ["T", "A"]]
    proteinList2 = [["A", "C", "G"], ["A", "C", "T"], ["G", "T"]]

    
    labels = makeAminoAcidLabels(proteinList1, proteinList2)
    freqList1 = setupChartData(labels, proteinList1)
    freqList2 = setupChartData(labels, proteinList2)
    differences = findAminoAcidDifferences(proteinList1, proteinList2, 0.005)
    edgeList = makeEdgeList(labels, differences)

   
    displayTextResults([], differences)
    createChart(labels, freqList1, "Gene 1", freqList2, "Gene 2", edgeList=edgeList)


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    test.week1Tests()
    print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    runWeek1()

    ## Uncomment these for Week 2 ##
  
   
   
   
   
   
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
   

    ## Uncomment these for Week 3 ##

    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
 
