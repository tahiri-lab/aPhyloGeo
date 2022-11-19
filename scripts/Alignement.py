import sys
import time  

import os
import Bio as Bio 
from Bio import pairwise2
from Bio.Seq import Seq
from MultiProcessor import Multi
import Params as p

class AlignSequences:

    def __init__(self,sequences):

        self.sequences = sequences
        self.centroidKey = self.getSequenceCentroid()[0]
        self.centroidSeq = self.sequences.pop(self.centroidKey)

        self.aligned = self.alignSequences()
        self.windowed = self.slidingWindow()
        self.msa = ""


    """
    Method that aligns two DNA sequences using an algorithm.

    ex.: alignSingle( "ACTTTCG" , "ACTACG" )
    Could output "ACT--ACG"

    Args:
        original    (String) The DNA sequence to compare to
        next        (String) The DNA sequence to modify using the original
        resultList  (List) The list containing the results.
    Return:
    """
    def alignSingle(self,args):
        scID = args[0]
        sc = args[1]
        seqBID = args[2]
        seqB = args[3]
        aligned = pairwise2.align.globalxx(str(self.centroidSeq), str(seqB), one_alignment_only = True)
        return [seqBID, aligned, scID]

    def ScoreSingle(self,seqA, seqB):
        score = pairwise2.align.globalxx( 
            str(seqA), str(seqB), 
            one_alignment_only = True, 
            score_only=True
            )
        return score

    def getSequenceCentroid(self):
        seqs = self.sequences
        resultKey= ""
        resultSum= sys.maxsize

        for seqID in seqs.keys():
            sum = 0
            for seqID2 in seqs.keys():
                sum += self.ScoreSingle(seqs[seqID], seqs[seqID2])
            if resultSum > sum:
                resultSum = sum; resultKey = seqID

        print("The centroid is \'", resultKey, "\' with a score of ", resultSum,"\n")

        return [resultKey,resultSum]

    """
    Method that aligns multiple DNA sequences.
    The first speciment of the dataset is used as the main pivot.
    This method uses parrallel computing.

    Args:
        sequences (Dictionary) 
            Key String) is the ID of the specimen
            Data (Seq(String)) is the specimen's DNS sequence

    Return:
        resultList (Dictionary) 
            Key is the ID of the specimen
            Data is the specimen's DNS sequence
    """
    def alignSequences(self):
        seqs = self.sequences

        list = []
        for seqXID in seqs.keys():
            list.append( [self.centroidKey, self.centroidSeq, seqXID, seqs[seqXID] ] )

        result = Multi(list,self.alignSingle).processingLargeData()

        resultDict = { result[0][2] : Seq(result[0][1][0].seqB) }

        for i in result:
            resultDict[i[0]] = Seq(i[1][0].seqB)

        #######JUST TO MAKE THE DEBUG FILES
        temp={}
        for i in result:
            temp2 = {}
            temp2[i[0]] = Seq((i[1][0].seqA))
            temp2[i[2]] = Seq((i[1][0].seqB))
            temp[str(i[0]+" vs "+i[2])]=temp2
        os.mkdir("./debug/1_alignSequences")
        time.sleep(1)
        for w in temp.keys():
            self.dictToFile(temp[w],str("1_alignSequences/"+w),".fasta")
        time.sleep(1)

        self.dictToFile(resultDict,"1_alignSequences_OLD",".fasta")
        #vielle version avec les sequences non allignee
        ######JUST TO MAKE THE DEBUG FILES

        return resultDict
    
    """
    Method that slices all the sequences in a dictionary to a specific window (substring)

    ex.:
        step_size=3
        window_size=5

        123 : CGGCTCAGCT  -->   123_3_7 : GCTCA
        456 : TAGCTTCAGT  -->   456_3_7 : GCTTC

    Args:
        alignedSequences (Dictionary)
            Key (String) is the ID of the specimen
            Data (Seq(String)) is the specimen's DNS sequence
        others* (var) see param.yaml

    Return:
        resultDict (Dictionary)
            Key is originalKey_i_j
                originalKey = the name of the key before the window
                i = The starting position of the window, relative to the original sequence
                j = The ending position of the window, relative to the original sequence
    """
    def slidingWindow(self):
        alignedSequences = self.aligned 
        before=time.time()

        windowsDict={}
        
        longKey = max(alignedSequences, key=alignedSequences.get)
        maxLength = len(alignedSequences[longKey])
        
        winSize = p.window_size #longueur
        stepSize = p.step_size #avance de x
        stepStart = 0
        stepEnd = winSize -1

        while stepStart < maxLength:
            if stepEnd > maxLength:
                stepEnd = maxLength
            windowsBySpecies={}
            for key in alignedSequences.keys():
                seq = alignedSequences[key]
                winSeq = seq[stepStart : stepEnd ]
                winKey = str(key) 
                windowsBySpecies[winKey]=Seq(winSeq)
            windowKey = str(stepStart) + "_" + str(stepEnd)
            windowsDict[windowKey] = windowsBySpecies
            stepStart += stepSize
            stepEnd += stepSize

        print(time.time()-before)

        ##############
        os.mkdir("./debug/2_slidingWindow")
        for w in windowsDict.keys():
            self.dictToFile(windowsDict[w],"2_slidingWindow/"+w,".fasta")
        ##############

        return windowsDict

    def dictToFile(self,dict,filename,ext):
        dir = "./debug"
        if not os.path.exists(dir):
            os.mkdir(dir)

        f=open(dir+"/"+filename+ext,"w")
        for key in dict.keys():
            f.write(">"+str(key)+"\n")
            f.write(str(dict[key]+"\n"))
        return dict