import sys
import time  

import os
import Bio as Bio 
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from MultiProcessor import Multi
import Params as p

class AlignSequences:

    def __init__(self):

        self.sequences = self.openFastaFile()
        self.centroidKey = self.getSequenceCentroid()[0]
        self.centroidSeq = self.sequences.pop(self.centroidKey)

        self.aligned = self.alignSequences()
        self.heuristicMSA = self.starAlignement()
        self.windowed = self.slidingWindow()
        
        self.msa = ""

    def openFastaFile(self):
        '''
        Reads the .fasta file and read every line to get the
        sequence to analyze.

        Args:
            reference_gene_file (the fasta file to read)

        Return:
            sequences (a dictionnary containing the data from fasta file)
        '''
        sequences = {}
        with open(p.reference_gene_file) as sequencesFile:
            for sequence in SeqIO.parse(sequencesFile,"fasta"):
                sequences[sequence.id] = sequence.seq
        return sequences

    def getSequenceCentroid(self):
        seqs = self.sequences
        resultKey= ""
        resultSum= sys.maxsize

        list = []
        for seqID in seqs.keys():
            for seqID2 in seqs.keys():
                list.append([seqs[seqID],seqID, seqs[seqID2],seqID2])

        results = Multi(list,self.ScoreSingle).processingSmallData()
        rDict = {}

        for tuple in results:
            rDict[tuple[0]] = 0

        for tuple in results:
            rDict[tuple[0]] = rDict[tuple[0]]+tuple[2]

        for k in rDict.keys():
            if rDict[k]<resultSum:
                resultSum = rDict[k]
                resultKey = k

        print("The centroid is \'", resultKey, "\' with a score of ", resultSum,"\n")

        return [resultKey,resultSum]

    def ScoreSingle(self, args):
        seqA= args[0]
        seqAID = args[1]
        seqB = args[2]
        seqBID = args[3]
        score = pairwise2.align.globalxx( 
            str(seqA), str(seqB), 
            one_alignment_only = True, 
            score_only=True
            )
        return (seqAID, seqBID, score)

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

        ####### JUST TO MAKE THE DEBUG FILES ####### 
        temp={}
        for i in result:
            temp2 = {}
            temp2[i[2]] = Seq((i[1][0].seqA))
            temp2[i[0]] = Seq((i[1][0].seqB))
            temp[str(i[0]+" vs "+i[2])]=temp2
        os.mkdir("./debug/1_alignSequences")
        time.sleep(1)
        for w in temp.keys():
            self.dictToFile(temp[w],str("1_alignSequences/"+w),".fasta")
        time.sleep(1)
        ####### JUST TO MAKE THE DEBUG FILES ####### 

        return temp
   
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

    def starAlignement(self):

        scKey = self.centroidKey
        starAlign = {}

        for k in self.aligned.keys():
            couple = self.aligned[k]    #a couple is SeqA and SeqB of a pairwise alignement

            a = list(couple.keys())
            a.remove(scKey)
            sNewKey = a[0]  #SeqB ID, *not* the reference

            starAlign[ scKey ] = couple[ scKey ] #SeqA, the reference
            starAlign[ sNewKey ] = couple[ sNewKey ] #SeqB, *not* the reference

            if len( starAlign ) > 2:
                starAlign = self.merge( starAlign, scKey, sNewKey )
                starAlign = self.equalizeLength(starAlign)

            starAlign[ "temp" ] = starAlign[ scKey ] #SeqA, the *old* reference
        starAlign.pop( "temp" )

        ####### JUST TO MAKE THE DEBUG FILES ####### 
        os.mkdir("./debug/2_starAlignement")
        self.dictToFile(starAlign,"2_starAlignement/starAligned",".fasta")
        ####### JUST TO MAKE THE DEBUG FILES ####### 

        return starAlign
        
        
    def merge(self, result, k1, k2):
        newRef= result[k1]
        tempRef = result["temp"]
        minLen = 1
        pos = 0
        while pos< minLen:
            #The sequence length could change at each iteration
            #these assignemnts needs to be done each loop to get to the true end
            #and not to skip the '-' we just added
            newRef= result[k1]
            tempRef = result["temp"]
            minLen = min( len(newRef), len(tempRef) ) 

            nChar = newRef[pos]
            tChar = tempRef[pos]

            if nChar != tChar:
                if nChar =='-': #- found in the new reference; change all but the 2 new elements
                    keyList = list(result.keys())
                    keyList.remove(k1)
                    keyList.remove(k2)
                elif tChar =='-':   #- found in the old reference; change the 2 new elements
                    keyList=[k1,k2]
                else:
                    errStr = str(
                    "Alignement error. Merge() found \""+str(nChar) + "\" and \"" + str(tChar) + "\" " 
                    "at position " + str(pos) + " of two versions of the centroid sequences\n"+
                    "Please check the previous methods and ensure the pairwise alignemnt is correct"+
                    "\nCentroid ID: "+str(self.centroidKey)+
                    "\nPairwise seq ID last inserted: "+str(k2)
                    )
                    raise Exception(errStr)

                result = self.insertDash(result, pos, keyList)

            pos+=1
        
        return result

    def insertDash(self, dict, pos, keyList):
        for k in keyList:
            char = '-'
            s = dict[k]
            s = s[:pos]+ char +s[pos:]
            dict[k] = s
        return dict

    def equalizeLength(self, unEqualSeqs):
        equalizedSeqs = {}
        maxLen = len( str( max( list( unEqualSeqs.values() ) ) ) ) #number of chars in the longest sequence

        for k in unEqualSeqs.keys():
            equalizedSeqs[k]= Seq(str(unEqualSeqs[k]).ljust(maxLen,'-'))

        return equalizedSeqs

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
        alignedSequences = self.heuristicMSA
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
        os.mkdir("./debug/3_slidingWindow")
        for w in windowsDict.keys():
            self.dictToFile(windowsDict[w],"3_slidingWindow/"+w,".fasta")
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

