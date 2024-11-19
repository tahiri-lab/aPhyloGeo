import glob
import os
import sys
import json
import statistics as st
import subprocess
from collections import defaultdict
from io import StringIO
from itertools import combinations
from pathlib import Path


import Bio.SeqIO
from Bio import AlignIO
from Bio.Align import PairwiseAligner
from Bio.Align.Applications import ClustalwCommandline, MafftCommandline
from Bio.Seq import Seq
import pandas as pd

from .params import Params
from .multiProcessor import Multi

import itertools
import textdistance as td


class Alignment:
    """
    Class that contains the data of a multiple sequence alignment.
    """

    def __init__(self, alignment_method: str, msa):
        """
        init method
        args:
            alignment_method (str): The method used to align the sequences
            msa (dict): The multiple sequence alignment data
        """
        self.type = "Alignment"
        self.alignment_method = alignment_method
        self.msa = msa

    def to_dict(self):
        """
        Method that converts the alignment data to a dictionary
        args:
            None
        """
        # Convert the msa data to string format
        msa_str_dict = {key: self.msa_to_string(self.msa[key]) for key in self.msa}
        return {"type": self.type, "alignment_method": self.alignment_method, "msa": msa_str_dict}

    @staticmethod
    def msa_to_string(msa_obj):
        """
        Method that converts the alignment data to a string

        Args:
            MSA (AlignIO): The multiple sequence alignment data
        """
        return "\n".join([f">{record.id}\n{str(record.seq)}" for record in msa_obj])

    @staticmethod
    def msa_from_string(msa_str):
        """
        Method that converts the alignment string to MSA

        Args:
            MSA string
        """
        return AlignIO.read(StringIO(msa_str), "fasta")

    @classmethod
    def from_dict(cls, d):
        # Convert msa data back to its original format
        msa_dict = {key: cls.msa_from_string(d["msa"][key]) for key in d["msa"]}
        return cls(d["alignment_method"], msa_dict)

    def save_to_json(self, filename):
        """
        Method that saves a json sequence file

        Args:
            filename (str): The name of the json file to save
        """
        with open(filename, "w") as f:
            json.dump(self.to_dict(), f)

    @classmethod
    def load_from_json(cls, filename):
        """
        Method that loads a json sequence file

        Args:
            filename (str): The name of the json file to load

        Returns:
            alignment class
        """
        with open(filename, "r") as f:
            data = json.load(f)
            return cls.from_dict(data)

    @classmethod
    def from_json_string(cls, json_string):
        """
        Method that loads an Alignment object from a JSON string.

        Args:
            json_string (str): The JSON string to load.

        Returns:
            An instance of the Alignment class.
        """
        data = json.loads(json_string)
        return cls.from_dict(data)

    @classmethod
    def from_fasta_file(cls, filename, alignment_method):
        """
        Method that loads a sequence from a fasta file

        Args:
            filename (str): The name of the fasta file to load
            alignment_method: The method used to align the sequences

        Returns:
            _type_: MSA
        """
        with open(filename, "r") as f:
            msa = AlignIO.read(f, "fasta")
        msa_dict = {msa.id: msa}
        return cls(alignment_method, msa_dict)


class AlignSequences:
    """
    Class that perform a heuristic Multiple Sequence Alignement and windi from a single fasta file.
    """

    def __init__(self, sequences, makeDebugFiles=False):
        """
        Constructor if the alignment object.
        All parts of the process are available as variables.

        Inputs:
            sequences (dict) key = Sequence ID, value = Seq()
            makeDebugFiles (Boolean) if True, will create a folder with all the intermediate files


        args:
            self.sequences (Dictionary) Read data from the fasta file
                key = Sequence ID
                value = Seq()

            self.centroidKey (String) Specimen ID of the centroid sequence

            self.centroidSeq (Seq()) The centroid sequence

            self.aligned (Dictionary) Sequences after pairwise alignement
                key = Sequence ID (String) the specimen's ID
                value = Sequence (Seq()) the specimen's DNA sequence

            self.heuristicMSA (Dictionary) Sequences after star style alignement
                key = Sequence ID (String) the specimen's ID
                value = Sequence (Seq()) the specimen's DNA sequence

            self.windowed (Dictionary of Dictionary) Sequences after being windowed
                key = Window ID (String) as "startPosition_endPosition"
                value = (Dictionary)
                    key = Sequence ID (String) the specimen's ID
                    value = Sequence (Seq()) the specimen's DNA sequence

            ##todo
            self.msa (AlignIO())
            self.rate_similarity
            self.method_similarity
        """

        self.sequences = sequences
        self.makeDebugFiles = makeDebugFiles

    def align(self) -> Alignment:
        """
        Method that align sequences

        Returns:
            Alignment: The alignment object
        """
        if Params.alignment_method == "0":
            heuristicMSA = self.sequences
            
        elif Params.alignment_method == "1":
            centroidKey = self.getSequenceCentroid()[0]
            centroidSeq = self.sequences.pop(centroidKey)
            aligned = self.alignSequencesWithPairwise(centroidKey, centroidSeq)
            if Params.fit_method == "1":
                heuristicMSA = self.starAlignement(centroidKey, aligned)
            elif Params.fit_method == "2":
                heuristicMSA = self.narrowFitPairwise(aligned)

        elif Params.alignment_method == "2":
            heuristicMSA = self.muscleAlign()

        elif Params.alignment_method == "3":
            heuristicMSA = self.clustalAlign()

        elif Params.alignment_method == "4":
            heuristicMSA = self.mafftAlign()

        else:
            raise ValueError("Invalid alignment method")
        [os.remove(file) for file in glob.glob("aphylogeo/bin/tmp/*.fasta")]  # Remove temp fasta files
        windowed = self.slidingWindow(heuristicMSA)
        self.msa = self.makeMSA(windowed)
        self.alignment = Alignment(Params.alignment_method, self.msa)
        return self.alignment

    def getSequenceCentroid(self):
        """
        Method that picks the centroid sequence in a dictionary of sequences

        variables:
            list (list) Needed variable format for using Multiprocessor
        return: a list of:
            resultKey (String)
            resultSum (int)
        """
        seqs = list((self.sequences).items())
        print("\nSearching for the centroid")

        # formats input as a list for the Multiprocess
        seq_pairs_swapped = [((i[0][1], i[0][0]), (i[1][1], i[1][0])) for i in combinations(seqs, r=2)]
        seq_pairs = [list(sum(i, ())) for i in seq_pairs_swapped]

        # starts all the processes
        align_scores = Multi(seq_pairs, self.ScoreSingle).processingLargeData()

        # formats the multiprocess output back in a dictionnary
        scores = defaultdict(list)
        for pair in align_scores:
            for i in range(0, 2):
                scores[pair[i]].append(pair[2])

        centroid_accession = max({key: st.median(val) for key, val in scores.items()})
        centroid_score = st.median(scores[centroid_accession])

        print(f"The centroid is '{centroid_accession}' with a median score of {centroid_score}\n")

        return [centroid_accession, centroid_score]

    def ScoreSingle(self, args):
        """
        Method the gets only the score of a couple of sequence regarding the
        pairwise alignement

        Args: a list:
            seqA    (Seq()) Sequence A; considered the refenrence
            seqAID  (String) Specimen A's ID
            seqB    (Seq()) Sequence B
            seqBID  (String) Speciment B's ID

        return:
            seqAID  see above
            seqBID  see above
            score   (float) the resulting score of this couple of alignement
        """
        seqA, seqB = args[0], args[2]
        seqAID, seqBID = args[1], args[3]
        aligner = PairwiseAligner()
        score = aligner.score(seqA, seqB)
        return (seqAID, seqBID, score)

    def alignSequencesWithPairwise(self, centroidKey, centroidSeq):
        """
        Method that aligns multiple DNA sequences.
        The first speciment of the dataset is used as the main pivot.
        This method uses parrallel computing.

        Variables:
            seqs (Dictionary) see self.sequences
            list (list) Needed variable format for using Multiprocessor
            result (list) output of all the processes

        Return:
            resultList (Dictionary) see self.aligned

        """
        print("\nStarting sequence alignement")
        seqs = self.sequences

        seq_pairs = []
        for seqXID in seqs.keys():
            seq_pairs.append([centroidKey, centroidSeq, seqXID, seqs[seqXID]])

        align_scores = Multi(seq_pairs, self.alignSingle).processingLargeData()
        aligned = {}

        # reformats the output in a  dictionnary
        for i in align_scores:
            temp = {}
            temp[i[2]] = Seq((i[1][0]))
            temp[i[0]] = Seq((i[1][1]))
            aligned[str(i[0] + " vs " + i[2])] = temp

        # JUST TO MAKE THE DEBUG FILES
        if self.makeDebugFiles:
            directory = os.path.abspath("./debug/1_alignSequences")
            os.makedirs(directory, exist_ok=True)
            # os.mkdir("./debug/1_alignSequences")
            for w in aligned.keys():
                self.dictToFile(aligned[w], str("1_alignSequences/" + w), ".fasta")
        # JUST TO MAKE THE DEBUG FILES

        return aligned

    def muscleAlign(self):
        """Method to perform a multiple DNA sequence alignment using Muscle Algorithm

        Returns (Dict):
            heuristicMSA
                - Keys: accession ID
                - Values: Aligned sequences
        """
        try:
            if sys.platform == "win32":
                muscle_exe = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "\\aphylogeo\\bin\\muscle5.1.win64.exe"
                out_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "\\aphylogeo\\bin\\tmp"
            elif (sys.platform == "linux1") | (sys.platform == "linux2") | (sys.platform == "linux") | (sys.platform == "darwin"):
                muscle_exe = r"aphylogeo/bin/muscle5.1.linux_intel64"
                out_dir = r"aphylogeo/bin/tmp/"
            in_file = Params.reference_gene_filepath
            out_file = os.path.splitext(os.path.basename(Params.reference_gene_filepath))[0]
            out_fullname = str(out_dir + out_file + "_muscle_aligned.fasta")
            process = subprocess.Popen([muscle_exe, "-align", in_file, "-output", out_fullname])
            process.wait()
            records = Bio.SeqIO.parse(out_fullname, "fasta")
            return {rec.id: str(rec.seq) for rec in records}
        # If the muscle executable is not found
        except FileNotFoundError:
            print("Muscle executable not found inside the bin folder. Please check the documentation for installation instructions.")
            sys.exit(1)            
            
    def clustalAlign(self):
        """Method to perform a multiple DNA sequence alignment using ClustalW2 Algorithm

        Returns (Dict):
            heuristicMSA
                - Keys: accession ID
                - Values: Aligned sequences
        """
        try:
            if sys.platform == "win32":
                clustal_exe = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "\\aphylogeo\\bin\\clustalw2.exe"
                fasta_out = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "\\aphylogeo\\bin\\tmp\\clustal_alignment.fasta"
            elif (sys.platform == "linux1") | (sys.platform == "linux2") | (sys.platform == "linux") | (sys.platform == "darwin"):
                clustal_exe = r"aphylogeo/bin/clustalw2"
                fasta_out = r"aphylogeo/bin/tmp/clustal_alignment.fasta"
            in_file = Params.reference_gene_filepath
            clustalw_cline = ClustalwCommandline(clustal_exe, infile=in_file, outfile=fasta_out, output="FASTA")
            out, err = clustalw_cline()
            records = Bio.SeqIO.parse(fasta_out, "fasta")
            return {rec.id: str(rec.seq) for rec in records}
        # If the clustalw executable is not found
        except FileNotFoundError:
            print("ClustalW2 executable not found inside the bin folder. Please check the documentation for installation instructions.")
            sys.exit(1)
            
    def mafftAlign(self):
        """Method to perform a multiple DNA sequence alignment using MAFFT Algorithm

        Returns (Dict):
            heuristicMSA
                - Keys: accession ID
                - Values: Aligned sequences
        """
        try:
            if sys.platform == "win32":
                mafft_exe = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "\\aphylogeo\\bin\\mafft-win\\mafft.bat"
            elif (sys.platform == "linux1") | (sys.platform == "linux2") | (sys.platform == "linux") | (sys.platform == "darwin"):
                mafft_exe = r"aphylogeo/bin/mafft-linux64/mafft.bat"
            in_file = Params.reference_gene_filepath
            mafft_cline = MafftCommandline(mafft_exe, input=in_file)
            out, err = mafft_cline()
            fasta_io = StringIO(out)
            records = Bio.SeqIO.parse(fasta_io, "fasta")
            return {rec.id: str(rec.seq).upper() for rec in records}
        # If the mafft executable is not found
        except FileNotFoundError:
            print("MAFFT executable not found inside the bin folder. Please check the documentation for installation instructions.")
            sys.exit(1)
            
    def alignSingle(self, args):
        """
        Method that aligns two DNA sequences using the pairwise2 algorithm.

        Args:
            args (list)
                scID    (String) The centroid sequence ID to compare to
                sc      (Seq()) The DNA sequence of scID
                seqBID  (String) The sequence ID to compare with
                seqB    (Seq()) The DNA sequence of seqBID

        Return: (list)
            seqBID see above
            aligned (List) The list containing the results.
            scID see above
        """
        scID, seqBID = args[0], args[2]
        sc, seqB = args[1], args[3]

        # Must return Alignment at index [0] of BioPairwiseAlignment
        aligner = PairwiseAligner()
        aligned = aligner.align(sc, seqB)[0]
        return [seqBID, aligned, scID]

    def narrowFitPairwise(self, aligned):
        """Fit length of a centroid sequence and its pairwise aligned sequences

        The length of each sequence from the pairwise alignment are set equal by
        inserting dash (-) in most appropriate location of a given sequence.

        args:
            alignment: dict of nested dict
                {accession couple #1 : {Centroid Acc:Centroid Aligned Seq, Non-centroid Acc #1: non-centroid Aligned Seq #1},
                ... ,
                {accession couple #n : {Centroid Acc:Centroid Aligned Seq, Non-centroid Acc #n: non-centroid Aligned Seq #n}}

        Returns:
            A dictionary of all accessions and their fitted aligned sequences.
        """
        seqs = self.getAlignSeqs(aligned)
        max_len = max(self.getAlignSeqLens(aligned))
        for nucleo_i in range(0, max_len):
            for seq_i in range(0, len(seqs)):
                if self.isCurrentCharDash(seqs, seq_i, nucleo_i):
                    seqs = self.insertDashToShorterSeq(seqs, nucleo_i, aligned)
        seqs = self.appendDashToShorterSeqs(seqs, max_len)
        return self.mergeFitPairwise(aligned, seqs)

    def getAlignSeqs(self, aligned):
        """Extract all sequences aligned using a pairwise alignment

        args:
            alignment: see fitPairwise(alignment) docstring

        Returns:
            List of sequences aligned through pairwise alignment
        """
        seqs = []
        for alignment in aligned:
            seqs.append([str(seq) for seq in aligned[alignment].values()])
        return list(sum(seqs, []))

    def getAlignSeqLens(self, aligned):
        """Get length of all sequences aligned using a pairwise alignment

        args:
            alignment: see fitPairwise(alignment) docstring

        Returns:
            List of the length of each aligned sequences
        """
        return [len(seq) for seq in self.getAlignSeqs(aligned)]

    def getAlignCouple(self, aligned):
        """Get nested couple accessions and their respective sequences

        args:
            alignment: see fitPairwise(alignment) docstring

        Returns:
            List of paired accessions and their aligned sequences
        """
        return [val for val in list(aligned.values())]

    def extractOneAlignAcc(self, aligned, nest_ord=0):
        """Extract the accession from a nested alignment couple

        args:
            alignment: see fitPairwise(alignment) docstring
            nest_ord (int) optional:
                The position of the nested accessions (Default = 0 (centroid), 1 (aligned sequence))

        Returns:
            The list of either centroid (nest_ord = 0 (Default)) or non-centroid (nest_ord = 1)
            accessions of a group of sequences aligned throug pairwise alignment.
        """
        try:
            return [list(i)[nest_ord] for i in self.getAlignCouple(aligned)]

        # Return the centroid sequence if an invalid position is queried
        except IndexError:
            return [list(i)[0] for i in self.getAlignCouple(aligned)]

    def isCurrentCharDash(self, seqs, seq_i, ch_i):
        """Assess whether the character at current cursor position is a dash

        args:
            seqs  (list): aligned sequences to fit
            seq_i (int):  index of the current sequence
            ch_i  (int):  index of the currenct character

        Returns:
            True if the current character assessed is a dash, False otherwise
        """
        try:
            return seqs[seq_i][ch_i] == "-"
        except IndexError:
            return False

    def insertDashToShorterSeq(self, seqs, ch_i, aligned):
        """Insert a dash at the current position of a sequence

        Insert a dash (-) character in a sequence if its length is shorter
        than the longest one in the group of aligned sequence.

        args:
            seqs  (list):    aligned sequences to fit
            seq_i (int):     index of the current sequence

        Returns (List):
            - The fitted sequences of a pairwise alignment
        """
        for seq_j in range(0, len(seqs)):
            try:
                if (len(seqs[seq_j]) < max(self.getAlignSeqLens(aligned))) & (seqs[seq_j][ch_i] != "-"):
                    seqs[seq_j] = seqs[seq_j][:ch_i] + "-" + seqs[seq_j][ch_i:]
            except IndexError:
                seqs[seq_j] = seqs[seq_j][:ch_i] + "-"
        return seqs

    def mergeFitPairwise(self, aligned, seqs):
        """Generate a dictionary of all accessions and their fitted sequences

        args:
            alignment: see fitPairwise(alignment) docstring
            seqs  (list):    aligned sequences to fit

        Returns (Dict):
             Group of accessions and their fitted sequences from a pairwise alignment
        """
        centroid = {list(set(self.extractOneAlignAcc(aligned)))[0]: seqs[0]}
        non_centroid = dict(zip(self.extractOneAlignAcc(aligned, 1), seqs[1::2]))
        return centroid | non_centroid

    def appendDashToShorterSeqs(self, seqs, max_len):
        """Append dash to all sequences shorter than the longest one from a list of sequences

        args:
            seqs, list:  List of fitted sequences post pairwise alignment
            max_len int: Length of the longest aligned sequence, including the blank/dash

        Returns:
            List of sequences with dash appended where applicable
        """
        return [f"{str(seq):-<{max_len}}" for seq in seqs]

    def starAlignement(self, centroidKey, aligned):
        """
        Method that combs through all the pairwise alignments couples and makes
        it so that every sequenced is aligned with every other sequences. If a
        "-" is found in the seqA of a pair, but not another, it is inserted
        into every other ones.

        Example:
            pair1:          pair2:

            seqA1: TACTAC   seqA2: TAC-TAC
            seqB1: TACTAC   seqB2: TACTTAC

            becomes:
            seqA1: TAC-TAC   seqA2: TAC-TAC
            seqB1: TAC-TAC   seqB2: TACTTAC

            and outputs:
            seqA : TAC-TAC   #now combines SeqA1 and SeqA2
            seqB1: TAC-TAC
            seqB2: TACTTAC

            then, we compare the aligned set with the next pair:

            SeqA : TAC-TAC  seqA3: TACTA-C
            seqB1: TAC-TAC  seqB3: TACTAAC
            seqB2: TACTTAC

            wich makes:
            SeqA : TAC-TA-C  seqA3: TAC-TA-C
            seqB1: TAC-TA-C  seqB3: TAC-TAAC
            seqB2: TACTTA-C

            and outputs:
            SeqA : TAC-TA-C     #now combines SeqA1, SeqA2 and SeqA3
            seqB1: TAC-TA-C
            seqB2: TACTTA-C
            seqB3: TAC-TAAC

            over and over again

        Return:
            starAlign (dict) see self.heuristicMSA
        """
        scKey = centroidKey
        starAlign = {}

        for k in aligned.keys():
            # couple is SeqA and SeqB of a pairwise alignement
            couple = aligned[k]

            a = list(couple.keys())
            a.remove(scKey)
            sNewKey = a[0]  # SeqB ID, *not* the reference

            starAlign[scKey] = str(couple[scKey])  # SeqA, the reference
            starAlign[sNewKey] = str(couple[sNewKey])  # SeqB, *not* the reference

            if len(starAlign) > 2:
                starAlign = self.merge(starAlign, scKey, sNewKey, centroidKey)

            starAlign["temp"] = starAlign[scKey]  # SeqA, the *old* reference
        starAlign.pop("temp")

        # JUST TO MAKE THE DEBUG FILES
        if self.makeDebugFiles:
            directory_path = "./debug/2_starAlignement"
            os.makedirs(directory_path, exist_ok=True)
            self.dictToFile(starAlign, "2_starAlignement/starAligned", ".fasta")
        # JUST TO MAKE THE DEBUG FILES

        return starAlign

    def merge(self, result, k1, k2, centroidKey):
        """
        Method that loops through each position of two strings ans compares the Chars.

        Arguments:
            result (dict) the dictionnary of objects to compare;
                contains only object that have already been aligned + a new pair to align
                can be refered to as "aligned set"
            k1 (String) The Key of the object we want compared
            k2 (String) The Key of the object we want compared
        Variables:
            minLen  (int)   The number of char in the smallest of the two strings
            pos     (int)   The char position at wich we are now; it loops
            nChar   (char)  The char from k1
            tChar   (char)  The char from K2
            keylist (list)  Ultimatly, contains all the keys of the object that need to change
        Return:
            result (dict)   The same object we started with, but with one more aligned pair inside.
        """
        newRef = result[k1]
        tempRef = result["temp"]
        minLen = 1
        pos = 0
        while pos < minLen:
            # The sequence length could change at each iteration
            # these assignemnts needs to be done each loop to get
            # to the true end and not to skip the '-' we just added
            newRef = result[k1]
            tempRef = result["temp"]
            minLen = min(len(newRef), len(tempRef))

            nChar = newRef[pos]
            tChar = tempRef[pos]

            if nChar != tChar:
                # - found in the new ref, change all but the 2 new elements
                if nChar == "-":
                    keyList = list(result.keys())
                    keyList.remove(k1)
                    keyList.remove(k2)

                # - found in the old ref; change the 2 new elements
                elif tChar == "-":
                    keyList = [k1, k2]
                else:
                    errStr = str(
                        'Alignement error. Merge() found "'
                        + str(nChar)
                        + '" and "'
                        + str(tChar)
                        + '" '
                        + "at position "
                        + str(pos)
                        + " of two versions of the centroid "
                        + "sequences\n"
                        + "Please check the previous methods"
                        + " and ensure the pairwise alignemnt is correct"
                        + "\nCentroid ID: "
                        + str(centroidKey)
                        + "\nPairwise seq ID"
                        + " last inserted: "
                        + str(k2)
                    )

                    raise Exception(errStr)

                result = self.insertDash(result, pos, keyList)

            pos += 1

        return result

    def insertDash(self, dict, pos, keyList):
        """
        Method that inserts a "-" at [pos] in a string at every Key in a dict

        Arguments:
            dict    (dict)  contains many objects as:
                - key = (string)
                - values = (string)
            pos     (int)   the char position at wich to insert
            keyList (list)  list of keys of objects to modify
        Variables:
            char    (char)  The char to insert
        Return:
            dict    (dict)  The same we started with, with the modifications
        """
        for k in keyList:
            char = "-"
            s = dict[k]
            s = s[:pos] + char + s[pos:]
            dict[k] = s
        return dict

    def slidingWindow(self, heuristicMSA, optimized=True):
        """
        Method that slices all the sequences in a dictionary to a specific window (substring)

        Example:
            step_size=3
            window_size=5

            123 : CGGCTCAGCT  -->   123_3_7 : GCTCA
            456 : TAGCTTCAGT  -->   456_3_7 : GCTTC

        Args:
            alignedSequences (Dictionary)
                - Key (String) is the ID of the specimen
                - Data (Seq(String)) is the specimen's DNS sequence
            others* (var) see param.yaml

        Return:
            resultDict (Dictionary)
                Key is originalKey_i_j
                    originalKey = the name of the key before the window
                    i = The starting position of the window, relative to the original sequence
                    j = The ending position of the window, relative to the original sequence
        """
        step = Params.window_size

        windowed_alignment = dict()
        seq_len = max([len(h) for h in heuristicMSA.values()])

        paddedMSA = {key: str(val).ljust(seq_len, "-") for key, val in heuristicMSA.items()}

        if optimized:
            for i in range(0, seq_len, step):
                if i + step < seq_len:
                    windowed_alignment[f"{i}_{i + step - 1}"] = {key: val[i : i + step - 1] for key, val in paddedMSA.items()}
                    combinations = itertools.combinations(windowed_alignment[f"{i}_{i + step - 1}"].values(), 2)
                    df = pd.DataFrame(list(combinations))
                    if Params.rate_similarity > self.similarity(df):
                        windowed_alignment.pop(f"{i}_{i + step - 1}")
                else:
                    windowed_alignment[f"{i}_{seq_len-1}"] = {key: val[i : i + seq_len - 1] for key, val in paddedMSA.items()}
                    combinations = itertools.combinations(windowed_alignment[f"{i}_{seq_len-1}"].values(), 2)
                    df = pd.DataFrame(list(combinations))
                    if Params.rate_similarity > self.similarity(df):
                        windowed_alignment.pop(f"{i}_{seq_len-1}")
        else:
            for i in range(0, seq_len, step):
                if i + step < seq_len:
                    windowed_alignment[f"{i}_{i + step - 1}"] = {key: val[i : i + step - 1] for key, val in paddedMSA.items()}
                else:
                    windowed_alignment[f"{i}_{seq_len-1}"] = {key: val[i : i + seq_len - 1] for key, val in paddedMSA.items()}

            # if self.rate_similarity[0] < self.similarity(df):
            #    windowed_alignment.pop(f"{i}_{seq_len-1}")
            #    print(self.rate_similarity)

        # JUST TO MAKE THE DEBUG FILES
        if self.makeDebugFiles:
            directory_path = "./debug/3_slidingWindow"
            os.makedirs(directory_path, exist_ok=True)
            for w in windowed_alignment.keys():
                self.dictToFile(windowed_alignment[w], "3_slidingWindow/" + w, ".fasta")
        # JUST TO MAKE THE DEBUG FILES

        return windowed_alignment

    def dictToFile(self, dict, filename, ext):
        """
        Debuging method that creates files from a dictonnary of sequences.
        File is put in the debug file of the cwd

        args:
            dict        (dict)      the objects to write in the file
                - key = (string)
                - values = (string)
            filename    (String)    the name of the future file
            ext         (String)    the file extension

        return:
            dict        (dict)       see dict from arguments
        """
        dir = "./debug"
        if not os.path.exists(dir):
            os.mkdir(dir)

        f = open(dir + "/" + filename + ext, "w")
        for key in dict.keys():
            f.write(">" + str(key) + "\n")
            f.write(str(dict[key] + "\n"))
        return dict

    def makeMSA(self, windowed):
        """
        Method that create a dictionnary of Multiple Sequence Alignment(MSA)
        objects from bioPython. Each entry in the dictionnary is a MSA object
        of a single sliding window.

        returns: 
            msaSet (dict)
                - key (String) - the window name
                - value (AlignIO) - the MSA object
        """
        msaSet = {}
        for windowSet in windowed.keys():
            data = ""
            window = windowed[windowSet]
            for seq in window.keys():
                data += str(">" + seq + "\n" + window[seq] + "\n")
            msaSet[windowSet] = AlignIO.read(StringIO(data), "fasta")

        return msaSet

    def fileToDict(filename, ext):
        """
        Method that reads a fasta file and returns a dictionnary of Seq objects

        arguments:
            filename    (String)    
                the name of the file
            ext         (String)    
                the file extension

        return:
            dict        (dict)
                - key = (string)
                - values = (string)
        """
        f = open(Path(filename + ext), "r")
        dict = {}
        for line in f:
            if line[0] == ">":
                key = line[1:-1]
                dict[key] = ""
            else:
                dict[key] += line[:-1]
        f.close()
        return dict

    def fileToAlignIO(filename, ext):
        """
        Method that reads a fasta file and returns a AlignIO object

        arguments:
            filename    (String)    the name of the file
            ext         (String)    the file extension

        return:
            alignIO     (AlignIO)   the AlignIO object
        """
        f = open(Path(filename + ext), "r")
        data = ""
        for line in f:
            data += line
        f.close()
        alignIO = AlignIO.read(StringIO(data), "fasta")
        return alignIO

    def similarity(self, df):
        """
        Method that compute similarity in all string in dataframe (df)
        Method source:
        https://yassineelkhal.medium.com/the-complete-guide-to-string-similarity-algorithms-1290ad07c6b7

        arguments:
            method (int) similarity method
                1: Hamming distance
                2: Levenshtein distance
                3: Damerau-Levenshtein distance
                4: Jaro similarity
                5: Jaro-Winkler similarity
                6: Smith–Waterman similarity
                7: Jaccard similarity
                8: Sørensen-Dice similarity
                --NOT IMPLEMENTED--
                9: Tversky similarity
                10: Overlap similarity
                11: Cosine similarity
                12: N-gram similarity
                13: Ratcliff-Obershelp similarity
                14: Longest common substring/subsequence similarity
            df    (dataframe)    dataframe of all string sequances

        return:
            percentage similarity     (float)   the similarity between the dataframe of all strings elements
        """
        rateSimilarity = 0.0
        if Params.method_similarity[0] == "1":
            rateSimilarity = df.apply(lambda x: td.hamming.normalized_similarity(x[0], x[1]), axis=1).mean() * 100
        elif Params.method_similarity[0] == "2":
            rateSimilarity = df.apply(lambda x: td.levenshtein.normalized_similarity(x[0], x[1]), axis=1).mean() * 100
        elif Params.method_similarity[0] == "3":
            rateSimilarity = df.apply(lambda x: td.damerau_levenshtein.normalized_similarity(x[0], x[1]), axis=1).mean() * 100
        elif Params.method_similarity[0] == "4":
            rateSimilarity = df.apply(lambda x: td.jaro(x[0], x[1]), axis=1).mean() * 100
        elif Params.method_similarity[0] == "5":
            rateSimilarity = df.apply(lambda x: td.jaro_winkler(x[0], x[1]), axis=1).mean() * 100
        elif Params.method_similarity[0] == "6":
            rateSimilarity = df.apply(lambda x: td.smith_waterman(x[0], x[1]), axis=1).mean() * 100
        elif Params.method_similarity[0] == "7":
            rateSimilarity = df.apply(lambda x: td.jaccard(x[0], x[1]), axis=1).mean() * 100
        elif Params.method_similarity[0] == "8":
            rateSimilarity = df.apply(lambda x: td.sorencen(x[0], x[1]), axis=1).mean() * 100

        return rateSimilarity
