#!/usr/bin/env python3.6

import random
import os
import argparse
import gzip

import igraph

import numpy
import pandas
from scipy.stats import pearsonr

from cod.util.log.Logger import Logger
from cod.util.time.Timer import Timer
from cod.util.stats.Stats import Stats

#===============================================================================
DESC_COMMENT = "Script to identify cis co-expressed gene domains (COD)."
SCRIPT_NAME = "CODer.py"
#===============================================================================

"""
#===============================================================================
@author: Diogo Ribeiro
@date: 18 Mar 2019
@copyright: Copyright 2019, University of Lausanne
Script to identify cis co-expressed gene domains (COD).
#===============================================================================
"""

class CODer(object):
    
    ## Set default arguments
    DEFAULT_PARAM_PERMUTATIONS = 1000
    DEFAULT_PARAM_WINDOW_SIZE = 1000000
    DEFAULT_PARAM_RANK_APPROACH = 0
    DEFAULT_PARAM_FDR_CUTOFF = -1
    DEFAULT_PARAM_MINIMUM_CORRELATION = float("-inf")
    DEFAULT_PARAM_NEGATIVE_CORRELATION = 2
    DEFAULT_PARAM_CONSECUTIVE_RANK = 0
    DEFAULT_PARAM_CORRECTION_PROCEDURE = 0
    DEFAULT_PARAM_WANTED_SEED = -1
    DEFAULT_PARAM_WRITE_LOG = ""
    DEFAULT_PARAM_LOW_MEM = 0
    DEFAULT_PARAM_NULL_MAX_DISTANCE = 100.0
    DEFAULT_PARAM_DETERMINE_TSS = 0
    DEFAULT_PARAM_RERUN = 0
    DEFAULT_PARAM_VARIBLE_MATCH = ""
    DEFAULT_PARAM_VARIBLE_MATCH_PROP = 0.1
    
    ## Other constants
    # Correlation method to be used
    CORRELATION_METHOD = "pearson"
    # Number of column from input matrix from which individuals/samples start
    FIRST_SAMPLE_COLUMN = 6
    # Flag for performing 'whole genome' pvalue correction globally or only using first rank
    GLOBAL_PVALUE_CORRECTION = 0
    
    ## Output files    
    OUTPUT_COP_RANDOM_RESULTS = ""
    OUTPUT_COP_RAW_RESULTS = "CODer_raw_results.bed"
    OUTPUT_COP_DISTANCE_CONTROL_NULL = "CODer_distance_controlled_null.bed"
    OUTPUT_COP_DISTANCE_VARIABLE_CONTROL_NULL = "CODer_distance_and_variable_controlled_null.bed"
    OUTPUT_PVALUE_FDR = "CODer_pvalue_fdr_rank0.bed"
    OUTPUT_COD_IDENTIFICATION_COPS = "CODer_cod_identification_cops.bed"
    OUTPUT_COD_IDENTIFICATION_PER_PHENOTYPE = "CODer_cod_identification_per_phenotype.bed"
    OUTPUT_COD_IDENTIFICATION_PER_COD = "CODer_cod_identification_per_cod.bed"
    OUTPUT_LOG = "CODer.log"
        

    def __init__(self, matrixFile, outputFolder, permutations, rankApproach, \
                  windowSize, fdrCutoff, minimumCorrelation, negativeCorrelation, 
                  consecutiveRank, correctionProcedure, nullMaxDistance, determineTSS, wantedSeed, 
                  verbosityLevel, writeLog, lowMem, rerun, variableMatch, variableMatchProp):
        
        Timer.get_instance().start_chrono()

        self.matrixFile = matrixFile
        self.outputFolder = outputFolder
        self.permutations = permutations
        self.rankApproach = rankApproach
        self.windowSize = windowSize
        self.fdrCutoff = fdrCutoff
        self.minimumCorrelation = minimumCorrelation
        self.negativeCorrelation = negativeCorrelation
        self.consecutiveRank = consecutiveRank
        self.correctionProcedure = correctionProcedure
        self.nullMaxDistance = nullMaxDistance
        self.determineTSS = determineTSS
        self.wantedSeed = wantedSeed        
        self.verbosityLevel = verbosityLevel
        self.writeLog = writeLog
        self.lowMem = lowMem
        self.rerun = rerun
        self.variableMatch = variableMatch
        self.variableMatchProp = variableMatchProp
        
        Logger.get_instance().set_level( verbosityLevel)

        if not os.path.isdir(self.outputFolder) :
            os.mkdir(self.outputFolder)
            
        if self.writeLog:
            self.outLog = open( self.writeLog, "w")

        # Decide if randomising seed or using a given seed
        if self.wantedSeed == CODer.DEFAULT_PARAM_WANTED_SEED:
            seed = numpy.random.randint(0,100000000)
        else:
            seed = self.wantedSeed

        # Set seed
        numpy.random.seed( seed)

        t = "Matrix file: %s\nPermutations: %s\nRank approach: %s\nWindow size: %s\
            \nFDR cutoff: %s\nMinimum correlation: %s\nNegative correlation: %s\nConsective Rank: %s\
            \nCorrection procedure: %s\nMaximum null distance: %s\nDetermine TSS: %s\nSeed: %s\nLow mem: %s\
            \nRerun: %s\nVariable match: %s\nVariable match proportion: %s\nOutput folder: %s\n" \
              % (self.matrixFile, self.permutations, self.rankApproach, self.windowSize,
                 self.fdrCutoff, self.minimumCorrelation, self.negativeCorrelation, self.consecutiveRank,
                 self.correctionProcedure, self.nullMaxDistance, self.determineTSS, seed, self.lowMem, 
                 self.rerun, self.variableMatch, self.variableMatchProp, self.outputFolder)
        if self.writeLog: self.outLog.write(t + "\n\n")
        Logger.get_instance().info(t)
        
        
    def q_value(self, x):
        """
        Returns qvalue p-value correction. 
        """                
        return Stats.qvalue(x)


    def bh_correct(self, x):
        """
        Returns Benjamini-Hochberg FDR procedure. 
        """                
        return Stats.multiple_test_correction(x, meth="fdr_bh")


    def read_matrix_file(self):
        """
        Read input matrix file using pandas.
        Transposes matrix and separates info from data values.
        Randomises values for each phenotype (e.g. equivalent to shuffling individuals IDs)
        Processes each chromosome separately.
        
        Example format:
            #chr    start   end     gene    info    strand  sample1  sample2
            GL000192.1      495564  495565  ENSG00000277655.1_5     L=451;T=unprocessed_pseudogene;R=GL000192.1:493155-495565;N=AC245407.1  -       0.4       0
            GL000193.1      81322   81323   ENSG00000280081.3_5     L=2485;T=lincRNA;R=GL000193.1:49232-81323;N=LINC01667   -       0       1.34

        Note: there is an alternative format in which the header labels are as follows:
            #chr    start   end     id      gid     strd    sample1    sample2
        """

        
        # Read matrix as dataframe
        if self.matrixFile.endswith(".gz"):
            fileHandler = gzip.open(self.matrixFile, "rt")
        else:
            fileHandler = open(self.matrixFile, "r")
        matrix = pandas.read_csv(fileHandler, sep="\t")

        # Get set of chromosomes in file
        chromosomes = set(matrix["#chr"])

        # Get list of individuals/samples (this should be the same for all chromosomes)
        individualsList = matrix.iloc[:,CODer.FIRST_SAMPLE_COLUMN:].columns
        
        # Process each chromosome
        matrixData = {} # key -> chromosome, val -> row:individuals, col:phenotype; dataframe regarding that chromosome containing phenotype values
        randMatrixData = {} # key -> chromosome, val -> dict. key -> permutation index, val -> row:individuals, col:phenotype; 
                                                                                        #dataframe contaning shuffled phenotype values
        
        matrixHeaders = {} # key -> chromosome, val -> row:info, col:phenotype; dataframe regarding that chromosome containing phenotype info
        phenotypeInfo = {} # key -> chromosome, val -> row:phenotype, col:info; dataframe regarding that chromosome containing phenotype info

        for chro in chromosomes:
            chroMatrix = matrix[matrix["#chr"] == chro].copy()

            # Get list of phenotypes (rows in input file, different depending on chromosome)
            try:
                # input file may have different header
                phenotypeList = chroMatrix["id"]
            except:
                phenotypeList = chroMatrix["gene"]

            # Get phenotype info (non-transposed)
            phenotypeInfo[chro] = chroMatrix.iloc[:,:CODer.FIRST_SAMPLE_COLUMN]
            phenotypeInfo[chro] = phenotypeInfo[chro].set_index(phenotypeList.copy())
            
            # Get the part with info about the phenotype, transpose matrix
            matrixHeaders[chro] = chroMatrix.iloc[:,:CODer.FIRST_SAMPLE_COLUMN].transpose()
            matrixHeaders[chro].columns = phenotypeList.copy()

            # Get the part of matrix with values, transpose
            matrixData[chro] = chroMatrix.iloc[:,CODer.FIRST_SAMPLE_COLUMN:].transpose()
            matrixData[chro].columns = phenotypeList.copy()      

            t = "read_matrix_file: Chromosome %s. Individuals: %s, Phenotypes: %s" % (chro, matrixData[chro].shape[0], matrixData[chro].shape[1])
            Logger.get_instance().info(t)
            if self.writeLog: self.outLog.write(t + "\n\n")

            if self.lowMem == 0:
                # Randomise matrix X times
                Logger.get_instance().debug( "read_matrix_file: Randomizing expression matrix %s times" % (self.permutations))
     
                randMatrixData[chro] = {} 
                for perm in range(self.permutations):                        
                    randMatrixData[chro][perm] = matrixData[chro].copy(deep=True)
    
                    # permute order of the data (expression values per individual), then attribute the same indexes as in the real case                    
                    randMatrixData[chro][perm] = randMatrixData[chro][perm].reindex(index=numpy.random.permutation(individualsList.copy(deep=True)))
                    randMatrixData[chro][perm].index = individualsList
                
        Logger.get_instance().info( "read_matrix_file: Input matrix with %s chromosomes." % (len(matrixData)))

        # Simple dictionary with phenotype info
        globalMatrixHeaders = {} # key -> phenotype name, val -> info
        for chro in matrixHeaders:
            for idx, row in matrixHeaders[ chro].iteritems():
                globalMatrixHeaders[ idx] = row
        
        self.matrixData = matrixData
        self.randMatrixData = randMatrixData
        self.matrixHeaders = matrixHeaders
        self.phenotypeInfo = phenotypeInfo
        self.chromosomes = chromosomes
        self.globalMatrixHeaders = globalMatrixHeaders
    

    def perform_correlation(self, x, y):
        """
        Performs correlation between two arrays

        @param x and y: arrays to be compared.
        
        @return correlation coefficient and the sign of the correlation (positive/negative)
            If --negativeCorrelation is used, return the square or absolute of the correlation coefficient.
        """
        
        corr = pearsonr(x, y)[0]
        if corr < 0: sign = "-"
        else: sign = "+"
        
        if self.negativeCorrelation == 1:
            # return square
            return numpy.power(corr, 2), sign
        elif self.negativeCorrelation == 2:
            # return absolute
            return abs(corr),sign
        else:
            return corr, sign


#     def calculate_coefficient_variation(self):
#         '''
#         Calculate coefficient of variation (aka relative standard deviation; standard deviation divided by the mean) 
#         for each gene in the gene expression matrix.
#         '''
# 
#         coefVar = {} # key -> pheno ID, val -> coefficient of variation
#         for chro in self.matrixData:
#             for pheno,_ in self.matrixData[chro].iteritems():
#                 vals = list(self.matrixData[chro][pheno])
#                 coefVar[pheno] = numpy.std(vals) / numpy.mean(vals)
# 
#         self.coefVar = coefVar


    def read_real_coordinates(self):
        '''
        Get the real start and end coordinates for each gene using the "info" field, without accounting for strand.
        
        Example format:
            #chr    start   end     gene    info    strand  sample1  sample2
            GL000192.1      495564  495565  ENSG00000277655.1_5     L=451;T=unprocessed_pseudogene;R=GL000192.1:493155-495565;N=AC245407.1  -       0.4       0
            GL000193.1      81322   81323   ENSG00000280081.3_5     L=2485;T=lincRNA;R=GL000193.1:49232-81323;N=LINC01667   -       0       1.34

        Note: there is an alternative format in which the header labels are as follows:
            #chr    start   end     id      gid     strd    sample1    sample2
        '''

        realCoordinates = {} # key -> phenotype name, val -> (start,end)
        phenoStart = {} # key -> phenotype name, val -> coordinate to be considered as the start (e.g. TSS of a gene)
        
        if self.determineTSS:
            # if the coordinates are not provided by the "info" column but rather coded in start, end and strand
            for chro in self.matrixHeaders:
                tsss = []
                for pheno in self.matrixHeaders[chro]:
                    if abs(self.matrixHeaders[chro][pheno]["end"] - self.matrixHeaders[chro][pheno]["start"]) <= 1:
                        Logger.get_instance().warning( "read_real_coordinates: the start and end coordinates do not differ. These do not resemble real coordinates: start: %s end: %s" % \
                                                      (self.matrixHeaders[chro][pheno]["start"], self.matrixHeaders[chro][pheno]["end"]) )
                    if self.matrixHeaders[chro][pheno]["strd"] == "+":
                        tss = self.matrixHeaders[chro][pheno]["start"]
                    elif self.matrixHeaders[chro][pheno]["strd"] == "-":
                        tss = self.matrixHeaders[chro][pheno]["end"]
                    else:
                        raise Exception("read_real_coordinates: strand should either be '+' or '-'.")
                    
                    realCoordinates[pheno] = (self.matrixHeaders[chro][pheno]["start"],self.matrixHeaders[chro][pheno]["end"])
                    phenoStart[pheno] = tss
                    tsss.append(tss)
                    
                # Replace the start with the tss, to be used when finding cis phenotypes
                self.phenotypeInfo[chro]["start"] = tsss
        else:
            # if the real coordinates are provided by the "info" column
            for chro in self.matrixHeaders:
                for pheno in self.matrixHeaders[chro]:
                    try:
                        info = self.matrixHeaders[chro][pheno]["gid"]
                    except:
                        info = self.matrixHeaders[chro][pheno]["info"]
                    #strand = self.matrixHeaders[chro][pheno]["strd"]
                    region = info.split(";")[2].split(":")[1]
                    spl = region.split("-")
                    start = int(spl[0])
                    end = int(spl[1])
                    realCoordinates[ pheno] = (start,end)

        self.realCoordinates = realCoordinates
        self.phenoStart = phenoStart


    def get_cis_phenotypes(self, chro, start):
        """
        Get indexes of all phenotypes within the cis window.
        Only based on start position
        """
        query = (self.phenotypeInfo[chro]["start"] > start-self.windowSize) & (self.phenotypeInfo[chro]["start"] < start+self.windowSize)                                
        indexes = list(query[query == True].index)
                
        return indexes


    def cis_window_rank_approach(self):
        """

        Calculates correlation between phenotypes in cis. Performs permutations to calculate an empirical p-value.
        Ranks correlations found in cis window.
        
        Two approaches:
        - Rank approach: for each rank, a null distribution is built based on permutations and ranking permutation results.
        Each phenotype correlation (rank) is compared to the null distribution of the same rank to withdraw an empirical p-value.
        - First rank approach (more conservative): a null distribution is built only for the first rank (best correlations per permutation).
        Each phenotype correlation is compared to this same null distribution to withdraw an empirical p-value.
        
        The physical distance between phenotypes is only used to set a cis window limit ('start' positions are used), 
        proximity between phenotypes does not affect identification, but is reported.
        
        Reports several output files.
        """

        # Store results for all phenotypes
        outResults = open(self.outputFolder + "/" + CODer.OUTPUT_COP_RAW_RESULTS, "w")
        outResults.write("#chr\tcentralStart\tcentralEnd\tcentralPhenotype\tcentralInfo\tcentralStrand\tcisStart\tcisEnd\tcisPheno\tcisInfo\tcisStrand\tpairID\tcorr\tcorrSign\trank\tdistance\tadjustedPval\n")

        # stores phenotypes not in CODs because they have no phenotype in cis window
        noCisPheno = set() 
        # counts all phenos
        phenoCount = 0
        # Identify CODs for each chromosome separately
        for chro in self.chromosomes:
            
            Logger.get_instance().info( "cis_window_rank_approach: Processing chromosome: %s" % (chro) )
            
            if self.verbosityLevel == "debug":
                # Create chromosome-specific folder
                if not os.path.isdir(self.outputFolder + "/chr" + str(chro)) :
                    os.mkdir(self.outputFolder + "/chr" + str(chro))

            # Loop correlation matrix for each phenotype in chromosome
            chrPhenoCount = 0
            for pheno,_ in self.matrixData[chro].iteritems():
                
                phenoCount+=1
                chrPhenoCount+=1
                if chrPhenoCount % 500 == 0: 
                    Timer.get_instance().step( "cis_window_rank_approach: chr%s: phenotypes processed: %s" % (chro, chrPhenoCount) )  

                if self.verbosityLevel == "debug":                    
                    # Write a file for each phenotype
                    outRand = open("%s/chr%s/%s%s" % (self.outputFolder, chro, CODer.OUTPUT_COP_RANDOM_RESULTS, pheno), "w")
                    outRand.write("#chr\tstart\tend\tpheno_in_cis\tpheno_in_cis_info\tstrand\tcorrelation\tcorrSign\trank\tdistanceAbs\tdataset\tadjusted_pval\n")
                                
                # Get central phenotype start site (note that the start is like TSS in input file, relative to strand)                                
                chroInfo,start,_,phenoID,phenoInfo,strand = list(self.matrixHeaders[chro].loc[:,pheno])
                assert (pheno == phenoID), "confirm that corrMatrix and phenotypeInfo phenotypes match."
                assert (chro == chroInfo), "confirm that corrMatrix and phenotypeInfo phenotypes match."
                
                ######################
                # Get all phenotypes within the window distance              
                ######################

                # get indexes of phenotypes within the window range
                if self.determineTSS:
                    # if real coordinates are given in start and end columns, determine TSS using strand
                    indexes = self.get_cis_phenotypes(chro, self.phenoStart[pheno])
                else:
                    # if the start column already represent TSS of gene
                    indexes = self.get_cis_phenotypes(chro, start)

                # exclude central phenotype from the cis window
                if pheno in indexes:
                    indexes.remove(pheno)   
                else:
                    Logger.get_instance().warning( "cis_window_rank_approach: central gene not in indexes: %s. Indexes: %s" % (pheno, indexes))                    
                
                if len(indexes) == 0:
                    Logger.get_instance().warning( "cis_window_rank_approach: no phenotypes in cis window %s:%s for %s" % 
                                                   (start-self.windowSize, start+self.windowSize, pheno))                    
                    noCisPheno.add( pheno)
                    continue
                    
                # Get real phenotype matrix values for phenotypes in cis
                cisPhenoData = self.matrixData[chro][indexes].copy(deep=True)

                ######################
                # Get randomisation values
                ######################
                # RANDOMISATION APPROACH: for each randomisation, get shuffled matrix values for the central phenotype,
                #                  and correlate these with real values of phenotypes in cis                   
                # RANK APPROACH: Rank correlations by highest correlation value. 
                #                Real case against random null distribution is then compared for each rank separately.
                # FIRST RANK APPROACH: Rank correlations by highest correlation value. 
                #                      Real case against random null distribution of the first rank is then compared.
                                            
                randomCorrelations = {} # key -> rank, val -> list of random correlation values
                
                for perm in range(self.permutations):
                    rank = 0 # rank counter for phenotypes sorted by highest correlation
                    
                    if self.lowMem:
                        randVals = random.sample(list(self.matrixData[chro][pheno]), len(self.matrixData[chro][pheno]))
                    else:
                        randVals = self.randMatrixData[chro][perm][pheno]

                    # perform correlation of random central pheno values against real cis pheno values
                    randCorrelations = { cisPheno: self.perform_correlation(randVals, cisVals) for cisPheno, cisVals in cisPhenoData.items()}

                    # sort by highest correlation value
                    sortedRandCorrelations = sorted(randCorrelations.items(), key=lambda kv: kv[1], reverse = 1)
                                        
                    # Loop central phenotype correlations to each cis phenotype
                    for cisPhenoRandom,corrRandom in sortedRandCorrelations:
                         
                        # store correlation values of randomisations by rank, to be compared to real case ranks
                        if rank not in randomCorrelations: randomCorrelations[rank] = []
                        randomCorrelations[rank].append( corrRandom[0])

                        if self.verbosityLevel == "debug":
                            # write entry on file for each phenotype
                            outRand.write("%s\tNA\tNA\t%s\tNA\tNA\t%.5f\t%s\t%s\tNA\tperm%s\tNA\n" % (chro, cisPhenoRandom, corrRandom[0], corrRandom[1], rank, perm) )
                        
                        rank+=1
                    
                        if self.rankApproach == 0:
                            # only need to calculate first rank
                            break
 
                if self.rankApproach:
                    assert(len(randomCorrelations) == len(indexes)), "number of ranks should be the same as number of cis phenotypes"
                else:
                    assert(len(randomCorrelations) == 1), "only first rank should be present"

                ######################
                # Get real case and compare to random
                ######################

                # Get real coordinates of central phenotype (not only TSS)
                realStart, realEnd = self.realCoordinates[ pheno]

                # perform correlation of real central pheno values against real cis pheno values
                realCorrelations = {cisPheno: self.perform_correlation(self.matrixData[chro][pheno], cisVals) for cisPheno, cisVals in cisPhenoData.items() }

                # sort by highest correlation value
                sortedCorrelations = sorted(realCorrelations.items(), key=lambda kv: kv[1], reverse = 1)
 
                # iterate real correlations sorted by highest correlation
                rank = 0 # rank counter for phenotypes sorted by highest correlation
                for cisPhenoReal,corrReal in sortedCorrelations:

                    # get empirical pvalue 
                    if self.rankApproach:
                        # Rank approach: compare real correlation of a phenotype in this rank to random correlations of same rank
                        pval = Stats.empirical_pvalue(randomCorrelations[ rank], corrReal[0])[0]
                    else:
                        # First rank approach: compare real correlation against random correlations of first rank (the best correlations)
                        pval = Stats.empirical_pvalue(randomCorrelations[ 0], corrReal[0])[0]

                    # calculate absolute distance between central pheno and cis pheno
                    distance = abs( self.matrixHeaders[chro][cisPhenoReal]["start"] - start)

                    ######################
                    # Report results
                    ######################
                    # get info about cis phenotype
                    cisChro, _, _, _, cisInfo, cisStrand = self.matrixHeaders[chro][cisPhenoReal]
                    # get real coodinates of cis phenotype
                    cisStart,cisEnd = self.realCoordinates[ cisPhenoReal]

                    # Process ID
                    gene1ID = pheno
                    if gene1ID.startswith("ENS"):
                        # special case for Ensembl IDs. e.g. replace ENSG00000187634.11_5 with ENSG00000187634
                        if "." in gene1ID:
                            gene1ID = gene1ID.split(".")[0]
                    gene2ID = cisPhenoReal
                    if gene2ID.startswith("ENS"):
                        if "." in gene2ID:
                            gene2ID = gene2ID.split(".")[0]
                    # Sort IDs and make a pair ID that will be the same for gene1_gene2 and gene2_gene1
                    pairID = "|".join(sorted([gene1ID, gene2ID]))

                    # Write on file for each phenotype
                    resultText = "%s\t%i\t%i\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\t%.5f\t%s\t%i\t%i\t%.3g\n" % (chro, realStart, realEnd, pheno, phenoInfo, strand, 
                                                                                                           cisStart, cisEnd, cisPhenoReal, cisInfo, cisStrand, 
                                                                                                           pairID, corrReal[0], corrReal[1], rank, distance, pval )
                    outResults.write(resultText)

                    if self.verbosityLevel == "debug":
                        text = "%s\t%i\t%i\t%s\t%s\t%s\t%.5f\t%s\t%i\t%i\treal\t%.3g\n" % (cisChro, cisStart, cisEnd, cisPhenoReal, cisInfo, cisStrand, 
                                                                                                corrReal[0], corrReal[1], rank, distance, pval)
                        outRand.write(text)
                     
                    rank+=1
                     
                # add central pheno
                outResults.write("%s\t%i\t%i\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tcentral\tNA\tNA\n" % (chro, realStart, realEnd, pheno, phenoInfo, strand))
                if self.verbosityLevel == "debug":
                    outRand.write("%s\t%i\t%i\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tcentral\tNA\n" % (chro, realStart, realEnd, pheno, phenoInfo, strand ))
                    outRand.close()
                 
        t = "cis_window_rank_approach: total phenotypes processed: %s" % (phenoCount)
        Logger.get_instance().info(t)
        if self.writeLog: self.outLog.write(t + "\n\n")
        t = "cis_window_rank_approach: total phenotypes with no other phenotype in cis window: %s" % (len(noCisPheno))
        if self.writeLog: self.outLog.write(t + "\n\n")
        Logger.get_instance().info(t)
        
        outResults.close()


    def estimate_adjusted_pval_cutoff(self, data):
        '''
        Estimate adjusted p-value (empirical p-value) for which FDR is still below --fdrCutoff
        
        Note that this is only an approximation, since it takes the adjusted p-value in the real data,
        for which we are guaranteed to have an FDR below the --fdrCutoff.
        If no FDR is below the --fdrCutoff, the highest adjusted p-value is returned.
        
        @param data: pandas Dataframe with columns "adjustedPval" and "fdr"
        @return pvalCutoff : float - adjusted pvalue for which FDR is guaranteed to be below --fdrCutoff
        @return correspondingFDR : float - the actual FDR value respective of the chose cutoff adjusted pvalue
        '''

        assert (max(data["adjustedPval"]) <= 1.0 and min(data["adjustedPval"]) >= 0.0)
        assert (max(data["fdr"]) <= 1.0 and min(data["fdr"]) >= 0.0)

        # sort highest to lowest FDR
        data = data.sort_values(by=["fdr"], ascending=False)
        
        data.to_csv(open(self.outputFolder + "/" + CODer.OUTPUT_PVALUE_FDR, "w"), sep = "\t", index = False)

        # loop through FDR and select last adjusted pvalue where FDR was still above --fdrCutoff 
        pvalCutoff = min(data["adjustedPval"])
        correspondingFDR = None                
        for _, row in data.iterrows():
            if row["fdr"] < self.fdrCutoff:
                correspondingFDR = row["fdr"]
                pvalCutoff = row["adjustedPval"]
                break

        if correspondingFDR == None:
            Logger.get_instance().warning( "estimate_adjusted_pval_cutoff: could not estimate adjusted pvalue cutoff \
            (FDR never goes below --fdrCutoff %s)" % (self.fdrCutoff))
        else:
            t = "estimate_adjusted_pval_cutoff: adjusted pvalue cutoff: %.3g (for FDR = %3.g)" % (pvalCutoff, correspondingFDR)
            if self.writeLog: self.outLog.write(t + "\n\n")
            Logger.get_instance().info(t)

        return pvalCutoff, correspondingFDR


    def cop_identification(self):
        '''
        Apply multiple test correction and filter by FDR and/or minimum correlation.
        Creates output file with pairwise correlated phenotypes.
        '''

        data = pandas.read_csv(self.outputFolder + "/" + CODer.OUTPUT_COP_RAW_RESULTS, sep="\t", index_col=False)

        # Remove entries with central from file
        data = data[data["rank"] != "central"]
        data["rank"] = data["rank"].astype('int64')

        t = "cop_identification: processing %s phenotype-cis pairs " % (data.shape[0])
        if self.writeLog: self.outLog.write(t + "\n\n")
        Logger.get_instance().info(t)

        ######################
        # Apply multiple test correction over all phenotypes
        ######################

        if self.fdrCutoff != CODer.DEFAULT_PARAM_FDR_CUTOFF:
                
            # Filter by FDR value 
            if CODer.GLOBAL_PVALUE_CORRECTION:
                # Globally for the whole dataset
                if self.correctionProcedure:  data["fdr"] = self.q_value( data["adjustedPval"])
                else: data["fdr"] = self.bh_correct( data["adjustedPval"])

                # Only consider cases passing pvalue cutoff
                data = data[data["fdr"] < self.fdrCutoff]                
            else:
                # Use first rank to calculate FDR
                rank0 = data[data["rank"] == 0].copy(deep=True)
                
                # choice of correction method
                if self.correctionProcedure: rank0["fdr"] = self.q_value( rank0)
                else: rank0["fdr"] = self.bh_correct( rank0["adjustedPval"])
        
                pvals = rank0.copy(deep=True)
                
                # Estimate adjusted p-value (empirical p-value) for which FDR is below --fdrCutoff
                pvalCutoff, correspondingFDR = self.estimate_adjusted_pval_cutoff(pvals)

                # Only consider cases passing pvalue cutoff                
                data = data[data["adjustedPval"] <= pvalCutoff]
                                
                # Add FDR column
                data["fdr"] = "<%.1f%%" % (float(self.fdrCutoff)*100.0)

            Logger.get_instance().info( "cop_identification: %s phenotype-cis pairs after FDR filtering" % (data.shape[0]))

        ######################
        # Apply minimumCorrelation filter
        ######################
        if self.minimumCorrelation != CODer.DEFAULT_PARAM_MINIMUM_CORRELATION:
            # Only consider cases passing minimum correlation value
            data = data[data["corr"] > self.minimumCorrelation]
            t = "cop_identification: %s phenotype-cis pairs after minimum correlation filtering" % (data.shape[0])
            if self.writeLog: self.outLog.write(t + "\n\n")
            Logger.get_instance().info(t)

        ######################
        # Apply consecutive rank filter
        ######################
        # filter out cases in which first rank(s) do not pass filters but subsequent ranks do pass the filters

        if self.consecutiveRank:
            indexesToKeep = [] # stores all the rows that are fine
            
            # get first phenotype
            fPheno = str( data[:1]["centralPhenotype"])            
            expectedRank = 0
            consecutiveFlag = 1
            # the results are ordered by phenotype and then by rank
            for idx, row in data.iterrows():
                                
                pheno = str(row["centralPhenotype"])
                
                if pheno != fPheno:
                    # i.e. if this is a new phenotype
                    # reset flags
                    fPheno = pheno
                    expectedRank = 0
                    consecutiveFlag = 1

                # if the counter and observed rank do not match, it means that the observed rank jumped,
                # which means that the entry was filtered out with the previous filters
                rank = int(row["rank"])
                if expectedRank != rank:
                    consecutiveFlag = 0

                if consecutiveFlag:
                    indexesToKeep.append( idx)
                    
                expectedRank += 1

            t = "cop_identification: %s phenotype-cis pairs after consecutive rank filtering" % (len( indexesToKeep))
            if self.writeLog: self.outLog.write(t + "\n\n")
            Logger.get_instance().info(t)
            data = data.loc[indexesToKeep].copy(deep=True)

        ######################
        # Writing to output
        ######################
        
        data.to_csv(open(self.outputFolder + "/" + CODer.OUTPUT_COD_IDENTIFICATION_COPS, "w"), sep = "\t", index = False)

        # total number of COPs
        totalCOPs = len(set(data["centralPhenotype"]))
        totalPhenoInCOPs = len( set(data["cisPheno"]).union(set(data["centralPhenotype"])) )

        t = "cop_identification: total COPs found: %s" % (totalCOPs)
        if self.writeLog: self.outLog.write(t + "\n\n")
        Logger.get_instance().info(t)
        t = "cop_identification: total number of phenotypes in COPs: %s" % (totalPhenoInCOPs)
        if self.writeLog: self.outLog.write(t + "\n\n")
        Logger.get_instance().info(t)

        self.outputData = data

        # return only for testing purposes
        return (totalCOPs)
                          

    def cod_identification_network(self):
        '''
        Use igraph to turn pairwise central_pheno-cis_pheno into a network and get all connected components.

        Writes these components to two files: 
            - file with a line per phenotype (stating their component ID) 
            - file with a line per COD.
        '''
                
        ######################
        # Build igraph network
        ######################
        # Note: igraph needs node IDs to be numeric, not strings
        
        edges = [] # stores edges
        nodeNametoID = {} # key -> phenotype name, val -> phenotype ID in graph
        nodeIDtoName = {} # key -> phenotype ID in graph, val -> phenotype name
        idCounter = 0 # incremental number to attribute numerical ID to phenotypes
        for node1,node2 in self.outputData[["cisPheno","centralPhenotype"]].values:
            
            if node1 not in nodeNametoID:
                nodeNametoID[node1] = idCounter
                nodeIDtoName[idCounter] = node1
                idCounter += 1
     
            if node2 not in nodeNametoID:
                nodeNametoID[node2] = idCounter
                nodeIDtoName[idCounter] = node2
                idCounter += 1
     
            edges.append( (nodeNametoID[node1], nodeNametoID[node2]))
     
        # Create graph, only edges need to be given
        graph = igraph.Graph(edges)

        # Get the components (CODs) of the graph
        components = graph.components()

        t = "cod_identification_network: total edges in network: %s" % (len( edges))
        if self.writeLog: self.outLog.write(t + "\n\n")
        Logger.get_instance().info(t)
        t = "cod_identification_network: total nodes in network: %s" % ( idCounter) 
        if self.writeLog: self.outLog.write(t + "\n\n")
        Logger.get_instance().info(t)
        t = "cod_identification_network: total components (CODs) in network: %s" % ( len( components)) 
        if self.writeLog: self.outLog.write(t + "\n\n")
        Logger.get_instance().info(t)
        t = "cod_identification_network: clustering coefficient: %s" % (graph.transitivity_undirected()) 
        if self.writeLog: self.outLog.write(t + "\n\n")
        Logger.get_instance().info(t)

        ######################
        # Make CODs / report results
        ######################
        # output per phenotype
        outPerPheno = open(self.outputFolder + "/" + CODer.OUTPUT_COD_IDENTIFICATION_PER_PHENOTYPE, "w")
        outPerPheno.write("#chr\tstart\tend\tname\tinfo\tstrand\tcod_id\n")
        # output per COD
        outPerCOD = open(self.outputFolder + "/" + CODer.OUTPUT_COD_IDENTIFICATION_PER_COD, "w")
        outPerCOD.write("#chr\tstart\tend\tcod_id\tinfo\tstrand_ratio\tcod_size\tcod_size_genomic\tedgeMeanCorr\tcodMeanCorr\n")

        # for each COD
        componentID = 0
        for component in components:
            componentID+=1
                        
            strands = []
            coordinates = set() # to calculate COD boundaries
            phenos = [] # store phenotype names in COD
            correlations = [] # store correlations
            
            # for each phenotype in the COD
            for phenoIndex in component:

                try:
                    chro,_,_,phenoName,info,strand = self.globalMatrixHeaders[ nodeIDtoName[ phenoIndex]].values 
                    # get real coordinates instead of TSS
                    start,end = self.realCoordinates[ phenoName]
                except KeyError:
                    # do not provide this information in output file
                    chro,start,end,phenoName,info,strand = ("NA",0,0,nodeIDtoName[ phenoIndex],"NA","NA")
                
                strands.append(strand)
                coordinates.add(int(start) )
                coordinates.add(int(end) )
                phenos.append( phenoName)
                
                # Get all correlation values for this phenotype as central phenotype
                cors = self.outputData[self.outputData["centralPhenotype"] == phenoName][["corr"]].values
                for cor in cors:
                    correlations.append( cor)

                assert phenoName == nodeIDtoName[ phenoIndex]

                outPerPheno.write("%s\t%i\t%i\t%s\t%s\t%s\t%i\n" % (chro,start,end,phenoName,info,strand,componentID))

            # Calculate the mean correlation among all edges of a COD
            edgeMeanCorr = numpy.mean(correlations) 

            # Calculate the mean of the mean correlation among all phenotypes in the COD straight from the expression matrix
            try:
                codMeanCorr = "%.3f" % numpy.mean(numpy.mean( self.matrixData[chro][phenos].corr() )) 
            except:
                codMeanCorr = "NA"
            
            # Calculate ratio between positive and negative stranded phenotypes
            positiveStrand = len([x for x in strands if x == "+"])
            negativeStrand = len( strands) - positiveStrand
            if negativeStrand > 0:
                strandRatio = positiveStrand / float(negativeStrand)
            else:
                strandRatio = 1

            #Note: the COD coordinates are the minimum/maximum among all coordinates of its constituent phenotypes            
            outPerCOD.write("%s\t%i\t%i\t%s\t%s\t%s\t%i\t%i\t%.3f\t%s\n" % 
                            (chro, min(coordinates), max(coordinates), componentID, ",".join( phenos), strandRatio, len(component), max(coordinates)-min(coordinates), edgeMeanCorr, codMeanCorr))
            
        outPerPheno.close()
        
        return len(edges),idCounter,len(components)
        
        
    def build_distance_matched_null(self):
        """
        Create a file containing results for COPs and a set of (gene pair TSS) distance-matched non-COPs (null distribution).
        
        The self.nullMaxDistance controls the level of matching between COPs and non-COP (null) cases.
        Null cases within +/- self.nullMaxDistance are picked randomly without replacement.
        
        Note: only outputs COPs for which there is a matched null control, others are ignored.
        Note: null cases can be picked from any chromosome and output file will not be sorted.
        
        """
        
        ######################
        # Read all (positive and negative) results
        ######################

        rawData = pandas.read_csv(self.outputFolder + "/" + CODer.OUTPUT_COP_RAW_RESULTS, sep="\t", index_col=False)        
        # Remove entries with central from file
        if "rank" in rawData:
            rawData = rawData[rawData["rank"] != "central"]
            rawData["rank"] = rawData["rank"].astype('int64')
            rawData["distance"] = rawData["distance"].astype('int64')

        distanceDist = rawData["distance"].values

        # Store all information in a dictionary
        # Note: this is memory hungry, but faster than using pandas
        fullLines = {} # key -> pairID, val -> string with all info about this pair 
        # Note: the two directions are stored (e.g. gene1-gene2, gene2-gene1)
        with open(self.outputFolder + "/" + CODer.OUTPUT_COP_RAW_RESULTS, "r") as inFile:
            header = inFile.readline()
            for line in inFile:
                line = line.strip()
                spl = line.split()
                pairID = spl[11]
                if pairID == "NA":
                    continue
                if pairID not in fullLines:
                    fullLines[pairID] = []
                fullLines[pairID].append(line)

        ######################
        # Read positive cases and randomly pick distance-matched null
        ######################
        # Note: self.outputData should hold the set of "positive/significant" COPs
        # Note: assuming the two directions are being stored here (e.g. gene1-gene2, gene2-gene1

        allPositives = set(self.outputData["pairID"].values)
        
        outFile = open(self.outputFolder + "/" + CODer.OUTPUT_COP_DISTANCE_CONTROL_NULL, "w")
        outFile.write(header.strip() + "\tsignificant\tnullId\n")
                
        missingNullCount = 0
        processingCount = 0
        writeCount = 0
        alreadyPicked = set() # stores set of negative cases that have already been used in randomization
        alreadyProcessedPositive = {} # stores mapping between positive and negative cases
        for _, positive in self.outputData.iterrows():
            processingCount += 1
            if processingCount % 5000 == 0:
                Timer.get_instance().step( "process_null: positives processed: %s.." % (processingCount) )  

            # check if "inverse" comparison has already been processed
            if positive["pairID"] in alreadyProcessedPositive:
                # if this pair appeared before, use previous randomly picked pair and write inverse comparison
                pair = alreadyProcessedPositive[positive["pairID"]]

                # Write positive
                outFile.write("\t".join([str(x) for x in positive.values[:-1]]) + "\t1\t" + str(writeCount) + "\n")
                # Write negative (second comparison)
                outFile.write(fullLines[pair][1] + "\t0\t" + str(writeCount) + "\n")
                writeCount+=1
                
            else:
                dist = positive["distance"]

                # Get all cases within max distance (exclusive boundaries)
                # Note: this includes the positive case itself
                l1 = distanceDist[distanceDist > dist - self.nullMaxDistance]
                possibleDistances = set(l1[l1 < dist + self.nullMaxDistance])
    
                ## Pick one of the cases randomly, without replacement
                possiblePairs = set(rawData[rawData["distance"].isin(possibleDistances)]["pairID"])
                # Exclude nulls that were already used
                possiblePairs = possiblePairs - alreadyPicked
                # Exclude any pair identified as a COP
                possiblePairs = possiblePairs - allPositives
                
                if len(possiblePairs) > 0:
                    pair = numpy.random.choice(list(possiblePairs))
                    alreadyPicked.add(pair)
                    # Store randomly chosen pair so that it can be reused for inverse comparison
                    alreadyProcessedPositive[positive["pairID"]] = pair
                    
                    # Write positive
                    outFile.write("\t".join([str(x) for x in positive.values[:-1]]) + "\t1\t" + str(writeCount) + "\n")
                    # Write negative (first comparison)
                    outFile.write(fullLines[pair][0] + "\t0\t" + str(writeCount) + "\n")
                    writeCount+=1
                    
                else:                    
                    missingNullCount+=1
                    continue
        
        t = "build_distance_matched_null: non-redundant COPs without null (and thus excluded): %s " % (missingNullCount)
        if self.writeLog: self.outLog.write(t + "\n\n")
        Logger.get_instance().info(t)

        t = "build_distance_matched_null: wrote %s total non-redundant COP-null cases" % (writeCount)
        if self.writeLog: self.outLog.write(t + "\n\n")
        Logger.get_instance().info(t)

        outFile.close()
       

    def build_distance_and_variable_matched_null(self):
        """
        Create a file containing results for COPs and a set of non-COP gene pairs matched for the chosen variable 
        in addition to a distance-matched non-COPs (null distribution).
        
        The self.variableMatchProp controls how spread (in proportion terms) is the value of the variable between a COP and a matching non-COP
        
        The self.nullMaxDistance is interpreted here as a proportion in the same way as self.variableMatchProp.
        
        Null cases are picked randomly without replacement.
        
        Note: only outputs COPs for which there is a matched null control, others are ignored.
        Note: null cases can be picked from any chromosome and output file will not be sorted.
        
        """
        
        variable = self.variableMatch
        
        ######################
        # Read all (positive and negative) results
        ######################

        rawData = pandas.read_csv(self.outputFolder + "/" + CODer.OUTPUT_COP_RAW_RESULTS, sep="\t", index_col=False)        
        # Remove entries with central from file
        if "rank" in rawData:
            rawData = rawData[rawData["rank"] != "central"]
            rawData["rank"] = rawData["rank"].astype('int64')
            rawData["distance"] = rawData["distance"].astype('int64')

        distanceDist = rawData["distance"].values
        
        if variable in rawData:
            varDist = rawData[variable].values
               
            # Mapping of gene pair and variable
            varPair = {key : val for key, val in rawData[["pairID",variable]].values}
            
            # Store all information in a dictionary
            # Note: this is memory hungry, but faster than using pandas
            fullLines = {} # key -> pairID, val -> string with all info about this pair 
            # Note: the two directions are stored (e.g. gene1-gene2, gene2-gene1)
            with open(self.outputFolder + "/" + CODer.OUTPUT_COP_RAW_RESULTS, "r") as inFile:
                header = inFile.readline()
                for line in inFile:
                    line = line.strip()
                    spl = line.split()
                    pairID = spl[11]
                    if pairID == "NA":
                        continue
                    if pairID not in fullLines:
                        fullLines[pairID] = []
                    fullLines[pairID].append(line)
    
            ######################
            # Read positive cases and randomly pick distance-matched null
            ######################
            # Note: self.outputData should hold the set of "positive/significant" COPs
            # Note: assuming the two directions are being stored here (e.g. gene1-gene2, gene2-gene1
    
            allPositives = set(self.outputData["pairID"].values)
            
            outFile = open(self.outputFolder + "/" + CODer.OUTPUT_COP_DISTANCE_VARIABLE_CONTROL_NULL, "w")
            outFile.write(header.strip() + "\tsignificant\tnullId\n")
                    
            missingNullCount = 0
            processingCount = 0
            writeCount = 0
            alreadyPicked = set() # stores set of negative cases that have already been used in randomization
            alreadyProcessedPositive = {} # stores mapping between positive and negative cases
            for _, positive in self.outputData.iterrows():
                            
                processingCount += 1
                if processingCount % 5000 == 0:
                    Timer.get_instance().step( "process_null: positives processed: %s.." % (processingCount) )  
    
                # check if "inverse" comparison has already been processed
                if positive["pairID"] in alreadyProcessedPositive:
                    # if this pair appeared before, use previous randomly picked pair and write inverse comparison
                    pair = alreadyProcessedPositive[positive["pairID"]]
    
                    # Write positive
                    outFile.write("\t".join([str(x) for x in positive.values[:-1]]) + "\t1\t" + str(writeCount) + "\n")
                    # Write negative (second comparison)
                    outFile.write(fullLines[pair][1] + "\t0\t" + str(writeCount) + "\n")
                    writeCount+=1
                    
                else:
                    dist = positive["distance"]
                    var = varPair[positive["pairID"]]
    
                    # Get all cases within max distance proportion (exclusive boundaries)
                    # Note: this includes the positive case itself
                    l1 = distanceDist[distanceDist > dist - (dist * self.nullMaxDistance)]
                    possibleDistances = set(l1[l1 < dist + (dist * self.nullMaxDistance)])
                    possibleDistancePairs = set(rawData[rawData["distance"].isin(possibleDistances)]["pairID"])
                    
                    # Get all cases within max variable proportion (exclusive boundaries)
                    # Note: this includes the positive case itself
                    l2 = varDist[varDist > var - (var * self.variableMatchProp)]
                    possibleVars = set(l2[l2 < var + (var * self.variableMatchProp)])
                    possibleVarPairs = set(rawData[rawData[variable].isin(possibleVars)]["pairID"])
                    
                    # intersect the two filters
                    possiblePairs = possibleDistancePairs.intersection(possibleVarPairs)
                    
                    # Exclude nulls that were already used
                    possiblePairs = possiblePairs - alreadyPicked
                    # Exclude any pair identified as a COP
                    possiblePairs = possiblePairs - allPositives
    
                    ## Pick one of the cases randomly, without replacement
                    if len(possiblePairs) > 0:
                        pair = numpy.random.choice(list(possiblePairs))
                        alreadyPicked.add(pair)
                        # Store randomly chosen pair so that it can be reused for inverse comparison
                        alreadyProcessedPositive[positive["pairID"]] = pair
                        
                        # Write positive
                        outFile.write("\t".join([str(x) for x in positive.values[:-1]]) + "\t1\t" + str(writeCount) + "\n")
                        # Write negative (first comparison)
                        outFile.write(fullLines[pair][0] + "\t0\t" + str(writeCount) + "\n")
                        writeCount+=1
                        
                    else:  
                        missingNullCount+=1
                        continue
            
            t = "build_distance_and_variable_matched_null: non-redundant COPs without null (and thus excluded): %s " % (missingNullCount)
            if self.writeLog: self.outLog.write(t + "\n\n")
            Logger.get_instance().info(t)
    
            t = "build_distance_and_variable_matched_null: wrote %s total non-redundant COP-null cases" % (writeCount)
            if self.writeLog: self.outLog.write(t + "\n\n")
            Logger.get_instance().info(t)
    
            outFile.close()
        else:
            t = "build_distance_and_variable_matched_null: variable %s not present. Could not calculate matched null." % (variable)
            if self.writeLog: self.outLog.write(t + "\n\n")
            Logger.get_instance().info(t)

       
    def run(self):
        """
        Run functions in order
        """

        Timer.get_instance().step( "Reading matrix file.." )        
        self.read_matrix_file()

        Timer.get_instance().step( "Processing phenotype coordinates.." )        
        self.read_real_coordinates()

        if self.rerun:
            # rerun from CODer_raw_results.bed
            Timer.get_instance().step( "Skipping COD identification, loading provided file instead.." )        
        else:
            Timer.get_instance().step( "Identifying CODs using Rank Approach.." )        
            self.cis_window_rank_approach()

        Timer.get_instance().step( "Filtering and reporting COPs.." )        
        self.cop_identification()
            
        Timer.get_instance().step( "Network analysis and reporting final CODs.." )        
        self.cod_identification_network()

        if self.rerun and len(self.variableMatch) > 0:
            Timer.get_instance().step( "Building distance and variable-matched null distribution.." )        
            self.build_distance_and_variable_matched_null()
        else:
            Timer.get_instance().step( "Building distance-matched null distribution.." )        
            self.build_distance_matched_null()

        if self.writeLog: self.outLog.close()


if __name__ == "__main__":

    try:
    
        # Start chrono
        print ("STARTING " + SCRIPT_NAME)

        #===============================================================================
        # Get input arguments
        #===============================================================================
        parser = argparse.ArgumentParser(description= DESC_COMMENT) 
    
        # positional args
        parser.add_argument('matrixFile', metavar='matrixFile', type=str,
                             help='Input phenotype matrix in BED format to be processed. Individuals in columns, phenotypes in rows, header should be present. \
                             BED columns expected: "#chr    start   end     id      gid     strd".')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Folder where all output files will be written.')
        parser.add_argument('permutations', metavar='permutations', type=int,
                             help='Number of permutations of expression values per individual used to identify CODs.')
        parser.add_argument('--rankApproach', metavar='rankApproach', type=int, default = CODer.DEFAULT_PARAM_RANK_APPROACH, 
                             help='Choice of the "rank approach" or "first rank approach" in which phenotypes are sorted by their correlation value.\
                             If 0, a single null distribution for the first rank (best correlations) is made. If 1, each rank has its own null distribution (default = 0)')
        parser.add_argument('--windowSize', metavar='windowSize', type=int, default = CODer.DEFAULT_PARAM_WINDOW_SIZE, 
                             help='Window size (bp) in which to find correlation, upstream and downstream of phenotype. (default = 1000000 (1Mb))')
        parser.add_argument('--fdrCutoff', metavar='fdrCutoff', type=float, default = CODer.DEFAULT_PARAM_FDR_CUTOFF, 
                             help='If different than default, apply multiple test correction for all phenotypes.\
                              Results with FDR above given value will be filtered out. Use values from 0 to 1 (not percentage). (default = -1 / "OFF")')
        parser.add_argument('--minimumCorrelation', metavar='minimumCorrelation', type=float, default = CODer.DEFAULT_PARAM_MINIMUM_CORRELATION, 
                             help='Minimum correlation level between phenotypes to be considered into a COD. (default = OFF)')
        parser.add_argument('--negativeCorrelation', metavar='negativeCorrelation', type=int, default = CODer.DEFAULT_PARAM_NEGATIVE_CORRELATION, 
                             help='If 0, use normal correlation values (positive/negative). If 1, square the correlation value. \
                             If 2, use absolute correlation value (default = 2)')
        parser.add_argument('--consecutiveRank', metavar='consecutiveRank', type=int, default = CODer.DEFAULT_PARAM_CONSECUTIVE_RANK, 
                             help='Whether to consider sequentiality of ranks when identifying CODs (if one rank fails, all subsequent fail too). Relevant for --rankApproach 1 only. (default = 1)')
        parser.add_argument('--correctionProcedure', metavar='correctionProcedure', type=int, default = CODer.DEFAULT_PARAM_CORRECTION_PROCEDURE, 
                             help='This is only used if --fdrCutoff is >0. If 0, use BH correction, if 1, use qvalue correction  (default = 0)')
        parser.add_argument('--nullMaxDistance', metavar='nullMaxDistance', type=float, default = CODer.DEFAULT_PARAM_NULL_MAX_DISTANCE, 
                             help='Regarding distance-matched null control. Determines the maximum distance for which to look for a null case. \
                             This value may need to be increased in number of possible null cases is low. (default = 100 (bp))')
        parser.add_argument('--determineTSS', metavar='determineTSS', type=int, default = CODer.DEFAULT_PARAM_DETERMINE_TSS, 
                             help='If TSS is not given as start coordinate, determine it using start, end and strand. \
                             Use this if there is no "info" column with coordinates. (default = 0)')
        parser.add_argument('--wantedSeed', metavar='wantedSeed', type=int, default = CODer.DEFAULT_PARAM_WANTED_SEED, 
                             help='Choice of a specific seed for the randomisation (to have reproducible results). (default = OFF)')
        parser.add_argument('--verbosityLevel', metavar='verbosityLevel', type=str, default = "info", 
                             choices = ["debug", "info", "warning", "error", "critical", "fatal"],
                             help='Level of verbosity. Choices: "debug", "info", "warning", "error", "critical", "fatal"')
        parser.add_argument('--writeLog', metavar='writeLog', type=str, default = CODer.DEFAULT_PARAM_WRITE_LOG,
                             help='Provide file name to write log to a file as well as standard output. (default = OFF)')
        parser.add_argument('--lowMem', metavar='lowMem', type=int, default = CODer.DEFAULT_PARAM_LOW_MEM,
                             help='If 1, reduce memory consumption at the expense of a longer runtime. (default = 0)')
        parser.add_argument('--rerun', metavar='rerun', type=int, default = CODer.DEFAULT_PARAM_RERUN,
                             help='Rerun CODer from CODer_raw_results.bed results onwards (possibly using new options). (default = 0)')
        parser.add_argument('--variableMatch', metavar='variableMatch', type=str, default = CODer.DEFAULT_PARAM_VARIBLE_MATCH,
                             help='This works with --rerun option only. If a variable is provided, in addition to distance, use another variable to match the null. \
                             This variable has to be a column in the CODer_raw_results.bed file. Use the --variableMatchProp option to set the match spread used for this variable. \
                             Moreover, --nullMaxDistance is interperted as a percentage instead of base-pairs. (default = OFF)')
        parser.add_argument('--variableMatchProp', metavar='variableMatchProp', type=float, default = CODer.DEFAULT_PARAM_VARIBLE_MATCH_PROP,
                             help='Option to set the match spread used for this variable. \
                             Moreover, --nullMaxDistance is interperted as a percentage instead of base-pairs. (default = 0.1)')
                
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        # Initialise class    
        run = CODer( args.matrixFile, args.outputFolder, args.permutations, args.rankApproach,
                     args.windowSize, args.fdrCutoff, args.minimumCorrelation, args.negativeCorrelation, 
                     args.consecutiveRank, args.correctionProcedure, args.nullMaxDistance, args.determineTSS, args.wantedSeed, 
                     args.verbosityLevel, args.writeLog, args.lowMem, args.rerun, args.variableMatch, args.variableMatchProp)

        run.run()

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except Exception as e:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + str(e))

