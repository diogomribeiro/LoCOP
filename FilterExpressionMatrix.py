#!/usr/bin/env python3.6

import argparse
import pandas
import numpy
import gzip

from cod.util.log.Logger import Logger
from cod.util.time.Timer import Timer

#===============================================================================
DESC_COMMENT = "Simple script to filter phenotype matrix by gene type, samples with zeroes or value threshold, and genomic regions."
SCRIPT_NAME = "FilterExpressionMatrix.py"
#===============================================================================

"""
#===============================================================================
Created on 27 Feb 2019
@author: Diogo Ribeiro
Simple script to filter phenotype matrix by gene type, samples with zeroes or value threshold, and genomic regions.
#===============================================================================

#===============================================================================
# Processing notes:
# 1) Using QTLtools quan to quantify expression, only the latest "quan_rewrite" QTLtools version displays gene types, not older versions.
# 2) Default options of filtering matrix by values assume values are positive.      
#===============================================================================


"""

class FilterExpressionMatrix(object):
    
    ## Set default arguments
    DEFAULT_GENE_TYPES = ""
    #DEFAULT_REGIONS_TO_EXCLUDE = ""
    DEFAULT_MINIMUM_SAMPLES = -1
    DEFAULT_MINIMUM_VALUE = 0.0
    DEFAULT_ADD_GTF_INFO = ""
    
    ## Other constants
    # number of column from input matrix from which individuals/samples start (from QTLtools quan)
    FIRST_SAMPLE_COLUMN = 6
    
    def __init__(self, matrixFile, outputFile, geneTypes, minimumSamples, minimumValue, addGtfInfo, verbosityLevel):
        
        self.matrixFile = matrixFile
        self.outputFile = outputFile
        self.geneTypes = geneTypes
#         self.regionsToExclude = regionsToExclude
        self.minimumSamples = minimumSamples
        self.minimumValue = minimumValue
        self.addGtfInfo = addGtfInfo

        Logger.get_instance().set_level( verbosityLevel)


    def read_matrix_file(self):
        """
        Read input matrix file using pandas.
        Example format:
            #chr    start   end     gene    info    strand  sample1  sample2
            GL000192.1      495564  495565  ENSG00000277655.1_5     L=451;T=unprocessed_pseudogene;R=GL000192.1:493155-495565;N=AC245407.1  -       0       0
            GL000193.1      81322   81323   ENSG00000280081.3_5     L=2485;T=lincRNA;R=GL000193.1:49232-81323;N=LINC01667   -       0       0

        """
        
        if self.matrixFile.endswith(".gz"):
            fileHandler = gzip.open(self.matrixFile, "rt")
        else:
            fileHandler = open(self.matrixFile, "r")
                    
        self.matrix = pandas.read_csv(fileHandler, sep="\t")
        

    def filter_by_gene_types(self):
        """
        Filter out rows by the gene type. 
        
        Note that this only works on updated QTLtools quan (branch: quan_rewrite).
        The gene type information is stored on the "info" column.        
        """
        
        # Read wanted gene types
        wantedGeneTypes = self.geneTypes.replace(",","|")

        # copy matrix
        filteredMatrix = self.matrix.copy()        

        # filter matrix by gene type, i.e. info column containing at least one of the wanted gene types
        filteredMatrix = filteredMatrix[filteredMatrix["info"].str.contains(wantedGeneTypes)]
        
#         # reindex matrix
#         filteredMatrix = filteredMatrix.reset_index( drop = True)
        
        Logger.get_instance().info( "Matrix size after filtering gene type: %s." % (len(filteredMatrix)))

        self.matrix = filteredMatrix


    def filter_by_minimum_samples_with_value(self):
        """
        Filter rows (e.g. genes/loci) in which the proportion of samples with a minimum value (--minimumValue) is equal or below the threshold (--minimumSamples).
        In other words, keep only rows in which at least X samples have values that pass the given cutoff.
        This can be used to remove rows containing only zeroes by setting --minimumValue at 0.0 (default value).  
        This method assumes values are positive.      
        """

        # Get number of samples from first matrix row
        totalSamples = self.matrix.shape[1] - FilterExpressionMatrix.FIRST_SAMPLE_COLUMN
        
        Logger.get_instance().info( "Filtering out rows if less than %s%% samples have value below or equal %s." % (self.minimumSamples*100.0, self.minimumValue))
        
        # iterate each matrix row
        rowsToExclude = set()
        for index, row in self.matrix.iterrows():
             
            # get values from all sample columns
            values = row[FilterExpressionMatrix.FIRST_SAMPLE_COLUMN:]
                        
            # calculate proportion of samples above minimum value
            proportion = (numpy.sum(values > self.minimumValue)) / float(totalSamples)
 
            # filter out row if number of samples with minimum value is below threshold
            if proportion <= self.minimumSamples:
                rowsToExclude.add( index)
                
        self.matrix = self.matrix[~self.matrix.index.isin(rowsToExclude)]     
        Logger.get_instance().info( "Matrix size after filtering minimum samples with value: %s." % (len(self.matrix)))


    def add_gtf_info(self):
        """
        Read optional gencode GTF file and add information about the feature and its strand.
        Example format (TSV, no header):
            1       .  gene    11869   14362   .       +       .       gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; \
                gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; \
                transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";

        Example output:
            ('info' column) L=3726;T=lincRNA;R=chr1:89295-133723;N=AL627309.1
            ('strd' column) +
        """

        ###############
        # Read GFT file and store in memory
        ###############

        featureDict = {} # key -> gene ID, val -> info (formatted)
        strandDict = {} # key -> gene ID, val -> strand

        with open(self.addGtfInfo, "r") as inFile:
            for line in inFile:
                if line.startswith("#"):
                    continue
                line = line.strip()
                chro,_,typ,start,end,_,strand,_,info = line.split("\t")
                
                if typ != "gene":
                    # skip entries that are not relative to a gene
                    continue
                
                #length = int(end)-int(start) ## the "L=" term on that this should correspond to is the sum of exon length. Not relevant here.
                #assert length >= 0
                
                # info example: gene_id "ENSG00000204481.6"; transcript_id "ENSG00000204481.6"; gene_type "protein_coding";
                geneID = ""
                geneType = ""
                geneName = ""
                for inf in info.split(";"):
                    if "gene_id" in inf:
                        geneID = inf.split('"')[1]
                    if "gene_type" in inf:
                        geneType = inf.split('"')[1]
                    if "gene_name" in inf:
                        geneName = inf.split('"')[1]

                if geneID == "":
                    raise Exception("add_gtf_info: could not identify the gene ID in info: " % info)
                if geneType == "":
                    Logger.get_instance().warning("add_gtf_info: gene_type not present in GTF file info: %s" % info)
                if geneName == "":
                    Logger.get_instance().warning("add_gtf_info: gene_name not present in GTF file info: %s" % info)
                    
                # Compose string like: L=3726;T=lincRNA;R=chr1:89295-133723;N=AL627309.1
                infoText = "L=%s;T=%s;R=chr%s:%s-%s;N=%s" %  ("NA", geneType, chro, start, end, geneName)
                
                if geneID not in featureDict:
                    featureDict[geneID] = ""
                    strandDict[geneID] = ""
                else:
                    raise Exception("add_gtf_info: duplicated gene ID: " % geneID)

                featureDict[geneID] = infoText
                strandDict[geneID] = strand

        Logger.get_instance().info("add_gtf_info: read %s gene entries" % len(featureDict) )

        ###############
        # Read info to each gene
        ###############

        # try different gene ID column names
        if "gene_id" in self.matrix.columns:
            geneCol = "gene_id"
        elif "gid" in self.matrix.columns:
            geneCol = "gid"
        else:
            raise Exception("add_gtf_info: cannot find column with gene IDs %s" % self.matrix.columns)

        
        infos = {} # key -> idx, val -> formatted feature info
        strds = {} # key -> idx, val -> strand 
        for idx, row in self.matrix.iterrows():
            gid = row[geneCol]
            
            if gid in featureDict:
                infos[idx] = featureDict[gid]
                strds[idx] = strandDict[gid]
            else:
                Logger.get_instance().warning("add_gtf_info: geneID not present in GTF file: %s" % gid)
                infos[idx] = "NA"
                strds[idx] = "NA"

        # Add columns to matrix
        self.matrix["info"] = pandas.Series(infos)
        self.matrix["strand"] = pandas.Series(strds)
        
        #reorder columns
        columns = list(self.matrix.columns)
        columns.insert(4,"info")
        columns.insert(5,"strand")
        columns = columns[:-2]
        assert len(columns) == len(list(self.matrix.columns))
        self.matrix = self.matrix[columns]

        Logger.get_instance().info( "Matrix size with GTF info: %s." % (len(self.matrix)))
        

    def write_filtered_matrix(self):
        """
        Write matrix after all filtering steps.
        """    
        self.matrix.to_csv( self.outputFile, sep="\t", index = False)


    def run(self):
        """
        Run functions in order
        """

        Timer.get_instance().step( "Reading matrix file.." )        
        self.read_matrix_file()
        Logger.get_instance().info( "Original matrix size: %s." % (len(self.matrix)))

        if self.addGtfInfo != FilterExpressionMatrix.DEFAULT_ADD_GTF_INFO:
            Timer.get_instance().step( "Adding information from GTF file.." )        
            self.add_gtf_info()

        if self.geneTypes != FilterExpressionMatrix.DEFAULT_GENE_TYPES:
            Timer.get_instance().step( "Filtering by gene type.." )        
            self.filter_by_gene_types()

        if self.minimumSamples != FilterExpressionMatrix.DEFAULT_MINIMUM_SAMPLES:
            Timer.get_instance().step( "Filtering by minimum samples with value.." )        
            self.filter_by_minimum_samples_with_value()

        Timer.get_instance().step( "Write filtered table.." )                    
        self.write_filtered_matrix()

        Logger.get_instance().info( "Final matrix size: %s." % (len(self.matrix)))


if __name__ == "__main__":

    try:
    
        # Start chrono
        print ("STARTING " + SCRIPT_NAME)
        Timer.get_instance().start_chrono()

        #===============================================================================
        # Get input arguments
        #===============================================================================
        parser = argparse.ArgumentParser(description= DESC_COMMENT) 
    
        # positional args
        parser.add_argument('matrixFile', metavar='matrixFile', type=str,
                             help='Input phenotype matrix to be filtered. Can be a gzip (.gz) file.')
        parser.add_argument('outputFile', metavar='outputFile', type=str,
                             help='Filtered phenotype matrix.')
        parser.add_argument('--geneTypes', metavar='geneTypes', type=str, default=FilterExpressionMatrix.DEFAULT_GENE_TYPES,
                             help='Comma-separated list of gene types to be included (Default: "").')
        parser.add_argument('--minimumSamples', metavar='minimumSamples', type=float, default=FilterExpressionMatrix.DEFAULT_MINIMUM_SAMPLES,
                             help='Minimum proportion (0.0 to 1.0) of samples (columns) meeting the "--minimumValue" in order to NOT exclude row. \
                             Use 0.0 to exclude rows with only zeroes. (Default: -1 [OFF]).')
        parser.add_argument('--minimumValue', metavar='minimumValue', type=float, default=FilterExpressionMatrix.DEFAULT_MINIMUM_VALUE,
                             help='Minimum value (row value, e.g. RPKM) required in order to NOT exclude row (this works together with "--minimumSamples") (Default: 0.0).')
        parser.add_argument('--addGtfInfo', metavar='addGtfInfo', type=str, default=FilterExpressionMatrix.DEFAULT_ADD_GTF_INFO,
                             help='If a gencode GTF file is given, add the columns relative to info (e.g. L=3726;T=lincRNA;R=chr1:89295-133723;N=AL627309.1) and strand to output expression matrix. (Default: OFF)')
        parser.add_argument('--verbosityLevel', metavar='verbosityLevel', type=str, default = "debug", 
                             choices = ["debug", "info", "warning", "error", "critical", "fatal"],
                             help='Level of verbosity. Choices: "debug", "info", "warning", "error", "critical", "fatal"')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        # Initialise class    
        run = FilterExpressionMatrix( args.matrixFile, args.outputFile, args.geneTypes, args.minimumSamples, args.minimumValue, args.addGtfInfo, args.verbosityLevel)

        run.run()

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except Exception as e:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + str(e))


