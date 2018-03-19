###########################
# convert Silicodart Data formatted variant table into vcf
##########################
library(data.table)

#
# enter filename of dart csv
#
fileName = "Report-DW17-2996/Report_DW17-2996_SNP_2.csv"


# read in csv
SNP2 = read.csv(fileName, skip = 5) #skip info rows
numInfoCols = 21 # the number of columns before we get to the SNP matrix

# get the number of samples and snps.
numRows = length(SNP2$AlleleID)
numSNPs = numRows/2
numSamples = dim(SNP2)[2]- numInfoCols  

# create indexes
# these indexs are used for a number of things.
# we are using 2 line dart so odd  index = ref, alt index = alt, etc
index = seq(1, numRows, 2) #list of uneven numbers
indexOp = seq(2, numRows, 2) #list of even numbers

chrom = SNP2$Chrom_Wheat_ChineseSpring04[index] # chromosome of SNP
pos = SNP2$ChromPos_Wheat_ChineseSpring04[index]+SNP2$SnpPosition[index] +1 # position of SNP - count starts at 0 so add 1
id = as.character(SNP2$AlleleID[index]) #SNP ID

# get ref and alt allele from SNP column. SNP = 
snp = as.character(SNP2$SNP)
snp = snp[indexOp]
snp = unlist(strsplit(snp, ":"))
alleles = snp[indexOp]
ref = unlist(strsplit(alleles, ">"))[index]
alt = unlist(strsplit(alleles, ">"))[indexOp]

#set up other vestors needed for output to vcf.
qual = rep(".", numSNPs)
filter = rep(".", numSNPs)
info = rep(".", numSNPs)
format = rep("GT", numSNPs)

#fill genotype matrix and get sample names
genos = SNP2[(numInfoCols+1):dim(SNP2)[2]]
genos[genos == "-"] = NA

samples = colnames(genos)
samples = paste(samples, collapse = "\t")

gt = matrix(nrow= numSNPs, ncol=numSamples)

gtI = 1
for(i in seq(1, numSNPs*2, 2))
{
  gt[gtI,] = paste(genos[i,],genos[i+1,], sep = "/")
  gtI = gtI +1
}

# R changes 0 and 1 to 2 and 3 due to a character conversion.
# change these back to 0 and 1 and make all NA = .
gt[gt == "2/2"] = "0/0"
gt[gt == "2/3"] = "0/1"
gt[gt == "3/2"] = "1/0"
gt[gt == "3/3"] = "1/1"

gt[gt == "NA/NA"] = "./."
gt[gt == "NA/0"] = "./0"
gt[gt == "NA/1"] = "./1"
gt[gt == "0/NA"] = "0/."
gt[gt == "1/NA"] = "1/."


############
# build up data for vcf
############
vcfData = data.frame(chrom, pos, id, ref, alt, qual, filter, info, format, gt)

## filter out SNPs with no chromosome
vcfData = vcfData[!vcfData$chrom =="",]

##order by chromosome and marker pos
vcfData = setorder(vcfData, chrom, pos)

## check for and remove duplicated variants, based on chrom and pos (not id)
checkId = paste(vcfData$chrom, vcfData$pos, sep = "/")
dupIds = !duplicated(checkId)
vcfData = vcfData[dupIds,]

############
#Header
############
date = Sys.Date()

header = sprintf("##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##fileDate=20180805
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s",samples)

#############
# Save to a new vcf
############
write.table(header,"output.vcf", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(vcfData,"output.vcf", quote = FALSE, append=TRUE, col.names = FALSE, row.names = FALSE, sep="\t")  


