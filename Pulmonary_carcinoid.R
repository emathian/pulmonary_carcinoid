######################## Pulmmonary carcinoid #################

# Figure 1 : 
# LNEN MOFA MAP 
# Attributes to include 
# S1 : Sample Overview -> And ML predictions -> and Survival 
# S4 : Mutation 
# Immune check point with VST Table
# retinoid and xenobiotic metabolism pathways : CYP and UGT
# Cancer revelant core genes
# IHC Data : LAMP3 et CD1 ? NO
# Homeobox genes : HNF1A and HNFA4
# Expresison of DLL3 and ASCL1 NOTCH1-3 -MKI67 - EIF1AX  Cf:  fig 5B , fig 5C and fig 6 (very important)

#######################
#     LIBRAIRIES      #
#######################
library(openxlsx)


#######################
# IMPORTATION OF DATA #
#######################
Sample_overvieqw <- read.xlsx("SupplementaryTables_R1_20190318.xlsx", sheet = 1, startRow = 41, colNames = TRUE,
          rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
          skipEmptyCols = TRUE, rows = NULL, cols = NULL, check.names = TRUE,
          namedRegion = NULL, na.strings = "NA")

Somatic_mutation  <- read.xlsx("SupplementaryTables_R1_20190318.xlsx", sheet = 4, startRow = 51, colNames = TRUE,
                              rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                              skipEmptyCols = TRUE, rows = NULL, cols = NULL, check.names = TRUE,
                              namedRegion = NULL, na.strings = "NA")


Data_vst_50 <- read.table("data/VST_nosex_50pc_TCACLCNECSCLC.txt",  sep = " ", dec="." , header = TRUE,   quote="")

Ref_gene <- read.table("data/VST_nosex_50pc_TCACLCNECSCLC_annot.txt",  sep = " ", dec="." , header =FALSE,   quote="")
