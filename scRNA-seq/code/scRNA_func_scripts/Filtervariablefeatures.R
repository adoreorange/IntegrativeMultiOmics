Filtervariablefeatures <- function(object,
                                   pattern=NULL){
  #filter non-essential variable features of lineages.
  if(!pattern %in% c("T_cell", "B_cell", "NBT")){
    cat("pattern must be one of T_cell, B_cell and NBT.\n")
  }
  
  #bcr.gene <- read.csv("/data/khaoxian/project/metagene/TCR_gene.csv")
  #bcr.gene <- bcr.gene$Approved.symbol
  
  if(pattern=="T_cell"){
    bcr.gene <- bcr.gene[-grep("^TRA$|^TRAC$|^TRBC1$|^TRBC2$|^TRD$|^TRDC$|^TRG$|^TRGC1$|^TRGC2$", bcr.gene)]
    #retain TRA, TRAC, TRBC1, TRBC2, TRD, TRDC, TRG, TRGC1, TRGC2 
    grepl <- "^IG[HJKL]|^RNA|^MT-|^RPS|^RPL"
  }
  if(pattern=="B_cell"){
    #grepl <- "^IG[JKL]|^IGH[VDJ]|^RNA|^MT-|^RPS|^RPL"
    grepl <- "^Ig[jkl]|^Igh[vdj]|^Rna|^Mt-|^Rps|^Rpl"
    grepl <- "^Ig[jkl]|^Igh[vdj]"
    #grepl <- "^Ig[jkl]|^Igh[vdj]|^Rna|^mt-|^Rps|^Rpl|^Hsp|^Xist|^Gm|*Rik$"
    #grepl <- "^Ig[jkl]|^Igh[vdj]|^Rna|^mt-|^Rps|^Rpl|^Hsp|^Gm|*Rik$"
  }
  if(pattern=="NBT"){
    grepl <- "^IG[HJKL]|^RNA|^MT-|^RPS|^RPL"
    #grepl <- "^Ig[hjkl]|^Rna|^Mt-|^Rps|^Rpl"
  }
  
  hvg <- VariableFeatures(object)
  hvg <- hvg[-grep(grepl,hvg)]
  print(hvg)
  #hvg <- hvg[!hvg %in% bcr.gene]
  VariableFeatures(object) <- hvg
  return(object)
}
