#install.packages("VennDiagram")
getwd()
setwd("./cop_sanmple/") 
dir.create('venn')
setwd('venn')

library(VennDiagram)               #引用包
outFile="intersectGenes.txt"       #输出交集基因文件
outPic="venn.pdf"                  #输出图片文件

files=dir()                        #获取目录下所有文件
files=grep("*.CSV",files,value = T)   #提取TXT结尾的文件
geneList=list()

#读取所有CSV文件中的基因信息，保存到GENELIST
for(i in 1:length(files)){
  inputFile=files[i]
  if(inputFile==outFile){next}
  rt=read.csv(inputFile,header=T,sep = ',')        #读取
  #rt2=subset(rt,avg_log2FC>0)
  rt2=subset(rt,avg_log2FC<0)
  geneNames=as.vector(rt2[,1])              #提取基因名
  uniqGene=unique(geneNames)               #基因取unique
  header=unlist(strsplit(inputFile,"\\.|\\-"))[1]
  geneList[[header[1]]]=uniqGene
  uniqLength=length(uniqGene)
  print(paste(header[1],uniqLength,sep=" "))
}

#绘制vennͼ
outPic="venn_up.pdf"   
outPic="venn_down.pdf"
library(VennDiagram)
venn.plot=venn.diagram(geneList,filename=NULL,fill=rainbow(length(geneList)) )
pdf(file=outPic, width=5, height=5)
grid.draw(venn.plot)
dev.off()

#保存交集基因
outFile="intersectGenes_up.txt" 
outFile="intersectGenes_down.txt" 
intersectGenes=Reduce(intersect,geneList)
write.table(file=outFile,intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)
