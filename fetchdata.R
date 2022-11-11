suppressWarnings({lapply(c("dplyr","Seurat","HGNChelper",'GEOquery','Matrix','gridExtra','tidyverse','ggplot2'), library, character.only = T)})

workdir="C:/Users/dell/Documents/Mega/LearnscRNAAnalysis/Retina"
setwd(workdir)
accession_no = "GSE188280"
# download rawfile to the directory
getGEOSuppFiles(accession_no)
# To open the tar file, access the folder with
setwd(paste0(workdir,"/",accession_no))
# Open the tar file 
untar(paste0(accession_no,"_RAW.tar"))


# sapply(file_list, gunzip)

# Select some samples to analysis
# samples=c("GSM5676873","GSM5676874","GSM5676875","GSM5676876","GSM5676877",
#           "GSM5676878","GSM5676879","GSM5676880","GSM5676881","GSM5676882",
#           "GSM5676883","GSM5676884")

samples=c("GSM5676874","GSM5676876","GSM5676878","GSM5676880")



## Create subfolder  for each sample 
paths='C:/Users/dell/Documents/Mega/LearnscRNAAnalysis/Retina/GSE188280'

#fileList=list.files(paths, pattern = '.gz')


for (j in 1:length(samples)){
  folder<-dir.create(paste0(paths,'/',samples[j]))
}
## Move each sample file to the folder 

library(filesstrings)
dir.list= list.dirs(paths,recursive = FALSE)


for (i in 1:length(samples)){
  file_list<-list.files(paths, pattern  = paste0(samples[i],'_'))
  file.move(file_list, dir.list[i], overwrite=TRUE)
  
}

## Rename files in directory to fit into 10X 
for (i in dir.list){
  x=list.files(i, pattern = '.gz')
  for (a in x){
    b=gsub(".*_", "", a)
    setwd(i)
    file.rename(a,b)
  }
  
}

setwd(paths)

### Merge all together to prepare seurat object 

# get data location
dirs <- list.dirs(path = paths, recursive = F, full.names = F)

for(x in dirs){
  name <- x
  
  cts <- Read10X(data.dir = x)
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts, min.cells = 3, min.genes =200,
                                  project = 'AMD'))
}


## Merge the dataset 
merge_obj = merge(GSM5676874, y=c(GSM5676876, GSM5676878, GSM5676880),
                  add.cell.ids=c("GSM5676874","GSM5676876","GSM5676878","GSM5676880"), 
                  project = 'AMD'
)


# Save seurat object in file 

save(merge_obj, 
     file=paste0(workdir,'/','merge_obj.RData'))


# merge_max = merge(GSM5676873, y=c(GSM5676874,GSM5676875,GSM5676876,GSM5676877,
#                                   GSM5676878,GSM5676879,GSM5676880,GSM5676881,GSM5676882,
#                                   GSM5676883,GSM5676884),
#                   add.cell.ids=c("GSM5676873","GSM5676874","GSM5676875","GSM5676876","GSM5676877",
#                                  "GSM5676878","GSM5676879","GSM5676880","GSM5676881","GSM5676882",
#                                  "GSM5676883","GSM5676884"), 
#                   project = 'AMD'
# )



save(merge_max, 
     file=paste0(workdir,'/','mergeMax_obj.RData'))

