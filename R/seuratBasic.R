rankCells<-function(disease,path,annotate=TRUE,scenario="Malacards",checkdrug=TRUE,userlabel="label",usercelltype="celltype"){
  if (missing(disease)) cat("Argument disease is missing") else cat(paste("Argument disease =", disease));
  cat ("\n");
  if (missing(path)) cat("Argument path is missing") else cat(paste("Argument path =", path));
  cat("\n\n");

  print(paste("Starting Analysis for",disease,sep=" "))
  setwd(path)

  if(annotate==TRUE){
    userlabel<-"label"
    usercelltype<-"celltype"
  }


  ##########################################Run Integrated analysis#############################################################################
  dirs <- list.dirs()

  loaded.dataSO.list<-c()

  plot.list<-vector("list", length(dirs)-1)
  index<-1
  # if(disease=="MM"){
  #   dirs<-dirs[c(1,2,13,20,27,28,29)] #changed this for mm
  # }

  #start from 2 to avoid home directory
  for(i in 2:length(dirs)){
    data_dir<-dirs[i]

    filesindir<-list.files(data_dir)
    #check if dir contains seurat raw data using different formats (rds, H5, ect)
    foundfalse<-which(c("barcodes.tsv","genes.tsv","matrix.mtx") %in% filesindir==FALSE)

    foundH5<-grep("h5",filesindir)
    foundrds<-grep("rds",filesindir)
    foundmeta<-"meta.txt"%in% filesindir
    if(length(foundfalse)==0){
      # Load the dataset
      loaded.data <- Read10X(data.dir = data_dir)
      project_name<-gsub("[\\.\\/]","",data_dir,)
      print(project_name)
      # Initialize the Seurat object with the raw (non-normalized data).
      if(foundmeta){
        meta <- read.table(paste(project_name,"/meta.txt",sep = ""), header=T, sep="\t", as.is=T, row.names=1)
        loaded.dataSO <- CreateSeuratObject(counts = loaded.data, project = project_name, min.cells = 3, min.features = 200,meta.data=meta)
      }else{
        loaded.dataSO <- CreateSeuratObject(counts = loaded.data, project = project_name, min.cells = 3, min.features = 200)
      }
      #QC and selecting cells for further analysis
      # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
      loaded.dataSO[["percent.mt"]] <- PercentageFeatureSet(loaded.dataSO, pattern = "^MT-")
      # Visualize QC metrics as a violin plot
      # VlnPlot(loaded.dataSO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
      # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

      #plot1 <- FeatureScatter(loaded.dataSO, feature1 = "nCount_RNA", feature2 = "percent.mt")
      plot2 <- FeatureScatter(loaded.dataSO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

      plot.list[[index]]<-plot2
      index<-index+1

      #Filtering data
      #loaded.dataSO <- subset(loaded.dataSO, subset = nFeature_RNA > 100 & nFeature_RNA < 6500 & percent.mt < 75)

      #Normalizing the data. Normalized values are stored in UnWounded1[["RNA"]]@data. *
      loaded.dataSO <- NormalizeData(loaded.dataSO)#, normalization.method = "LogNormalize", scale.factor = 10000)
      #Identification of highly variable features (feature selection) *
      loaded.dataSO <- FindVariableFeatures(loaded.dataSO, selection.method = "vst", nfeatures = 2000)

      loaded.dataSO.list<-c(loaded.dataSO.list,loaded.dataSO)
    }else if(length(foundH5) !=0){
      loaded.data <- Read10X_h5(paste(data_dir,"/",filesindir,sep=""), use.names = TRUE, unique.features = TRUE)
      project_name<-gsub("[\\.\\/]","",data_dir,)
      print(project_name)
      # Initialize the Seurat object with the raw (non-normalized data).
      if(foundmeta){
        meta <- read.table(paste(project_name,"/meta.txt",sep = ""), header=T, sep="\t", as.is=T, row.names=1)
        loaded.dataSO <- CreateSeuratObject(counts = loaded.data, project = project_name, min.cells = 3, min.features = 200,meta.data=meta)
      }else{
        loaded.dataSO <- CreateSeuratObject(counts = loaded.data, project = project_name, min.cells = 3, min.features = 200)
      }
      #QC and selecting cells for further analysis
      # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
      loaded.dataSO[["percent.mt"]] <- PercentageFeatureSet(loaded.dataSO, pattern = "^MT-")
      # Visualize QC metrics as a violin plot
      # VlnPlot(loaded.dataSO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
      # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

      #plot1 <- FeatureScatter(loaded.dataSO, feature1 = "nCount_RNA", feature2 = "percent.mt")
      plot2 <- FeatureScatter(loaded.dataSO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

      plot.list[[index]]<-plot2
      index<-index+1

      #Filtering data
      #loaded.dataSO <- subset(loaded.dataSO, subset = nFeature_RNA > 100 & nFeature_RNA < 6500 & percent.mt < 75) #changed this

      #Normalizing the data. Normalized values are stored in UnWounded1[["RNA"]]@data. *
      loaded.dataSO <- NormalizeData(loaded.dataSO)#, normalization.method = "LogNormalize", scale.factor = 10000)
      #Identification of highly variable features (feature selection) *
      loaded.dataSO <- FindVariableFeatures(loaded.dataSO, selection.method = "vst", nfeatures = 2000)
      loaded.dataSO.list<-c(loaded.dataSO.list,loaded.dataSO)
    }else if(length(foundrds) !=0){

      loaded.dataSO<- readRDS(paste(data_dir,"/",filesindir[foundrds],sep=""))
      project_name<-gsub("[\\.\\/]","",data_dir,)
      print(project_name)
      # Initialize the Seurat object with the raw (non-normalized data).
      # if(foundmeta){
      #   meta <- read.table(paste(project_name,"/meta.txt",sep = ""), header=T, sep="\t", as.is=T, row.names=1)
      #   loaded.dataSO <- CreateSeuratObject(counts = loaded.data, project = project_name, min.cells = 3, min.features = 200,meta.data=meta)
      # }else{
      #   loaded.dataSO <- CreateSeuratObject(counts = loaded.data, project = project_name, min.cells = 3, min.features = 200)
      # }
      #QC and selecting cells for further analysis
      # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
      # loaded.dataSO[["percent.mt"]] <- PercentageFeatureSet(loaded.dataSO, pattern = "^MT-")
      # Visualize QC metrics as a violin plot
      # VlnPlot(loaded.dataSO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
      # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

      #plot1 <- FeatureScatter(loaded.dataSO, feature1 = "nCount_RNA", feature2 = "percent.mt")
      plot2 <- FeatureScatter(loaded.dataSO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

      plot.list[[index]]<-plot2
      index<-index+1

      #Filtering data
      #loaded.dataSO <- subset(loaded.dataSO, subset = nFeature_RNA > 100 & nFeature_RNA < 6500 & percent.mt < 75) #changed this

      DefaultAssay(loaded.dataSO)<-"RNA"

      #Normalizing the data. *
      #loaded.dataSO <- NormalizeData(loaded.dataSO)#, normalization.method = "LogNormalize", scale.factor = 10000)
      #Identification of highly variable features (feature selection) *
      loaded.dataSO <- FindVariableFeatures(loaded.dataSO, selection.method = "vst", nfeatures = 2000)
      loaded.dataSO.list<-c(loaded.dataSO.list,loaded.dataSO)
    }
  }

  #increase memory size to hold larger objects
  options(future.globals.maxSize = 8000 * 1024^2)
  if(index==2){
    loaded.dataSO.combined <-loaded.dataSO
    plot2
  }else{
    do.call("grid.arrange", c(plot.list, ncol = (index-1)))
    features <- SelectIntegrationFeatures(object.list = loaded.dataSO.list)
    loaded.dataSO.anchors <- FindIntegrationAnchors(object.list = loaded.dataSO.list, dims = 1:20)#,reference = 1,k.filter = NA,anchor.features = features,verbose = TRUE
    loaded.dataSO.combined <- IntegrateData(anchorset = loaded.dataSO.anchors, dims = 1:20)
    DefaultAssay(loaded.dataSO.combined) <- "integrated"
  }

  # Run the standard workflow for visualization and clustering
  loaded.dataSO.combined <- ScaleData(loaded.dataSO.combined, verbose = FALSE)
  loaded.dataSO.combined <- RunPCA(loaded.dataSO.combined, npcs = 30, verbose = FALSE)
  # t-SNE and Clustering
  loaded.dataSO.combined <- RunUMAP(loaded.dataSO.combined, reduction = "pca", dims = 1:20)
  #loaded.dataSO.combined <- RunTSNE(loaded.dataSO.combined,reduction = "pca",dims = 1:10)
  #Louvein method or improved Leiden
  loaded.dataSO.combined <- FindNeighbors(loaded.dataSO.combined, reduction = "pca", dims = 1:20)
  loaded.dataSO.combined <- FindClusters(loaded.dataSO.combined, resolution  = 0.3)#was at 0.3

  #Can view object
  #loaded.dataSO.combined[[]]

  if(annotate==FALSE){
    p1combined <- DimPlot(loaded.dataSO.combined, reduction = "umap", group.by = userlabel)
  }
  p2combined <- DimPlot(loaded.dataSO.combined, reduction = "umap", label = TRUE)
  if(annotate==FALSE){
    p1combined+p2combined
  }else{
    p2combined
  }

  if(annotate==FALSE){
    DimPlot(loaded.dataSO.combined, reduction = "umap", split.by = userlabel,group.by = usercelltype,raster=FALSE)
  }

  if(annotate==TRUE){
    # Identify conserved cell type markers
    # For performing differential expression after integration, we switch back to the original
    # data
    DefaultAssay(loaded.dataSO.combined) <- "RNA"

    loaded.dataSO.combined.markers <- FindAllMarkers(loaded.dataSO.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

    loaded.dataSO.combined.markerstop1<-loaded.dataSO.combined.markers %>%
      group_by(cluster) %>%
      slice_max(n = 1, order_by = avg_log2FC) #100 gives best results

    loaded.dataSO.combined.markerstop2<-loaded.dataSO.combined.markers %>%
      group_by(cluster) %>%
      slice_max(n = 2, order_by = avg_log2FC) #100 gives best results

    loaded.dataSO.combined.markerstop100<-loaded.dataSO.combined.markers %>%
      group_by(cluster) %>%
      slice_max(n = 100, order_by = avg_log2FC) #100  gives best results

    loaded.dataSO.combined.markerstop200<-loaded.dataSO.combined.markers %>%
      group_by(cluster) %>%
      slice_max(n = 200, order_by = avg_log2FC) #200 gives second best results

    loaded.dataSO.combined.markerstop6<-loaded.dataSO.combined.markers %>%
      group_by(cluster) %>%
      slice_max(n = 6, order_by = avg_log2FC) #100 gives best results

    genesAll<-as.data.frame(loaded.dataSO.combined.markerstop6$gene)
    unique(loaded.dataSO.combined$seurat_clusters)

    #Visualize top 4 markers
    FeaturePlot(loaded.dataSO.combined, features = loaded.dataSO.combined.markerstop1$gene[c(1,2,3,4)])

    #Visualize volcano plot top 1 markers
    VlnPlot(loaded.dataSO.combined, features = loaded.dataSO.combined.markerstop1$gene[1])

    #FeaturePlot(loaded.dataSO.combined, features = loaded.dataSO.combined.markerstop1$gene[1],label = T)& theme(legend.position = c(0.1,0.2))

    loaded.dataSO.combined <- ScaleData(loaded.dataSO.combined, verbose = FALSE)

    DoHeatmap(loaded.dataSO.combined, features = loaded.dataSO.combined.markerstop2$gene) + NoLegend()

    #Run Enrich R on top 100 Markers
    listEnrichrSites()
    setEnrichrSite("Enrichr") # Human genes
    websiteLive <- TRUE
    #dbs <- listEnrichrDbs()
    #if (is.null(dbs)) websiteLive <- FALSE
    #f (websiteLive) head(dbs)
    annotatedclusters<-c()
    listofresults<-c()
    dbstouse <- c("PanglaoDB_Augmented_2021","CellMarker_Augmented_2021","Tabula_Sapiens")
    for(m in as.numeric(as.character(unique(loaded.dataSO.combined.markerstop100$cluster)))){
      markers_to_search<-loaded.dataSO.combined.markerstop100$gene[loaded.dataSO.combined.markerstop100$cluster==m]
      if (websiteLive) {
        enriched <- enrichr(markers_to_search, dbstouse)
      }
      celltype<-c(enriched[[1]]$Term[1],enriched[[2]]$Term[1],enriched[[3]]$Term[1])
      print(paste("Cluster",m,celltype,sep = " "))
      annotatedclusters<-c(annotatedclusters,toString(celltype[3]))#changed to Tabula Sapiens changed to Panglao for MM
      names(enriched)<-paste(names(enriched),m,sep = "_")
      listofresults<-c(listofresults,enriched)

      if (websiteLive) plot(plotEnrich(enriched[[paste("PanglaoDB_Augmented_2021_",m,sep="")]],  numChar = 40, y = "Count", orderBy = "P.value",title = paste("PanglaoDB",m))+plotEnrich(enriched[[paste("CellMarker_Augmented_2021_",m,sep="")]], numChar = 40, y = "Count", orderBy = "P.value",title = paste("CellMarker",m))+plotEnrich(enriched[[paste("Tabula_Sapiens_",m,sep="")]], numChar = 40, y = "Count", orderBy = "P.value",title = paste("Tabula_Sapiens",m)))
    }

    #Rename indents
    new.cluster.ids <- make.names(annotatedclusters,unique = T)
    names(new.cluster.ids) <- levels(loaded.dataSO.combined)
    loaded.dataSO.combined <- RenameIdents(loaded.dataSO.combined, new.cluster.ids)

    DimPlot(loaded.dataSO.combined, reduction = "umap",label = TRUE)

    #Change the Labels
    loaded.dataSO.combined$celltype.label <- paste(Idents(loaded.dataSO.combined), gsub("[0-9]+","",loaded.dataSO.combined$orig.ident), sep = "_")
    loaded.dataSO.combined$label <- gsub("[0-9]+","",loaded.dataSO.combined$orig.ident)
    loaded.dataSO.combined$celltype <- Idents(loaded.dataSO.combined)

    DimPlot(loaded.dataSO.combined, reduction = "umap", split.by = "label",label = TRUE)

    markers.to.plot <- unique(loaded.dataSO.combined.markerstop2$gene)
    DotPlot(loaded.dataSO.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "label") +
      RotatedAxis()
    DotPlot(loaded.dataSO.combined, features = c(loaded.dataSO.combined.markerstop1$gene[1]), cols = c("blue", "red"), dot.scale = 8, split.by = "label") +
      RotatedAxis()
    VlnPlot(loaded.dataSO.combined, features = c(loaded.dataSO.combined.markerstop1$gene[1]),split.by = "label")

  }

  xx <- as.list(GOTERM)

  if(scenario=="Malacards"){
    #Decided to include all the pathways in the KEGG to increase chances of finding KEGG and MSIG pathways
    KEGG<-list.files(pattern = "PathwaysKEGG.txt")
    GOs<-list.files(pattern = "PathwaysGO.txt")
    React<-list.files(pattern = "PathwaysReactome.txt")
    Wiki<-list.files(pattern = "PathwaysWiki.txt")

    file.size(Wiki) == 0L
    keywordsWiki<-read.table(Wiki, header=F, sep="\t")
    keywordsWiki<-as.array(keywordsWiki[,1])
    keywordsOther<-read.table(KEGG, header=F, sep="\t")
    keywordsOther<-as.array(keywordsOther[,1])
    keywordsGO<-read.table(GOs, header=F, sep="\t")
    keywordsGO<-as.array(keywordsGO[,1])
    keywordsMSIG<-keywordsOther
    keywordsKEGG<-keywordsOther
    termsWiki<-keywordsWiki
    termsWiki<-unique(tolower(termsWiki))
    keywordsReact<-c("")
    tryCatch(
      #try to do this
      {
        keywordsReact<-read.table(React, header=F, sep="\t")
        keywordsReact<-as.array(keywordsReact[,1])
      },
      #if an error occurs, print the error
      error=function(e) {
        message('An Error Occurred... empty file')
        print(e)
      },
      #if a warning occurs, print the warning
      warning=function(w) {
        message('A Warning Occurred')
        print(w)
        return(NA)
      }
    )

  }else{
    if(disease=="LAM"){
      keywordsWiki<-c("MTOR","PI3K","MAPK","apoptosis","NF-k","TNF")#,"estrogen" #LAM
      keywordsKEGG<-c("MTOR","PI3K","MAPK","apoptosis","NF-k","TNF")
      keywordsGO<-c("MTOR","PI3K","MAPK","apoptosis","NF-k","TNF")#"JAK","STAT","Cytokine","Inflammation","Th17","Th1","IL-"
      keywordsMSIG<-c("MTOR","PI3K","MAPK","apoptosis","NF-k","TNF")
      keywordsReact<-c("MTOR signalling","PI3K","MAPK","apoptosis","NF-k","TNF")
    }else if(disease=="autism"){
      keywordsWiki<-c("autism","Thyroid hormone")
      keywordsKEGG<-c("Pathways of neurodegeneration","Prion","Alzheimer","Parkinson","ALS","axon")
      keywordsGO<-c("axon","synapse","neurotransmitter","neuron","PI3K","MAPK","apoptosis","NF-k","TNF","JAK","STAT","Cytokine","Inflammation","Th17","Th1","IL-")
      keywordsMSIG<-c("Oxidative","PI3K","MAPK","apoptosis","NF-k","TNF","JAK","STAT","Cytokine","Inflammation","Th17","Th1","IL-")
      keywordsReact<-c("Pathways of neurodegeneration","Prion","Alzheimer","Parkinson","ALS","axon")
    }else{
      keywordsWiki<-c("SARS-CoV-2","Infectious disease","Innate Immune","Alzheimer","Huntington","autism")
      keywordsKEGG<-c("Corona disease","viral","Prion","Alzheimer","Parkinson","ALS","axon")
      keywordsGO<-c("axon","synapse","neurotransmitter","neuron","PI3K","MAPK","apoptosis","NF-k","TNF","JAK","STAT","Cytokine","Inflammation","Th17","Th1","IL-")
      keywordsMSIG<-c("Oxidative","PI3K","MAPK","apoptosis","NF-k","TNF","JAK","STAT","Cytokine","Inflammation","Th17","Th1","IL-")
      keywordsReact<-c("SARS-CoV-2","Infectious disease","Innate Immune","Alzheimer","Huntington","autism")
    }

    termsWiki<-c()
    indexWiki<-1
    for(WI in keywordsWiki){
      print(WI)
      termsWikiFirstSearch<-findPathwayNamesByText(WI)

      if(indexWiki !=1){
        foundinsearch<-grepl(WI, termsWikiFirstSearch,ignore.case = T)#not sure here double check changed from previous runs
        termsWiki<-c(termsWiki,termsWikiFirstSearch[foundinsearch])
      }else{
        termsWiki<-c(termsWiki,termsWikiFirstSearch)
      }
      indexWiki<-indexWiki+1
    }
    termsWiki<-termsWiki[!is.na(termsWiki)]
    termsWiki<-unique(tolower(termsWiki))
  }


  if(scenario=="Hypothesis"){
    termsGO<-c()
    for(KI in 1:length(keywordsGO)){
      indexsGO<-grep(keywordsGO[KI], xx,ignore.case = T)
      if(length(indexsGO)!=0){
        for(GOI in 1:length(indexsGO)){
          print(Term(xx[[indexsGO[GOI]]]))
          print(Ontology(xx[[indexsGO[GOI]]]))
          print(GOID(xx[[indexsGO[GOI]]]))
          if(Ontology(xx[[indexsGO[GOI]]])=="BP"){
            termsGO<-c(termsGO,Term(xx[[indexsGO[GOI]]]))
          }
        }
      }
    }
    termsGO<-unique(termsGO)
  }else{
    termsGO<-keywordsGO
  }
  termsKEGG<-c()
  for(KI in 1:length(keywordsKEGG)){
    tryCatch(
      #try to do this
      {
        termsKEGG<-c(termsKEGG,as.character(keggFind("pathway", c(keywordsKEGG[KI])))[1])
        termsKEGG<-c(termsKEGG,as.character(keggFind("disease", c(keywordsKEGG[KI])))[1])
      },
      #if an error occurs, print the error
      error=function(e) {
        message('An Error Occurred... term was not found in KEGG')
        print(e)
      },
      #if a warning occurs, print the warning
      warning=function(w) {
        message('A Warning Occurred')
        print(w)
        return(NA)
      }
    )
  }

  termsKEGG<-termsKEGG[!is.na(termsKEGG)]
  termsKEGG<-unique(termsKEGG)

  misig<-getMsigdb( org = c("hs"), id = c("SYM"), version = getMsigdbVersions() )
  hallmark<-subsetCollection(misig, 'h')
  termsMSIG<-c()

  for(KI in 1:length(keywordsMSIG)){
    indexsMSIG<-grep(keywordsMSIG[KI], hallmark,ignore.case = T)
    if(length(indexsMSIG)!=0){
      for(MSIGI in 1:length(indexsMSIG)){
        print(setName(hallmark[[indexsMSIG[MSIGI]]]))
        MSIGPath<-setName(hallmark[[indexsMSIG[MSIGI]]])
        MSIGPath<-gsub("HALLMARK_","",MSIGPath)
        termsMSIG<-c(termsMSIG,MSIGPath)
      }
    }
  }
  termsMSIG<-unique(termsMSIG)

  termsReact<-c()
  for(KI in 1:length(keywordsReact)){
    bdd.search <- searchQuery(query = keywordsReact[KI], species = "human",types = "Pathway")
    bdd.search$results$entries[[1]]$name[1]
    id<-query(id = paste("R-HSA-",bdd.search$results$entries[[1]]$dbId[1],sep=""))
    print(id$displayName)
    termsReact<-c(termsReact,id$displayName)
  }

  if(checkdrug==TRUE){
    Drugs<-list.files(pattern = "DrugsSorted.txt")
    keywordsMOA<-read.table(Drugs, header=F, sep="\t")
    keywordsMOA<-as.array(keywordsMOA[,1])
    termsMOA<-keywordsMOA
  }else{
    drugInfo<-read.delim("../drug_repurposing_hub.txt")
    #LAM
    if(disease=="LAM"){
      keywordsMOA<-c("CDK inhibitor","MTOR inhibitor","MEK inhibitor")
    }else if(disease=="autism"){
      #ASD
      keywordsMOA<-c("dopamine receptor antagonist","choline","histamine","serotonin")
    }else{
      #COVID
      keywordsMOA<-c("dopamine receptor antagonist","choline","histamine","serotonin")
    }


    termsMOA<-c()
    for(KI in 1:length(keywordsMOA)){
      indexsMOA<-grep(keywordsMOA[KI], drugInfo$moa,ignore.case = T)
      if(length(indexsMOA)!=0){
        for(MMOA in indexsMOA){
          if(grepl(keywordsMOA[KI], drugInfo$moa[MMOA],ignore.case = T)){
            termsMOA<-c(termsMOA,drugInfo$moa[MMOA])
          }
        }
      }
    }
    termsMOA<-unique(termsMOA)
  }
  #Perform mapping and ranking steps
  listofCellRanks<-scanCells(loaded.dataSO.combined,"Cell",termsKEGG,termsGO,termsMSIG,termsWiki,termsReact,termsMOA,userlabel,usercelltype,12,12,12,checkdrug)


  Ranks<-c()
  EucValues<-c()
  for(li in c(1,4,5,6,7,8)){
    print(li)
    fromlist<-as.data.frame(listofCellRanks[[li]])
    fromlist[,2]<-as.numeric(fromlist[,2])
    fromlist <- arrange(fromlist, fromlist[,2])
    fromlist2 <- arrange(fromlist, fromlist[,1])
    Ranks<-cbind(Ranks,rank(fromlist2[,2]))
    EucValues<-cbind(EucValues,fromlist2[,2])
  }
  colnames(Ranks)<-c("DRUG","KEGG","GOBP","MSIG","WIKI","REACT")
  rownames(Ranks)<-fromlist2$cellID

  colnames(EucValues)<-c("DRUG","KEGG","GOBP","MSIG","WIKI","REACT")
  rownames(EucValues)<-fromlist2$cellID
  write.table(Ranks,paste("Ranks",scenario,"PlusWikiReactFinalPac.txt",sep=""),quote = F,row.names = T,sep = "\t")
  write.table(EucValues,paste("EucValues",scenario,"PlusWikiReactFinalPac.txt",sep=""),quote = F,row.names = T,sep = "\t")

  MeanRanks<-as.data.frame(meanranks(t(Ranks))$mean.ranks)
  colnames(MeanRanks)<-c("Mean Ranking")
  MeanRanks<-arrange(MeanRanks,MeanRanks$`Mean Ranking`)

  write.table(MeanRanks,paste("MeanRanks",scenario,"PlusWikiReactFinalPac.txt",sep=""),quote = F,row.names = T,sep = "\t")


  Allavelog2FC<-as.data.frame(listofCellRanks[[11]])
  write.table(Allavelog2FC,"Allavelog2FCPac.txt",quote = F,row.names = T,sep = "\t")

  TotalNumberDEGs<-as.data.frame(listofCellRanks[[12]])
  TotalNumberDEGs$Var1<-rownames(TotalNumberDEGs)
  TotalNumberDEGs<-as.data.frame(TotalNumberDEGs)
  TotalNumberDEGs <-TotalNumberDEGs[order(TotalNumberDEGs$Var1), ]
  write.table(TotalNumberDEGs,"TotalNumberDEGsPac.txt",quote = F,row.names = T,sep = "\t")

  foldchangeInterMat<-runCellChat(loaded.dataSO.combined,userlabel,usercelltype)
  write.table(foldchangeInterMat,"CellChatPac.txt",quote = F,row.names = F,sep = "\t")

  # # How many cells are in each cell type or condition?
  Nocellspercelltype<-as.data.frame(table(loaded.dataSO.combined[[usercelltype]][,1]))
  Nocellspercelltype<-arrange(Nocellspercelltype, Nocellspercelltype$Var1)

  TotalNumberDEGsNorm<-as.data.frame(TotalNumberDEGs$TotalNumberDEGs/Nocellspercelltype$Freq)
  rownames(TotalNumberDEGsNorm)<-Nocellspercelltype$Var1
  colnames(TotalNumberDEGsNorm)[1]<-c("TotalNumberDEGsNorm")
  write.table(TotalNumberDEGsNorm,"TotalNumberDEGsNormPac.txt",quote = F,row.names = T,sep = "\t")

  # What proportion of cells are in each cluster?
  pie(prop.table(table(loaded.dataSO.combined$celltype)))
  #prop.table(table(Idents(loaded.dataSO.combined))) #same result

  # How does cluster membership vary by condition?
  tableofcellconts<-table(loaded.dataSO.combined$celltype, loaded.dataSO.combined$label)

  # What proportion of cells are in each cluster by condition
  proportions<-as.data.frame(prop.table(table(loaded.dataSO.combined[[usercelltype]][,1], loaded.dataSO.combined[[userlabel]][,1]), margin = 2))
  #proportions1<-as.data.frame(prop.table(table(loaded.dataSO.combined$celltype, loaded.dataSO.combined$label), margin = 2))

  write.table(proportions,"proportionsFinalPac.txt",quote = F,row.names = T,sep = "\t")
}
