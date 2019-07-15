setwd("/home/lauterbur/Documents/orthoCapture/evaluation")
## orthoCapture evaluation metrics
library(ggplot2)
library(dplyr)
library(ape)
library(phytools)
library(gridBase)
library(ggtextures)
library(grid)
library(stringr)
library(viridis)

# load data
cfam<-read.csv("Cfam/Cfam_metrics.csv")
    cfam$GeneID<-gsub("d","",cfam$GeneID)
mmus<-read.csv("Mmus/Mmus_metrics.csv")
    mmus$GeneID<-gsub("d","",mmus$GeneID)
ppan<-read.csv("Ppan/Ppan_metrics.csv")
    ppan$GeneID<-gsub("d","",ppan$GeneID)
rnor<-read.csv("Rnor/Rnor_metrics.csv")
    rnor$GeneID<-gsub("d","",rnor$GeneID)
ptro<-read.csv("Ptro/Ptro_metrics.csv")
    ptro$GeneID<-gsub("d","",ptro$GeneID)

# make dataframe of % recovered for each gene, each species

data<-data.frame(rbind(cfam,mmus,ppan,rnor,ptro),
                 Species=c(rep("Canis familiaris",length(cfam[,1])),rep("Mus musculus",length(mmus[,1])),
                           rep("Pan paniscus",length(ppan[,1])),rep("Rattus norvegicus",length(rnor[,1])),
                           rep("Pan troglodytes",length(ptro[,1]))))

recovered<-data.frame(GeneID=data$GeneID,Species=data$Species,
                      perc_recovered=
                          (100*data$Query.Aligned.Length)/data$Target.Length)
    recovered_only<-na.omit(recovered)

    recovered_only %>% filter(Species=="Canis familiaris") %>% summarize(mean(perc_recovered))
## boxplot of exon coverage by oC (of those correctly retrieved and also annotated) for each species
    cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    levels(recovered_only$Species)<-c("Canis familiaris","Mus musculus","Rattus norvegicus","Pan paniscus","Pan troglodytes")
    p<-ggplot(recovered_only, aes(Species,perc_recovered)) +
        geom_boxplot(aes(fill=Species,color=Species)) +
        theme_bw() +
        #scale_fill_grey(start=0.4,end=1) +
        #scale_color_grey(start = 0,end = 0.6) +
        scale_fill_viridis(discrete=TRUE,alpha=0.75) +
        scale_color_viridis(discrete=TRUE) +
        theme(axis.text.x = element_text(face = "italic",size=6),
              legend.text = element_text(face = "italic",size=6), 
              axis.title = element_text(size=10),
              legend.title = element_text(size=10)) +
        ylab("Percent of exon\nrecovered by orthoCapture") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
    p    
    ggsave(p,filename = "percent_recovered.png",dpi = 300)
    ggsave(p,filename = "percent_recovered.eps",dpi = 300, device=cairo_ps)
    
# exons recovered only by orthoCapture or annotation - none
# oC_only<-data %>% filter(Alignment..ID=="noAnn")
#     oC_only<-select(oC_only,GeneID,Species)
# ann_only<-data %>% filter(Alignment..ID=="noOC")
#     ann_only %>% group_by(Species) %>% count()
    ## consider it to have the correct exon if >50% ID match to something with that name in BLAST or to design sequence
#    oc_correct<-read.csv("unique_to_orthocap.csv",stringsAsFactors = FALSE)
    # oc_correct<-
    # #oc_correct<-select(oc_correct,Species,Gene,correct) %>% filter(correct=="yes")
    #     oc_correct$Species[which(oc_correct$Species=="Cfam")]<-"Canis familiaris"
    #     oc_correct$Species[which(oc_correct$Species=="Mmus")]<-"Mus musculus"
    #     oc_correct$Species[which(oc_correct$Species=="rnor")]<-"Rattus norvegicus"
    #     oc_correct$Species[which(oc_correct$Species=="Ppan")]<-"Pan paniscus"
    #     oc_correct$Species[which(oc_correct$Species=="Ptro")]<-"Pan troglodytes"
        
# recovered method
oC_vs_ann<-data.frame(GeneID=c(), Species=c(),method=c(),correct=c())
    oC_vs_ann<-rbind(cbind(recovered_only[,1:2],
                           method=rep("annotation",length(recovered_only[,1])),
                           correct=rep("yes",length(recovered_only[,1]))),
                  cbind(recovered_only[,1:2],
                        method=rep("orthoCapture",length(recovered_only[,1])),
                        correct=rep("yes",length(recovered_only[,1]))))
    levels(oC_vs_ann$correct)<-c(levels(oC_vs_ann$correct),"no")
    for (i in oC_vs_ann$GeneID){
        #print(i)
        #print(data$Alignment..ID[which(data$GeneID==i)])
        if (min(data$Alignment..ID[which(data$GeneID==i)]) < 80) {
            print(oC_vs_ann$correct[which(oC_vs_ann$GeneID==i & oC_vs_ann$method=="orthoCapture")])
            oC_vs_ann$correct[which(oC_vs_ann$GeneID==i & oC_vs_ann$method=="orthoCapture")]<-"no"
            print(oC_vs_ann$correct[which(oC_vs_ann$GeneID==i & oC_vs_ann$method=="orthoCapture")])
        }
    }
    # oC_vs_ann<-rbind(oC_vs_ann,cbind(GeneID=as.character(oc_correct$Gene),
    #                                  Species=as.character(oc_correct$Species),
    #                                  method=rep("orthoCapture",length(oc_correct$correct)),
    #                                  correct=oc_correct$correct))
    # oC_vs_ann<-rbind(oC_vs_ann,cbind(GeneID=as.character(ann_only$GeneID),
    #                                  Species=as.character(ann_only$Species),
    #                                  method=rep("annotation",length(ann_only$GeneID)),
    #                                  correct=rep("yes",length(ann_only$GeneID))))
         
        oC_vs_ann_count<-oC_vs_ann %>% group_by(method,correct,Species) %>% count(method)
        oC_vs_ann_count$method<-factor(oC_vs_ann_count$method,levels(oC_vs_ann_count$method)[c(2,1)])
        oC_vs_ann_count$correct<-factor(oC_vs_ann_count$correct,levels(oC_vs_ann_count$correct)[c(2,1)])
        oC_vs_ann_count$Species<-factor(oC_vs_ann_count$Species,levels(oC_vs_ann_count$Species)[c(1,2,5,3,4)])
        
        
        
## plot by number recovered
      # images<-c("https://chrissinerantzi.files.wordpress.com/2018/05/future_earth.jpg",
          "https://images.fabric.com/images/693/693/RFR-005.jpg",
          "https://i.stack.imgur.com/WLsci.png")

    levels(oC_vs_ann_count$Species)<-c("Canis familiaris","Mus musculus","Rattus norvegicus","Pan paniscus","Pan troglodytes")

    # q<-ggplot(data = oC_vs_ann_count, aes(x = method, y = n, image= interaction(correct,method))) +
       q<-ggplot(data = oC_vs_ann_count, aes(x = method, y = n, fill=interaction(correct,method))) +
       # geom_textured_bar(stat = "identity", width = 1, position = "stack") +
        geom_col(position = "stack") +
       # scale_fill_manual(name = "legend",
       # scale_image_manual(name="legend",
         scale_fill_viridis(name="legend", discrete = TRUE,
         #             values = images,
                      labels = c(" orthoCapture,\n incorrect", " orthoCapture,\n correct"," annotation\n ")
                      ) +
        facet_grid(. ~ Species, labeller = label_wrap_gen(width=10)) +
        ylab("Number of genes\nretrieved") +
       # labs(tag="Homo sapiens",face="italic") +
        theme_bw() +
        theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.y = element_text(size=15), axis.title.y = element_text(size=20), 
              legend.title = element_blank(),
              legend.text = element_text(size=15),
              strip.text = element_text(size=17, face="italic"),
              plot.tag.position = c(.9,.01)) +
        coord_cartesian(ylim=c(500,710))    
       q

# add tree
    tree<-read.newick("species-list.nwk")
        #tree<-drop.tip(tree,"Pan_troglodytes")
        tree$tip.label[which(tree$tip.label=="Canis_lupus")]<-"Canis_familiaris"
       # tree$tip.label[which(tree$tip.label=="Rattus_norvegicus")]<-"Rattus_rattus"
        tree$tip.label[which(tree$tip.label=="Homo_sapiens")]<-"Homo\nsapiens"
        
#png("oC_vs_ann_tree.png",width=1200,height=1000)
postscript("oC_vs_ann_tree.eps",width=12,height=10)
#pdf("oC_vs_ann_tree.pdf",width=12,height=10)        
    par(mfrow=c(2,1)) 
    plot.new()
    vps<-baseViewports()
        pushViewport(vps$figure)        
        vp1 <-plotViewport(c(0,0,0,10))
        print(q,vp = vp1)
    par(mar=c(1,10,0,10),new=FALSE)
    plot(tree,direction="upwards",edge.width=10,show.tip.label=FALSE)
    text(x=6.1,y=50,labels = "million years",srt=90,cex=2)
    #text(x=6.1,y=100,labels = "Homo sapiens",cex=2)
    mtext("Homo \nsapiens",side=3,cex=1.5,adj=1,font=3)
    # corners<-par("usr")
    # par(xpd=TRUE)
    # text(x=)
        #,cex=c(.01,.01,.01,.01,1),srt=-90,adj=0.5,align.tip.label=TRUE,label.offset=-10
    axisPhylo(side=4,cex.axis=1.5)
dev.off()

        
# make heat map
### ARGH lost the new code for this
# genes<-read.table("gene-list.csv")
#     colnames(genes)<-"GeneID"
#     genes<-data.frame(genes,annotation=rep(0,length(genes[,1])),orthoCapture=rep(0,length(genes[,1])))
#     
# cfam<-genes
#     cfam$GeneID[which(cfam$GeneID %in% (oC_vs_ann %>% filter(Species=="Canis familiaris", method=="annotation"))$GeneID)]
#     cfam$annotation[which(cfam$GeneID %in% (oC_vs_ann %>% filter(Species=="Canis familiaris", method=="annotation"))$GeneID)]<-"present"
#     cfam$orthoCapture[which(cfam$GeneID %in% (oC_vs_ann %>% filter(Species=="Canis familiaris", method=="orthoCapture", correct=="yes"))$GeneID)]<-"present"
#     cfam$orthoCapture[which(cfam$GeneID %in% (oC_vs_ann %>% filter(Species=="Canis familiaris", method=="orthoCapture", correct=="no"))$GeneID)]<-"present but\nincorrect"
# 
#     cfam_melt<-melt(cfam,id.vars="GeneID",measure.vars=c("annotation","orthoCapture"))
#         colnames(cfam_melt)<-c("GeneID","method","presence")
#         cfam_melt$presence[which(cfam_melt$presence==0)] <- "absent"
# cfam_plot<-ggplot(data=cfam_melt, aes(x=method, y=GeneID)) +
#     geom_tile(aes(fill = presence),color="grey") +
#     scale_fill_manual(values=c("white","black","grey")) + 
#     ggtitle("Canis familiaris") +
#     theme(axis.text.x=element_text(size=20), axis.title.x=element_blank(), axis.title.y = element_text(size=20),
#           legend.title=element_text(size=20),legend.text=element_text(size=17),
#           plot.title=element_text(size=25, face="italic",hjust=0.5))
# cfam_plot
# ggsave(cfam_plot,filename = "cfam_heatmap.png")
# 
# mmus<-genes
#     mmus$GeneID[which(mmus$GeneID %in% (oC_vs_ann %>% filter(Species=="Mus musculus", method=="annotation"))$GeneID)]
#     mmus$annotation[which(mmus$GeneID %in% (oC_vs_ann %>% filter(Species=="Mus musculus", method=="annotation"))$GeneID)]<-"present"
#     mmus$orthoCapture[which(mmus$GeneID %in% (oC_vs_ann %>% filter(Species=="Mus musculus", method=="orthoCapture", correct=="yes"))$GeneID)]<-"present"
#     mmus$orthoCapture[which(mmus$GeneID %in% (oC_vs_ann %>% filter(Species=="Mus musculus", method=="orthoCapture", correct=="no"))$GeneID)]<-"present but\nincorrect"
# 
#     mmus_melt<-melt(mmus,id.vars="GeneID",measure.vars=c("annotation","orthoCapture"))
#         colnames(mmus_melt)<-c("GeneID","method","presence")
#         mmus_melt$presence[which(mmus_melt$presence==0)] <- "absent"
# mmus_plot<-ggplot(data=mmus_melt, aes(x=method, y=GeneID)) +
#     geom_tile(aes(fill = presence),color="grey") +
#     scale_fill_manual(values=c("white","black","grey")) + 
#     ggtitle("Mus musculus") +
#     theme(axis.text.x=element_text(size=20), axis.title.x=element_blank(), axis.title.y = element_text(size=20),
#           legend.title=element_text(size=20),legend.text=element_text(size=17),
#           plot.title=element_text(size=25, face="italic",hjust=0.5))
# mmus_plot
# ggsave(mmus_plot,filename = "mmus_heatmap.png")
# 
# rnor<-genes
#     rnor$GeneID[which(rnor$GeneID %in% (oC_vs_ann %>% filter(Species=="Rattus norvegicus", method=="annotation"))$GeneID)]
#     rnor$annotation[which(rnor$GeneID %in% (oC_vs_ann %>% filter(Species=="Rattus norvegicus", method=="annotation"))$GeneID)]<-"present"
#     rnor$orthoCapture[which(rnor$GeneID %in% (oC_vs_ann %>% filter(Species=="Rattus norvegicus", method=="orthoCapture", correct=="yes"))$GeneID)]<-"present"
#     rnor$orthoCapture[which(rnor$GeneID %in% (oC_vs_ann %>% filter(Species=="Rattus norvegicus", method=="orthoCapture", correct=="no"))$GeneID)]<-"present but\nincorrect"
# 
#     rnor_melt<-melt(rnor,id.vars="GeneID",measure.vars=c("annotation","orthoCapture"))
#         colnames(rnor_melt)<-c("GeneID","method","presence")
#         rnor_melt$presence[which(rnor_melt$presence==0)] <- "absent"
# rnor_plot<-ggplot(data=rnor_melt, aes(x=method, y=GeneID)) +
#     geom_tile(aes(fill = presence),color="grey") +
#     scale_fill_manual(values=c("white","black","grey")) + 
#     ggtitle("Rattus norvegicus") +
#     theme(axis.text.x=element_text(size=20), axis.title.x=element_blank(), axis.title.y = element_text(size=20),
#           legend.title=element_text(size=20),legend.text=element_text(size=17),
#           plot.title=element_text(size=25, face="italic",hjust=0.5))
# rnor_plot
# ggsave(rnor_plot,filename = "rnor_heatmap.png")
# 
# cfam<-genes
#     cfam$GeneID[which(cfam$GeneID %in% (oC_vs_ann %>% filter(Species=="Canis familiaris", method=="annotation"))$GeneID)]
#     cfam$annotation[which(cfam$GeneID %in% (oC_vs_ann %>% filter(Species=="Canis familiaris", method=="annotation"))$GeneID)]<-"present"
#     cfam$orthoCapture[which(cfam$GeneID %in% (oC_vs_ann %>% filter(Species=="Canis familiaris", method=="orthoCapture", correct=="yes"))$GeneID)]<-"present"
#     cfam$orthoCapture[which(cfam$GeneID %in% (oC_vs_ann %>% filter(Species=="Canis familiaris", method=="orthoCapture", correct=="no"))$GeneID)]<-"present but\nincorrect"
# 
#     cfam_melt<-melt(cfam,id.vars="GeneID",measure.vars=c("annotation","orthoCapture"))
#         colnames(cfam_melt)<-c("GeneID","method","presence")
#         cfam_melt$presence[which(cfam_melt$presence==0)] <- "absent"
# cfam_plot<-ggplot(data=cfam_melt, aes(x=method, y=GeneID)) +
#     geom_tile(aes(fill = presence),color="grey") +
#     scale_fill_manual(values=c("white","black","grey")) + 
#     ggtitle("Canis familiaris") +
#     theme(axis.text.x=element_text(size=20), axis.title.x=element_blank(), axis.title.y = element_text(size=20),
#           legend.title=element_text(size=20),legend.text=element_text(size=17),
#           plot.title=element_text(size=25, face="italic",hjust=0.5))
# cfam_plot
# ggsave(cfam_plot,filename = "Cfam_heatmap.png")