#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html

nw <- 1

library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

base_path <- "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_outputs/CellChat/"

print("Subset data for microglia and inhibitory parvalbumin neurons in females")

harmonized_object <- readRDS(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")
harmonized_subset <- subset(harmonized_object, Cluster %in% c("Mic1", "InN1_PV", "InN9_PV") & Sex == "Female")

harmonized_subset$Cluster <- recode(harmonized_subset$Cluster, InN1_PV = "InN_PV", InN9_PV = "InN_PV")
harmonized_subset$Cluster <- as.character(harmonized_subset$Cluster)

print("Create CellChat objects for cases and controls")

harmonized_subset_control <- subset(harmonized_subset, Condition == "Control")
harmonized_subset_case <- subset(harmonized_subset, Condition == "Case")

print(table(harmonized_subset_control$Cluster))
print(table(harmonized_subset_case$Cluster))

cellchat_control <- createCellChat(object = harmonized_subset_control, assay = "RNA",  group.by = "Cluster")
cellchat_case <- createCellChat(object = harmonized_subset_case, assay = "RNA",  group.by = "Cluster")

CellChatDB <- CellChatDB.human

future::plan("multicore", workers = nw)

print("Process case and control CellChat objects")

pdf(file = paste0(base_path, "12_individual_plots.pdf"))
for(cc_obj_name in c("cellchat_control", "cellchat_case")) { 
	print(cc_obj_name)
	cc_obj <- get(cc_obj_name)
	cc_obj@DB <- CellChatDB
	cc_obj <- subsetData(cc_obj)
	print("Identify over expressed signaling genes and L-R pairs")
	cc_obj <- identifyOverExpressedGenes(cc_obj)
	cc_obj <- identifyOverExpressedInteractions(cc_obj)
	#Using the PPI gives a lot of significantly different interactions. Using without the PPI for now, as this seems more conservative. 
	#cc_obj <- projectData(cc_obj, PPI.human)
	print("Calculate the communication probability between cell-types at L-R pair and signaling pathway level")
	cc_obj <- computeCommunProb(cc_obj, nboot = 1000)
		#This needs to be uncommented us using PPI.
		#raw.use = FALSE)
	cc_obj <- filterCommunication(cc_obj, min.cells = 10)
	cc_obj <- computeCommunProbPathway(cc_obj)
	print("Build aggregated network and compute centrality, that is the importance of different senders and receivers")
	cc_obj <- aggregateNet(cc_obj)
	cc_obj <- netAnalysis_computeCentrality(cc_obj, slot.name = "netP")
	groupSize <- as.numeric(table(cc_obj@idents))
	par(mfrow = c(1,2), xpd=TRUE)
	netVisual_circle(cc_obj@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, 
		title.name = "Number of interactions", vertex.label.cex = 2, edge.label.cex = 2) %>% print
	netVisual_circle(cc_obj@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, 
		title.name = "Interaction weights/strength", vertex.label.cex = 2, edge.label.cex = 2) %>% print
	#pathways.show <- c("GAS", "L1CAM", "CX3C", "EPHB", "TGFb", "NRXN", "COLLAGEN", "TENASCIN", "JAM", "MPZ", "NECTIN", "NGL", "SPP1") 
	#for(pt in pathways.show) {
	#	if(pt %in% cc_obj@netP$pathways) {
	#		par(mfrow=c(1,1))
	#		netVisual_heatmap(cc_obj, signaling = pt, color.heatmap = "Reds") %>% print
	#	}
	#}
	#netAnalysis_contribution(cc_obj, signaling = cc_obj@netP$pathways) %>% print
	#netAnalysis_signalingRole_heatmap(cc_obj, pattern = "outgoing") %>% print
	#netAnalysis_signalingRole_heatmap(cc_obj, pattern = "incoming") %>% print
	(netVisual_bubble(cc_obj, sources.use = c("Mic1"), targets.use = c("InN_PV"), font.size = 20) + theme(axis.text.x = element_text(angle = 45))) %>% print
	(netVisual_bubble(cc_obj, sources.use = c("InN_PV"), targets.use = c("Mic1"), font.size = 20) + theme(axis.text.x = element_text(angle = 45))) %>% print
	assign(cc_obj_name, cc_obj)
	write.csv(subsetCommunication(cc_obj), paste0(base_path, "12_", cc_obj_name,"_df.csv"))	
}
dev.off()

print("Combine CellChat objects")

object.list <- list(control = cellchat_control, case = cellchat_case)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

print("Compute the similarity of the networks between cases and controls")
#Need python and UMAP for some of the network steps
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
i <- 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

rm(list = ls(pattern = "^harmonized"))

save.image(file = paste0(base_path, "12_CellChat.RData"))

pdf(file = paste0(base_path, "12_combined_plots.pdf"), width = 14, height = 8)
print("Summary comparison of interactions in cases and controls")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), size.text = 25)   
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight", size.text = 25)  
gg1 + gg2
write.csv(gg1$data, file = paste0(base_path, "12_source_data_interactions_number.csv"))
write.csv(gg2$data, file = paste0(base_path, "12_source_data_interactions_weights.csv"))
#gg1 <- netVisual_heatmap(cellchat)
#gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#gg1 + gg2
print("Plot the change in sginaling network for each cell-type on a per pathway basis between cases and controls") 
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mic1")  + 
	theme(axis.text.x = element_text(angle = 45, size = 16), text = element_text(size = 16))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "InN_PV") +
	 theme(axis.text.x = element_text(angle = 45, size = 16), text = element_text(size = 16))
gg1 + gg2
print("Plot the differential number of interactions or interaction strength between Mic-InNPV in a circle plot")  
gg1 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", label.edge = T,
	vertex.label.cex = 2, edge.label.cex = 2)  
gg2 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", label.edge = T,
	vertex.label.cex = 2, edge.label.cex = 2) 
gg1
gg2
#netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#rankSimilarity(cellchat, type = "functional")
#print("Statistically assess whether certain signaling pathways have more weight in cases versus controls")
print("Assess whether certain signaling pathways have more weight in cases versus controls")
res <- rankNet(cellchat, mode = "comparison", stacked = F, 
	#do.stat = TRUE, 
	tol = 0.3, return.data = TRUE) 
write.csv(res$signaling.contribution, file = paste0(base_path, "12_rankNet_stats.csv"))
gg1 <- res$gg.obj + theme(text = element_text(size = 22), axis.text = element_text(angle = 45, size = 22))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = T, 
	#do.stat = TRUE, 
	tol = 0.3)  + 
	theme(text = element_text(size = 22), axis.text = element_text(angle = 45, size = 22))
gg1 + gg2
write.csv(gg2$data, file = paste0(base_path, "12_source_data_rankNet.csv"))
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c("Mic1"), targets.use = c("InN_PV"), tol = 0.1)+
	title("Mic1 to InN_PV") +  theme(text = element_text(size = 22), axis.text = element_text(angle = 45, size = 22))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c("InN_PV"), targets.use = c("Mic1"), tol = 0.1)+
	title("InN_PV to Mic1") +  theme(text = element_text(size = 22), axis.text = element_text(angle = 45, size = 22))
gg1 + gg2
write.csv(gg1$data, file = paste0(base_path, "12_source_data_rankNet_Mic1InPV.csv"))
write.csv(gg2$data, file = paste0(base_path, "12_source_data_rankNEt_InPVMic1.csv"))
print("Plot L-R pairs different communication probability between cases and controls. This does not use differential expression.")
res <- netVisual_bubble(cellchat, sources.use = c("Mic1"), targets.use = c("InN_PV"), comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Case",
        angle.x = 45, remove.isolate = T, return.data = TRUE)
write.csv(res$communication, file = paste0(base_path, "12_netvisual_bubble_MicInPV_inc_stats.csv"))
gg1 <- res$gg.obj + theme(text = element_text(size = 22), axis.text.x = element_text(angle = 45, size = 22))
gg1 + plot_spacer() # + gg2
write.csv(gg1$data, file = paste0(base_path, "12_source_data_netvisual_bubble_Mic1InPV_inc.csv"))
res <- netVisual_bubble(cellchat, sources.use = c("InN_PV"), targets.use = c("Mic1"), comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Case",
        angle.x = 45, remove.isolate = T, return.data = TRUE)
gg1 <- res$gg.obj+theme(text = element_text(size = 22), axis.text.x = element_text(angle = 45, size = 22))
write.csv(res$communication, file = paste0(base_path, "12_netvisual_bubble_InPVMic_inc_stats.csv"))
write.csv(gg1$data, file = paste0(base_path, "12_source_data_netvisual_bubble_InPVMic1_inc.csv"))
res <- netVisual_bubble(cellchat, sources.use = c("InN_PV"), targets.use = c("Mic1"), comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Case",
        angle.x = 45, remove.isolate = T, return.data = TRUE)
write.csv(res$communication, file = paste0(base_path, "12_netvisual_bubble_InPVMic_dec_stats.csv"))
gg2 <- res$gg.obj + theme(text = element_text(size = 22), axis.text.x = element_text(angle = 45, size = 22))
write.csv(gg2$data, file = paste0(base_path, "12_source_data_netvisual_bubble_InPVMic1_dec.csv"))
gg1 + gg2
res <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Case",
        angle.x = 45, remove.isolate = T, return.data = TRUE)
write.csv(res$communication, file = paste0(base_path, "12_netvisual_bubble_all_inc_stats.csv"))
res$gg.obj
res <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Case",
        angle.x = 45, remove.isolate = T, return.data = TRUE)
write.csv(res$communication, file = paste0(base_path, "12_netvisual_bubble_all_dec_stats.csv"))
res$gg.obj
print("For all signaling pathways detected in cases or controls plot the strength for cases and controls separately (all, outgoing, incoming)")
netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd") +  
	 #theme_classic(base_size = 25) #ComplexHeatmap::draw()
netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd") +  
	 #theme_classic(base_size = 25) #ComplexHeatmap::draw()
netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6) + 
	 #theme_classic(base_size = 25) #ComplexHeatmap::draw()
netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6) + 
	 #theme_classic(base_size = 25) #ComplexHeatmap::draw()
netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, 
	color.heatmap = "GnBu") +  
	#theme_classic(base_size = 25)
	#ComplexHeatmap::draw()
#netVisual_bubble(cellchat,  sources.use = c("Mic1"), targets.use = c("InN_PV"), comparison = c(1, 2), angle.x = 45)
#netVisual_bubble(cellchat,  sources.use = c("InN_PV"), targets.use = c("Mic1"), comparison = c(1, 2), angle.x = 45)
dev.off()

pdf(file = paste0(base_path, "12_combined_plots_additional.pdf"), width = 14, height = 8)
#Apparently no values returned for the line below.
tryCatch({netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6,
       color.heatmap = "GnBu")}, error = function(e) {print(e)})
       #%>%
       #ComplexHeatmap::draw()
#No interactions detected with the below line
tryCatch({netVisual_bubble(cellchat, sources.use = c("Mic1"), targets.use = c("InN_PV"), comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Case",
        angle.x = 45, remove.isolate = T)}, error = function(e) {print(e)})
dev.off()

print(summary(warnings()))
sessionInfo()
