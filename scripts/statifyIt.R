#'Post-merge function to quantify compound emission rates relative to internal or external standard(s).
#'
#'@description If using internal standard (IS), will use the 'standardify()' function with user-specified 
#'inputs. If using external standard (ES), will require an input "matrix" from which standard curves will 
#'be derived.
#'
#'@param input uafR 'theMerger()' output. 
#'

data_in = unknowns_standardized

statifyIt = function(data_in){
 data_clms = data_in[,-c(1,2,3)]
 meta_clms = data_in[,c(1,2,3)]
 
 
 stat_format = t(data_clms)
 colnames(stat_format) = meta_clms$Chemical
 STD_missing = is.na(matrix(rowSums(stat_format)))
 
 stat_ready = data.frame(stat_format[!STD_missing,])
 Treatment = rownames(stat_ready)
 Treatment_split = strsplit(Treatment, "[\\_ \\-]")
 
 Treatment_split_df = do.call(rbind.data.frame, Treatment_split)
 row.names(Treatment_split_df) = NULL
 
 Treatment_clms = paste0("Treatment", 1:ncol(Treatment_split_df))
 colnames(Treatment_split_df) = Treatment_clms
 
 # as.matrix(unlist(Treatment_split))
 # stat_ready = data.frame(cbind(Treatment_split_df, stat_ready))
 
 dat_pca = prcomp(stat_ready, center = T)
 pca_dat_x = data.frame(dat_pca$x)
 pca_dat_rotated = dat_pca$rotation
 
 row.names(pca_dat_x) = NULL
 row.names(pca_dat_rotated) = NULL
 
 dat_analyzed = as.data.frame(cbind(Treatment_split_df,
                                 pca_dat_x))
 
 dat_manova = manova(data = dat_analyzed, cbind(PC1,PC2,PC3) ~ Treatment2-1)
 manova_summary = summary(dat_manova)
 print(manova_summary)
 return(dat_manova)
 # dat_plot = ggbiplot(dat_pca, choices = 1:2, 
 #                     ellipse = T, 
 #                     obs.scale = T, 
 #                     # var.scale = 8, 
 #                     var.axes = T, 
 #                     # varname.size = 14, 
 #                     # varname.adjust = 1.2, 
 #                     groups = dat_analyzed$Treatment, 
 #                     alpha=0)
 #  # scale_color_manual(name="Species", 
 #  #                    values=c(1:10))+
 #  # scale_shape_manual(name="Treatment", 
 #  #                    values=c(16:17)) +
 #  # geom_point(aes(color=gc_PCA$PlantSampled, 
 #  #                shape=gc_PCA$Treatment),
 #  #            size=8)+
 #  # labs(title = "Interchem Bioassay GC-MS PCA") +
 #  # theme_classic(base_size = 80) +
 #  # theme(legend.position = 'right',
 #  #       axis.text = element_text(size = 70),
 #  #       axis.text.x = element_text(size = 75),
 #  #       axis.title.x = element_text(size = 80),
 #  #       axis.title.y = element_text(size = 80),
 #  #       axis.text.y = element_text(size = 75),
 #  #       legend.text = element_text(size = 75))
 
 # dat_pca_summary = summary(dat_pca)
}

# morrison_test_standardized = unknowns_standardized
# morrison_necromones_standardized = unknowns_standardized
# hansen_silphium_AMF_standardized = unknowns_standardized
hansen_kernza_standardized = unknowns_standardized

dat_standardized = t(morrison_test_standardized[,-c(1,2,3)])
no_standard = is.na(matrix(rowSums(dat_standardized)))
dat_standardized2 = dat_standardized[!no_standard,]

dat_standardized = t(morrison_necromones_standardized[,-c(1,2,3)])
no_standard = is.na(matrix(rowSums(dat_standardized)))
dat_standardized2 = dat_standardized[!no_standard,]

dat_standardized = t(hansen_silphium_AMF_standardized[,-c(1,2,3)])
no_standard = is.na(matrix(rowSums(dat_standardized)))
dat_standardized2 = dat_standardized[!no_standard,]

dat_standardized = t(hansen_kernza_standardized[,-c(1,2,3)])
no_standard = is.na(matrix(rowSums(dat_standardized)))
dat_standardized2 = dat_standardized[!no_standard,]

dat_pca = prcomp(dat_standardized2)
summary(dat_pca)
dat_meta = unknowns_standardized[,c(1,2,3)]

Treatments = colnames(unknowns_standardized[,-c(1,2,3)])

dat_analyzed = cbind(Treatments, dat_meta, dat_pca$x)

dat_manova = manova(data = dat_analyzed, cbind(PC1, PC2)~Treatments-1)
summary.aov(dat_manova)

dat_plot = ggbiplot(dat_analyzed, choices = 1:2, 
                    ellipse = T, 
                    obs.scale = T, 
                    var.scale = 8, 
                    var.axes = T, 
                    varname.size = 14, 
                    varname.adjust = 1.2, 
                    groups = dat_analyzed$Chemical, 
                    alpha=0)+
 # scale_color_manual(name="Species", 
 #                    values=c(1:10))+
 # scale_shape_manual(name="Treatment", 
 #                    values=c(16:17)) +
 # geom_point(aes(color=gc_PCA$PlantSampled, 
 #                shape=gc_PCA$Treatment),
 #            size=8)+
 # labs(title = "Interchem Bioassay GC-MS PCA") +
 theme_classic(base_size = 80) +
 theme(legend.position = 'right',
       axis.text = element_text(size = 70),
       axis.text.x = element_text(size = 75),
       axis.title.x = element_text(size = 80),
       axis.title.y = element_text(size = 80),
       axis.text.y = element_text(size = 75),
       legend.text = element_text(size = 75))
