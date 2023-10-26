select = dplyr::select

clean_fun = function(){
  rm(list = setdiff(ls(), lsf.str()))
}
colorstudy = function(name){
  colors_sets = ''
  if(name == 'oldTissue'){
    colors_sets = c("Adipose" = "#ba7a8d", 
                    "Brain" = "#6876a4",
                    "Liver" = "#f0b64d",
                    "Muscle" = "#754817")
  }
  if(name == 'Tissue'){
    colors_sets = c("vWAT" = "#ba7a8d",  
                    "Brain" = "#6876a4",
                    "Liver" = "#f0b64d",
                    "SkM" = "#754817")
  }
  if(name == 'oldGroup'){
    colors_sets = c("CD" = "#bbc6ca",
                    "HSD" ="#e5aa06",
                    "HSD+JNK_D1" = "#9c9629", 
                    "HSD+JNK_D2" = "#6674b0")
  }
  if(name == 'newGroup'){
    colors_sets = c("Control" = "#bbc6ca",
                    "Sucrose" ="#e5aa06",
                    "Sucrose+JNK_D1" = "#9c9629", 
                    "Sucrose+JNK_D2" = "#6674b0")
  }
  return(colors_sets)
}

reorder_col = function(data,name){
  temp = ''
  if(name == 'Tissue'){
    temp = data %>% 
      mutate(Tissue = case_when(Tissue == 'Brain' ~ 'Brain',
                                Tissue == 'Adipose' ~ 'vWAT',
                                Tissue == 'Muscle' ~ 'SkM',
                                Tissue == 'Liver' ~ 'Liver')) %>% 
      mutate(Tissue = factor(Tissue, levels = c("Liver","vWAT",
                                                "SkM","Brain")))
  }
  if(name == 'Group'){
    temp = data %>% mutate(Group = case_when(Group == "CD" ~ "Control",
                                             Group == "HSD" ~ "Sucrose",
                                             Group == "HSD+JNK_D1" ~ "Sucrose+JNK_D1",
                                             Group == "HSD+JNK_D2" ~ "Sucrose+JNK_D2",
                                             TRUE ~ "Check")) %>% 
      mutate(Group = factor(Group, levels= c('Control',
                                             'Sucrose',
                                             'Sucrose+JNK_D1',
                                             'Sucrose+JNK_D2')))
  }
  return(temp)
}

plot_theme = function(name){
  p_theme = ''
  if(name == 'grid'){
    p_theme = theme(
      axis.text.y = element_text(color = 'black', size = 8),
      legend.title = element_text(color = 'black', size = 8),
      legend.text = element_text(color = 'black', size = 8),
      legend.key.height = unit(.4, 'cm'),
      legend.key.width = unit(.4, 'cm'),
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(color = 'black', size = 8))
  }
  if(name == 'nogrid'){
    p_theme = theme(
      axis.text.y = element_text(color = 'black', size = 8),
      legend.title = element_text(color = 'black', size = 8),
      legend.text = element_text(color = 'black', size = 8),
      legend.key.height = unit(.4, 'cm'),
      legend.key.width = unit(.4, 'cm'),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(color = 'black', size = 8))
  }
  return(p_theme)
}
gene2FC = function(tissue, genes){
  file1 = paste0('../03.Results/DESeq.',tissue,'.HSDvsCD.txt')
  DEG_res1 = read.table(file1, sep = '\t', header = T) %>%
    mutate(Group = 'Sucrose')
  
  file2 = paste0('../03.Results/DESeq.',tissue,'.HSDJNK_D1vsHSD.txt')
  DEG_res2 = read.table(file2, sep = '\t', header = T) %>%
    mutate(Group = 'Sucrose+JNK_D1')
  
  file3 = paste0('../03.Results/DESeq.',tissue,'.HSDJNK_D2vsHSD.txt')
  DEG_res3 = read.table(file3, sep = '\t', header = T) %>%
    mutate(Group = 'Sucrose+JNK_D2')
  
  gene2fcs = rbind(DEG_res1, DEG_res2,DEG_res3) %>% 
    #filter(gene_names %in% genes) %>% 
    mutate(Tissue = tissue)
  return(gene2fcs)
}
removeindex = function(name){
  index = c("Cancer","Infection","Disease","Glioma","Influenza",
            "Arrhythmogenic Right Ventricular","Amyotrophic lateral sclerosis",
            "Asthma","Leishmaniasis","Viral myocarditis","Toxoplasmosis",
            "Diabetic cardiomyopathy","syndrome","diabetic complications",
            "Cocaine addiction","Pertussis","Malaria","Bacterial invasion",
            "Hepatitis","Carcinoma","Leukemia","Dilated Cardiomyopathy",
            "Hypertrophic Cardiomyopathy Hcm","atherosclerosis","Depression",
            "Melanoma","Diabetes Mellitus","Allograft rejection",
            "Chemical carcinogenesis - reactive oxygen species",
            "Rheumatoid arthritis","Spinocerebellar ataxia","Tuberculosis",
            "Hypertrophic cardiomyopathy","Viral life cycle - HIV-1","addiction")
  
  index = index %>% paste(., collapse = '|')
  return(index)
}

sig_label = function(data){
  data %>%  mutate(p.signif = case_when(padj <= 0.0001 ~ "****",
                            padj < 0.001 ~ "***",
                            padj <= 0.01 ~ "**",
                            padj <= 0.05 ~ "*",
                            TRUE ~ ""))
}

plotGseaTable_modify <- function(pathways, stats, fgseaRes,
                                 gseaParam=1,
                                 colwidths=c(5, 3, 0.8, 1.2, 1.2),
                                 render=TRUE){
  rnk <- rank(-stats)
  ord <- order(rnk)
  
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  
  pathways <- lapply(pathways, function(p) {
    unname(as.vector(na.omit(match(p, names(statsAdj)))))
  })
  
  # fixes #40
  pathways <- pathways[sapply(pathways, length) > 0]
  
  ps <- lapply(names(pathways), function(pn) {
    p <- pathways[[pn]]
    annotation <- fgseaRes[match(pn, fgseaRes$pathway), ]
    list(
      textGrob(pn, just="right", x=unit(0.95, "npc"),
               gp = gpar(col = "black", fontsize = 8)),
      ggplot() +
        geom_segment(aes(x=p, xend=p,
                         y=0, yend=statsAdj[p]),
                     size=0.2) +
        scale_x_continuous(limits=c(0, length(statsAdj)),
                           expand=c(0, 0)) +
        scale_y_continuous(limits=c(-1, 1),
                           expand=c(0, 0)) +
        xlab(NULL) + ylab(NULL) +
        theme(panel.background = element_blank(),
              axis.line=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              panel.grid = element_blank(),
              axis.title=element_blank(),
              plot.margin = rep(unit(0,"null"),4),
              panel.spacing = rep(unit(0,"null"),4)
        ),
      textGrob(sprintf("%.2f", annotation$NES),
               gp = gpar(col = "black", fontsize = 8)),
      textGrob(sprintf("%.1e", annotation$pval),
               gp = gpar(col = "black", fontsize = 8)),
      textGrob(sprintf("%.1e", annotation$padj),
               gp = gpar(col = "black", fontsize = 8))
    )
  })
  
  rankPlot <-
    ggplot() +
    geom_blank() +
    scale_x_continuous(limits=c(0, length(statsAdj)),
                       expand=c(0, 0)) +
    scale_y_continuous(limits=c(-1, 1),
                       expand=c(0, 0)) +
    xlab(NULL) + ylab(NULL) +
    theme(panel.background = element_blank(),
          axis.line=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid = element_blank(),
          axis.title=element_blank(),
          plot.margin = unit(c(0,0,0.5,0), "npc"),
          panel.spacing = unit(c(0,0,0,0), "npc")
    )
  
  grobs <- c(
    lapply(c("Pathway", "Gene ranks", "NES", "pval", "padj"), textGrob),
    unlist(ps, recursive = FALSE),
    list(nullGrob(),
         rankPlot,
         nullGrob(),
         nullGrob(),
         nullGrob()))
  
  # not drawing column if corresponding colwidth is set to zero
  grobsToDraw <- rep(colwidths != 0, length(grobs)/length(colwidths))
  
  
  p <- arrangeGrob(grobs=grobs[grobsToDraw],
                   ncol=sum(colwidths != 0),
                   widths=colwidths[colwidths != 0])
  
  if (render) {
    grid.draw(p)
  } else {
    p
  }
}