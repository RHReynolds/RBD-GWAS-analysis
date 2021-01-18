#' Plot colocalisations
#'
#' @param coloc_hits Character. Vector of coloc hits as ensembl ids.
#' @param eQTL_GWAS_to_plot Dataframe. Dataframe of p-values for SNPs
#'   overlapping between GWAS and eQTL. Columns to include: chr, chromosome;
#'   pos, position; gene, ensembl id; hgnc_symbol, gene symbol; Dataset, whether
#'   p value is derived from GWAS or eQTL, denoted as 'pvalue_gwas' or
#'   'pvalue_eqtl'; p.value; pos_mb, position divided by 1000000 to give Mb
#'   position; log_pval, -log10(p.value).
#' @param facet_labels Named vector. Named character vector that that maps
#'   original names (Dataset i.e. pvalue_gwas or pvalue_eqtl) to new names of
#'   the datasets.
#' @param figure_labels Character. Vector with numbering/lettering that should
#'   be used to label each figure.
#' @param theme_base_size Integer. Argument for plot theme; base font size,
#'   given in pts. Default is 10.
#'
#' @return Plot of colocalisations.
#'   

plot_coloc_hits <- function(coloc_hits, eQTL_GWAS_to_plot, facet_labels = NULL, figure_labels = NULL, theme_base_size = 10){
  
  plot_list <- vector(mode = "list", length = length(coloc_hits))
  
  # Setting facet labels
  if(is.null(facet_labels)){
    
    labels <- c(pvalue_gwas = "pvalue_gwas", pvalue_eqtl = "pvalue_eqtl")
    
  } else{
    
    labels <- facet_labels 
    
  }
  
  # Setting figure labels
  if(is.null(figure_labels)){
    
    figure_labels <- NULL
    print(str_c("No figure labels provided, so no figure labels assigned."))
    
  } else{
    
    figure_labels <- figure_labels 
    
  }
  
  # Figure loop
  for(i in 1:length(coloc_hits)){
    
    chr <- 
      eQTL_GWAS_to_plot %>% 
      dplyr::filter(gene %in% coloc_hits[i]) %>% 
      .[["chr"]] %>% 
      unique
    
    plot_list[[i]] <- 
      eQTL_GWAS_to_plot %>% 
      dplyr::filter(gene %in% coloc_hits[i]) %>% 
      ggplot(aes(x = pos_mb, y = log_pval)) +
      geom_point(size = 0.7, alpha = 0.3) +
      facet_wrap(vars(Dataset, hgnc_symbol), ncol = 1, scale = "free", labeller=labeller(Dataset = labels)) +
      labs(x = str_c("Chromosome ", chr, " position (Mb)"), y = "-log10(p-value)") +
      theme_bw(base_size = theme_base_size) +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank())
    
  }
  
  ggarrange(plotlist = plot_list, 
            # ncol = 2,
            labels = figure_labels, 
            align = "hv",
            common.legend = TRUE, 
            legend = "none")
  
}

#' Plot tissue specificity of human genes.
#'
#' @param specificity Dataframe. Derived from specificity values calculated for
#'   the GTEx dataset (GTEx_v8.Rds) from
#'   \url{https://github.com/RHReynolds/MarkerGenes}.
#' @param theme_base_size Integer. Argument for plot theme; base font size,
#'   given in pts. Default is 10.
#'
#' @return Plot of tissue specificity.
#' 

plot_gtex_specificity <- function(specificity, theme_base_size = 10){
  
  data_to_plot <- specificity %>% 
    dplyr::mutate(`Brain tissue?` = case_when(str_detect(Organ, "Brain") ~ "TRUE",
                                              TRUE ~ "FALSE") %>% 
                    fct_relevel(., levels = c("TRUE", "FALSE")),
                  Organ = Organ %>% 
                    str_replace(., "\\(.*", "") %>% 
                    str_replace(., "Brain - ", "") %>% 
                    stringr::str_wrap(., width = 30)) 
  
  
  plot <- data_to_plot %>% 
    ggplot(aes(x = MarkerGenes::reorder_within(x = Organ,
                                               by = specificity,
                                               within = Description,
                                               fun = median,
                                               desc = FALSE),
               y = specificity,
               fill = `Brain tissue?`)
    ) +
    geom_col() +
    MarkerGenes::scale_x_reordered() +
    facet_wrap(vars(Description), scales = "free_y", ncol = 3) +
    labs(x = "Tissue", y = expression("Low specificity" %->% "High specificity"), title = "") +
    scale_y_continuous(limits = c(0,1)) +
    scale_fill_manual(values = c("#00BFC4", "#888888")) +
    theme_bw(base_size = theme_base_size) + 
    theme(legend.position = "top") +
    coord_flip() 
  
  return(plot)
  
}

#' Plot cell-type specificity of human genes.
#'
#' @param specificity Dataframe. Output of \code{MarkerGenes::query_gene_ctd()}
#'   using the AIBS dataset (AIBS2018_MTG.rda) from \url{https://github.com/RHReynolds/MarkerGenes}.
#' @param theme_base_size Integer. Argument for plot theme; base font size,
#'   given in pts. Default is 10.
#'
#' @return Plot of cell-type specificity.
#' 

plot_aibs_specificity <- function(specificity, theme_base_size = 10){
  
  plot <- specificity %>% 
    ggplot(aes(x = MarkerGenes::reorder_within(x = CellType,
                                               by = Specificity,
                                               within = Gene,
                                               fun = median,
                                               desc = FALSE),
               y = Specificity)
    ) +
    geom_col(fill = "#888888") +
    MarkerGenes::scale_x_reordered() +
    facet_wrap(vars(Gene), scales = "free_y", nrow = 1) +
    labs(x = "Cell type", y = expression("Low specificity" %->% "High specificity"), title = "") +
    scale_y_continuous(limits = c(0,1)) +
    theme_bw(base_size = theme_base_size) +
    coord_flip()
  
  return(plot)
  
}
