#' Tidy GWAS
#'
#' @param GWAS df. RBD GWAS dataframe.
#'
#' @return Tidied LBD GWAS dataframe with appropriately named columns.
#' 

tidy_GWAS <- function(GWAS){
  
  GWAS_tidy <-
    GWAS %>% 
    as_tibble() %>% 
    dplyr::mutate(GWAS = "RBD",
                  Al1 = stringr::str_to_upper(A1),
                  Al2 = stringr::str_to_upper(A2)) %>% 
    dplyr::select(GWAS, SNP = MarkerName, beta, se, p.value = `P-value`, Al1, Al2, maf = MAF)
  
  return(GWAS_tidy)
  
}
