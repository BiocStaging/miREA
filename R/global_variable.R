#' Supported enrichment methods
#'
#' Character vector of the nine node-based and edge-based enrichment methods
#' supported by [miREA()] and [get_all_input_data()].
#'
#' @format A character vector of length 9.
#' @return A character vector of length 9.
#' @examples
#' default.method_list
#' @export
default.method_list <- c("TG_ORA", "TG_Score", "MiR_ORA", "MiR_Score", "Edge_ORA", "Edge_Score1D", "Edge_Score2D", "Edge_Topology", "Edge_Network")

#' Supported score transformations
#'
#' Character vector of the score transformations supported by the
#' `scoreFun` argument of [get_all_input_data()].
#'
#' @format A character vector of length 4.
#' @return A character vector of length 4.
#' @examples
#' default.scoreFun
#' @export
default.scoreFun <- c("none", "rank", "sigmoid", "normalize")

#' Supported p-value adjustment methods
#'
#' Character vector of the multiple-testing correction methods supported by
#' the `pAdjMethod` argument of [miREA()] and the individual enrichment
#' method functions (see [stats::p.adjust()]).
#'
#' @format A character vector of length 8.
#' @return A character vector of length 8.
#' @examples
#' default.pAdjMethod
#' @export
default.pAdjMethod <- c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none")

#' Supported predefined pathway sets
#'
#' Character vector naming the pathway sets shipped with the miREA
#' manuscript's companion data deposit.
#'
#' @format A character vector of length 3.
#' @return A character vector of length 3.
#' @examples
#' default.pathway_list
#' @export
default.pathway_list <- c("Reactome", "hallmark", "TP")

#' Cancer types used for the TP pathway set
#'
#' Character vector of TCGA cancer-type abbreviations used to build the
#' `"TP"` (cancer-specific) pathway set in [default.pathway_list].
#'
#' @format A character vector of length 16.
#' @return A character vector of length 16.
#' @examples
#' default.TP_cancer
#' @export
default.TP_cancer <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "KIRC",
                      "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD",
                      "READ", "STAD", "THCA", "UCEC")

#' Supported GSEA score types
#'
#' Character vector of the score types supported by the `pvalueType`
#' argument of the scoring-based enrichment methods (`TG_Score`, `MiR_Score`,
#' `Edge_Score1D`, `Edge_Score2D`, `Edge_Topology`).
#'
#' @format A character vector of length 3.
#' @return A character vector of length 3.
#' @examples
#' default.pvalueType
#' @export
default.pvalueType <- c("pos", "std", "neg")

#' Supported gene-gene interaction sources
#'
#' Character vector of the background gene-gene interaction (GGI) sources
#' supported by the `GGI_source` argument of [get_all_input_data()].
#'
#' @format A character vector of length 2.
#' @return A character vector of length 2.
#' @examples
#' default.GGI_source
#' @export
default.GGI_source <- c("Omnipath", "Reactome")

#' @noRd
snapshot_rng_state <- function() {
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  function() {
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old, envir = .GlobalEnv)
    }
  }
}

#' @noRd
combine_parallel <- function(res) {
  obj <- data.frame()
  for (i in seq_along(res)) {
    obj <- rbind(obj, as.data.frame(res[[i]]))
  }
  obj
}
