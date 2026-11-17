# Functions reused from Practical 5. Add documentation in your new package.
#' Encode nucleotide contexts as categorical features
#' 
#' Internal helper preserving the course model's position-specific encoding.
#' @param dna_strings Nonempty character vector of equal-length uppercase DNA
#'   strings containing A, T, C and G. For the supplied model use five bases
#'   with A at the centre. Inputs are not normalized or validated here.
#' @return A data.frame with one row per sequence and columns nt_pos1 through
#'   nt_posN, each a factor with levels A, T, C, G. Unknown bases become NA.
#' @keywords internal
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A status for multiple candidate sites
#' 
#' Apply the supplied fitted random forest without retraining or transforming
#' numeric features. Scores are the model's Positive-class probabilities.
#' @param ml_fit Fitted randomForest classification model, normally loaded from
#'   the bundled rf_fit.rds using system.file().
#' @param feature_df Nonempty data.frame containing gc_content (numeric),
#'   RNA_type (mRNA, lincRNA, lncRNA or pseudogene), RNA_region (CDS, intron,
#'   3'UTR or 5'UTR), exon_length and distance_to_junction (already log2
#'   transformed), evolutionary_conservation (numeric), and DNA_5mer
#'   (uppercase five-base DNA context with central A). Use the same feature
#'   definitions and preprocessing as training. Missing or unsupported inputs
#'   may cause prediction errors; no imputation is performed.
#' @param positive_threshold Numeric probability threshold, normally between
#'   zero and one. Only scores strictly greater than this value are Positive.
#' @return Input data.frame with RNA_type and RNA_region converted to factors,
#'   plus numeric predicted_m6A_prob and character predicted_m6A_status columns.
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                             package = "m6APrediction", mustWork = TRUE))
#' feature_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                   package = "m6APrediction", mustWork = TRUE))
#' prediction_multiple(ml_fit, feature_df[1:3, ], positive_threshold = 0.6)
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}

#' Predict m6A status for one candidate site
#' 
#' Construct a single feature row and apply prediction_multiple(). Numeric
#' inputs must already use the training feature scales.
#' @param ml_fit Fitted randomForest classification model from rf_fit.rds.
#' @param gc_content Numeric GC-content feature on the training scale.
#' @param RNA_type One of mRNA, lincRNA, lncRNA or pseudogene.
#' @param RNA_region One of CDS, intron, 3'UTR or 5'UTR.
#' @param exon_length Numeric exon length already log2-transformed.
#' @param distance_to_junction Numeric distance already log2-transformed.
#' @param evolutionary_conservation Numeric evolutionary conservation feature
#'   on the training scale.
#' @param DNA_5mer Uppercase five-base DNA context with central A, e.g. GGACA.
#' @param positive_threshold Numeric probability threshold (default 0.5).
#'   Equality is Negative; only a strictly greater score is Positive.
#' @return Named character vector of length two: predicted_m6A_prob (the score
#'   coerced to text) and predicted_m6A_status (Positive or Negative).
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                             package = "m6APrediction", mustWork = TRUE))
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA",
#'                   RNA_region = "CDS", exon_length = 12,
#'                   distance_to_junction = 5, evolutionary_conservation = 0.8,
#'                   DNA_5mer = "GGACA", positive_threshold = 0.5)
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
