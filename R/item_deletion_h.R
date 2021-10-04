#' Effect size of item deletion on selection accuracy
#' 
#' \code{item_deletion_h} computes the effect size of the impact of item bias on  
#' selection accuracy indices under strict vs. partial invariance using Cohen's  
#' h (1988).
#' 
#' @param propsel: proportion of selection. If missing, computed using `cut_z`
#' @param cut_z: pre-specified cutoff score on the observed composite. This 
#'        argument is ignored when `propsel` has input
#' @param weights_item: a vector of item weights
#' @param n_dim: number of dimensions, 1 by default. If the user does not supply 
#'        a different value, proceeds with the assumption that the scale is 
#'        unidimensional
#' @param n_i_per_dim: a vector containing the number of items in each 
#'        dimension; NULL by default. If the user provides a value for n_dim 
#'        that is > 1 but leaves n_i_per_dim = NULL, assumes that the subscales 
#'        have an equal number of items. 
#' @param weights_latent: a  vector of latent factor weights
#' @param alpha_r: a vector of latent factor mean for the reference group.
#' @param alpha_f: (optional) a vector of latent factor mean for the focal  
#'        group;if no input, set equal to alpha_r.
#' @param psi_r: a matrix of latent factor variance for the reference group.
#' @param psi_f: (optional) a matrix of latent factor variance for the focal  
#'        group; if no input, set equal to psi_r.
#' @param lambda_r_p: a matrix of factor loadings for the reference group under
#'        the partial invariance condition.
#' @param lambda_f_p: (optional) a matrix of factor loadings for the focal group 
#'        under the partial invariance condition; if no input, set equal to 
#'        lambda_r.
#' @param nu_r_p: a matrix of measurement intercepts for the reference group 
#'        under the partial invariance condition.
#' @param nu_f_p: (optional) a matrix of measurement intercepts for the focal 
#'        group under the partial invariance condition; if no input, set equal 
#'        to nu_r.
#' @param Theta_r_p: a matrix of the unique factor variances and covariances 
#'        for the reference group under the partial invariance condition.
#' @param Theta_f_p: (optional) a matrix of the unique factor variances and 
#'        covariances for the focal group under the partial invariance
#'        condition; if no input, set equal to Theta_r.
#' @param pmix_ref: Proportion of the reference group; 
#'        default to 0.5 (i.e., two populations have equal size)
#' @param plot_contour: logical; whether the contour of the two populations 
#'        should be plotted; default to TRUE
#' @param return_results: logical; whether the outputs from each call 
#'        of \code{PartInv_Multi_we} should also be returned as part of the
#'        returned object; default to FALSE
#' @return a list of the following two elements (or a list of the following 
#'        four elements if return_results == TRUE):
#'        - delta_table: a (8 x number of items) dataframe that stores Cohen's 
#'          h values for the effect size of the impact of item bias for each i
#'        - store_h: a list of length (number of items + 1) containing outputs 
#'          from \code{ref_acc_indices_h}. Each item in the list is a 
#'          restructured data frame comparing the various selection accuracy 
#'          indices for the reference group under the strict and partial 
#'          invariance conditions, and the corresponding h for each index. The
#'          first item in the list, 'full', contains the dataframe output for 
#'          when all items are included, and the remaining dataframes
#'          can be accessed by specifying i in '...$store_h$deleteitem_i'
#'  
#'          The following two lists of length (number of items + 1) are also 
#'          returned if the user indicates return_results == TRUE:
#'        
#'       - strict_results: contains outputs from \code{PartInvMulti_we} 
#'         under strict invariance.
#'       - partial_results: contains outputs from \code{PartInvMulti_we} 
#'         under partial invariance. 
#'         The first item in either list can be accessed through 
#'         '...$(partial/strict)_results$full' and contains the 
#'         \code{PartInvMulti_we} output when all items are included. Remaining 
#'         items in the list are of the 'deleteitem_i' form, and contain the 
#'         \code{PartInvMulti_we} output when a given item i is deleted
#'
#' @examples
#' # Multidimensional example
#' lambda_matrix <- matrix(0, nrow = 5, ncol = 2)
#' lambda_matrix[1:2, 1] <- c(.322, .655)
#' lambda_matrix[3:5, 2] <- c(.398, .745, .543)
#' multi_dim <- item_deletion_h(propsel = .05, n_dim = 5,
#'                             weights_item = c(1/4, 1/4, 1/6, 1/6, 1/6),
#'                             weights_latent = c(0.5, 0.5),
#'                             alpha_r = c(0, 0),
#'                             alpha_f = c(-0.3, 0.1),
#'                             psi_r = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
#'                             lambda_r_p = lambda_matrix,
#'                             nu_r_p = c(.225, .025, .010, .240, .125),
#'                             nu_f_p = c(.225, -.05, .240, -.025, .125),
#'                             theta_r_p = diag(1, 5),
#'                             theta_f_p = c(1, .95, .80, .75, 1), 
#'                             plot_contour = FALSE,
#'                             return_results = TRUE)
#' multi_dim$delta_table
#' multi_dim$store_h$full
#' multi_dim$store_h$deleteitem_1
#' multi_dim$strict_results$full
#' multi_dim$strict_results$deleteitem_1
#' multi_dim$partial_results$full
#' multi_dim$strict_results$deleteitem_1
#' 
#' # Single dimension example
#' single_dim <- item_deletion_h(propsel = .10,
#'                              weights_item = c(1, 0.9, 0.8, 1),
#'                              weights_latent = 0.9,
#'                              alpha_r = 0.5,
#'                              alpha_f = 0,
#'                              psi_r = 1,
#'                              lambda_r_p = c(.3, .5, .9, .7),
#'                              nu_r_p = c(.225, .025, .010, .240),
#'                              nu_f_p = c(.225, -.05, .240, -.025),
#'                              theta_r_p = diag(.96, 4),
#'                              n_dim = 1, plot_contour = FALSE, 
#'                              return_results = FALSE)
#' single_dim$delta_table
#' single_dim$store_h$full
#' single_dim$store_h$deleteitem_1
#' single_dim$strict_results$full
#' single_dim$strict_results$deleteitem_1
#' single_dim$partial_results$full
#' single_dim$strict_results$deleteitem_1  
#'
#' # If we specify return_results = FALSE
#' single_dim2 <- item_deletion_h(propsel = .10,
#'                                weights_item = c(1, 0.9, 0.8, 1),
#'                                weights_latent = 0.9,
#'                                alpha_r = 0.5,
#'                                alpha_f = 0,
#'                                psi_r = 1,
#'                                lambda_r_p = c(.3, .5, .9, .7),
#'                                nu_r_p = c(.225, .025, .010, .240),
#'                                nu_f_p = c(.225, -.05, .240, -.025),
#'                                theta_r_p = diag(.96, 4),
#'                                n_dim = 1, plot_contour = FALSE,
#'                                return_results = FALSE)
#' single_dim2$delta_table
#' single_dim2$store_h$full
#' single_dim2$store_h$deleteitem_1                     

item_deletion_h <- function(propsel, cut_z = NULL, 
                            weights_item, 
                            n_dim,
                            weights_latent,
                            alpha_r, alpha_f = alpha_r,
                            psi_r, psi_f = psi_r,
                            lambda_r_p, lambda_f_p = lambda_r_p, 
                            nu_r_p, nu_f_p = nu_r_p,
                            theta_r_p, theta_f_p = theta_r_p,
                            pmix_ref = 0.5, 
                            plot_contour = TRUE,
                            return_results = FALSE,
                            ...) {

  n_items <- length(weights_item)
  # Start with pre-allocating space:
  store_str <- store_par <- store_h <- vector(mode = "list", n_items + 1)
  delta_table <- as.data.frame(matrix(nrow = 8, ncol = n_items))
  
  # Include all items first
  # Under strict invariance:
  store_str[[1]] <- PartInvMulti_we(propsel, cut_z = cut_z, 
                                    weights_item, weights_latent,
                                    alpha_r = alpha_r,
                                    alpha_f = alpha_f,
                                    psi_r = psi_r,
                                    psi_f = psi_f,
                                    lambda_r = lambda_f_p * (1 - pmix_ref) + 
                                               lambda_r_p * pmix_ref,
                                    nu_r = nu_f_p * (1 - pmix_ref) + 
                                           nu_r_p * pmix_ref,
                                    Theta_r = theta_f_p * (1 - pmix_ref) + 
                                              theta_r_p * pmix_ref,
                                    pmix_ref = pmix_ref, 
                                    plot_contour = plot_contour)
  # Under partial invariance:
  store_par[[1]] <- PartInvMulti_we(propsel, cut_z = cut_z, 
                                    weights_item, weights_latent,
                                    alpha_r = alpha_r,
                                    alpha_f = alpha_f,
                                    psi_r = psi_r,
                                    psi_f = psi_f,
                                    lambda_r = lambda_r_p,
                                    nu_r = nu_r_p,
                                    nu_f = nu_f_p,
                                    Theta_r = theta_r_p,
                                    Theta_f = theta_f_p,
                                    pmix_ref = pmix_ref, 
                                    plot_contour = plot_contour) 
  
  # Compare accuracy indices for reference group under strict vs. partial
  # invariance conditions and compute h for full item set:
  store_h[[1]] <- ref_acc_indices_h(store_str[[1]], store_par[[1]]) 
  # Delete one item at a time by assigning 0 to index i and redistributing
  # weights for the subscale; populate store_str and store_par
  for (i in seq_len(length(weights_item) + 1)[-1]) {
    # Strict invariance
    take_one_out <- del_redist_weights(weights_item, n_dim = n_dim, del_i = i-1)
    store_str[[i]] <- PartInvMulti_we(propsel, cut_z = cut_z, 
                                      take_one_out, weights_latent,
                                      alpha_r = alpha_r,
                                      alpha_f = alpha_f,
                                      psi_r = psi_r,
                                      psi_f = psi_f,
                                      lambda_r = lambda_f_p * (1 - pmix_ref) + 
                                                 lambda_r_p * pmix_ref,
                                      nu_r = nu_f_p * (1 - pmix_ref) +
                                             nu_r_p * pmix_ref,
                                      Theta_r = theta_f_p * (1 - pmix_ref) + 
                                                theta_r_p * pmix_ref,
                                      pmix_ref = pmix_ref, 
                                      plot_contour = plot_contour)
    # Partial invariance
    store_par[[i]] <- PartInvMulti_we(propsel, cut_z = cut_z,
                                      take_one_out, weights_latent,
                                      alpha_r = alpha_r,
                                      alpha_f = alpha_f,
                                      psi_r = psi_r,
                                      psi_f = psi_f,
                                      lambda_r = lambda_r_p,
                                      nu_r = nu_r_p,
                                      nu_f = nu_f_p,
                                      Theta_r = theta_r_p,
                                      Theta_f = theta_f_p,
                                      pmix_ref = pmix_ref, 
                                      plot_contour = plot_contour)
                    
    # Compare accuracy indices for reference group under strict vs. partial
    # invariance conditions and compute h for each index i
    store_h[[i]] <- ref_acc_indices_h(store_str[[i]], store_par[[i]]) 
    
    # Compute Cohen's h for the effect of item bias for each i
    delta_table[i - 1] <- delta_h(store_h[[1]]$h, store_h[[i]]$h)
  }
  
  # Format stored variables
  colnames <- paste0(c("full", rep("deleteitem_", n_items)), c("", 1:n_items))
  names(store_str) <- names(store_par) <- names(store_h) <- colnames
  names(delta_table) <- paste0(rep(paste0('\u0394', "h(-")),
                               c(1:n_items), c(rep(")")))
  rownames(delta_table) <- c("TP", "FP", "TN", "FN", "PS", "SR", "SE", "SP")
  
  if (return_results == TRUE) {
    return(list("delta_table" = delta_table, 
                "store_h" = store_h,
                "strict_results" = store_str, 
                "partial_results" = store_par))
  } else {
    return(list("delta_table" = delta_table, 
                "store_h" = store_h))
  }
}


# Delete item i and redistribute its weight within subscale

# \code{del_redist_weights} replaces the item weight with 0 for the item to be
# deleted, and redistributes this item's weight across the remaining items.

#' @param weights_item: a vector of item weights
#' @param n_dim: number of dimensions, 1 by default. If the user does not supply 
#'        a different value, proceeds with the assumption that the scale is 
#'        unidimensional
#' @param n_i_per_dim: a vector containing the number of items in each 
#'        dimension; NULL by default. If the user provides a value for n_dim 
#'        that is > 1 but leaves n_i_per_dim = NULL, assumes that the subscales 
#'        have an equal number of items. 
#' @param del_i: index of the item to be deleted
#' 
#' @return take_one_out: weights vector with redistributed weights
#' @examples
#' one_dim_weights <- c(1:7)
#' del_redist_weights(one_dim_weights, del_i = 2)
#' one_dim_weights2 <- c(1:7)
#' del_redist_weights(one_dim_weights2, n_dim = 1, n_i_per_dim = 7, del_i = 2)
#' multi_equal_len_weights <- c(1:9)
#' del_redist_weights(multi_equal_len_weights, n_dim = 3, del_i = 2)
#' multi_equal_len_weights2 <- c(1:9)
#' del_redist_weights(multi_equal_len_weights2, n_dim = 3, 
#'                    n_i_per_dim = c(3, 3, 3), del_i = 2)
#' multi_unequal_len_weights <- c(1:12)
#' del_redist_weights(multi_unequal_len_weights, n_dim = 3, 
#'                    n_i_per_dim = c(3, 6, 3), del_i = 2)
#' error_ex <- c(1:12)
#' del_redist_weights(error_ex, n_dim = -3, 
#'                    n_i_per_dim = c(3, 6, 3), del_i = 2)

del_redist_weights <- function(weights_item, n_dim, n_i_per_dim,  n_i_per_dim = NULL,
                               del_i) {
  n_items <- length(weights_item)
  take_one_out <- weights_item; take_one_out[del_i] <- 0
  
  # Unidimensional
  if ((n_dim == 1) & (is.null(n_i_per_dim) | length(n_i_per_dim) == 1)) {
    take_one_out <- take_one_out / (n_items - 1) * n_items
    
  # Multidimensional, equal n  
  } else if ((n_dim > 1) & is.null(n_i_per_dim)) { 
    subscale_len <- n_items / n_dim # subscale length
    # Split indices into dimensions
    i_by_dim <- split(1:n_items, cut(seq_along(1:n_items), n_dim, 
                                     labels = FALSE))
    for(k in 1:n_dim) { # Loop through the number of dimensions
      if(del_i %in% i_by_dim[[k]]) { # If del_i is in dimension k
        # Create temporary vector to store the remaining indices in dimension k
        temp_i <- i_by_dim[[k]][i_by_dim[[k]] != del_i] 
        for(j in temp_i) { # Re-weight the remaining indices in the subscale
          take_one_out[j] <- take_one_out[j] / (subscale_len - 1) * subscale_len
        }
      }
    }
  # Multidimensional, unequal n
  } else if ((n_dim > 1) & !is.null(n_i_per_dim)) { 
    # Split indices into dimensions
    i_by_dim <- split(1:n_items, cut(seq_along(1:n_items), 
                                       breaks = cumsum(c(0, n_i_per_dim)), 
                                       labels = FALSE))
    for(k in 1:n_dim) { #Loop through the number of dimensions
      if(del_i %in% i_by_dim[[k]]){ # If del_i is in dimension k
        subscale_len <- n_i_per_dim[k]
        # Create temporary vector to store the remaining indices in dimension k
        temp_i <- i_by_dim[[k]][i_by_dim[[k]] != del_i] 
        for(j in temp_i) { # Re-weight the remaining indices in the subscale
          take_one_out[j] <- take_one_out[j] / (subscale_len - 1) * subscale_len
        }
      }
    }
  } else {
    stop('Check n_dim and n_i_per_dim')
  }
  return(take_one_out)
}


#' Computes Cohen's h (Cohen, 1988) for the difference in two proportions:
#'  $h = 2arcsin(\sqrt{p1}) - 2arcsin(\sqrt{p2})$ 
#' 
#' @param p1 The first proportion.
#' @param y The second proportion.
#' @return The Cohen's h value.
#' @examples
#' cohens_h(0.7, 0.75)
#' cohens_h(0.3, 0.4)

cohens_h <- function(p1, p2) {
  h <-  2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))
  return(h)
}

#' Compute effect size for the impact of item deletion

#' Uses the formula below to compute the effect size for impact of item bias by
#' comparing Cohen's h values for a given selection accuracy index when an item
#' is deleted vs. included e.g. improvement in SE when item i is deleted:
#' $\Delta h_{SE^{(-i)}}=\text{sign}(h_{SE^{(R)}}-h_{SE^{(-i)}})
#'           \left|\left|h_{SE^{(-i)}}\right|\left|h_{SE^{(-i)}}\right|\right|$

#' @param h_R: h effect sizes for when the item is included
#' @param h_i_del: h effect sizes for when the item is deleted

#' @return Cohen's h for the difference in the selection accuracy index when the
#' item is deleted
delta_h <- function(h_R, h_i_del) {
  sign(h_R - h_i_del) * abs(abs(h_i_del) - abs(h_R))
}


#' Selection accuracy indices and Cohen's h for the reference group under strict
#' and partial invariance
#' 
#' Takes in the outputs from the PartInv_Multi_we function for the strict 
#' invariance and partial invariance conditions, and returns a restructured 
#' data frame comparing the various indices (P(PS),..., SE, SP etc.) for 
#' the reference group under the strict and partial invariance conditions, and 
#' the corresponding h for each of the indices.

#' @param strict_output: \code{PartInv_Multi_we} output (a list) under strict 
#'        invariance 
#' @param partial_output: \code{PartInv_Multi_we} output (a list) under partial
#'        invariance
#' @return A 8x3 dataframe, columns for strict invariance, partial invariance,
#'        and h.
ref_acc_indices_h <- function(strict_output, partial_output) {
  ref_par_strict <- partial_output[4]$summary[1][, 1]
  ref_strict <- strict_output[4]$summary[1][, 1]
  r_names <- c("A (true positive)", "B (false positive)", "C (true negative)", 
               "D (false negative)", "Proportion selected", "Success ratio", 
               "Sensitivity", "Specificity")
  df <- data.frame(strict_invariance =  ref_strict, 
                   partial_invariance = ref_par_strict, row.names = r_names)
  df["h"] <- round(cohens_h(df$strict_invariance, df$partial_invariance), 3)
  
  return(df)
}
